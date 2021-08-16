# # Single molecules centered at region of interest
# ## Annie Maslan, modified by Kousik Sundararajan

# Input: bed file of coordinates where single molecules should be centered, mod_mappings.bam, mod_mappings.bam.bai

# 1. Input bed file with windows over which to extract and center reads (e.g. CTCF sites +/- 1kb)
# 2. Extract reads within the windows
# 3. Plot modified bases within the windows colored by probability of modification


import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import colors

import scipy as sp
from scipy import signal
from scipy.optimize import curve_fit

import pysam

import multiprocessing
from joblib import Parallel, delayed
from pybedtools import BedTool


# ## Extracting methylation data


class Region(object):
	def __init__(self, region):
		self.chromosome = region[1][0]
		self.begin = region[1][1]
		self.end = region[1][2]
		self.size = self.end - self.begin
		self.string = f"{self.chromosome}_{self.begin}_{self.end}"

class Methylation(object):
	def __init__(self, table, data_type, name, called_sites):
		self.table = table
		self.data_type = data_type
		self.name = name
		self.called_sites = called_sites

def make_windows(bed):
	reg = []
	for row in bed.iterrows():
		reg.append(Region(row))
	return reg

def get_data(methylation_file, name, windows, thresh, num_cores, window_size):
	"""
	Import methylation data from all files in the list methylation_files
	Data in bam format
	data is extracted within the window
	"""
	# data is for each window the methylation data and the position data
	meth_data = Parallel(n_jobs=num_cores)(delayed(parse_ont_bam)(methylation_file, name, w, thresh, window_size) for w in windows)
	print('methyltation data collection done')
	base_pos_data = Parallel(n_jobs=num_cores)(delayed(parse_ont_bam_base_pos)(methylation_file, name, w, thresh, window_size) for w in windows)
	print('base abundance data collection done')
	return meth_data, base_pos_data

#############
def parse_ont_bam_base_pos(filename, name, window, thresh, window_size):
	bam = pysam.AlignmentFile(filename, "rb")
	base_pos_data = []
	for read in bam.fetch(reference=window.chromosome, start=window.begin, end=window.end):
		#if not read.is_supplementary and not read.is_secondary:
		[(mod, base_positions), (mod2, base_positions2)] = get_reference_positions(read, window)
		for base_pos in base_positions:
			if base_pos is not None:
				if abs(base_pos) <= window_size: # to decrease memory, only store bases within the window
					base_pos_data.append((read.query_name,
						'-' if read.is_reverse else '+',
						base_pos,
						mod))
		for base_pos in base_positions2:
			if base_pos is not None:
				if abs(base_pos) <= window_size: # to decrease memory, only store bases within the window
					base_pos_data.append((read.query_name,
						'-' if read.is_reverse else '+',
						base_pos,
						mod2))
	pos_data_return = pd.DataFrame(base_pos_data, columns=['read_name', 'strand', 'base_pos', 'mod']).astype(dtype={'mod': 'category', 'base_pos': 'int32'})
	return pos_data_return

def get_reference_positions(read, window):
	mod1 = read.get_tag('Mm').split(';')[0].split(',', 1)[0] # all reads have these tags
	mod2 = read.get_tag('Mm').split(';')[1].split(',', 1)[0] # all reads have these tags
	# don't require there be modified bases
	mod1_return = get_pos_base(read, mod1, 0, window)
	mod2_return = get_pos_base(read, mod2, 1, window)
	return (mod1_return, mod2_return)

	# if (read.has_tag('Mm')) & (';' in read.get_tag('Mm')): # all read should have this
	# 	mod1 = read.get_tag('Mm').split(';')[0].split(',', 1)[0]
	# 	mod2 = read.get_tag('Mm').split(';')[1].split(',', 1)[0]
	# 	if mod1 is not None: # don't require there be modified bases
	# 		mod1_return = get_pos_base(read, mod1, 0, window)
	# 	else:
	# 		print('mod1 none')
	# 		mod1_return = (None, [None])
	# 	if mod2 is not None: # don't require there be modified bases
	# 		mod2_return = get_pos_base(read, mod2, 1, window)
	# 		return (mod1_return, mod2_return)
	# 	else:
	# 		print('mod2 none')
	# 		return (mod1_return, (None, [None]))
	# else:
	# 	print('missing tag')
	# 	return ((None, [None]), (None, [None]))

def get_pos_base(read, basemod, index, window):
	if '-' in basemod:
		sys.exit("ERROR: modifications on negative strand currently unsupported.")
	if 'A' not in basemod: #if 'A+a' not in basemod:
		if 'C' not in basemod: #if 'C+m' not in basemod
			return (None, [None])
	seq = read.get_forward_sequence()
	base, mod = basemod.split('+')
	base_index = np.array([i for i, letter in enumerate(seq) if letter == base])
	refpos = np.array(read.get_reference_positions(full_length=True))
	if read.is_reverse:
		refpos = np.flipud(refpos)
	base_keep = [] # to track A or CG abundance
	# deal with None for refpos from soft clipped / unaligned bases
	if 'C' in basemod: #if 'C+m' in basemod:
		for b in base_index:
			if b < len(seq) - 1:
				if (refpos[b] is not None) & (refpos[b+1] is not None):
					if seq[b + 1] == 'G':
						if abs(refpos[b+1] - refpos[b]) == 1: # ensure there isn't a gap
							base_keep.append(b)
	# for m6A no need to look at neighboring base; do need to remove refpos that are None
	else:
		for b in base_index:
			if refpos[b] is not None:
				base_keep.append(b)
	return (basemod, np.array(refpos[base_keep]) - round(((window.end-window.begin)/2 + window.begin)))

#############
def parse_ont_bam(filename, name, window, thresh, window_size):
	'''
	parse mod_mappings.bam file to create methylation object with read_name, strand, pos, quality, and mod 
	in window above threshold specified
	'''
	bam = pysam.AlignmentFile(filename, "rb")
	data = []
	for read in bam.fetch(reference=window.chromosome, start=window.begin, end=window.end):
		#if not read.is_supplementary and not read.is_secondary:
		[(mod, positions, quals), (mod2, positions2, quals2)] = get_modified_reference_positions(read, window)
		for pos, qual in zip(positions, quals):
			if pos is not None:
				if abs(pos) <= window_size: # to decrease memory, only store bases within the window
					if qual >= 0: # decide whether to include quality filter here; include for less memory use
						data.append((read.query_name,
							'-' if read.is_reverse else '+',
							pos,
							qual,
							mod))
		for pos, qual in zip(positions2, quals2):
			if pos is not None:
				if abs(pos) <= window_size: # to decrease memory, only store bases within the window
					if qual >= 0: # decide whether to include quality filter here; include for less memory use
						data.append((read.query_name,
							'-' if read.is_reverse else '+',
							pos,
							qual,
							mod2))
		

	data_return = Methylation(
		table=pd.DataFrame(data, columns=['read_name', 'strand', 'pos', 'quality', 'mod'])
				.astype(dtype={'mod': 'category', 'pos': 'int32', 'quality': 'int32'})# .astype(dtype={'mod': 'category', 'quality': 'float'})
				.sort_values(['read_name', 'pos']),
		data_type="ont-bam",
		name=name,
		called_sites=len(data))
	
	
	return data_return


def get_modified_reference_positions(read, window):
	'''
	extract mA and mC pos & prob information for the read
	'''
	if (read.has_tag('Mm')) & (';' in read.get_tag('Mm')):
		mod1 = read.get_tag('Mm').split(';')[0].split(',', 1)[0]
		mod2 = read.get_tag('Mm').split(';')[1].split(',', 1)[0]
		mod1_list = read.get_tag('Mm').split(';')[0].split(',', 1)
		mod2_list = read.get_tag('Mm').split(';')[1].split(',', 1)
		if len(mod1_list) > 1:
			mod1_return = get_pos_prob(read, mod1, 0, window)
		else:
			mod1_return = (None, [None], [None])
		if len(mod2_list) > 1:
			mod2_return = get_pos_prob(read, mod2, 1, window)
			return (mod1_return, mod2_return)
		else:
			return (mod1_return, (None, [None], [None]))
	else:
		return ((None, [None], [None]), (None, [None], [None]))


def get_pos_prob(read, basemod, index, window):
	'''
	get (modified base type, position of modified base, probability of modification)
	'''
	if '-' in basemod:
		sys.exit("ERROR: modifications on negative strand currently unsupported.")
	if 'A' not in basemod: #if 'A+a' not in basemod:
		if 'C' not in basemod: #if 'C+m' not in basemod
			return (None, [None], [None])
	base, mod = basemod.split('+')
	deltas = [int(i) for i in read.get_tag('Mm').split(';')[index].split(',')[1:]]
	num_base = len(read.get_tag('Mm').split(';')[index].split(','))-1
	Ml = read.get_tag('Ml')
	if index == 0:
		probabilities = np.array(Ml[0:num_base],dtype=int)
	if index == 1:
		probabilities = np.array(Ml[0-num_base:],dtype=int)
	base_index = np.array([i for i, letter in enumerate(read.get_forward_sequence()) if letter == base])
	# determine locations of the modified bases, where index_adj is the adjustment of the base_index
	# based on the cumulative sum of the deltas
	locations = np.cumsum(deltas)
	# loop through locations and increment index_adj by the difference between the next location and current one + 1
	# if the difference is zero, therefore, the index adjustment just gets incremented by one because no base should be skipped
	index_adj = []
	index_adj.append(locations[0])
	i = 0
	for i in range(len(locations) - 1):
		diff = locations[i+1] - locations[i]
		index_adj.append(index_adj[i] + diff + 1)
	# get the indices of the modified bases
	index_adj = index_adj[0:len(base_index)]
	modified_bases = base_index[index_adj]
	refpos = np.array(read.get_reference_positions(full_length=True))
	if read.is_reverse:
		refpos = np.flipud(refpos)
		probabilities = probabilities[::-1]	
	# extract CpG sites only rather than all mC
	keep = []
	prob_keep = []
	i = 0
	seq = read.get_forward_sequence()
	# deal with None for refpos from soft clipped / unaligned bases
	if 'C' in basemod: #if 'C+m' in basemod:
		for m in modified_bases:
			if m < len(seq) - 1: # if modified C is not the last base in the read	
				if (refpos[m] is not None) & (refpos[m+1] is not None):
					if seq[m + 1] == 'G':
						if abs(refpos[m+1] - refpos[m]) == 1: # ensure there isn't a gap
							keep.append(m)
							prob_keep.append(i)
			i=i+1
	# for m6A no need to look at neighboring base; do need to remove refpos that are None
	else:
		for m in modified_bases:
			if refpos[m] is not None:
				keep.append(m)
				prob_keep.append(i)
			i=i+1
	# adjust position to be centered at 0 at the center of the window; round in case is at 0.5
	# add returning base_index for plotting mod/base_abundance
	return (basemod, np.array(refpos[keep]) - round(((window.end-window.begin)/2 + window.begin)), probabilities[prob_keep])


def remove_overlapping(bed, window_size, out):
	'''
	remove regions of interest that are overlapping in the window that will be plotted
	'''
	b = BedTool(bed)
	b.sort().saveas(out + '/tmp.sorted.bed');
	bed = pd.read_csv(out + '/tmp.sorted.bed', sep='\t', header=None)
	# find middle
	bed['middle'] = (bed[1] + bed[2]) / 2
	# assign left end of window
	bed['window_left'] = bed['middle'] - window_size
	# assign right end of window
	bed['window_right'] = bed['middle'] + window_size
	# rows to keep or discard
	keep = [0] # keep first entry
	discard = []
	# if window overlaps, remove that row
	for i in range(0, len(bed) - 1):
		if bed.iloc[i+1].window_left < bed.iloc[i].window_right:
			discard.append(i+1)
		else:
			keep.append(i+1)
	return bed.iloc[keep]


def make_plots(all_data, window_size, out, name, smooth, all_base_pos_data, thresh):
	'''
	input: dataframe with all the data for plotting, window size, output directory, name for plot
	calls make_mod_plot to make plot for each modification type
	'''

	cmapA = colors.ListedColormap(['white', '#053C5E'])
	cmapC = colors.ListedColormap(['white', '#BB4430'])

	mean_rolling_mA = make_mod_plot(all_data, 'A', cmapA, window_size, out, name, smooth, all_base_pos_data, thresh) #make_mod_plot(all_data, 'A+a', 'Blues', window_size, out, name, smooth)
	#mean_rolling_mC = make_mod_plot(all_data, 'C', cmapC, window_size, out, name, smooth, all_base_pos_data, thresh) #make_mod_plot(all_data, 'C+m', 'Reds', window_size, out, name, smooth)

	return mean_rolling_mA, mean_rolling_mC


def make_mod_plot(all_data, mod, color, window_size, out, name, smooth, all_base_pos_data, thresh):
	'''
	input: dataframe with all the data for plotting, modification type, cmap for heatmap, window size, output directory, name for plot
	output: heatmap saved as png
	'''
	# plot mA for single molecules
	all_data_mod = all_data[all_data['mod'].str.contains(mod)] # all_data[all_data['mod'] == mod]
	all_data_mod = all_data_mod[all_data_mod['pos'] >= - window_size]
	all_data_mod = all_data_mod[all_data_mod['pos'] <= window_size]
	#all_data_mod = all_data_mod[all_data_mod['quality'] > thresh]
	# make binary for calls above threshold
	#all_data_mod['quality'] = 1
	all_data_pivoted_mod = pd.pivot_table(all_data_mod, values = 'quality', columns = 'pos', index='read_name')

	all_base_pos_data_mod = all_base_pos_data[all_base_pos_data['mod'].str.contains(mod)]
	all_base_pos_data_mod = all_base_pos_data_mod[all_base_pos_data_mod['base_pos'] >= - window_size]
	all_base_pos_data_mod = all_base_pos_data_mod[all_base_pos_data_mod['base_pos'] <= window_size]
	# assign 1 for all positions where base is there
	all_base_pos_data_mod['presence'] = 1

	all_base_pos_data_pivoted_mod = pd.pivot_table(all_base_pos_data_mod, values = 'presence', columns='base_pos', index='read_name')

	# fill in missing positions with 0 for heatmap because x-axis is non-numeric so need position at each bp
	# sort by columns so is from -window_size to +window_size
	r = range(-window_size, window_size+1, 1)
	for bp in r:
		if bp not in all_data_pivoted_mod.columns:
			all_data_pivoted_mod[bp] = np.nan
		if bp not in all_base_pos_data_pivoted_mod.columns:
			all_base_pos_data_pivoted_mod[bp] = np.nan
	all_data_pivoted_mod = all_data_pivoted_mod.sort_index(axis=1)
	all_base_pos_data_pivoted_mod = all_base_pos_data_pivoted_mod.sort_index(axis=1)

	#write df for clustering
	outfile_prob = open(out + '/' + name + '_'+ mod + '_prob.csv', 'w')
	print('writing file ' + out + '/' + name + '_'+ mod + '_prob.csv')
	all_data_pivoted_mod.to_csv(outfile_prob)
	outfile_prob.close()

	outfile_pos = open(out + '/' + name + '_'+ mod + '_pos.csv', 'w')
	print('writing file ' + out + '/' + name + '_'+ mod + '_pos.csv')
	all_base_pos_data_pivoted_mod.to_csv(outfile_pos)
	outfile_pos.close()

	return

def main():
	#insert filepath to input hybrid bam file (from Guppy + Winnowmap), alternately, use Megalodon bam file
	#bams = [<filepath_1>,<filepath_2>, etc]
	bams = ["/scratch/groups/astraigh/minion_seq/guppy_winnow_merge/DiMeLo_cen_enrich_merge_2021.06.19/CA/CA_merge.sorted_q10.bam"]

	# insert output file directory
	out = <output_directory_path>

	# Threshold of p(mA) for including 'A' in output csv, 0-255
	#set to 0 to read out all A's. Will need to do this if using A position for coverage information
	thresh = 0

	# will plot non-overlapping features in -window_size to +window_size
	window_size = 150000

	# name for plots; 4 per target
	#names = ["free", "CA","IgG","unt"]
	names = ["CA"]

	#bed file with chromosome region of interest
	bed = "/home/groups/astraigh/ont_minion/single_molecule_roi/beds/HG002_ChrX_dip_200kb_5745_5775.bed"
	## This example bed file looks like this:
	'''
	chrX	57450000	57750000
	'''

	#### end of parameters and run-specific file paths ####

	bed_no_overlap = remove_overlapping(bed, window_size, out)
	quarts = [bed_no_overlap]

	num_cores = multiprocessing.cpu_count()
	print('executing with ' + str(num_cores) + ' cores')

	i = 0
	for bam in bams:

		mod_mean_rolling_list_mA = []
		mod_mean_rolling_list_mC = []

		for q in quarts:
			# make windows 
			windows = make_windows(q)

			# get methylation data within all windows in parallelized way
			meth_data, base_pos_data = get_data(bam, names[i], windows, thresh, num_cores, window_size)

			# combine all methylation tables into a single dataframe
			list_tables = []
			for m in meth_data:
				list_tables.append(m.table)
			all_data = pd.concat(list_tables)

			base_list_tables = []
			for b in base_pos_data:
				base_list_tables.append(b)
			all_base_pos_data = pd.concat(base_list_tables)


			print('processing ' + str(len(all_data['read_name'].unique())) + ' methylated reads for ' + names[i] + ' for bam: ' + bam)
			print('processing ' + str(len(all_base_pos_data['read_name'].unique())) + ' total reads for ' + names[i] + ' for bam: ' + bam)

			# make heatmaps for mA and CpG methylation
			make_plots(all_data, window_size, out, names[i], smooth, all_base_pos_data, thresh)

		i += 1



if __name__ == '__main__':
	main()




