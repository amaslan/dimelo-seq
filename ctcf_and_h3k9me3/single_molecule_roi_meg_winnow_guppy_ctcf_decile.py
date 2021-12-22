# # Single molecules centered at region of interest
# # Here analysis is by decile
# ## Annie Maslan
# ## 06.04.21

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

from mpl_toolkits.axes_grid1 import make_axes_locatable


# ## Extracting methylation data


class Region(object):
	def __init__(self, region):
		self.chromosome = region[1][0]
		self.begin = region[1][1]
		self.end = region[1][2]
		self.size = self.end - self.begin
		self.string = f"{self.chromosome}_{self.begin}_{self.end}"
		self.strand = region[1][3] # store strand of the motif here

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
	pos_data_return = pd.DataFrame(base_pos_data, columns=['read_name', 'strand', 'base_pos', 'mod']).astype(dtype={'mod': 'category', 'base_pos': 'int16'})
	return pos_data_return

def get_reference_positions(read, window):
	mod1 = read.get_tag('Mm').split(';')[0].split(',', 1)[0] # all reads have these tags
	mod2 = read.get_tag('Mm').split(';')[1].split(',', 1)[0] # all reads have these tags
	# don't require there be modified bases
	mod1_return = get_pos_base(read, mod1, 0, window)
	mod2_return = get_pos_base(read, mod2, 1, window)
	return (mod1_return, mod2_return)

def get_pos_base(read, basemod, index, window):
	if '-' in basemod:
		sys.exit("ERROR: modifications on negative strand currently unsupported.")
	if 'A' not in basemod: 
		if 'C' not in basemod: 
			return (None, [None])
	seq = read.get_forward_sequence()
	base, mod = basemod.split('+')
	base_index = np.array([i for i, letter in enumerate(seq) if letter == base])
	refpos = np.array(read.get_reference_positions(full_length=True))
	if read.is_reverse:
		refpos = np.flipud(refpos)
	base_keep = [] # to track A or CG abundance
	# deal with None for refpos from soft clipped / unaligned bases
	if 'C' in basemod: 
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
	# perform strand adjustment for -
	# center at center of the motif
	if window.strand == '+':
		refpos_adjusted = np.array(refpos[base_keep]) - round(((window.end-window.begin)/2 + window.begin))
	if window.strand == '-':
		refpos_adjusted = -1*(np.array(refpos[base_keep]) - round(((window.end-window.begin)/2 + window.begin)))
	return (basemod, refpos_adjusted)

#############
def parse_ont_bam(filename, name, window, thresh, window_size):
	'''
	parse mod_mappings.bam file to create methylation object with read_name, strand, pos, quality, and mod 
	in window above threshold specified
	'''
	bam = pysam.AlignmentFile(filename, "rb")
	data = []
	for read in bam.fetch(reference=window.chromosome, start=window.begin, end=window.end):
		[(mod, positions, quals), (mod2, positions2, quals2)] = get_modified_reference_positions(read, window)
		for pos, qual in zip(positions, quals):
			if pos is not None:
				if abs(pos) <= window_size: # to decrease memory, only store bases within the window
					if qual >= thresh: # decide whether to include quality filter here; include for less memory use
						data.append((read.query_name,
							'-' if read.is_reverse else '+',
							pos,
							qual,
							mod))
		for pos, qual in zip(positions2, quals2):
			if pos is not None:
				if abs(pos) <= window_size: # to decrease memory, only store bases within the window
					if qual >= thresh: # decide whether to include quality filter here; include for less memory use
						data.append((read.query_name,
							'-' if read.is_reverse else '+',
							pos,
							qual,
							mod2))
		

	data_return = Methylation(
		table=pd.DataFrame(data, columns=['read_name', 'strand', 'pos', 'quality', 'mod'])
				.astype(dtype={'mod': 'category', 'pos': 'int16', 'quality': 'int16'})
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
	if 'A' not in basemod: 
		if 'C' not in basemod: 
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
	modified_bases = base_index[index_adj]
	refpos = np.array(read.get_reference_positions(full_length=True))
	if read.is_reverse:
		refpos = np.flipud(refpos)
		# probabilities = probabilities[::-1]	
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
	# adjust position to be centered at 0 at the center of the motif; round in case is at 0.5
	# add returning base_index for plotting mod/base_abundance
	if window.strand == '+':
		refpos_adjusted = np.array(refpos[keep]) - round(((window.end-window.begin)/2 + window.begin))
	if window.strand == '-':
		refpos_adjusted = -1*(np.array(refpos[keep]) - round(((window.end-window.begin)/2 + window.begin)))
	return (basemod, refpos_adjusted, probabilities[prob_keep])


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
	

def resolution_analysis(all_data, all_base_pos_data, window_size, out, name):
	'''
	fit exponential decay curve to signal in q4 peaks; also confidence in mapping peak center
	instead print to csv to analyze in jupyter notebook
	'''
	all_data_mA = all_data[all_data['mod'].str.contains('A')] #all_data[all_data['mod'] == 'A+a']
	all_data_mA = all_data_mA[all_data_mA['pos'] >= - window_size]
	all_data_mA = all_data_mA[all_data_mA['pos'] <= window_size]

	all_data_mA.to_csv(out + '/' + name + '_all_data_mA.csv', index=False) # save as csv for resolution analysis


def main():
	#### start of parameters and run-specific file paths ####
	bams = ["prod_ctcf_mod_mappings_merge.sorted.bam",
	"prod_IgG_mod_mappings.sorted.bam",
	"prod_untreated_mod_mappings.sorted.bam"]

	# output file directory
	out = "results"

	# 0-255
	thresh = 0 

	# will plot non-overlapping features in -window_size to +window_size
	window_size = 1000

	# number of bases to smooth over for moving average curve
	smooth = 50 

	# name for plots; 10 per target
	# always have q10-1 in the name
	names = ["meg_CTCF_CTCF_q10", "meg_CTCF_CTCF_q9", "meg_CTCF_CTCF_q8", "meg_CTCF_CTCF_q7",
	"meg_CTCF_CTCF_q6", "meg_CTCF_CTCF_q5", "meg_CTCF_CTCF_q4", "meg_CTCF_CTCF_q3",
	"meg_CTCF_CTCF_q2", "meg_CTCF_CTCF_1",
	"meg_IgG_CTCF_q10", "meg_IgG_CTCF_q9", "meg_IgG_CTCF_q8", "meg_IgG_CTCF_q7",
	"meg_IgG_CTCF_q6", "meg_IgG_CTCF_q5", "meg_IgG_CTCF_q4", "meg_IgG_CTCF_q3",
	"meg_IgG_CTCF_q2", "meg_IgG_CTCF_q1",
	"meg_untreated_CTCF_q10", "meg_untreated_CTCF_q9", "meg_untreated_CTCF_q8", "meg_untreated_CTCF_q7",
	"meg_untreated_CTCF_q6", "meg_untreated_CTCF_q5", "meg_untreated_CTCF_q4", "meg_untreated_CTCF_q3",
	"meg_untreated_CTCF_q2", "meg_untreated_CTCF_q1"]

	bed = "intersection.motifs.chip.formatted.chm13.bed" 

	#### end of parameters and run-specific file paths ####

	bed_no_overlap = remove_overlapping(bed, window_size, out)
	# signalValue is in column 5 (index 4) of the motif-based bed file
	bed_no_overlap_sorted = bed_no_overlap.sort_values(by=bed_no_overlap.columns[4], ascending=False)

	# make quartiles by signalValue
	# signalValue is in column 5 (index 4) of the motif-based bed file
	quants = bed_no_overlap_sorted[4].quantile(q=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
	q1 = bed_no_overlap_sorted[bed_no_overlap_sorted[4] <= quants[0.1]]
	q2 = bed_no_overlap_sorted[(bed_no_overlap_sorted[4] > quants[0.1]) & (bed_no_overlap_sorted[4] <= quants[0.2])]
	q3 = bed_no_overlap_sorted[(bed_no_overlap_sorted[4] > quants[0.2]) & (bed_no_overlap_sorted[4] <= quants[0.3])]
	q4 = bed_no_overlap_sorted[(bed_no_overlap_sorted[4] > quants[0.3]) & (bed_no_overlap_sorted[4] <= quants[0.4])]
	q5 = bed_no_overlap_sorted[(bed_no_overlap_sorted[4] > quants[0.4]) & (bed_no_overlap_sorted[4] <= quants[0.5])]
	q6 = bed_no_overlap_sorted[(bed_no_overlap_sorted[4] > quants[0.5]) & (bed_no_overlap_sorted[4] <= quants[0.6])]
	q7 = bed_no_overlap_sorted[(bed_no_overlap_sorted[4] > quants[0.6]) & (bed_no_overlap_sorted[4] <= quants[0.7])]
	q8 = bed_no_overlap_sorted[(bed_no_overlap_sorted[4] > quants[0.7]) & (bed_no_overlap_sorted[4] <= quants[0.8])]
	q9 = bed_no_overlap_sorted[(bed_no_overlap_sorted[4] > quants[0.8]) & (bed_no_overlap_sorted[4] <= quants[0.9])]
	q10 = bed_no_overlap_sorted[bed_no_overlap_sorted[4] > quants[0.9]]

	print('q10 peak number: ' + str(q10.shape[0]))
	print('q9 peak number: ' + str(q9.shape[0]))
	print('q8 peak number: ' + str(q8.shape[0]))
	print('q7 peak number: ' + str(q7.shape[0]))
	print('q6 peak number: ' + str(q6.shape[0]))
	print('q5 peak number: ' + str(q5.shape[0]))
	print('q4 peak number: ' + str(q4.shape[0]))
	print('q3 peak number: ' + str(q3.shape[0]))
	print('q2 peak number: ' + str(q2.shape[0]))
	print('q1 peak number: ' + str(q1.shape[0]))

	quarts = [q10, q9, q8, q7, q6, q5, q4, q3, q2, q1]

	num_cores = multiprocessing.cpu_count()
	print('executing with ' + str(num_cores) + ' cores')

	i = 0
	for bam in bams:

		plt.close('all')

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


			resolution_analysis(all_data, all_base_pos_data, window_size, out, names[i])

			i=i+1

			plt.close('all')





if __name__ == '__main__':
	main()




