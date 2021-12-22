# # Single molecules spanning HOR boundaries aggregate
# ## Annie Maslan
# ## 06.04.21

# Input: bed file of coordinates where single molecules should be centered, guppy_winnow_hybrid.bam, guppy_winnow_hybrid.bam.bai

# 1. Input bed file with windows over which to extract reads
# 2. Extract reads
# 3. Plot modified bases and aggregate across boundaries


import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import colors

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
		self.side = region[1][3] # add this for whether L or R boundary of cenpA region

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
	meth_data = Parallel(n_jobs=num_cores)(delayed(parse_ont_bam)(methylation_file, name, w, thresh) for w in windows)
	print('methyltation data collection done')
	base_pos_data = Parallel(n_jobs=num_cores)(delayed(parse_ont_bam_base_pos)(methylation_file, name, w, thresh, window_size) for w in windows)
	print('base abundance data collection done')
	return meth_data, base_pos_data


######
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
	pos_data_return = pd.DataFrame(base_pos_data, columns=['read_name', 'strand', 'base_pos', 'mod']).astype(dtype={'mod': 'category'}) #, 'base_pos': 'int16'})
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
	# perform adjustment for side of the boundary
	# center at center of the motif

	if window.side == 'R':
		return (basemod, np.array(refpos[base_keep]) - round(((window.end-window.begin)/2 + window.begin)))
	if window.side == 'L':
		return (basemod, -1*(np.array(refpos[base_keep]) - round(((window.end-window.begin)/2 + window.begin))))

######



def parse_ont_bam(filename, name, window, thresh):
	'''
	parse hybrid.bam file to create methylation object with read_name, strand, pos, quality, and mod 
	in window above threshold specified
	'''
	bam = pysam.AlignmentFile(filename, "rb")
	data = []
	for read in bam.fetch(reference=window.chromosome, start=window.begin, end=window.end):
		[(mod, positions, quals), (mod2, positions2, quals2)] = get_modified_reference_positions(read, window)
		for pos, qual in zip(positions, quals):
			if pos is not None:
				if qual >= thresh: 
					data.append((read.query_name,
								 '-' if read.is_reverse else '+',
								 pos,
								 qual,
								 mod))
		for pos, qual in zip(positions2, quals2):
			if pos is not None:
				if qual >= thresh:
					data.append((read.query_name,
								 '-' if read.is_reverse else '+',
								 pos,
								 qual,
								 mod2))
	return Methylation(
		table=pd.DataFrame(data, columns=['read_name', 'strand', 'pos', 'quality', 'mod'])
				.astype(dtype={'mod': 'category', 'quality': 'float'})
				.sort_values(['read_name', 'pos']),
		data_type="ont-bam",
		name=name,
		called_sites=len(data))


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
	if 'A+a' not in basemod:
		if 'C+m' not in basemod:
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
	if 'C+m' in basemod:
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
	# account for whether is left or right boundary 
	if window.side == 'R':
		return (basemod, np.array(refpos[keep]) - round(((window.end-window.begin)/2 + window.begin)), probabilities[prob_keep])
	if window.side == 'L':
		return (basemod, -1*(np.array(refpos[keep]) - round(((window.end-window.begin)/2 + window.begin))), probabilities[prob_keep])


def make_plots(all_data, all_base_pos_data, window_size, out, name, smooth):
	'''
	input: dataframe with all the data for plotting, window size, output directory, name for plot
	calls make_mod_plot to make plot for each modification type
	'''

	make_mod_plot(all_data, all_base_pos_data, 'A+a', '#053C5E', window_size, out, name, smooth)
	make_mod_plot(all_data, all_base_pos_data, 'C+m', '#BB4430', window_size, out, name, smooth)


def make_mod_plot(all_data, all_base_pos_data, mod, color, window_size, out, name, smooth):
	'''
	input: dataframe with all the data for plotting, modification type, cmap for heatmap, window size, output directory, name for plot
	output: heatmap saved as png
	'''
	# plot mA for single molecules
	all_data_mod = all_data[all_data['mod'] == mod]
	all_data_mod = all_data_mod[all_data_mod['pos'] >= - window_size]
	all_data_mod = all_data_mod[all_data_mod['pos'] <= window_size]
	all_data_mod['quality'] = 1
	all_data_pivoted_mod = pd.pivot_table(all_data_mod, values = 'quality', columns = 'pos', index='read_name')

	print('number of reads with a ' + mod + ' above thresh: ' + str(len(all_data_mod['read_name'].unique())))


	all_data_pivoted_mod = all_data_pivoted_mod.sort_index(axis=1)

	all_base_pos_data_mod = all_base_pos_data[all_base_pos_data['mod'] == mod]
	all_base_pos_data_mod = all_base_pos_data_mod[all_base_pos_data_mod['base_pos'] >= - window_size]
	all_base_pos_data_mod = all_base_pos_data_mod[all_base_pos_data_mod['base_pos'] <= window_size]
	all_base_pos_data_mod['presence'] = 1
	all_base_pos_data_pivoted_mod = pd.pivot_table(all_base_pos_data_mod, values = 'presence', columns='base_pos', index='read_name')

	# scatterplot
	fig=plt.figure()
	sns.scatterplot(data=all_data_mod, x='pos', y='read_name', color=color, s=0.5, marker='s', linewidth=0, legend=None)
	plt.yticks([])
	plt.title(mod)
	plt.show()
	fig.savefig(out + '/' + name + '_' + mod + '_sm_scatter.png', dpi=300)

	for bp in all_base_pos_data_pivoted_mod.columns:
		if bp not in all_data_pivoted_mod.columns:
			all_data_pivoted_mod[bp] = np.nan
	all_data_pivoted_mod = all_data_pivoted_mod.sort_index(axis=1)
	all_base_pos_data_pivoted_mod = all_base_pos_data_pivoted_mod.sort_index(axis=1)

	all_data_pivoted_mod_0 = all_data_pivoted_mod.fillna(0)
	all_data_pivoted_base_0 = all_base_pos_data_pivoted_mod.fillna(0)

	mod_count = all_data_pivoted_mod_0.sum(axis=0)
	base_count = all_data_pivoted_base_0.sum(axis=0)

	mod_frac = mod_count / base_count
	mod_frac_rolling = mod_frac.rolling(window=smooth).mean() #, min_periods=10 #no center #mean

	pos = [int(i) for i in all_data_pivoted_mod_0.columns]
	d = {'pos': pos, 'frac': mod_frac_rolling.tolist()}
	mod_frac_rolling_df = pd.DataFrame.from_dict(d)

	fig=plt.figure()
	sns.lineplot(data=mod_frac_rolling_df, x='pos', y='frac', color=color);
	sns.lineplot(data=mod_frac_rolling_df, x='pos', y='frac', color=color);
	plt.title(mod)
	plt.show()
	fig.savefig(out + '/' + name + '_' + mod + '_sm_rolling_avg_frac.pdf')


def main():
	#### start of parameters and run-specific file paths ####
	bams = ["prod_H3K9me3_winnowmap_guppy_merge.sorted.q10.bam",
	"prod_Hia5_HG002_winnowmap_guppy_merge.sorted.q10.bam",
	"prod_IgG_HG002_winnowmap_guppy_merge.sorted.q10.bam"]

	# output file directory
	out = "out"

	# 0-255
	thresh = 230 #thresh90

	# will plot non-overlapping features in -window_size to +window_size
	window_size = 300000

	# number of bases to smooth over for moving average curve
	smooth = 60000

	beds = ["sharp_HOR_boundaries.300kb.bed"]

	# name for plots
	names = ["H3K9me3_HOR_boundaries_300kb_60", "Hia5_HOR_boundaries_300kb_60", "IgG_HOR_boundaries_300kb_60"]
	#### end of parameters and run-specific file paths ####




	# perform read extraction in parallelized way over windows
	num_cores = multiprocessing.cpu_count()

	i = 0;
	for bam in bams:
		for bed in beds:
			# remove regions of interest that overlap in plot window
			b = pd.read_csv(bed, sep='\t', header=None)

			print('processing ' + str(b.shape[0]) + ' boundaries')

			# make windows with non-overlapping bedfile
			windows = make_windows(b)

			# get methylation data within all windows in parallelized way
			meth_data, base_pos_data = get_data(bam, names[i], windows, thresh, num_cores, window_size)


			# combine all methylation tables into a single dataframe
			list_tables = []
			for m in meth_data:
				list_tables.append(m.table)
			all_data = pd.concat(list_tables)

			# combine all base data into a single dataframe
			base_list_tables = []
			for b in base_pos_data:
				base_list_tables.append(b)
			all_base_pos_data = pd.concat(base_list_tables)


			print('processing ' + str(len(all_data['read_name'].unique())) + ' modified reads for ' + names[i])
			print('processing ' + str(len(all_base_pos_data['read_name'].unique())) + ' reads for ' + names[i])

			# make heatmaps for mA and CpG methylation
			make_plots(all_data, all_base_pos_data, window_size, out, names[i], smooth)

			i=i+1


if __name__ == '__main__':
	main()




