# create plot as in Figure 4d to study joint occupancy at neighboring CTCF sites
# ## Annie Maslan
# ## 06.30.21
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pysam
import multiprocessing
from joblib import Parallel, delayed
from pybedtools import BedTool
from matplotlib import colors
from sklearn.cluster import KMeans
from mpl_toolkits.axes_grid1 import make_axes_locatable


class Region(object):
	def __init__(self, region):
		self.chromosome = region[1][0]
		# now being and end correspond to motif begin and end
		self.begin = region[1][1]
		self.end = region[1][2]
		self.size = self.end - self.begin
		self.string = f"{self.chromosome}_{self.begin}_{self.end}"
		# additional info to store about CTCF regions
		self.strand = region[1][3]
		self.peak_strength = region[1][4]

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

def extract_peak_pairs(bed, min_dist, max_dist, out):
	b = BedTool(bed)
	b.sort().saveas(out + '/tmp.sorted.bed');
	bed = pd.read_csv(out + '/tmp.sorted.bed', sep='\t', header=None)
	# find middle
	left_keep = []
	right_keep = []
	bed['middle'] = (bed[1] + bed[2]) / 2
	for i in range(0, len(bed) - 1):
		if (bed.iloc[i+1].middle - bed.iloc[i].middle >= min_dist) & (bed.iloc[i+1].middle - bed.iloc[i].middle <= max_dist):
			left_keep.append(i)
			right_keep.append(i+1)
	return bed.iloc[left_keep], bed.iloc[right_keep]

def get_data(bam, peak_left, peak_right, num_cores):
	meth_data = Parallel(n_jobs=num_cores)(delayed(parse_ont_bam)(bam, w1, w2) for w1, w2 in zip(peak_left, peak_right))
	return meth_data

def parse_ont_bam(filename, w1, w2):
	bam = pysam.AlignmentFile(filename, "rb")
	data = []
	reads1 = bam.fetch(reference=w1.chromosome, start=w1.begin, end=w1.end)
	reads2 = bam.fetch(reference=w2.chromosome, start=w2.begin, end=w2.end)
	for read in reads1:
		# only add reads that span both sites
		# includeunique identifier that is w1-w2-read_name
		if read in reads2:
			[(mod, positionsL, positionsR, quals), (mod2, positions2L, positions2R, quals2)] = get_modified_reference_positions(read, w1, w2)
			for posL, posR, qual in zip(positionsL, positionsR, quals):
				if posL is not None:
					if abs(posL) <= 1000:
						data.append((read.query_name,
							'-' if read.is_reverse else '+',
							posL,
							qual,
							mod,
							w1.peak_strength,
							w1.string + '-' + w2.string + '-' + read.query_name))
					if abs(posR) <= 1000:
						posR_adj = posR + 3000
						data.append((read.query_name,
							'-' if read.is_reverse else '+',
							posR_adj,
							qual,
							mod,
							w2.peak_strength,
							w1.string + '-' + w2.string + '-' + read.query_name))
			for posL, posR, qual in zip(positions2L, positions2R, quals2):
				if posL is not None:
					if abs(posL) <= 1000: 
						data.append((read.query_name,
							'-' if read.is_reverse else '+',
							posL,
							qual,
							mod2,
							w1.peak_strength,
							w1.string + '-' + w2.string + '-' + read.query_name))
					if abs(posR) <= 1000:
						posR_adj = posR + 3000
						data.append((read.query_name,
							'-' if read.is_reverse else '+',
							posR_adj,
							qual,
							mod2,
							w2.peak_strength,
							w1.string + '-' + w2.string + '-' + read.query_name))
	data_return = Methylation(
		table=pd.DataFrame(data, columns=['read_name', 'strand', 'pos', 'quality', 'mod', 'peak_strength', 'id'])
				.astype(dtype={'mod': 'category', 'pos': 'int16', 'quality': 'int16'})# .astype(dtype={'mod': 'category', 'quality': 'float'})
				.sort_values(['read_name', 'pos']),
		data_type="ont-bam",
		name='double peaks',
		called_sites=len(data))
	return data_return

def get_modified_reference_positions(read, w1, w2):
	if (read.has_tag('Mm')) & (';' in read.get_tag('Mm')):
		mod1 = read.get_tag('Mm').split(';')[0].split(',', 1)[0]
		mod2 = read.get_tag('Mm').split(';')[1].split(',', 1)[0]
		mod1_list = read.get_tag('Mm').split(';')[0].split(',', 1)
		mod2_list = read.get_tag('Mm').split(';')[1].split(',', 1)
		if len(mod1_list) > 1:
			mod1_return = get_pos_prob(read, mod1, 0, w1, w2)
		else:
			mod1_return = (None, [None], [None], [None])
		if len(mod2_list) > 1:
			mod2_return = get_pos_prob(read, mod2, 1, w1, w2)
			return (mod1_return, mod2_return)
		else:
			return (mod1_return, (None, [None], [None], [None]))
	else:
		return ((None, [None], [None], [None]), (None, [None], [None], [None]))

def get_pos_prob(read, basemod, index, w1, w2):
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
	# adjust position and return both left anchored and right anchored positions
	# account for strand of motif and flip pos if on - strand
	if (w1.strand == '+'):
		w1_pos = np.array(refpos[keep]) - round(((w1.end-w1.begin)/2 + w1.begin))
	if (w1.strand == '-'):
		w1_pos = -1*(np.array(refpos[keep]) - round(((w1.end-w1.begin)/2 + w1.begin)))
	if (w2.strand == '+'):
		w2_pos = np.array(refpos[keep]) - round(((w2.end-w2.begin)/2 + w2.begin))
	if (w2.strand == '-'):
		w2_pos = -1*(np.array(refpos[keep]) - round(((w2.end-w2.begin)/2 + w2.begin)))

	return (basemod, w1_pos, w2_pos, probabilities[prob_keep])

	# instead return with respec to motif center

def bin_qualities(all_data, mod):
	''' 
	bin data for 10 A's as slide across the read
	calculate the probability at least one base is methylated
	'''
	all_data_mod = all_data[all_data['mod'].str.contains(mod)] 
	all_data_mod.loc[:, 'quality'] = all_data_mod['quality'] / 255
	read_ids = all_data_mod['id'].unique()
	for r in read_ids:
		subset = all_data_mod[all_data_mod['id'] == r] 
		qualities = subset['quality'] 
		binned_qualities = qualities.rolling(window=20, center=True).apply(lambda b: prob_bin(b))

		all_data_mod.loc[all_data_mod['id'] == r,'quality'] = binned_qualities 
	return all_data_mod


def prob_bin(bin):
	# probability a base in the window is methylated by:
	# calculating probability that no base in the window is methylated and then taking the complement
	# treat p=1 as 254/255 for prevent log(0)
	probs = [np.log(1-p) for p in bin if ((p < 1) and (p >= 0.5))] # only consider probabilities > 0.5 and handle 1 later
	probs1 = [np.log(1-254/255) for p in bin if p == 1]
	probsAll = probs + probs1
	prob = 1 - np.exp(sum(probsAll)) 
	return prob


def make_cluster_plot(all_data, all_data_C, thresh, out, name):
	print(name + ' A: ' + str(all_data.shape[0]))
	all_data_t = all_data[all_data['quality'] > thresh] 
	# require that quality > thresh be within 100 bp of the peak center on either side
	peak = all_data_t[abs(all_data_t['pos'] <= 100)]
	peak2 = all_data_t[abs(all_data_t['pos'] > 2900) & abs(all_data_t['pos'] < 3100) ]
	peak_ids = peak['id'].unique()
	peak_ids2 = peak2['id'].unique()
	boolean_keep_series = (all_data_t.id.isin(peak_ids) | all_data_t.id.isin(peak_ids2)) #reads_keep
	all_data_t_p = all_data_t[boolean_keep_series]
	print(name + ' mA: ' + str(all_data_t.shape[0]))
	print(name + ' mA with peak: ' + str(all_data_t_p.shape[0]))
	all_data_pivoted = pd.pivot_table(all_data_t_p, values = 'quality', columns = 'pos', index='id') # index was read_name #p
	r = range(-1000, 4000+1, 1)
	for bp in r:
		if bp not in all_data_pivoted.columns:
			all_data_pivoted[bp] = np.nan
	all_data_pivoted = all_data_pivoted.sort_index(axis=1)
	all_data_pivoted_0 = all_data_pivoted.fillna(0)
	cmapA = colors.LinearSegmentedColormap.from_list('custom A', ['white', '#053C5E'], N=255)

	all_data_pivoted_mod_0_rolling = pd.DataFrame()
	for i in range(0, all_data_pivoted_0.shape[0]):
		all_data_pivoted_mod_0_rolling = all_data_pivoted_mod_0_rolling.append(all_data_pivoted_0.iloc[i,:].rolling(window=5).mean()) #33
	all_data_pivoted_mod_0_rolling_0 = all_data_pivoted_mod_0_rolling.fillna(0)
	k = KMeans(n_clusters=4, random_state=1) # try different numbers of clusters #2
	k.fit(all_data_pivoted_mod_0_rolling_0)

	# sort by left peak signal strength after labels
	all_data_pivoted_0['left_sum'] = 0
	subset_left_peak = all_data_pivoted_0.iloc[:,900:1100]
	for idx, row in subset_left_peak.iterrows():
		left_sum = row.sum() 
		all_data_pivoted_0.loc[idx, 'left_sum'] = left_sum
	all_data_pivoted_0['labels'] = k.labels_
	all_data_pivoted_0 = all_data_pivoted_0.sort_values(by=['labels', 'left_sum'], axis=0, ascending=False) 
	to_plot = all_data_pivoted_0 
	to_plot_2 = to_plot.loc[:, (to_plot.columns != 'labels') & (to_plot.columns != 'left_sum')]

	fig=plt.figure()
	g = sns.heatmap(to_plot_2, cmap=cmapA, xticklabels=False, yticklabels=False,
		cbar_kws = dict(use_gridspec=False,location="top")) 
	plt.show()
	fig.savefig(out + '/' + name + '_cluster_double_peak.4.rolling.05.thresh90.w20.peak.noSmooth.png', dpi=500)

	# also plot 1D heatmap with the peak signal strength in order of reads shown
	# updated to purple
	cmapYellow = colors.LinearSegmentedColormap.from_list('custom yellow', ['white', '#610345'], N=200)
	
	ordered_read_ids = to_plot.index.values
	print('all clusters: ' + str(len(all_data_pivoted_mod_0_rolling_0.index.values)))
	w1_signal_strength = []
	w2_signal_strength = []
	for r in ordered_read_ids:
		read_data =  all_data[all_data['id'] == r]
		w1_peak = read_data[read_data['pos'] <= 1000]['peak_strength'].unique()
		w2_peak = read_data[read_data['pos'] >= 2000]['peak_strength'].unique()
		w1_signal_strength.append(w1_peak[0])
		w2_signal_strength.append(w2_peak[0])

	fig = plt.figure()
	x = np.linspace(0, len(ordered_read_ids), num=len(ordered_read_ids))
	y = pd.Series(w1_signal_strength) 
	fig, (ax,ax2) = plt.subplots(nrows=2, sharex=True)
	extent = [x[0]-(x[1]-x[0])/2., x[-1]+(x[1]-x[0])/2.,0,1]
	im = ax.imshow(y[np.newaxis,:], cmap=cmapYellow, aspect="auto", extent=extent) 
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('top', size='5%', pad=0.25)
	fig.colorbar(im, cax=cax, orientation='horizontal')
	ax.set_yticks([])
	ax.set_xlim(extent[0], extent[1])
	ax2.plot(x,y,'o',ms=0.5,color='#2D1E2F')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.get_xaxis().set_ticks([])
	plt.tight_layout()
	plt.show()
	fig.savefig(out + '/' + name + '_peak_signal_w1_t90.png', dpi=300)

	# window 2 
	# from positions 2000-4000
	fig = plt.figure()
	x = np.linspace(0, len(ordered_read_ids), num=len(ordered_read_ids))
	y = pd.Series(w2_signal_strength) 
	fig, (ax,ax2) = plt.subplots(nrows=2, sharex=True)
	extent = [x[0]-(x[1]-x[0])/2., x[-1]+(x[1]-x[0])/2.,0,1]
	im = ax.imshow(y[np.newaxis,:], cmap=cmapYellow, aspect="auto", extent=extent) 
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('top', size='5%', pad=0.25)
	fig.colorbar(im, cax=cax, orientation='horizontal')
	ax.set_yticks([])
	ax.set_xlim(extent[0], extent[1])
	ax2.plot(x,y,'o',ms=0.5,color='#2D1E2F')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.get_xaxis().set_ticks([])
	plt.tight_layout()
	plt.show()
	fig.savefig(out + '/' + name + '_peak_signal_w2_t90.png', dpi=300)

	peak_ids = peak['id'].unique()
	boolean_keep_series = all_data_t.id.isin(peak_ids)

	C_data_plot = pd.DataFrame()
	all_data_C_t = all_data_C[all_data_C['quality'] > thresh] 
	id_list = []
	quality_list = []
	pos_list = []
	for r in ordered_read_ids:
		sub = all_data_C_t[all_data_C_t['id'] == r]
		if sub is None:
			print('here')
		else:
			id_list.append(sub['id'].values.tolist())
			quality_list.append(sub['quality'].values.tolist())
			pos_list.append(sub['pos'].values.tolist())

	id_flat_list = [item for sublist in id_list for item in sublist]
	quality_flat_list = [item for sublist in quality_list for item in sublist]
	pos_flat_list = [item for sublist in pos_list for item in sublist]

	C_data_plot = pd.DataFrame({'id': id_flat_list, 'quality': quality_flat_list, 'pos': pos_flat_list})
	C_data_plot_pivoted = pd.pivot_table(C_data_plot, values = 'quality', columns = 'pos', index='id')

	r = range(-1000, 4000+1, 1)
	for bp in r:
		if bp not in C_data_plot_pivoted.columns:
			C_data_plot_pivoted[bp] = np.nan
	C_data_plot_pivoted = C_data_plot_pivoted.sort_index(axis=1)
	C_data_plot_pivoted_0 = C_data_plot_pivoted.fillna(0)

	cmapC = colors.LinearSegmentedColormap.from_list('custom A', ['white', '#BB4430'], N=255)

	fig=plt.figure()
	g = sns.heatmap(C_data_plot_pivoted_0, cmap=cmapC, xticklabels=False, yticklabels=False,
		cbar_kws = dict(use_gridspec=False,location="top")) 
	plt.show()
	fig.savefig(out + '/' + name + '_cluster_double_peak.4.rolling.05.thresh90.w20.peak.noSmooth.C.png', dpi=500)



def main():

	bams = ["prod_ctcf_mod_mappings_merge.sorted.bam",
	"prod_free_Hia5_mod_mappings.sorted.bam",
	"prod_IgG_mod_mappings.sorted.bam",
	"prod_untreated_mod_mappings.sorted.bam"]

	names = ['CTCF', 'Hia5', 'IgG', 'untreated']

	bed = "intersection.motifs.chip.formatted.chm13.bed" 
	top_bed_path = "top_ctcf_peaks_motif.chm13.bed"

	out = "results"

	thresh = 0.9 # this is threshold for probability at least 1 A is methylated - only include reads with a bin above this threshold

	# extract peaks that are 2kb - 10kb apart; don't allow any peaks < 2kb apart so not in the same visualization window
	# create pairs of peaks
	min_dist = 2000
	max_dist = 10000

	# just look at top 50% ChIP seq peaks sorted by signal strength
	b = pd.read_csv(bed, sep='\t', header=None)
	bed_sorted = b.sort_values(by=b.columns[4], ascending=False)
	mid = bed_sorted[4].quantile(q=0.5)
	top = bed_sorted[bed_sorted[4] >= mid]
	top.to_csv(top_bed_path, index=False, header=False, sep='\t')

	peak_left, peak_right = extract_peak_pairs(top_bed_path, min_dist, max_dist, out)

	# 625 peak pairs
	print('number of peaks: ' + str(peak_left.shape[0]))
	# make windows for the peak files
	peak_left_windows = make_windows(peak_left)
	peak_right_windows = make_windows(peak_right)

	# extract reads that overlap both peaks in the pair; only keep read ids that overlap both
	# adjust left overlap pos to be -1000-1000 and right overlap pos to be 2000-4000 so all in the same matrix
	num_cores = multiprocessing.cpu_count()

	i=0
	for bam in bams:
		meth_data = get_data(bam, peak_left_windows, peak_right_windows, num_cores)

		# combine all methylation tables into a single dataframe
		list_tables = []
		for m in meth_data:
			list_tables.append(m.table)
		all_data = pd.concat(list_tables)

		all_data_A_binned = bin_qualities(all_data, 'A')
		all_data_C_binned = bin_qualities(all_data, 'C')

		print('processing ' + str(len(all_data_A_binned['read_name'].unique())) + ' reads for for bam: ' + bam)

		# cluster matrix and plot
		make_cluster_plot(all_data_A_binned, all_data_C_binned, thresh, out, names[i])
		i=i+1


if __name__ == '__main__':
	main()