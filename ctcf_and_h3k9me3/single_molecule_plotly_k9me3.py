# plots mA and mC probabilities above threshold for single molecules
# modified from methplotlib
# input from guppy&winnow merge or megalodon
#
# use:
# python single_molecule.py -i <bam_file(s)> -w <window> -t <thresh> -b <bed_file> -n <name(s)> 
#
# example use:
# python single_molecule_plotly_meg_winnow_guppy_FINAL.py -i \
# /clusterfs/rosalind/groups/streetslab/amaslan/nanopore/paper/h3k9me3/bams/guppywinnow/prod_H3K9me3_winnowmap_guppy_merge.sorted.q10.bam \
# /clusterfs/rosalind/groups/streetslab/amaslan/nanopore/paper/h3k9me3/bams/guppywinnow/prod_Hia5_HG002_winnowmap_guppy_merge.sorted.q10.bam \
# /clusterfs/rosalind/groups/streetslab/amaslan/nanopore/paper/h3k9me3/bams/guppywinnow/prod_IgG_HG002_winnowmap_guppy_merge.sorted.q10.bam \
# -t 230 -w chr8:43000000-45000000 -n H3K9me3 Hia5 IgG

# python single_molecule_plotly_meg_winnow_guppy.py -i /clusterfs/rosalind/groups/streetslab/amaslan/nanopore/20210506/h3k9me3/sre_guppy_winnow/20210524_1_SRE_winnnowmap_guppy_merge.sorted.bam \
# -w chr6:58000000-62000000 -t 128 -b /clusterfs/rosalind/groups/streetslab/amaslan/nanopore/paper/h3k9me3/sm/t2t_cenAnnotation.v2.021921.cen6.all.bed \
# -n H3K9me3
#
# threshold is 0-255 --> probability of methylation is displayed when hover over a dot and also indicated by the intensity of the dot color
# window can be a range as in the example above 'chr14:95000000-105000000' or can be an entire chromosome, e.g. 'chr14'
# can specify multiple input files, just also specify the same number of names 
#
#
# Known plotly bug where for certain windows the x-axis labels don't show up. If this happens, change the window or uncomment line: showticklabels=True,
#
# Annie Maslan
# 06.08.21



import pandas as pd
import pyranges as pr
import numpy as np
import sys
import logging
import pysam
import plotly.io as pio
import plotly.graph_objs as go
import plotly
from math import ceil
from argparse import ArgumentParser
from matplotlib import colors
import matplotlib.pyplot as plt
import seaborn as sns

def main():
	args = get_args()
	windows = make_windows(args.window, fasta=args.fasta)

	for w in windows:
		meth_data = get_data(args.input, args.names, w, args.thresh)
		for m in meth_data:
			if m.table.empty:
				print('empty window')
			else:
				all_data = m.table
				all_data_mA = all_data[all_data['mod'].str.contains('A')] 
				all_data_mC = all_data[all_data['mod'].str.contains('C')] 

				meth_browser(meth_data=meth_data,
					window=w,
					gtf=False,
					bed=args.bed,
					simplify=False,
					split=False,
					outfile='static_supp_all_mA.pdf', #name for static #pdf
					dotsize=2,
					static=True, # True for static png!
					binary=False,
					thresh=args.thresh
					)

class Region(object):
	def __init__(self, region, fasta=None):
		if ':' in region:
			try:
				self.chromosome, interval = region.replace(',', '').split(':')
				try:
					# see if just integer chromosomes are used
					self.chromosome = int(self.chromosome)
				except ValueError:
					pass
				self.begin, self.end = [int(i) for i in interval.split('-')]
			except ValueError:
				sys.exit("\n\nERROR: Window (-w/--window) inproperly formatted, "
						 "examples of accepted formats are:\n"
						 "'chr5:150200605-150423790'\n\n")
			self.size = self.end - self.begin
			self.string = f"{self.chromosome}_{self.begin}_{self.end}"
		else:  # When region is an entire chromosome, contig or transcript
			if fasta is None:
				sys.exit("A fasta reference file is required if --window "
						 "is an entire chromosome, contig or transcript")
			else:
				from pyfaidx import Fasta
				self.chromosome = region
				self.begin = 0
				self.string = region
				self.end = len(Fasta(fasta)[region])
				self.size = self.end

class Methylation(object):
	def __init__(self, table, data_type, name, called_sites):
		self.table = table
		self.data_type = data_type
		self.name = name
		self.called_sites = called_sites

class DataTraces(object):
	def __init__(self, traces, types, names, split):
		self.traces = traces
		self.types = types
		self.names = names
		self.split = split
		self.index = 0
	def __iter__(self):
		return self
	def __next__(self):
		if self.index == len(self.traces):
			raise StopIteration
		else:
			self.index += 1
			return self.traces[self.index - 1], self.types[self.index - 1]

def get_args():
	parser = ArgumentParser(description="plotting megalodon methylation calls")
	parser.add_argument("-i", "--input",
						nargs='+',
						help="data in ont-cram format",
						required=True)
	parser.add_argument("-n", "--names",
						nargs='+',
						help="names of datasets in --input",
						required=True)
	parser.add_argument("-w", "--window",
						help="window (region) to which the visualisation has to be restricted",
						required=True)
	parser.add_argument("-b", "--bed",
						help="add annotation based on a bed file")
	parser.add_argument("-f", "--fasta",
						help="required when --window is an entire chromosome, contig or transcript")
	parser.add_argument("-t", "--thresh",
						help="The minimal phred quality to show",
						type=int,
						default=128) # default thresh50
	args = parser.parse_args()
	if not len(args.names) == len(args.input):
		sys.exit("INPUT ERROR: Expecting the same number of names as datasets!")
	return args

def make_windows(full_window, max_size=20e6, fasta=None): #max size 20M
	'''
	if window is > 1Mb split into 1 Mb files
	'''
	reg = Region(full_window, fasta)
	if reg.size > max_size:
		chunks = ceil(reg.size / max_size)
		chsize = ceil(reg.size / chunks)
		return [
			Region(f"{reg.chromosome}:{reg.begin + i * chsize}-{reg.begin + (i + 1) * chsize}")
			for i in range(chunks)]
	else:
		return [reg]

def get_data(methylation_files, names, window, thresh):
	"""
	Import methylation data from all files in the list methylation_files
	Data in cram format
	data is extracted within the window args.window
	Frequencies are smoothened using a sliding window
	"""
	return [parse_ont_bam(f, n, window, thresh) for f, n in zip(methylation_files, names)]

def parse_ont_bam(filename, name, window, thresh):
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
				data.append((read.query_name,
					'-' if read.is_reverse else '+',
					pos,
					qual,
					mod))
		for pos, qual in zip(positions2, quals2):
			if pos is not None:
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
	if 'C' in basemod: 
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
	return (basemod, np.array(refpos[keep]), probabilities[prob_keep])
	
####### below this point are all plotting functions #######
####### adapted from methplotlib & modified to visualize both mA and mC #######


def create_subplots(num_methrows, split, names=None, annotation=True):
	'''
	Prepare the panels (rows * 1 column) for the subplots.
	If splitting: one row for each dataset, taking 90%/len(datasets) for heights
	If not: one row spanning 4 rows and taking 90% of the heights
	if annotation is True (bed or gtf) then add a row with height 10%
	'''
	if split:
		return plotly.subplots.make_subplots(
			rows=num_methrows + annotation,
			cols=1,
			shared_xaxes=True,
			specs=[[{}] for i in range(num_methrows + annotation)],
			print_grid=False,
			subplot_titles=names,
			vertical_spacing=0.1 if num_methrows < 10 else 0.01,
			row_heights=[0.9 / num_methrows] * num_methrows + [0.1] * annotation

		)
	else:
		return plotly.subplots.make_subplots(
			rows=num_methrows + annotation,
			cols=1,
			shared_xaxes=True,
			specs=[[{'rowspan': num_methrows}], [None], [None], [None]] + [[{}]] * annotation,
			print_grid=False,
			vertical_spacing=0.1 if num_methrows < 10 else 0.01,
			row_heights=[0.9, 0, 0, 0] + [0.1] * annotation
		)

def create_browser_output(fig, outfile, window):
	outfile = f"methylation_browser_{window.string}.html"
	if outfile.endswith(".html"):
		write_html_output(fig, outfile)
	else:
		try:
			fig.write_image(outfile, scale=10)#, width=1000, height=500)
		except ValueError as e:
			sys.stderr.write("\n\nERROR: creating the image in this file format failed.\n")
			sys.stderr.write("ERROR: creating in default html format instead.\n")
			sys.stderr.write("ERROR: additional packages required. Detailed error:\n")
			sys.stderr.write(str(e))
			write_html_output(fig, outfile)


def write_html_output(fig, outfile):
	with open(outfile, "w+") as output:
		output.write(plotly.offline.plot(fig,
										 output_type="div",
										 show_link=False,
										 include_plotlyjs='cdn'))

def methylation(meth_data, dotsize=4, binary=False, thresh=128):
	"""
	Plot methylation traces from various data types
	"""
	traces = []
	types = []
	names = []
	split = False
	for meth in meth_data:
		traces.append(
			make_per_read_meth_traces_phred(table=meth.table,
				dotsize=dotsize, thresh=thresh)
			)
		split = True
		types.append(meth.data_type)
		names.append(meth.name)
	return DataTraces(traces=traces,
		types=types,
		names=names,
		split=split)

def make_per_read_meth_traces_phred(table, max_cov=100, dotsize=4, thresh=128):
	"""Make traces for each read"""
	minmax_table = find_min_and_max_pos_per_read(table)
	df_heights = assign_y_height_per_read(minmax_table, max_coverage=max_cov)
	table = pd.merge(table, df_heights, left_on="read_name", right_on="read")
	traces = []
	hidden = 0
	for read in table["read_name"].unique():
		strand = table.loc[table["read_name"] == read, "strand"].values[0]
		try:
			traces.append(
				make_per_read_line_trace(read_range=minmax_table.loc[read],
										 y_pos=df_heights.loc[read, 'height'],
										 strand=strand)
			)
		except KeyError:
			hidden += 1
			continue
	if hidden:
		sys.stderr.write(f"Warning: hiding {hidden} reads because coverage above {max_cov}x.\n")
	read_table_mC = table[table['mod'].str.contains('C')]
	read_table_mA = table[table['mod'].str.contains('A')]
	cmapA = ['white', '#053C5E']
	cmapC = ['white', '#BB4430']
	traces.append(
		make_per_position_phred_scatter(read_table=read_table_mC[read_table_mC['quality'] > thresh], mod="mC", dotsize=dotsize, colorscale=cmapC, offset=0.05)
	) # comment out for single mod plot
	traces.append(
		make_per_position_phred_scatter(read_table=read_table_mA[read_table_mA['quality'] > thresh], mod="mA", dotsize=dotsize, colorscale=cmapA, offset=0.15)
	) # comment out for single mod plot
	return traces

def make_per_position_phred_scatter(read_table, mod, dotsize=4, colorscale='Reds', offset=0):
	"""Make scatter plot per modified base per read"""
	return go.Scatter(x=read_table['pos'],
					  y=read_table['height'],
					  mode='markers',
					  showlegend=False,
					  text=round(read_table['quality'] / 255, 2),
					  hoverinfo="text",
					  marker=dict(size=dotsize,
								  color=read_table['quality'],
								  colorscale=colorscale,
								  colorbar=dict(title=mod + " modification probability",
												titleside="right",
												tickvals=[read_table['quality'].min(),
														  read_table['quality'].max()],
												ticktext=[str(round(read_table['quality'].min() / 255, 2)),
														  str(round(read_table['quality'].max() / 255, 3))],
												ticks="outside",
												x=offset+1)
								  ))


def find_min_and_max_pos_per_read(table, phased=False):
	"""Return a table with for every read the minimum and maximum position"""
	mm_table = table.loc[:, ["read_name", "pos"]] \
		.groupby('read_name') \
		.min() \
		.join(table.loc[:, ["read_name", "pos"]]
			  .groupby('read_name')
			  .max(),
			  lsuffix="min",
			  rsuffix="max")
	if phased:
		return mm_table.join(table.loc[:, ["read_name", "HP"]]
							 .drop_duplicates(subset='read_name')
							 .set_index('read_name'))
	else:
		return mm_table


def assign_y_height_per_read(df, phased=False, max_coverage=1000):
	"""Assign height of the read in the per read traces
	Gets a dataframe of read_name, posmin and posmax.
	Sorting by position, and optionally by phase block.
	Determines optimal height (y coordinate) for this read
	Returns a dictionary mapping read_name to y_coord
	"""
	if phased:
		dfs = df.sort_values(by=['HP', 'posmin', 'posmax'],
							 ascending=[True, True, False])
	else:
		dfs = df.sort_values(by=['posmin', 'posmax'],
							 ascending=[True, False])
	heights = [[] for i in range(max_coverage)]
	y_pos = dict()
	for read in dfs.itertuples():
		for y, layer in enumerate(heights, start=1):
			if len(layer) == 0:
				layer.append(read.posmax)
				y_pos[read.Index] = y
				break
			if read.posmin > layer[-1]:
				layer.append(read.posmax)
				y_pos[read.Index] = y
				break
	return pd.DataFrame({'read': list(y_pos.keys()),
						 'height': list(y_pos.values())}) \
		.set_index('read')

def make_per_read_line_trace(read_range, y_pos, strand, phase=None):
	"""Make a grey line trace for a single read,
	with black arrow symbols on the edges indicating strand"""
	symbol = "triangle-right" if strand == "+" else "triangle-left"
	if phase:
		if phase == 1:
			color = 'lightgreen'
		elif phase == 2:
			color = 'yellow'
		else:  # phase is np.nan
			color = '#ffffff' #black
	else:
		color = '#ffffff' #black
	return go.Scatter(x=[read_range['posmin'], read_range['posmax']],
					  y=[y_pos, y_pos],
					  mode='lines', #'lines+markers','
					  line=dict(width=1, color='lightgrey'),
					  showlegend=False)

def meth_browser(meth_data, window, gtf=False, bed=False, simplify=False,
				 split=False, outfile=None, dotsize=4, static=False, binary=False, thresh=128):
	"""
	meth_Data is a list of Methylation objects from the import_methylation submodule
	annotation is optional and is a gtf or bed file
	if the traces are to be --split per sample (or include raw data as flagged in data.split),
	 then show one line per sample and one for the annotation, with methrows = number of datasets
	if no splitting is needed,
	 then 4/5 of the browser is used for overlayed samples and one for gtf annotation
	the trace to be used for annotation is thus always num_methrows + 1
	"""
	meth_traces = methylation(meth_data, dotsize=dotsize, binary=binary, thresh=thresh)
	logging.info("Prepared methylation traces.")
	if split or meth_traces.split:
		num_methrows = len(meth_data)
		logging.info(f'Making browser in split mode, with {num_methrows} modification rows.')
		annot_row = num_methrows + 1
		annot_axis = f'yaxis{annot_row}'
		fig = create_subplots(num_methrows,
									split=True,
									names=meth_traces.names,
									annotation=bool(bed or gtf))
		for y, (sample_traces, sample_type) in enumerate(meth_traces, start=1):
			logging.info(f"Adding traces of type {sample_type} at height {y}")
			for meth_trace in sample_traces:
				fig.append_trace(trace=meth_trace, row=y, col=1)
			fig["layout"][f"yaxis{y}"].update(title="Reads")
	else:
		logging.info('Making browser in overlaying mode.')
		num_methrows = 4
		annot_row = 5
		annot_axis = 'yaxis2'
		fig = create_subplots(num_methrows, split=False, annotation=bool(bed or gtf))
		for meth_trace in meth_traces.traces:
			for trace in meth_trace:
				fig.append_trace(trace=trace, row=1, col=1)
		fig["layout"].update(legend=dict(orientation='h'))
	logging.info("Prepared modification plots.")
	if bed:
		for annot_trace in bed_annotation(bed, window):
			fig.append_trace(trace=annot_trace, row=annot_row, col=1)
		y_max = -2
	if bed:
		fig["layout"][annot_axis].update(range=[-2, y_max + 1],
										 showgrid=False,
										 zeroline=False,
										 showline=False,
										 ticks='',
										 showticklabels=False)
		logging.info("Prepared annotation plots.")
	fig["layout"]["xaxis"].update(tickformat='g',
								  separatethousands=True,
								  #showticklabels=True,
								  range=[window.begin, window.end])
	fig["layout"].update(barmode='overlay',
						 title=window.chromosome,
						 hovermode='closest',
						 plot_bgcolor='rgba(0,0,0,0)')
	if num_methrows > 10: 
		for i in fig['layout']['annotations']:
			i['font']['size'] = 10
	create_browser_output(fig, outfile, window)
	if static:
		pio.write_image(fig, outfile, format='pdf', scale=10)

def bed_annotation(bed, window):
	return [go.Scatter(x=[begin, end],
					   y=[-2, -2],
					   mode='lines',
					   line=dict(width=16, color='grey'),
					   text=name,
					   hoverinfo='text',
					   showlegend=False)
			for (begin, end, name) in parse_bed(bed, window)]

def parse_bed(bed, window):
	logging.info("Parsing BED file")
	gr = pr.read_bed(bed)[window.chromosome, window.begin:window.end]
	df = gr.unstrand().df
	df = df.drop(columns=["Chromosome", "Score", "Strand"], errors='ignore')
	if "Name" not in df.columns:
		df["Name"] = "noname"
	df_short = df[df.columns[0:3]]
	return df_short.itertuples(index=False, name=None)

if __name__ == '__main__':
	main()

