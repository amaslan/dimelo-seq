### README ###

# #Use this script for analyses of individual reads from DiMeLo-seq in vitro, for plotting mA, clustering reads based on mA positions, 
# for calling nucleosomes on individual reads, and for plotting mA probability with respect to the called nucleosome position
# for plotting histogram of called nucleosome center from 601 dyad axis for all reads in a library
# for plotting average mA probability score (out of 1) for each A position along the template, centered at the called nucleosome dyad center
#Input Megaladon or Guppy with modified base calling bam files - mod_mappings.bam (filter by read length > 3.6 kb (18x601) to keep file size small)
#Template - 18x601.fa or fasta file containing sequence of the chromatin or naked DNA template

# The following filepaths and variables need to be updated for running the script

outfolder = <output_director_path>

#enter name of folder that contains bam
folder_name = <input_directory_path_containing_bam>

Barcode_name = {}

#List Barcode titles to be used for graphs here:
#For example:
Barcode_name[1] = 'CA chr + free PA-Hia5'
Barcode_name[2] = 'CA chr + CA directed PA-Hia5'

#Select barcodes to iterate over for analyses
Barcodes_subset = [1,2]

#reading in template into string - 18x601.fa
array_seq = SeqIO.read(<path_to_template_fasta>, "fasta")
array_str = str(array_seq.seq)

#In this example code, threshold after smoothing is 0.6 x 255 = 150
thr = 150

#number of reads to consider for each library
numreadsmax = 2000

###

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import kde
get_ipython().magic('matplotlib inline')

import math

import pysam
from Bio import SeqIO

from scipy import sparse

from sklearn import cluster, datasets, mixture, decomposition
import scipy.spatial.distance as ssd
from scipy.cluster import hierarchy
from sklearn.metrics.pairwise import pairwise_distances

import time


#identifying indices of A's and T's for assigning mod probs to right positions
A_pos = []
T_pos = []
A_pos_color = ['White' for x in range(len(array_str))]
T_pos_color = ['White' for x in range(len(array_str))]
AT_pos_color = ['White' for x in range(len(array_str))]
for basepos in range(0,len(array_str)-1):
    if array_str[basepos]=='A':
        A_pos.append(basepos)
    elif array_str[basepos]=='T':
        T_pos.append(basepos)

# extracting m6A prob from Ml, counting A
# input Mm Ml from bam and thr
# output - Mm first value, truncated m6A prob, total A count

def extract_Ml_Aa_fast(Mm, Ml):
         
    first=Mm.split(';')[0]
    if 'C+m' in first:
        A_n=len(Mm.split(';')[1].split(','))-1
        Mm_Aa = Mm.split(';')[1].split(',')[1:]
        Mm_Aa_start = int(Mm_Aa[0])
        Ml_Aa = np.array(Ml[0-A_n:],dtype=int)

    elif 'A+a' in first:
        A_n=len(Mm.split(';')[0].split(','))-1        
        Mm_Aa = Mm.split(';')[0].split(',')[1:]
        Mm_Aa_start = int(Mm_Aa[0])
        Ml_Aa=np.array(Ml[0:A_n],dtype=int)
    
    return Mm_Aa_start, Ml_Aa, A_n

# m6A_seqmap_prob_fast_NaN - takes in A or T pos array and readmap start, returns mod prob values assigned to A's along the template
# input arguments - A pos array, T pos array, pos of read along template, if strand=rev, Mm_start, truncated Ml prob values, read length, template length) thresholding done later
# output - m6A prob score along template, coverage along template - for each read
# (fills NaN in non-A positions or positions outside read coverage on template)

def m6A_seqmap_prob_fast_NaN(A_pos_array, T_pos_array, Readmap_start, rev_strand, Mm_start, Ml_Aa, readlength, templatelength = 6302):
    
    if rev_strand == False:
      
    #finding position of first A
        if (Readmap_start < A_pos_array[0]):
                first_A_index = A_pos_array[0]
                
        for i in range(0,len(A_pos_array)-1):
            if ((A_pos_array[i]) <= Readmap_start) and (Readmap_start < (A_pos_array[i+1])):
                if A_pos_array[i] == Readmap_start:
                    first_A_index = i + Mm_start
                else:
                    first_A_index = i+1+Mm_start
                
        #sequence pos of all A's in read
        AT_prob_read = np.empty(templatelength)
        AT_prob_read[:] = np.NaN
        for j in range(0,(len(Ml_Aa)-1)):
            A_index = A_pos_array[first_A_index+j]
            AT_prob_read[A_index] = Ml_Aa[j]
            
    elif rev_strand == True:

        if (Readmap_start < T_pos_array[0]):
                first_T_index = T_pos_array[0]

    #finding position of first T
        for i in range(0,len(T_pos_array)-1):
            if ((T_pos_array[i]) <= Readmap_start) and (Readmap_start < (T_pos_array[i+1])):
                if T_pos_array[i] == Readmap_start:
                    first_T_index = i + Mm_start
                else:
                    first_T_index = i+1+Mm_start
                
                #print(T_pos_array[i], Readmap_start, T_pos_array[first_T_index], T_pos_array[i+1], rev_strand)
                
        #sequence pos of all A's in read
        AT_prob_read = np.empty(templatelength)
        AT_prob_read[:] = np.NaN
        for j in range(0,(len(Ml_Aa)-1)):
            T_index = T_pos_array[first_T_index+j]
            AT_prob_read[T_index] = Ml_Aa[j]
           
    #sequence coverage
    Readmap_end = Readmap_start + readlength - 1
    AT_pos_read = np.zeros(templatelength,dtype=int)
    for k in range(Readmap_start, Readmap_end):
        AT_pos_read[k] = 1
    
    return AT_prob_read, AT_pos_read





#Function for reading individual reads in library, top 'numreads' reads considered
#input - bam file from megalodon or guppy
#outputs - AT_prob_arr_bar, coverage_Fwd_bar, coverage_Rev_bar, A_n_bar, readids_bar, rev_read_list_bar for each barcode(library)
##AT_prob_array - array of mA prob for each position along template for all reads 
# (fills NaN in non-A positions or positions outside read coverage on template)
##coverage_Fwd_bar, coverage_Rev_bar - arrays of coverage(sum of reads aligned to template at each position) for forward and reverse reads in library
##A_n_bar - number of A's per read, used to calculate mA/A
##readids_bar - read ids for all reads in given library
##rev_read_list_bar - boolean array indicating strand of all reads, shape = 1x(number of reads in libary), False - forward read, True - reverse read


def parse_bam6_fast_NaN(filename, numreads=5000, templatelength = 6302):
    reader = pysam.AlignmentFile(filename, check_sq=False)
    bamiter= reader.fetch(until_eof=True)

    A_n=[]
    readids=[]
    rev_read_list=[]
    
    AT_prob_list = []
    coverage_Fwd = np.zeros(templatelength,dtype=int)
    coverage_Rev = np.zeros(templatelength,dtype=int)
    
    count = 0
    
    while True:
        try:
            r=bamiter.__next__()
            Mm=(r.get_tag('Mm') if r.has_tag('Mm') else None)
            Ml=(r.get_tag('Ml') if r.has_tag('Ml') else None)
            readid=r.query_name
            pos=r.reference_start
            #pos in bam is 1 indexed, python is 0 indexed
            pos=pos-1
            length=r.query_length
            rev_strand=r.is_reverse
                        
            if count < numreads:
            
                if not Mm is None:
                
                    Mm_Aa_start, Ml_Aa, read_A_n = extract_Ml_Aa_fast(Mm, Ml)
                    if rev_strand ==True:
                        Mm_Aa_start = 0
                        Ml_Aa = Ml_Aa[::-1]                
                    A_n+=[read_A_n]
                    readids+=[readid]
                    rev_read_list+=[rev_strand]
                
                    AT_prob_read, AT_pos_read = m6A_seqmap_prob_fast_NaN(A_pos, T_pos, pos, rev_strand, Mm_Aa_start, Ml_Aa, length)
                
                    AT_prob_list.append(AT_prob_read)
                
                    if rev_strand == False:
                        coverage_Fwd += AT_pos_read
                    elif rev_strand == True:
                        coverage_Rev += AT_pos_read
            
                    count += 1
                
        except StopIteration:
            reader.close()
            
            AT_prob_arr = np.array(AT_prob_list)
            
            return AT_prob_arr, coverage_Fwd, coverage_Rev, A_n, readids, rev_read_list


# BarCode All through for loop

def by_barcode_NaN(folder_path, barcodes_list = Barcode_list, numreads = 5000):
    AT_prob_arr_bar = {new_list: [] for new_list in barcodes_list} 
    coverage_Fwd_bar = {new_list: [] for new_list in barcodes_list} 
    coverage_Rev_bar = {new_list: [] for new_list in barcodes_list} 
    A_n_bar = {new_list: [] for new_list in barcodes_list} 
    readids_bar = {new_list: [] for new_list in barcodes_list} 
    rev_read_list_bar = {new_list: [] for new_list in barcodes_list} 
    filename = {}
    for i in barcodes_list:
        start_time = time.time()
        #enter filename format such that str(i) fills in the barcode number and the rest of the string matches path to bam file
        #for example:
        filename[i] = folder_path + 'bar' + str(i) + '.1/mod_mappings_3.6kb.bam'
        AT_prob_arr_bar[i], coverage_Fwd_bar[i], coverage_Rev_bar[i], A_n_bar[i], readids_bar[i], rev_read_list_bar[i] = parse_bam6_fast_NaN(filename[i], numreads)
        
        print('Finished barcode ', i, 'in ', "--- %s seconds ---" % (time.time() - start_time))
        
    return AT_prob_arr_bar, coverage_Fwd_bar, coverage_Rev_bar, A_n_bar, readids_bar, rev_read_list_bar


AT_prob_arr_bar_NaN, coverage_Fwd_bar, coverage_Rev_bar, A_n_bar, readids_bar, rev_read_list_bar= by_barcode_NaN(folder_name, Barcodes_subset, numreadsmax)

### The following lines of code are for smoothing (by rolling window averaging) of A methylation probability scores and then converting to binary (based on a given threshold)
# To avoid long processing times, writing smoothed_AT_prob_array to csv files (by barcode) is suggested

#function fr averaging mA probability scores over a rolling window. In this example, 33 bp rolling window is set as default

def smoothing_prob_edge_NaN(AT_prob_array, window_size = 33):
    sliding_avg_edge_NaN = np.zeros(len(AT_prob_array), dtype=float)

    for i in range(0,len(AT_prob_array)):
        sliding_window_adj = []
        window_start = i - math.floor(window_size/2)
        window_start_adj= max(0, window_start)
        window_end = i + math.floor(window_size/2) 
        window_end_adj = min(window_end, len(AT_prob_array))
        sliding_window_adj = list(range(window_start_adj,window_end_adj))
        sliding_avg_edge_NaN[i] = np.nanmean(AT_prob_array[sliding_window_adj])
        
    return sliding_avg_edge_NaN

smoothed_AT_prob_array = {newlist:[] for newlist in Barcodes_subset}
for i in Barcodes_subset:
    smoothed_AT_prob_array[i] = smoothing_prob_edge_NaN(AT_prob_arr_bar_NaN[i])

'''#if out files already exist run the following:
smoothed_AT_prob_array = {newlist:[] for newlist in Barcodes_subset}
for i in Barcodes_subset:
    filename = <filepath_to_smoothed_prob_score_csv>
    pd_smooth_prob = pd.read_csv(filename, header=None)
    smoothed_AT_prob_array[i] = pd_smooth_prob.to_numpy()'''

## Thresholding after smoothing

thr = 150

AT_mod_arr_postsmooth = {new_list: [] for new_list in Barcodes_subset}
for i in Barcodes_subset:
    AT_mod_arr_postsmooth[i] = smoothed_AT_prob_array[i] > thr

'''#if out files already exist run the following:
AT_mod_arr_postsmooth = {new_list: [] for new_list in Barcodes_subset}
for i in Barcodes_subset:
    filename_mod = <filepath_to_thresholded_smoothed_prob_score_csv>
    pd_smooth_mod = pd.read_csv(filename_mod, header=None)
    AT_mod_arr_postsmooth[i] = pd_smooth_mod.to_numpy()'''

### Following code outputs heatmaps of individual reads clustered based on Jaccard metric
#Input - smoothed and thresholded binary arrays of methyl A (NaN - not A or not covered, 1 - methylated A)

#Clustering reads by thresholded smoothed prob scores (for entire 18x601 array region of template)

m6ABlue = sns.light_palette("#053C5E", 1000)
m6ABlue[0] = (1,1,1)

d_rows={}
row_linkage_mod = {}

for i in Barcodes_subset:

    d_rows[i]=pairwise_distances(np.asarray(AT_mod_arr_postsmooth[i])[0:numreadsmax,444:4044],metric='jaccard')

    row_linkage_mod[i] = hierarchy.linkage(ssd.squareform(d_rows[i]), method='average',optimal_ordering=True)

for i in Barcodes_subset:
    
    clstmap=sns.clustermap(np.asarray(AT_mod_arr_postsmooth[i])[0:numreadsmax,444:4044], pivot_kws=None, method='single', metric='jaccard', z_score=None, standard_scale=None, figsize=(4.5,6), cbar_kws=None, cbar_pos=None, row_cluster=True, col_cluster=False, row_linkage=row_linkage_mod[i], col_linkage=None, cmap = m6ABlue, dendrogram_ratio=0.05)
    
    clstmap.fig.suptitle(str(Barcode_name[i]), fontsize=16)
    clstmap.fig.show()
    clstmap.savefig(outfolder + 'm6A_smoothed_prob_thr'+ str(thr) + '_clustered_' + str(i) + '.png', dpi=300)
    

#Clustering reads by thresholded smoothed prob scores (for an example 4x601 array region of template)


d_rows_4x = {}
row_linkage_mod_4x = {}

for i in Barcodes_subset:
    
    d_rows_4x[i]=pairwise_distances(np.asarray(AT_mod_arr_postsmooth[i])[0:numreadsmax,1044:1844],metric='jaccard')

    row_linkage_mod_4x[i] = hierarchy.linkage(ssd.squareform(d_rows_4x[i]), method='average',optimal_ordering=True)

for i in Barcodes_subset:
    clstmap=sns.clustermap(np.asarray(AT_mod_arr_postsmooth[i])[0:numreadsmax,1044:1844], pivot_kws=None, method='single', metric='jaccard', z_score=None, standard_scale=None, figsize=(4.5,6), cbar_kws=None, cbar_pos=None, row_cluster=True, col_cluster=False, row_linkage=row_linkage_mod_4x[i], col_linkage=None, cmap = m6ABlue, dendrogram_ratio=0.05)
    
    clstmap.fig.suptitle(str(Barcode_name[i]), fontsize=16)
    
    #clstmap.fig.suptitle(str(Barcode_name[i] + ', m6A smoothed prob > '+ str(thr) + ', clustered'), fontsize=16)
    clstmap.savefig(outfolder + 'm6A_smoothed_prob_thr'+ str(thr) + '_clustered_4x_' + str(i) + '.png', dpi=300)
    
    clstmap.fig.show()


### Nucleosome caller
# The following function gaps_in_signal_map allows one to:
# find gaps in methylation (AT_mod_arr_postsmooth_gap), 
# filter out gaps within a given bp range (nucleosome_limits) to call them as nucleosomes,
# find positions of called nucleosomes along template (AT_mod_arr_postsmooth_nuc_map - start to end of nucleosome = 1)
# find dyad position of called nucleosome along template (AT_mod_arr_postsmooth_nuc_dyad)
# count number of called nucleosomes per read (AT_mod_arr_postsmooth_nuc_count), 
# find size of called nucleosomes (AT_mod_arr_postsmooth_nuc_size)

nucleosome_limits = (100,180)

def gaps_in_signal_map(mod_array, nucleosome_size = nucleosome_limits):
    nuc_count = 0
    nuc_map = np.full(len(mod_array),False)
    nuc_dyad_pos = []
    nuc_size = []
    switch = np.subtract(mod_array[1:len(mod_array)], mod_array[0:len(mod_array)-1])
    pos1 = np.append(np.where(switch == 1),[6301])
    pos2 = np.append([0],np.where(switch == -1))
    if np.shape(pos1) == np.shape(pos2):
        gap = np.subtract(pos1,pos2).tolist()
    else:
        print('ignoring read ' + str(j) + ' from barcode ' + str(i))
        gap = np.subtract(pos2,pos2).tolist()
    for k in range(0,len(gap)):
        space = gap[k]
        if space > nucleosome_size[0] and space < nucleosome_size[1]:
            nuc_start = pos2[k]
            nuc_end = pos1[k]
            nuc_dyad = nuc_start + math.floor(space/2)
            nuc_map[nuc_start:nuc_end].fill(True)
            nuc_dyad_pos.append(nuc_dyad)
            nuc_count+= 1
            nuc_size.append(space)
    return gap, nuc_map, nuc_dyad_pos, nuc_count, nuc_size

AT_mod_arr_postsmooth_gap = {newlist:[] for newlist in Barcodes_subset}
AT_mod_arr_postsmooth_nuc_map = {newlist:[] for newlist in Barcodes_subset}
AT_mod_arr_postsmooth_nuc_count = {newlist:[] for newlist in Barcodes_subset}
AT_mod_arr_postsmooth_nuc_dyad = {newlist:[] for newlist in Barcodes_subset}
AT_mod_arr_postsmooth_nuc_size = {newlist:[] for newlist in Barcodes_subset}
for i in Barcodes_subset:
    for j in range(0,numreadsmax):
        gap, nuc_map, nuc_dyad_pos, nuc_count, nuc_size = gaps_in_signal_map(AT_mod_arr_postsmooth[i][j], nucleosome_limits)
        AT_mod_arr_postsmooth_gap[i].append(np.asarray(gap).reshape(1,len(gap)))
        AT_mod_arr_postsmooth_nuc_map[i].append(nuc_map)
        AT_mod_arr_postsmooth_nuc_dyad[i].append(nuc_dyad_pos)
        AT_mod_arr_postsmooth_nuc_count[i].append(nuc_count)
        AT_mod_arr_postsmooth_nuc_size[i].append(nuc_size)

## Plotting the offset (in bp) between estimated dyad from DiMeLo-seq and theoretical dyad axis of 601 positioning sequence

called_dyad_list_offset = {x : [] for x in Barcodes_subset}
for bar in Barcodes_subset:
    fig = plt.figure(figsize=(3,3))
    called_dyad_list = np.ndarray.tolist(np.concatenate(AT_mod_arr_postsmooth_nuc_dyad[bar], axis=0))
    called_dyad_list_reset = [a - dyad_list[0]-100 for a in called_dyad_list]
    called_dyad_list_offset[bar] = [a%200-100 for a in called_dyad_list_reset]
    plt.hist(called_dyad_list_offset[bar], bins=200, density=True, color = '#053C5E', histtype='stepfilled')
    plt.axvline(0, 0,1, color='r', alpha=0.7, linestyle='--', linewidth=3)
    plt.title(Barcode_name[bar], fontsize = 16)
    plt.xlabel('distance from 601 dyad (bp)', fontsize = 14)
    plt.ylabel('nucleosome density', fontsize = 14)
    plt.xlim(-100,100)
    plt.ylim(0,0.08)
    plt.show()
    fig.savefig(outfolder + 'called_dyad_offset_' + str(bar) + '.png', dpi=300, bbox_inches='tight')


dyad_list = np.arange(549,4049,200)

# The following lines of code allows to
# 1. identify/filter nucleosomes that do not have adjacent nucleosomes (lonesome nucleosome)
# 2. calculate and plot average probability score of methylation centered along each lonesome nucleosome 

## function for calculating sum of called nucleosome dyad centered signal around each identified dyad axis
def dyad_centered_sum_NaN(seq_array, dyad_pos = dyad_list, half_range=200):
    window_size = half_range*2+1
    collapsed_array_sum_NaN = np.zeros(window_size)
    for pos in dyad_pos:
        repeat_start = pos-half_range
        repeat_end = pos+half_range+1
        for i in range(repeat_start,repeat_end):
            j = (i-repeat_start)%(window_size)
            if math.isnan(seq_array[i]) == False:
                collapsed_array_sum_NaN[j] += seq_array[i]
    
    return collapsed_array_sum_NaN

## Filtering called nucleosomes with no immediately adjacent nucleosomes

def lonesome_nuc_filter(list_of_dyad):
    lonesome_nuc_list = []
    list_of_dyad = [249] + list_of_dyad + [4249]
    for i in range(1,len(list_of_dyad)-1):
        left_dyad_separation = list_of_dyad[i]-list_of_dyad[i-1]
        right_dyad_separation = list_of_dyad[i+1]-list_of_dyad[i]
        if left_dyad_separation >= 400:
            if right_dyad_separation >= 400:
                lonesome_nuc_list.append(list_of_dyad[i])
    return lonesome_nuc_list

AT_mod_arr_postsmooth_nuc_dyad_lonesome = {bar:[] for bar in Barcodes_subset}
for i in Barcodes_subset:
    for read in AT_mod_arr_postsmooth_nuc_dyad[i]:
        AT_mod_arr_postsmooth_nuc_dyad_lonesome[i].append(lonesome_nuc_filter(read))


# counting number of lonesome nucleosomes (for calcualting average methylation probability)
AT_mod_arr_postsmooth_nuc_count_lonesome = {bar:[] for bar in Barcodes_subset}
for i in Barcodes_subset:
    AT_mod_arr_postsmooth_nuc_count_lonesome[i] = len(np.concatenate(AT_mod_arr_postsmooth_nuc_dyad_lonesome[i]))
    print(AT_mod_arr_postsmooth_nuc_count_lonesome[i], Barcode_name[i])


# Calculating average methylation probability at each position centered at the called (lonesome) nucleosome dyad axis

AT_mod_arr_postsmooth_nuc_dyad_lonesome_collapsed_NaN = {newlist:[] for newlist in Barcodes_subset}
half_window = 300
for i in Barcodes_subset:
    dyad_collapsed_sum = np.zeros(2*half_window+1)
    num_nuc = 0
    for j in range(0,2000):
        if len(AT_mod_arr_postsmooth_nuc_dyad_lonesome[i][j]) > 0:
            dyad_collapsed_read = dyad_centered_sum_NaN(smoothed_AT_prob_array[i][j],AT_mod_arr_postsmooth_nuc_dyad_lonesome[i][j],half_range=half_window)
            dyad_collapsed_sum = np.add(dyad_collapsed_sum,dyad_collapsed_read)
    AT_mod_arr_postsmooth_nuc_dyad_lonesome_collapsed_NaN[i] = dyad_collapsed_sum/AT_mod_arr_postsmooth_nuc_count_lonesome[i]

# Plotting average methylation probability along each called (lonesome) nucleosome

x = np.arange(-half_window, half_window+1)
fig = plt.figure(figsize=(4,2))
Barcodes_str = ''
a = {1:'#559CAD', 2:'#053C5E'}
for i in [1,2]:
    Barcodes_str += '_' + str(i)
    plt.plot(x,AT_mod_arr_postsmooth_nuc_dyad_lonesome_collapsed_NaN[i]/256, label=Barcode_name[i], color = a[i])
plt.ylabel('probability of methylation', fontsize = 14)
plt.axvline(0,0,1, color='red', linestyle='--')    
legend = plt.legend(fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel('bp from dyad axis', fontsize = 14)
plt.show()
fig.savefig(outfolder + 'reach_lonesome'+ Barcodes_str + '.png', dpi=300, bbox_inches='tight')

def export_legend(legend, filename="legend.png"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi=300, bbox_inches=bbox)

export_legend(legend, str(outfolder + 'reach_lonesome_legend'+ Barcodes_str + '.png'))


