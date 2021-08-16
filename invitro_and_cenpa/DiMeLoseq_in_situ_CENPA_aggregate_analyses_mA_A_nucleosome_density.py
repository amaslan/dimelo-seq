### README ###

# #Use this script for analyses of all reads that align to a region of interest from DiMeLo-seq in vivo or in situ, for plotting aggregate mA/A,
# also outputs nucleosome density within rolling window for region of reads that map within the window.
# Input - hybrid Guppy + winnowmap bams processed to csv files (that are output from DiMeLoseq_in_situ_roi_single_molecule_v7_q10.py) containing the following:
# ReadID | array of NaN's or probability score of methylation where, NaN ==> non-A basepair position or not covered, probability score 0 - 255

# First part of the code reads in csv file, calculates average mA/A (with binning for ex. by averaging over 1kb rolling window for better visualization),
# calculates coverage (number of reads aligned to each base position inferred from A positions i.e. non-NaNs in input csv), plots the results

# Second part of the code uses the following formula to bin methylation probability scores
# p(at least 1 mA) = 1 - product(1-p(mA)) where product is calculated over all elements of rolling bin
# then, bins that satisfy the condition p(at least 1 mA) > 0.9 is considered as a methylated bin
# nucleosome density within a rolling 5kb window is then calcuated as (total number of methylated bins)/200 (assuming 200 bp maximum footprint of a nucleosome)
# For the nucleosome density measurements, only reads that map to > 1kb within a given rolling window is considered
# This part of the code plots mean nucleosome density at each rolling 5kb window


### The following filepaths and variables need to be updated for running the script

#Insert barcode name (to iterate over different barcodes)
#Insert modification - 'A' and/or 'C'
#for example:
mod = ['A']
Barcode_name = {18:'CA'}
filename_prob = <filepath> + str(Barcode_name[i]) + '_' + str('A') + '_prob.csv'
i = 18 #alternately, iterate over multiple barcodes where required
j = 'A' #alternately, iterate over 'A' and 'C' where required

outfolder = <output_directory_path>

#Variables to update for part 1
thr_raw = 230 #Threshold of p(mA) for aggregate mA/A plots, 0-255

#Variables to update for part 2.
# change bin size and threshold for p(at least 1 mA within bin) for nucleosome density measurements
prob_bin_size = 130 #number of bp to be considered for rolling window while binning
thr = 0.9 #threshold of p(at least 1 mA)


#Update the following variables for calcualting and plotting nucleosome density for a region of interest (0 - size of chromosome region of interest)
roi_start = 0 
roi_end = window[1]-window[0]

#update the following variables for different rolling window size and step sizes, as well as minimum read alignment 
#for calculating and plotting nucleosome density
rolling_window_size = 5000 
step_size = 2500 #Make sure this is smaller than window size
minimum_coverage = 1000 #Make sure this is smaller than window size

## Enter the correct bed start, bed end and windowsize (same as in the script used to generate the input csv file) info below
#for example:
Chromosome_Region = 'ChrX_CDR_300kb_5745_5775'
bedstart = 57450000
bedend = 57750000
windowsize = 150000

import numpy as np
import math

windowstart = math.floor((bedstart+bedend)/2-windowsize)
windowend = math.floor((bedstart+bedend)/2+windowsize)
#window = (int(windowstart), int(windowend))
#Alternately, provide window as a subset of the total window:
window = (int(windowstart),int(windowend)-50000)


###

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#from scipy.stats import kde
get_ipython().magic('matplotlib inline')
from matplotlib import colors


from scipy import sparse

mod_csv = {i:{} for i in Barcode_name.keys()}
read_prob = {i:{} for i in Barcode_name.keys()}
read_pos = {i:{} for i in Barcode_name.keys()}

#reading in prob score csv
#input format: ReadID | array of NaN's or probability score of methylation where, NaN ==> non-A basepair position or not covered, probability score 0 - 255
#csv shape = (num of reads, size of chromosomal region + 1(read id))

read_prob[i][j] = np.genfromtxt(filename_prob, dtype = float, delimiter = ',')[1:,1:]

#reading all A positions in reads with respect to chromosome region of interest as any base position with p(mA) >= 0
read_pos[i][j] = np.empty(np.shape(read_prob[i][j]))
read_pos[i][j][(read_prob[i][j])>=0] = 1

#determining coverage from A positions. First and last A positions on each read (with respect to chromosome region of interest) is determined
read_start = {i:{} for i in Barcode_name.keys()}
read_end = {i:{} for i in Barcode_name.keys()}
for j in mod:
    for i in Barcode_name.keys():
        read_start[i][j] = [np.nonzero(read)[0][0] for read in np.nan_to_num(read_pos[i][j])]
        read_end[i][j] = [np.nonzero(read)[0][-1] for read in np.nan_to_num(read_pos[i][j])]

#Filling in 1's between first and last A position on each read (with respect to chromosome region of interest)
#total_coverage is calculated by summing over all reads for each position in chromosome region of interest
read_coverage = {i:{} for i in Barcode_name.keys()}
total_coverage = {i:{} for i in Barcode_name.keys()}
for j in mod:
    for i in Barcode_name.keys():
        windowlen = np.shape(read_pos[i][j])[1]
        read_coverage[i][j] = []
        for k in range(0,np.shape(read_pos[i][j])[0]):
            read_coverage[i][j].append(np.zeros(windowlen))
            read_coverage[i][j][k][read_start[i][j][k]:read_end[i][j][k]+1] = 1
        total_coverage[i][j] = np.sum(read_coverage[i][j], axis=0)


### Part 1: Plotting average mA/A after binning over 1 kb
read_mod = {i:{} for i in Barcode_name.keys()}
binned_average_read_mod_thr_1kb= {i:{} for i in Barcode_name.keys()}
for j in mod:
    for i in Barcode_name.keys():
        read_mod[i][j] = np.empty(np.shape(read_prob[i][j]))
        read_mod[i][j][:] = np.NaN
        read_mod[i][j][read_prob[i][j] > thr_raw] = 1
        binned_average_read_mod_thr_1kb[i][j] = pd.DataFrame.to_numpy(pd.DataFrame(np.nan_to_num(np.nanmean(read_mod[i][j], axis = 0))).rolling(window = 1000).mean())

x = np.arange(windowstart,windowend+1)

for j in mod:
    if j == 'A':
        colorpal = '#053C5E'
    elif j =='C':
        colorpal = '#BB4430'
    for i in Barcode_name.keys():
        a = 0
        fig,ax = plt.subplots(2,1, figsize = (18,1.5), gridspec_kw={'height_ratios': [2, 1]})
        ax[a].plot(x,binned_average_read_mod_thr_1kb[i][j],colorpal)
        #ax[a].set_title(Barcode_name[i] + ' m6A', loc='left', pad=1, fontsize = 12)
        if j == 'A': 
            ax[a].set_ylim(0,0.075)
        elif j == 'C':
            ax[a].set_ylim(0,0.03)
        ax[a].set_xlim(window[0],window[1])
        a += 1
        #ax[a].plot(x,total_coverage[i][j], 'grey')
        ax[a].fill_between(x,np.nansum(read_pos[i][j],axis = 0), color = 'grey')
        #ax[a].set_title(Barcode_name[i] + ' coverage', loc='left', pad=1, fontsize = 12)
        #ax[a].axhline(5, color = 'r', linestyle = '--', alpha=0.5, linewidth = 1)
        ax[a].set_ylim(0,20)
        ax[a].set_xlim(window[0],window[1])
        fig.savefig(outfolder + Chromosome_Region + '_250kb_bin1kb_' + Barcode_name[i] + '_' + j + '_' + str(window[0]) + '_' + str(window[1]) + '.png', dpi=300, bbox_inches='tight')
        plt.show()


### Part 2: finding rolling 'prob_bin_size' bp bins with probability of at least 1 mA > 'thr', calculating nucleosome density per read

#binning probability and thresholding p(at least 1 mA) > 0.9 (for example)

thrq = 1- thr
logthrq = np.log(thrq)
read_logq = np.log(1-read_prob[i][j][:,:window[1]-window[0]]/256)
read_logq_nancumsum = np.pad(np.nancumsum(read_logq, axis = 1), ((0,0),(int(np.floor(prob_bin_size/2)),int(np.floor(prob_bin_size/2)))), constant_values=np.NaN)
read_sumlogq_rolling = read_logq_nancumsum[:,prob_bin_size:window[1]-window[0]+prob_bin_size] - read_logq_nancumsum[:,0:window[1]-window[0]]
read_prob_binned_array = read_sumlogq_rolling < logthrq

#function for counting nucleosomes per read within region of interest, maximum nucleosome footprint assumed to be 200 bp
def find_nuc_dist(mod_array_input):
    mod_array = 1-np.pad(mod_array_input, 1)
    mod_pos = np.nonzero(mod_array)[0]
    mod_switch = np.diff(mod_pos)-1
    mod_switch = mod_switch[np.nonzero(mod_switch)]
    num_nuc_sized = np.nansum(mod_switch)/200
    return num_nuc_sized

##This portion calculates the nucleosome density witin specific regions of interest with minimum coverage per read within each rolling window

read_count = []
nuc_count_perread = []
window_mid = []
nuc_list_mean = []
nuc_list_median = []
nuc_list_std = []
num_reads = np.shape(read_coverage[i][j])[0]
sized_nuc_count_perread = []
for roll_start in range(roi_start, roi_end, step_size):
    rolling_window_start = int(roll_start + 0)
    rolling_window_end = int(roll_start + rolling_window_size)
    rolling_window_mid = int(rolling_window_start + math.floor(rolling_window_size/2))
    covered = np.array(read_coverage[i][j])[:,rolling_window_start:rolling_window_end].sum(axis = 1)>minimum_coverage
    readlength = []
    gap_list = []
    per_read_total_meth = []
    nuc_list = []
    for row in range(0,num_reads):
        if covered[row] == True:
            readlength.append(np.array(read_coverage[i][j])[row,rolling_window_start:rolling_window_end].sum())
            per_read_total_meth.append(np.nansum(read_prob_binned_array[row, rolling_window_start:rolling_window_end]))
            nuc_list.append(find_nuc_dist(read_prob_binned_array[row, rolling_window_start:rolling_window_end]))
             
    sized_nuc_count_perread.append(np.nanmean(per_read_total_meth)/200)
    window_mid.append(rolling_window_mid)
    nuc_list_mean.append(np.nanmean(nuc_list))    


#Plotting nucleosomes counted as (total methyl bin size / 200) for each read
#for example, 5kb window, 1kb step size, read coverage > 1000 within each region; Binning threshold 0.9, binsize = 130
fig = plt.figure(figsize = (18,1.5))
plt.fill_between(np.asarray(window_mid)-np.floor(rolling_window_size/2)+windowstart,nuc_list_mean, color = '#053C5E')
plt.xlim(0+windowstart,window[1]-window[0]+windowstart)

fig.savefig(outfolder + Chromosome_Region + 'nucdensity' + Barcode_name[i] + '.png', dpi=300, bbox_inches='tight')

plt.show()

