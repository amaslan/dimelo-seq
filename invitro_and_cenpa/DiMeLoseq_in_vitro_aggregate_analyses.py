### README ###

###This file is obsolete and corresponds to the initial submission of the DiMeLoseq manuscript and BioRxiv preprint.###

# #Use this script for plotting in vitro DiMeLoseq data to show aggregate mA/A information
# plotting histogram of mA/A per read for all libraries
# plotting average mA/read for each A position along the template
#Input Megaladon or Guppy with modified base calling bam files - mod_mappings.bam (filter by read length > 3.6 kb (18x601) to keep file size small)
#Template - 18x601.fa or fasta file containing sequence of the chromatin or naked DNA template
#Threshold of probability score for calling methylated A can be adjusted in parse_bam6_fast function. Currently set to 0.9

### Script ###
# The following filepaths and variables need to be updated for running the script

outfolder = '<output_directory_path>'

#enter name of folder that contains guppy or megalodon bam file.
folder_name = <input_directory_path>


#List Barcode titles to be used for graphs here:
#For example:
Barcode_name = {}
Barcode_name[1] = 'CA chr + free PA-Hia5'
Barcode_name[2] = 'CA chr + CA directed PA-Hia5'
Barcode_name[7] = 'H3 chr + free PA-Hia5'
Barcode_name[8] = 'H3 chr + CA directed PA-Hia5'

#Select barcodes to iterate over for analyses
Barcode_list = [1,2,7,8]


#reading in template into string - For in vitro analyses, we use 18x601.fa
#for example, array_seq = SeqIO.read("ASP696_18x601.fa", "fasta")
array_seq = SeqIO.read(<template_sequence.fa>, "fasta")
array_str = str(array_seq.seq)

###

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

get_ipython().magic('matplotlib inline')

import math

import pysam
from Bio import SeqIO
import time


#identifying indices of A's and T's for assigning mod probs to right positions

A_pos = []
T_pos = []
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

# m6A_seqmap_prob_fast - takes in A or T pos array and readmap start, returns mod prob values assigned to A's along the template
# input arguments - A pos array, T pos array, pos of read along template, if strand=rev, Mm_start, truncated Ml prob values, read length, template length) thresholding done later
# output - m6A prob score along template, coverage along template - for each read

def m6A_seqmap_prob_fast(A_pos_array, T_pos_array, Readmap_start, rev_strand, Mm_start, Ml_Aa, readlength, templatelength = 6302):
    
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
        AT_prob_read = np.zeros(templatelength,dtype=int)
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
        AT_prob_read = np.zeros(templatelength,dtype=int)
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
##coverage_Fwd_bar, coverage_Rev_bar - arrays of coverage(sum of reads aligned to template at each position) for forward and reverse reads in library
##A_n_bar - number of A's per read, used to calculate mA/A
##readids_bar - read ids for all reads in given library
##rev_read_list_bar - boolean array indicating strand of all reads, shape = 1x(number of reads in libary), False - forward read, True - reverse read

def parse_bam6_fast(filename, numreads = 5000, templatelength = 6302):
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
                
                    AT_prob_read, AT_pos_read = m6A_seqmap_prob_fast(A_pos, T_pos, pos, rev_strand, Mm_Aa_start, Ml_Aa, length)
                
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

def by_barcode(folder_path, barcodes_list = Barcode_list, numreads = 5000):
    AT_prob_arr_bar = {new_list: [] for new_list in barcodes_list} 
    coverage_Fwd_bar = {new_list: [] for new_list in barcodes_list} 
    coverage_Rev_bar = {new_list: [] for new_list in barcodes_list} 
    A_n_bar = {new_list: [] for new_list in barcodes_list} 
    readids_bar = {new_list: [] for new_list in barcodes_list} 
    rev_read_list_bar = {new_list: [] for new_list in barcodes_list} 
    filename = {}
    for i in barcodes_list:
        start_time = time.time()
        filename[i] = folder_path + 'bar' + str(i) + '.1/mod_mappings_3.6kb.bam'
        AT_prob_arr_bar[i], coverage_Fwd_bar[i], coverage_Rev_bar[i], A_n_bar[i], readids_bar[i], rev_read_list_bar[i] = parse_bam6_fast(filename[i], numreads)
        
        print('Finished barcode ', i, 'in ', "--- %s seconds ---" % (time.time() - start_time))
        
    return AT_prob_arr_bar, coverage_Fwd_bar, coverage_Rev_bar, A_n_bar, readids_bar, rev_read_list_bar

AT_prob_arr_bar, coverage_Fwd_bar, coverage_Rev_bar, A_n_bar, readids_bar, rev_read_list_bar = by_barcode(folder_name, Barcode_list, 10000)


#COVERAGE - for calculating coverage (no. of reads that map to each position along the template) for all the barcodes
def coverage(barcodes_list, coverage_Fwd_bar, coverage_Rev_bar, templatelength = 6302, label_name = Barcode_name):
    coverage_bar = {new_list: [] for new_list in barcodes_list}
    x = np.arange(templatelength)    
    for i in barcodes_list:
        coverage_bar[i] = coverage_Fwd_bar[i] + coverage_Rev_bar[i]
    return coverage_bar


coverage_bar = coverage(Barcode_list, coverage_Fwd_bar, coverage_Rev_bar)

#Thresholding methylation probability scores, and converting mA information to binary for each read
def binary_mod(barcodes_list, AT_prob_arr_bar, coverage_bar, AT_mod_arr_bar={}, thr=90):
    thr256 = thr/100*256
    try: AT_mod_arr_bar
    except NameError: 
        AT_mod_arr_bar = {new_list: [] for new_list in barcodes_list}
    for i in barcodes_list:
        AT_mod_arr_bar[i] = np.ndarray.astype(AT_prob_arr_bar[i] >= thr256, int)
    return AT_mod_arr_bar

AT_mod_arr_bar = binary_mod(Barcode_list, AT_prob_arr_bar, coverage_bar)

#calculating fraction of A's with > threshold probability of methylation
def frac_m6A(barcodes_list, AT_mod_arr_bar, A_n_bar, m6A_over_A_bar = {}):
    N_mod_bar = {new_list: [] for new_list in barcodes_list}
    for i in barcodes_list:
        N_mod_bar[i] = AT_mod_arr_bar[i].sum(axis = 1)
        m6A_over_A_bar[i] = np.ndarray.astype(N_mod_bar[i]/A_n_bar[i], float)
    return m6A_over_A_bar


m6A_over_A_bar = frac_m6A(Barcode_list, AT_mod_arr_bar, A_n_bar)

#This function plots histogram of m6A/A for all reads by library
#cumul = False ==> probability distribution plot
#cumul = True ==> cumulative distribution plot

def plot_m6A_over_A(barcodes_list, m6A_over_A_bar, ylimit=(0,22), xlimit = (0,1), num_bins = 100, label_name = Barcode_name, cumul = False, colorpal = ['#FFBC0A','#053C5E','#559CAD','#610345','#A9E5BB','#2D1E2F','#BB4420','#5E747F']):

    def export_legend(legend, filename="legend.png"):
        fig  = legend.figure
        fig.canvas.draw()
        bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(filename, dpi=300, bbox_inches=bbox)
        return
    fig = plt.figure(figsize = (3.6,2.4))
    barcodes_list_str = ''
    a = 0
    for i in barcodes_list:
        barcodes_list_str += '_' + str(i)
        if cumul == True:
            ax = plt.hist(m6A_over_A_bar[i], bins = num_bins, density = True, alpha = 1, label = label_name[i], cumulative=cumul, histtype='step', rwidth = 1, color = colorpal[a])
            plt.ylabel('probability', fontsize = 14)
    
        else:
            if i == 2 or i==1:
                alp = 0.9
            else:
                alp = 0.7
            ax = plt.hist(m6A_over_A_bar[i], bins = num_bins, density = True, alpha = alp, label = label_name[i], color = colorpal[a])
            plt.ylabel('density', fontsize = 14)
        a += 1
    
    plt.xlabel('m6A/A (per read)', fontsize = 14)
    #plt.legend(fontsize=14)
    legend = plt.legend(fontsize=12, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.ylim(ylimit)
    plt.xlim(xlimit)
    plt.show()
    if cumul == True:
        fig.savefig(outfolder+'plot_m6A_over_A_Cumulative_'+barcodes_list_str, dpi = 300, bbox_inches='tight')
        export_legend(legend, str(outfolder + 'plot_m6A_overA_Cumulative_legend_'+ barcodes_list_str + '.png'))
    else:
        fig.savefig(outfolder+'plot_m6A_over_A_'+barcodes_list_str, dpi = 300, bbox_inches='tight')
        export_legend(legend, str(outfolder + 'plot_m6A_overA_legend_'+ barcodes_list_str + '.png'))
    
    return

#ax = plot_m6A_over_A([1,7,2,8],m6A_over_A_bar, ylimit = [0,1.05], xlimit = [0,0.55], num_bins=1000, cumul=True)
ax = plot_m6A_over_A([7,2,1,8],m6A_over_A_bar, ylimit = [0,22], xlimit = [0,0.55], num_bins=100)


#Function to calculate average mA/A for all reads at each A position along the template
def binary_mod_sum_norm(barcodes_list, AT_mod_arr_bar, AT_mod_sum_norm_bar={}):
    for i in barcodes_list:
        AT_mod_sum_norm_bar[i] = np.ndarray.astype((AT_mod_arr_bar[i].sum(axis=0))/(coverage_bar[i]), float)
    return AT_mod_sum_norm_bar

AT_mod_sum_norm_bar = binary_mod_sum_norm([1,2,7,8], AT_mod_arr_bar)

#Function to plot average mA/A over all reads for each A position in template
def plot_mod_sum_norm_squished(barcodes_list, AT_mod_sum_norm_bar, ylimit = (0,1), xlimit = (0,6302), label_name = Barcode_name, colorpal = ['#053C5E','#559CAD','#FFBC0A','#610345','#A9E5BB','#2D1E2F','#BB4420','#5E747F']):
    x = np.arange(6302)
    barcodes_list_str = ''
    a = 0
    fig = plt.figure(figsize=(4,0.5))
        
    for i in barcodes_list:
        barcodes_list_str += str(i) + '_'
        plt.plot(x, AT_mod_sum_norm_bar[i], alpha=1, color = colorpal[a])
        #plt.xlabel('bp of ASP 696 18x601', fontsize=16)
        plt.ylabel('m6A/read', fontsize=16)
        #plt.title(label_name[i], fontsize=18, loc='left')
        plt.ylim(ylimit)
        plt.xlim(xlimit)
        for mon in range(0,18):
            plt.axvline(549+200*mon, 0,0.75, color='red', linestyle='--', alpha = 0.5)
        a+=1
    plt.show()    
    fig.savefig(outfolder+'plot_mod_sum_norm_squished_'+barcodes_list_str, dpi = 300, bbox_inches='tight')
    
    return


plot_mod_sum_norm_squished([1,2], AT_mod_sum_norm_bar, (0,0.4),(445,4074))

plot_mod_sum_norm_squished([1], AT_mod_sum_norm_bar, (0,0.5),(445,4074))
plot_mod_sum_norm_squished([2], AT_mod_sum_norm_bar, (0,0.1),(445,4074))







