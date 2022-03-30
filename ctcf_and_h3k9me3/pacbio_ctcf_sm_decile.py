# parse pacbio bam file for raw IPD and IPD ratio and number of passes
# Significant portions of code for bam parsing adapted from methplotlib (Copyright (c) 2018 Wouter De Coster)
# https://github.com/wdecoster/methplotlib

import pandas as pd
import pyranges as pr
import numpy as np
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing
from joblib import Parallel, delayed

class Region(object):
    def __init__(self, region):
        self.chromosome = region[1][0]
        # now being and end correspond to motif begin and end
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

def get_data(methylation_file, name, windows):
    """
    Import methylation data from all files in the list methylation_files
    Data in bam format
    data is extracted within the window 
    """
    num_cores = multiprocessing.cpu_count()
    print('executing with ' + str(num_cores) + ' cores')
    data = Parallel(n_jobs=num_cores)(delayed(parse_pb_bam)(methylation_file, name, w) for w in windows)
    return data


def parse_pb_bam(filename, name, window):
    '''
    parse mod_mappings.bam file to create methylation object with read_name, strand, pos, quality, and mod 
    in window 
    '''
    bam = pysam.AlignmentFile(filename, "rb")
    data = []
    for read in bam.fetch(reference=window.chromosome, start=window.begin, end=window.end):
        positions, iprs, ipds, passes = get_ipr_reference_positions(read, window)
        for pos, ipr, ipd in zip(positions, iprs, ipds):
            if pos is not None:
                if abs(pos) < 1000:
                    data.append((read.query_name,
                        '-' if read.is_reverse else '+',
                        pos,
                        ipr,
                        ipd,
                        passes))
    return Methylation(
        table=pd.DataFrame(data, columns=['read_name', 'strand', 'pos', 'ipr', 'ipd', 'np'])
                .sort_values(['read_name', 'pos']),
        data_type="pb-bam",
        name=name,
        called_sites=len(data))


def get_ipr_reference_positions(read, window):
    '''
    extract inter-pulse durations
    '''
    if read.has_tag('ir'):
        return get_pos_prob(read, window)
    else:
        return ([None], [None], [None], [None])


def get_pos_prob(read, window):
    '''
    get (position of modified base, IPD ratio)
    '''
    base = 'T'
    iprs = np.array(read.get_tag('ir'), dtype=float)
    ipds = np.array(read.get_tag('ip'), dtype=int)
    passes = read.get_tag('np')
    # return original read sequence; reads mapping to reverse strand will be reverse complemented
    seq = read.get_forward_sequence()
    base_index = np.array([i for i, letter in enumerate(seq) if letter == base])
    refpos = np.array(read.get_reference_positions(full_length=True))
    if read.is_reverse:
        refpos = np.flipud(refpos)
        # iprs = iprs[::-1] # native orientation
        # ipds = ipds[::-1] # native orientation
    keep = []
    # deal with None for refpos from soft clipped / unaligned bases
    # for m6A no need to look at neighboring base; do need to remove refpos that are None
    for b in base_index:
        if refpos[b] is not None:
            keep.append(b)
    # center reads (for now just assume all + strand because sensitivity analysis will not matter oriented to motif strand)
    refpos_mod_adjusted = np.array(refpos[keep]) - round(((window.end - window.begin) / 2 + window.begin))
    return (np.array(refpos_mod_adjusted), iprs[keep], ipds[keep], passes)

def main():
    out = '/out'
    bams = ['PB4.strandify.ipr.aligned.bam', 'PB8.strandify.ipr.aligned.bam']
    names = ['CTCF', 'untreated']
    bed = pd.read_csv('intersection.motifs.chip.formatted.chm13.q10.bed', sep='\t', header=None) # top decile of peaks
    windows = make_windows(bed)
    i = 0
    for bam in bams:
        data = get_data(bam, names[i], windows) 
        # combine all methylation tables into a single dataframe
        list_tables = []
        for d in data:
           list_tables.append(d.table)
        all_data = pd.concat(list_tables)
        #all_data = data.table
        all_data.to_csv(out + '/' + names[i] + '_pb_top_decile.csv', index=False)
        i = i + 1

if __name__ == '__main__':
    main()


    



    