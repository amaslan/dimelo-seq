from pybedtools import BedTool
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt


BED_DIR = "/directory/with/bedfiles/by/decile"

beds = [BED_DIR + '/intersection.motifs.chip.formatted.chm13.q1.bed',
BED_DIR + '/intersection.motifs.chip.formatted.chm13.q2.bed',
BED_DIR + '/intersection.motifs.chip.formatted.chm13.q3.bed',
BED_DIR + '/intersection.motifs.chip.formatted.chm13.q4.bed',
BED_DIR + '/intersection.motifs.chip.formatted.chm13.q5.bed',
BED_DIR + '/intersection.motifs.chip.formatted.chm13.q6.bed',
BED_DIR + '/intersection.motifs.chip.formatted.chm13.q7.bed',
BED_DIR + '/intersection.motifs.chip.formatted.chm13.q8.bed',
BED_DIR + '/intersection.motifs.chip.formatted.chm13.q9.bed',
BED_DIR + '/intersection.motifs.chip.formatted.chm13.q10.bed']

names = ['ctcf_q1', 'ctcf_q2', 'ctcf_q3', 'ctcf_q4', 'ctcf_q5', 'ctcf_q6', 'ctcf_q7', 'ctcf_q8', 'ctcf_q9', 'ctcf_q10']

bams = ["deep_ctcf_mod_mappings_merge.sorted.bam"]

out = "out/dir"

def calc_frac_peak(data):
    data_peak = data[(data['pos'] >= -150) & (data['pos'] <= 150)] 
    data_peak_t75 = data_peak[data_peak['prob'] >= 190]
    cnt = data_peak_t75['read_name'].value_counts()
    pos = len(cnt.values[cnt.values >= 6]) # require at least 6mA
    return len(data_peak['read_name'].unique()), pos

# only looking +/- 100 on each side of the motif so only need to extract bases in that range
all_data = []
for b, n in zip(beds[::-1], names[::-1]):
    all_data.append(pd.read_csv("qX.csv")) # do for each decile output

reads = []
peak_reads = []
frac = []
deciles = [1,2,3,4,5,6,7,8,9,10]
for d in all_data:
    r, pr = calc_frac_peak(d)
    reads.append(r)
    reads.append(pr)
    frac.append(pr/r)
    print(pr)

print(frac)

# repeat for igg
bams = ["/clusterfs/rosalind/groups/streetslab/amaslan/nanopore/paper/ctcf/bams/megalodon/prod_IgG_mod_mappings.sorted.bam"]
names = ['igg_q1', 'igg_q2', 'igg_q3', 'igg_q4', 'igg_q5', 'igg_q6', 'igg_q7', 'igg_q8', 'igg_q9', 'igg_q10']
all_data_igg = []
for b, n in zip(beds[::-1], names[::-1]):
    all_data_igg.append(pd.read_csv("qX.csv"))

reads = []
peak_reads = []
frac_igg = []
for d in all_data_igg:
    r, pr = calc_frac_peak(d)
    reads.append(r)
    reads.append(pr)
    frac_igg.append(pr/r)
    print(pr)

print(frac_igg)

fig = plt.figure()
plt.plot(deciles[::-1], frac, 'o', color='#053C5E')
plt.plot(deciles[::-1], frac_igg, 'o', color='#5E747F');
plt.ylabel('Fraction of reads with mA in peaks')
plt.xlabel('Decile of ChIP-seq peak signal')
plt.legend(['CTCF-targeted', 'IgG control'])
fig.savefig(out + '/6mA/decile_analysis.pdf')
