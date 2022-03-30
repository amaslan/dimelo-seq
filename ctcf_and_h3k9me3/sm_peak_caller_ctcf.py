# ### Single-molecule peak calling to determine for reads with a peak, what is the distance of the
# ### predicted peak from the motif center
# #### Annie Maslan
# #### 11.27.21

import os
import numpy as np
import pandas as pd

# plotting & display imports
import matplotlib.pyplot as plt
import seaborn as sns

import scipy as sp
from scipy import signal
from scipy.optimize import curve_fit
import scipy.stats
import itertools



out = 'out'

ctcf_mA = pd.read_csv('top_decile_mA.csv') # similarly for untreated for control; output from single_molecule_roi_meg_winnow_guppy_ctcf_decile.py

print(len(ctcf_mA['read_name'].unique()))
print('ctcf A: ' + str(ctcf_mA.shape[0]))


def bin_qualities(all_data, bin_size):
    ''' 
    bin data for bin_size A's as slide across the read
    calculate the probability at least one base is methylated
    '''
    all_data.loc[:, 'prob'] = all_data['prob'] / 255
    read_names = all_data['read_name'].unique()
    for r in read_names:
        subset = all_data[all_data['read_name'] == r] 
        qualities = subset['prob'] 
        binned_qualities = qualities.rolling(window=bin_size, min_periods=1, center=True).apply(lambda b: prob_bin(b))
        all_data.loc[all_data['read_name'] == r,'prob'] = binned_qualities 
    return all_data


def prob_bin(bin):
    # probability a base in the window is methylated by:
    # calculating probability that no base in the window is methylated and then taking the complement
    # treat p=1 as 254/255 for prevent log(0)
    probs = [np.log(1-p) for p in bin if ((p < 1) and (p > 0.5))] # go back to 0.5
    probs1 = [np.log(1-254/255) for p in bin if p == 1]
    probsAll = probs + probs1
    prob = 1 - np.exp(sum(probsAll)) 
    return prob


# only keep reads that do span the motif center with at least 100 bp on each side
ctcf_mA_left = ctcf_mA[(ctcf_mA['pos'] >= -100) & (ctcf_mA['pos'] <= 0)]
ctcf_mA_right = ctcf_mA[(ctcf_mA['pos'] <= 100) & (ctcf_mA['pos'] >= 0)]
read_L = ctcf_mA_left['read_name'].unique()
read_R = ctcf_mA_right['read_name'].unique()


# only keep reads with at least one mA call about threshold
ctcf_thresh75 = ctcf_mA[ctcf_mA['prob'] >= 230] # 0.9 very confident
reads_t75 = ctcf_thresh75['read_name'].unique()
print(len(reads_t75))


reads_keep = []
for r in read_L:
    if r in read_R:
        if r in reads_t75:
            reads_keep.append(r)
print(len(reads_keep))


# filter full dataframe to only keep reads that span the motif center and have at least one mA call above threshold
boolean_keep_series = ctcf_mA.read_name.isin(reads_keep)
ctcf_mA_filtered = ctcf_mA[boolean_keep_series]


bin_size_1 = 20 # bin size of As to start



binned = bin_qualities(ctcf_mA_filtered, bin_size_1)



# only keep reads that did call a peak above threshold
reads_peak = binned[binned['prob'] >= 0.8]['read_name'].unique() # was 0.9

binned_t = binned[binned['prob'] >= 0.8]

print('ctcf mA: ' + str(binned_t.shape[0]))
print(len(reads_peak))

# also only keep reads that had peak called in a dense neighborhood of mA
# take longest continuous stretch > thresh 95 for a given read
# then take midpoint of that stretch

positions = []
read_names = binned['read_name'].unique()
for r in read_names: 
    if r in reads_peak:
        sub = binned[binned['read_name'] == r]
        sub = sub.reset_index()
        s = sub['prob'] >= 0.8
        l = []
        for i, (j, k) in enumerate(itertools.groupby(s)):
            l.append([i, j, len(list(k))])
        trues = [x[2] if x[1] else 0 for x in l]
        locs = [x[2] for x in l]
        max_trues_len = max(trues)
        max_trues_start = np.sum(locs[:trues.index(max_trues_len)])
        center = int(max_trues_start + round(max_trues_len/2))
        
        assert s[center], 'Error. False.'

        p = sub.iloc[center,:]['pos']
        positions.append(p)


fig = plt.figure()
sns.distplot(positions, color='#053C5E');
plt.show()
plt.savefig(out + '/diff_distribution.pdf', dpi=600)


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a) 
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

print(mean_confidence_interval(positions, confidence=0.95))
print(mean_confidence_interval(positions, confidence=0.90))

series_pos = pd.Series(positions)
quants = series_pos.quantile([0.05, 0.95])
print(quants)

series_pos = pd.Series(positions)
quants = series_pos.quantile([0.025, 0.975])
print(quants)

series_pos = pd.Series(positions)
quants = series_pos.quantile([0.15, 0.85])
print(quants)

series_pos.to_csv(out + '/positions.csv')

# fit and make plot with t

tdof, tloc, tscale = scipy.stats.t.fit(positions)
print(tdof, tloc, tscale)

fig, ax = plt.subplots(1, 1)
x = np.linspace(scipy.stats.t.ppf(0.01, tdof, loc=tloc, scale=tscale), scipy.stats.t.ppf(0.99, tdof, loc=tloc, scale=tscale), 2000)
ax.hist(positions, color='#053C5E', bins=100, density=True);
ax.plot(x, scipy.stats.t.pdf(x, tdof, loc=tloc, scale=tscale), lw=1, label='t pdf', color='#FFBC0A');
plt.xlim(-1000, 1000);
plt.savefig(out + '/diff_distribution_fit.pdf', dpi=600)
print(scipy.stats.t.interval(0.7, tdof, loc=tloc, scale=tscale))

# new analysis just using out data
positions = pd.read_csv(out + '/positions.csv', index_col=0)
positions.columns = ['pos']
quants = positions['pos'].quantile([0.2, 0.8])
print(quants)

# fraction of called peaks that are within 100 bp of motif center
# positions size: 27407
positions[abs(positions['pos']) <= 100].shape
# (11447, 1)
positions[abs(positions['pos']) <= 200].shape
# (13897, 1) # 50% of predicted peaks are within 200 bp of motif center
positions[abs(positions['pos']) <= 300].shape
# (16685, 1) # 61%

print(mean_confidence_interval(positions['pos'], confidence=0.95))
# (-4.6338891524063195, -9.720810510628882, 0.4530322058162426)

fig, ax = plt.subplots(1, 1)
ax.hist(positions['pos'], color='#053C5E', bins=100, density=True);
plt.xlim(-1000, 1000);
plt.savefig(out + '/diff_distribution_nofit.pdf', dpi=600)



