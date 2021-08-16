#!/usr/bin/env python
# coding: utf-8

# ### Analyze peak center resolution and methylation decay (separate L and R)
# ### At the strongest peaks, what's the achievable resolution
# #### Annie Maslan
# #### 06.27.21
# ##### Call peaks in quartile 4. For reads with peak, what is the distance from the motif center?

# In[26]:


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



out = 'results'

ctcf_mA = pd.read_csv('ctcf_allquart4_all_data_mA.csv')

print(len(ctcf_mA['read_name'].unique()))
print('ctcf A: ' + str(ctcf_mA.shape[0]))


def bin_qualities(all_data, bin_size):
    ''' 
    bin data for bin_size A's as slide across the read
    calculate the probability at least one base is methylated
    '''
    all_data.loc[:, 'quality'] = all_data['quality'] / 255
    read_names = all_data['read_name'].unique()
    for r in read_names:
        subset = all_data[all_data['read_name'] == r] 
        qualities = subset['quality'] 
        binned_qualities = qualities.rolling(window=bin_size, min_periods=1, center=True).apply(lambda b: prob_bin(b))
        all_data.loc[all_data['read_name'] == r,'quality'] = binned_qualities 
    return all_data


def prob_bin(bin):
    # probability a base in the window is methylated by:
    # calculating probability that no base in the window is methylated and then taking the complement
    # treat p=1 as 254/255 for prevent log(0)
    probs = [np.log(1-p) for p in bin if ((p < 1) and (p >= 0.5))] # only consider probabilities > 0.5 and handle 1 on next line
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
ctcf_thresh90 = ctcf_mA[ctcf_mA['quality'] >= 230]
reads_t90 = ctcf_thresh90['read_name'].unique()
print(len(reads_t90))


reads_keep = []
for r in read_L:
    if r in read_R:
        if r in reads_t90:
            reads_keep.append(r)
print(len(reads_keep))


# filter full dataframe to only keep reads that span the motif center and have at least one mA call above threshold
boolean_keep_series = ctcf_mA.read_name.isin(reads_keep)
ctcf_mA_filtered = ctcf_mA[boolean_keep_series]


bin_size_1 = 20 # bin size of As to start



binned = bin_qualities(ctcf_mA_filtered, bin_size_1)



# only keep reads that did call a peak above threshold
reads_peak = binned[binned['quality'] >= 0.9]['read_name'].unique()

binned_t = binned[binned['quality'] >= 0.9]

print('ctcf mA: ' + str(binned_t.shape[0]))

# also only keep reads that had peak called in a dense neighborhood of mA
# take longest continuous stretch > thresh 95 for a given read
# then take midpoint of that stretch

positions = []
read_names = binned['read_name'].unique()
for r in read_names: 
    if r in reads_peak:
        sub = binned[binned['read_name'] == r]
        sub = sub.reset_index()
        s = sub['quality'] >= 0.9
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
