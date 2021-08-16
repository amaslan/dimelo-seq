# Relevant scripts and files for CTCF and H3K9me3 analysis for creating Figures 4 & 5 and Supplementary Figures 4, 5, 6


## CTCF - Figure 4 and Supplementary Figures 4&5 

### Fig 4a, 4c, and Supplementary Figure 5.
#### **Script:**
single_molecule_roi_meg_winnow_guppy_ctcf_center.py
#### **Required files:**
mod_mappings.sorted.bam - from megalodon basecalling
beds/intersection.motifs.chip.formatted.chm13.bed

### Fig 4b.
#### **Script:**
ctcf_resolution.py
#### **Required files:**
After running script single_molecule_roi_meg_winnow_guppy_ctcf_center.py, use the output ctcf_allquart4_all_data_mA.csv as input to this script

### Fig 4d.
#### **Script:**
- single_molecule_joint_ctcf.py

#### **Required files:**
- mod_mappings.sorted.bam - from megalodon basecalling
- beds/intersection.motifs.chip.formatted.chm13.bed 


### Supplementary 4a
#### **Scripts:**
- ParseModMappings.hash&#46;mA&#46;pl - extract a list of all mA mod probabilities from the mod_mappings bam file output by megalodon
- run_mA_A_peaks_calcs_ctcf.sh - create bed files for computing enrichment for 4a
- enrichment.ipynb 

#### **Required files:**
- mod_mappings.sorted.bam
- beds/ENCFF797SDL.chm13.bed
- data/ctcf_enrich.csv - created from run_mA_A_peaks_calcs_ctcf.sh

### Supplementary 4b
#### **Scripts:**
ctcf_analysis.ipynb

#### **Required Files:**
After running script single_molecule_roi_meg_winnow_guppy_ctcf_center.py, use ctcf_allquart4_all_data_mA.csv as input


### Supplementary 4c

#### **Scripts:**
- single_molecule_roi_meg_winnow_guppy_ctcf_decile.py - analysis by decile
- ctcf_analysis.ipynb

#### **Required files:**
- mod_mappings.sorted.bam
- beds/intersection.motifs.chip.formatted.chm13.bed
- meg_CTCF_CTCF_qX_all_data_mA.csv, meg_IgG_CTCF_qX_all_data_mA.csv - use these files for q1 to q10 output from single_molecule_roi_meg_winnow_guppy_ctcf_decile.py for ctcf_analysis.ipynb



## H3K9me3 - Figure 5 and Supplementary Figure 6
### Preprocessing for bam files that are a hybrid from guppy and winnowmap. This guppy-winnowmap hybrid bam is used for H3K9me3 analysis.
See instructions here: hybrid_guppy_winnnowmap_bam_creation

### Figure 5a, 5b

#### **Scripts:**
- run_mA_A_peaks_calcs_h3k9me3.sh
- ParseModMappings.hash&#46;mA&#46;pl - extract a list of all mA mod probabilities from the hybrid guppy & winnowmap bam
- enrichment.ipynb

#### **Required files:**
- data/k9me3_enrich.csv - created from run_mA_A_peaks_calcs_h3k9me3.sh
- beds/HG2_K9Me3_ADGCCR09-H_S103_L002_2.chm13_HG002X_HG38Y_peaks.broadPeak - CUT&RUN
- beds/t2t_cenAnnotation.v2.021921.hg002_X.4.bed - centromeres
- beds/t2t_cenAnnotation.v2.021921.liveHORonly.HG002.4.sorted.bed - live HOR arrays
- franken.chromsizes.txt - chromosome sizes for franken reference

### Figure 5c, Supplementary Figure 6a
#### **Scripts:**
- single_molecule_roi_meg_winnow_guppy_hor_k9me3.py

#### **Required files:**
- beds/sharp_HOR_boundaries.300kb.bed

### Figure 5d, Supplementary Figure 6b & 6c
#### **Scripts:**
- single_molecule_plotly_k9me3.py
- For showing just a signle modification (mA or mCG in Supplementary Figure 6b, comment out traces.append() --> see "comment out for single mod plot")

#### **Required files:**
- guppy_winnow_hybrid.bam

#### **Other tools**
For aggregate curves:
- Plotted using the WashU epigenome browser - https://epigenomegateway.wustl.edu
- Bedgraph files for plotted created by outputting all_data_mA and all_data_mC dataframes from single_molecule_plotly_k9me3.py as csv's and then running the below for data_A and data_C analagously.
```
data_A = pd.read_csv('all_data_mA.csv')
starts = []
ends = []
totals = []
mAs = []
for r in data_A['read_name'].unique():
    sub = data_A[data_A['read_name'] == r]
    start = min(sub['pos'])
    end = max(sub['pos'])
    total = sub.shape[0]
    mA = sub[sub['quality'] > 230].shape[0]
    starts.append(start)
    ends.append(end)
    totals.append(total)
    mAs.append(mA)
dictionary = {'start': starts, 'end': ends, 'A': totals, 'mA': mAs}
read_bed = pd.DataFrame(dictionary)
read_bed['chr'] = 'chr7'
read_bed = read_bed[['chr', 'start', 'end', 'A', 'mA']]
read_bed.to_csv('reads_mA_cen7_h3k9me3.bed', sep='\t', header=False, index=False)
```

