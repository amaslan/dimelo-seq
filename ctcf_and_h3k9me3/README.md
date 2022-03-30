# Relevant scripts and files for CTCF and H3K9me3 analysis for creating Figures 4 & 5 and Extended Data Figures 5, 6, 7, 8, 9

Data used to generate figures are a merge of all bam files for a given target (CTCF, H3K9me3, IgG, Hia5) in the given cell type (GM12878, HG002).

## CTCF - Figure 4 and Extended Data Figures 5,6,7,8 

### Fig 4a, 4b; Extended Data Figure 5a, 5e; Extended Data Figure 6
#### **Script:**
- single_molecule_roi_meg_winnow_guppy_ctcf_center.py 
#### **Required files:**
- mod_mappings.sorted.bam - bam output from megalodon basecalling
- beds/intersection.motifs.chip.formatted.chm13.bed - bed file specifying CTCF motifs within CTCF ChIP-seq peaks (Fig 4a: quartiles 1,2,3,4; Fig 4b: top decile only; Extended Data Fig 5a,e: quartile 4)


### Fig 4c.
#### **Script:**
- single_molecule_joint_ctcf.py

#### **Required files:**
- mod_mappings.sorted.bam - bam output from megalodon basecalling
- beds/intersection.motifs.chip.formatted.chm13.bed - bed file specifying CTCF motifs within CTCF ChIP-seq peaks

### Fig 4d., Extended Data Fig 7
#### **Scripts:**
- HaplotypePhasing.txt - instructions to phase reads to create HP1 and HP2 bams
- single_molecule_browser.py - to create browser views (Fig 4d, Extended Data 7a,b,c)
- single_molecule_roi_meg_winnow_guppy_ctcf_center.py - overlay graph from this function as in Fig 4b but filter bed file (beds/intersection.motifs.chip.formatted.chm13.bed) to CTCF sites on just the X chromosome

#### **Required files:**
- bam output from guppy & winnowmap merge --> see instructions here: hybrid_guppy_winnnowmap_bam_creation
- beds/intersection.motifs.chip.formatted.chm13.bed - bed file specifying CTCF motifs within CTCF ChIP-seq peaks


### Extended Data Figures 5b,c
#### **Scripts:**
- ParseModMappings.hash&#46;mA&#46;pl - extract a list of all mA mod probabilities from the mod_mappings.bam file output by megalodon --> outputs bed file with mA calls by position
- ctcf_atac.sh
- CTCF_ATAC_resolution_H3K9me3enrich.ipynb - create barplots

#### **Required files:**
- mod_mappings.sorted.bam - bam output from megalodon basecalling
- CTCF ChIP-seq peaks, CTCF motifs, ATAC-seq peaks bed files (see Data Availability)




### Extended Data Figure 5d
#### **Scripts:**
- CTCF_H3K9me3.ipynb
#### **Required files:**
- top_decile_mA.csv from single_molecule_roi_meg_winnow_guppy_ctcf_decile.py --> csv output when considering top decile with mod_mappings.sorted.bam input


### Extended Data Figure 5f 
#### **Scripts:**
- CTCF_H3K9me3.ipynb
#### **Required files:**
- peaks_seq_saturation.csv

### Extended Data Figure 5g
#### **Scripts:**
- sm_peak_caller_ctcf.py

#### **Required files:**
- top_decile_mA.csv from single_molecule_roi_meg_winnow_guppy_ctcf_decile.py --> csv output when considering top decile with mod_mappings.sorted.bam input

### Extended Data Figure 5h
#### **Scripts:**
- single_molecule_roi_meg_winnow_guppy_ctcf_decile.ipynb - with window of 2kb specified to extract read spanning top decile peaks with 2kb window
- CTCF_H3K9me3.ipynb - to generate plots
#### **Required files:**
- mod_mappings.sorted.bam - bam output from megalodon basecalling

### Extended Data Figure 5i 
#### **Scripts:**
- single_molecule_roi_meg_winnow_guppy_ctcf_decile.py - analysis by decile
- ctcf_sensitivity_decile_6mA.py

#### **Required files:**
- mod_mappings.sorted.bam - bam output from megalodon basecalling
- beds/intersection.motifs.chip.formatted.chm13.bed - bed file specifying CTCF motifs within CTCF ChIP-seq peaks
- meg_CTCF_CTCF_qX_all_data_mA.csv, meg_IgG_CTCF_qX_all_data_mA.csv - use these files for q1 to q10 output from single_molecule_roi_meg_winnow_guppy_ctcf_decile.py for ctcf_sensitivity_decile_6mA.ipynb


### Extended Data Figure 8
#### **Scripts:**
- pacbio_ctcf_ipr.sh - processing steps with smrtlink commands and calling PerMoleculeIPDRatio.py described below
- PerMoleculeIPDRatio.py - provided from PacBio to compute IPD ratio on a per molecule basis; stored in bam output
- pacbio_ctcf_sm_decile.py - generate ctcf_pb_top_decile.csv, untreated_pb_top_decile.csv
- CTCF_H3K9me3.ipynb - take in CSVs and generate plots
- single_molecule_roi_meg_winnow_guppy_ctcf_decile.py - for nanopore analysis to generate CSVs for top decile

#### **Required files:**
Files generated from top decile of CTCF ChIP-seq peaks as from single_molecule_roi_meg_winnow_guppy_ctcf_decile.py
- ctcf_nano_top_decile.csv
- untreated_nano_top_decile.csv
- ctcf_pb_top_decile.csv
- untreated_pb_top_decile.csv


## H3K9me3 - Figure 5 and Extended Data Figure 9
### Preprocessing for bam files that are a hybrid from guppy and winnowmap. This guppy-winnowmap hybrid bam is used for H3K9me3 analysis.
See instructions here: hybrid_guppy_winnnowmap_bam_creation

### Figure 5a, 5b

#### **Scripts:**
- fig5ab_153.sh - create bed files for computing enrichment
- ParseModMappings.hash&#46;mA&#46;pl - extract a list of all mA mod probabilities from the hybrid guppy & winnowmap bam
- CTCF_ATAC_resolution_H3K9me3enrich.ipynb - create barplots

#### **Required files:**
- k9me3_enrich.csv - calculated from bed files produced from run_mA_A_peaks_calcs_h3k9me3.sh
- beds/HG2_K9Me3_ADGCCR09-H_S103_L002_2.chm13_HG002X_HG38Y_peaks.broadPeak - CUT&RUN peaks
- beds/t2t_cenAnnotation.v2.021921.hg002_X.4.bed - centromeres
- beds/t2t_cenAnnotation.v2.021921.liveHORonly.HG002.4.sorted.bed - live HOR arrays
- franken.chromsizes.txt - chromosome sizes for franken reference

### Figure 5c, Supplementary Figure 9a
#### **Scripts:**
- single_molecule_roi_meg_winnow_guppy_hor_k9me3.py

#### **Required files:**
- beds/sharp_HOR_boundaries.300kb.bed

### Figure 5d, Supplementary Figure 9b & 9c
#### **Scripts:**
- single_molecule_browser.py --> with coordinates chr7:55414368-68714496
- For showing just a single modification (mA or mCG in Supplementary Figure 6b, comment out traces.append() --> see "comment out for single mod plot")

#### **Required files:**
- guppy_winnow_hybrid.bam - created from instructions above

#### **Other tools**
For aggregate curves:
- Plotted usings the WashU epigenome browser - https://epigenomegateway.wustl.edu
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

