# Relevant scripts and files for CTCF and H3K9me3 analysis for creating Figures 4 & 5 and Extended Data Figures 5, 6, 7, 8, 9

Data used to generate figures are a merge of all bam files for a given target (CTCF, H3K9me3, IgG, Hia5) in the given cell type (GM12878, HG002).

## CTCF - Figure 4 and Extended Data Figures 5,6,7,8 

### Fig 4a, 4b; Extended Data Fig 5a, 5e; Extended Data Fig 6
#### **Script:**
- single_molecule_roi_meg_winnow_guppy_ctcf_center.py 
#### **Required files:**
- mod_mappings.sorted.bam - bam output from megalodon basecalling
- beds/intersection.motifs.chip.formatted.chm13.bed - bed file specifying CTCF motifs within CTCF ChIP-seq peaks (Fig 4a & Extended Data Fig 6: quartiles 1,2,3,4; Fig 4b: top decile only; Extended Data Fig 5a,e: quartile 4)


### Fig 4c.
#### **Script:**
- single_molecule_joint_ctcf.py

#### **Required files:**
- mod_mappings.sorted.bam - bam output from megalodon basecalling
- beds/intersection.motifs.chip.formatted.chm13.bed - bed file specifying CTCF motifs within CTCF ChIP-seq peaks

### Fig 4d., Extended Data Fig 7
#### **Scripts:**
- HaplotypePhasing.txt - instructions to phase reads to create HP1 and HP2 bams using nanomethphase
- single_molecule_browser.py - to create browser views (Fig 4d, Extended Data 7a,b,c)
- single_molecule_roi_meg_winnow_guppy_ctcf_center.py - overlay graph from this function as in Fig 4b but filter bed file (beds/intersection.motifs.chip.formatted.chm13.bed) to CTCF sites on just the X chromosome

#### **Required files:**
- bam output from guppy & winnowmap merge --> see instructions here: hybrid_guppy_winnnowmap_bam_creation
- beds/intersection.motifs.chip.formatted.chm13.bed - bed file specifying CTCF motifs within CTCF ChIP-seq peaks


### Extended Data Fig 5b,c
#### **Scripts:**
- ParseModMappings.hash&#46;mA&#46;pl - extract a list of all mA mod probabilities from the mod_mappings.bam file output by megalodon --> outputs bed file with mA calls by position
- ctcf_atac.sh --> calculate mA and A counts in regions of interest defined in bed files
- CTCF_ATAC_resolution_H3K9me3enrich.ipynb - create barplots

#### **Required files:**
- mod_mappings.sorted.bam - bam output from megalodon basecalling
- CTCF ChIP-seq peaks, CTCF motifs, ATAC-seq peaks bed files (see Data Availability)


### Extended Data Fig 5d
#### **Scripts:**
- CTCF_H3K9me3.ipynb --> fit exponential decay to mA/A vs. distance from motif center and generate plot
#### **Required files:**
- top_decile_mA.csv from single_molecule_roi_meg_winnow_guppy_ctcf_decile.py --> csv output when considering top decile with mod_mappings.sorted.bam input


### Extended Data Fig 5f 
#### **Scripts:**
- ParseModMappings.hash&#46;mA&#46;pl - extract a list of all mA mod probabilities from the mod_mappings.bam file output by megalodon --> outputs bed file with mA calls by position
- aggregate_peak_calling.sh - calculate TP, FN, FP, TN mA calls at various thresholds and sequencing depths for aggregate peak calling
- CTCF_H3K9me3.ipynb - generate saturation and AUC curves as a function of sequencing depth
#### **Required files:**
- mod_mappings.sorted.bam - bam output from megalodon basecalling
- ENCFF797SDL.chm13.bed - CTCF ChIP-seq bed file for ground truth
- peaks_seq_saturation.csv - manually create from output of aggregate_peak_calling.sh

### Extended Data Fig 5g
#### **Scripts:**
- sm_peak_caller_ctcf.py - perform single-molecule peak calling and determine distance of predicted peak center to generate Fig 5g

#### **Required files:**
- top_decile_mA.csv from single_molecule_roi_meg_winnow_guppy_ctcf_decile.py --> csv output when considering top decile with mod_mappings.sorted.bam input

### Extended Data Fig 5h
#### **Scripts:**
- single_molecule_roi_meg_winnow_guppy_ctcf_decile.py - with window of 2kb specified to extract reads spanning top decile peaks within 2kb window
- CTCF_H3K9me3.ipynb - perform TPR and FPR calculations and generate Extended Data Fig 5h plot
#### **Required files:**
- mod_mappings.sorted.bam - bam output from megalodon basecalling

### Extended Data Fig 5i 
#### **Scripts:**
- single_molecule_roi_meg_winnow_guppy_ctcf_decile.py - analysis by decile
- ctcf_sensitivity_decile_6mA.py - generate Extended Data Fig 5i

#### **Required files:**
- mod_mappings.sorted.bam - bam output from megalodon basecalling
- beds/intersection.motifs.chip.formatted.chm13.bed - bed file specifying CTCF motifs within CTCF ChIP-seq peaks
- meg_CTCF_CTCF_qX_all_data_mA.csv, meg_IgG_CTCF_qX_all_data_mA.csv - use these files for q1 to q10 output from single_molecule_roi_meg_winnow_guppy_ctcf_decile.py for ctcf_sensitivity_decile_6mA.py


### Extended Data Fig 8
#### **Scripts:**
- pacbio_ctcf_ipr.sh - processing steps with smrtlink commands and calling PerMoleculeIPDRatio.py described below
- PerMoleculeIPDRatio.py - provided from PacBio to compute IPD ratio on a per molecule basis; stored in bam output
- pacbio_ctcf_sm_decile.py - generate ctcf_pb_top_decile.csv, untreated_pb_top_decile.csv with IPD ratio and number of passes extracted from PacBio bam files
- single_molecule_roi_meg_winnow_guppy_ctcf_decile.py - for nanopore analysis to generate CSVs for top decile
- CTCF_H3K9me3.ipynb - take in CSVs from pacbio_ctcf_sm_decile.py and generate plots


#### **Required files:**
- mod_mappings.sorted.bam - bam output from megalodon basecalling
- bam with HiFi kinetics information output for PacBio


## H3K9me3 - Figure 5 and Extended Data Figure 9
### Preprocessing for bam files that are a hybrid from guppy and winnowmap. This guppy-winnowmap hybrid bam is used for H3K9me3 analysis.
See instructions here: hybrid_guppy_winnnowmap_bam_creation

### Figure 5a, 5b

#### **Scripts:**
- fig5ab_153.sh - create bed files for computing enrichment
- ParseModMappings.hash&#46;mA&#46;pl - extract a list of all mA mod probabilities from the hybrid guppy & winnowmap bam
- CTCF_ATAC_resolution_H3K9me3enrich.ipynb - create barplots

#### **Required files:**
- guppy_winnow_hybrid.bam - created from instructions above
- k9me3_enrich.csv - calculated from bed files produced from fig5ab_153.sh that using guppy_winnow_hybrid.bam as input
- beds/HG2_K9Me3_ADGCCR09-H_S103_L002_2.chm13_HG002X_HG38Y_peaks.broadPeak - CUT&RUN peaks
- beds/t2t_cenAnnotation.v2.021921.hg002_X.4.bed - centromeres
- beds/t2t_cenAnnotation.v2.021921.liveHORonly.HG002.4.sorted.bed - live HOR arrays
- franken.chromsizes.txt - chromosome sizes for franken reference

### Figure 5c, Extended Data Fig 9a
#### **Scripts:**
- single_molecule_roi_meg_winnow_guppy_hor_k9me3.py - create boundary plots

#### **Required files:**
- guppy_winnow_hybrid.bam - created from instructions above
- beds/sharp_HOR_boundaries.300kb.bed

### Figure 5d, Extended Data Fig 9b & 9c
#### **Scripts:**
- single_molecule_browser.py --> with coordinates chr7:55414368-68714496
- For showing just a single modification (mA or mCG in Extended Data Figure 9b, comment out traces.append() --> see "comment out for single mod plot" in single_molecule_browser.py)

#### **Required files:**
- guppy_winnow_hybrid.bam - created from instructions above

#### **Other tools**
For aggregate curves:
- Plotted usings the WashU epigenome browser - https://epigenomegateway.wustl.edu
- Bedgraph files for plot created by outputting all_data_mA and all_data_mC dataframes from single_molecule_browser.py as csv's and then running the below for both data_A and data_C.
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

