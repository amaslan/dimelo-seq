# create alignmed bam files with IPD ratio tag included for PacBio sequencing data


PB4_RAW=PB4.hifi.bam # CTCF-targeted sample - bam from PacBio sequencer with kinetics info
PB8_RAW=PB8.hifi.bam # untreated sample - bam from PacBio sequencer with kinetics info


BED=intersection.motifs.chip.formatted.chm13.q10.bed # top decile of peaks

OUT=/path/to/out/dir
REF_FASTA=chm13.draft_v1.0.fasta
REF_MMI=chm13.draft_v1.0.fasta.mmi

SMRT_BIN=/smrtlink/smrtcmds/bin

$SMRT_BIN/ccs-kinetics-bystrandify $PB4_RAW $OUT/PB4.strandify.bam
$SMRT_BIN/ccs-kinetics-bystrandify $PB8_RAW $OUT/PB8.strandify.bam

PB4=PB4.strandify.bam # CTCF-targeted sample
PB8=PB8.strandify.bam # untreated sample

# extract subset of reads that overlap top decile CTCF sites

$SMRT_BIN/pbmm2 align --preset CCS --sort $REF_MMI $PB4 $OUT/PB4.strandify.aligned.bam --log-file $OUT/log_align_PB4
$SMRT_BIN/pbmm2 align --preset CCS --sort $REF_MMI $PB8 $OUT/PB8.strandify.aligned.bam --log-file $OUT/log_align_PB8

# extract subset of reads overlapping top decile of ChIP-seq peaks
bedtools intersect -u -a $OUT/PB4.strandify.aligned.bam -b $BED > $OUT/PB4.strandify.aligned.q10.bam
bedtools intersect -u -a $OUT/PB8.strandify.aligned.bam -b $BED > $OUT/PB8.strandify.aligned.q10.bam

# run PerMoleculeIPDRatio 
python $OUT/PerMoleculeIPDRatio.py $OUT/PB4.strandify.aligned.q10.bam $OUT/PB4.strandify.ipr.bam
python $OUT/PerMoleculeIPDRatio.py $OUT/PB8.strandify.aligned.q10.bam $OUT/PB8.strandify.ipr.bam

# align output that contains ipr
$SMRT_BIN/pbmm2 align --preset CCS --sort $REF_MMI $OUT/PB4.strandify.ipr.bam $OUT/PB4.strandify.ipr.aligned.bam
$SMRT_BIN/pbmm2 align --preset CCS --sort $REF_MMI $OUT/PB8.strandify.ipr.bam $OUT/PB8.strandify.ipr.aligned.bam

samtools index $OUT/PB4.strandify.ipr.aligned.bam
samtools index $OUT/PB8.strandify.ipr.aligned.bam

