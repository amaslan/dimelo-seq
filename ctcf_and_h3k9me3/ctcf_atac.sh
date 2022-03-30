# calculations for overall mA/A in regions of interest to produce Extended Data Fig 5b,c

BEDS='/path/to/bed/files'
OUT='/path/to/out/dir'

CTCF_ChIP='ENCFF797SDL.chm13.bed'
# with bedtools slop +/-5kb around CTCF ChIP-seq peaks
CTCF_ChIP_5kb_buffer='ENCFF797SDL.chm13.bed.5kb'
# all CTCF motifs
CTCF_MOTIF='ctcf.motifs.known.chm13.bed'

ATAC='ENCFF748UZH.chm13.bed'
ATAC_noCTCF_buffer='ATAC_noCTCF_buffer.chm13.bed' 
ATAC_noCTCF_buffer_noMotif='ATAC_noCTCF_buffer_noMotif.chm13.bed'
ATAC_yesCTCF='ATAC_yesCTCF.chm13.bed'

bedtools intersect -v -a $ATAC -b $CTCF_ChIP_5kb_buffer > $ATAC_noCTCF_buffer
bedtools intersect -v -a $ATAC_noCTCF_buffer -b $CTCF_MOTIF > $ATAC_noCTCF_buffer_noMotif
bedtools intersect -u -a $ATAC -b $CTCF_ChIP > $ATAC_yesCTCF

# deep_ctcf is merge of all CTCF bam files
for NAME in deep_ctcf prod_Hia5 prod_IgG;
do

# calculate #A and #mA sequenced at each position and make new bed file rather than summarizing by max
# wrt to what we sequence rather than wrt to reference so we don't need to take max
# megalodon threshold is 190
awk '$1=$1' FS="," OFS="\t" $BEDS'/'$NAME'.mA.bed' | awk '{meth=0;for(i=4;i<=NF;i++)if($i>=190)meth=meth+1;print $1, $2, $3, NF-3, meth}' | awk -v OFS="\t" '{print $1, $2, $3, $4, $5}' > $BEDS'/'$NAME'.mA.t190.all.bed'


# A positions where mA is intersecting ATAC-seq peaks outside of CTCF peaks
bedtools intersect -sorted -a $BEDS'/'$NAME'.mA.t190.all.bed' -b $ATAC_noCTCF_buffer_noMotif'.sorted'  > $OUT'/'$NAME'.mA.t190.all.overlap.openNoCTCF.bed'
bedtools intersect -sorted -a $BEDS'/'$NAME'.mA.t190.all.bed' -b $ATAC_yesCTCF'.sorted'  > $OUT'/'$NAME'.mA.t190.all.overlap.openYesCTCF.bed'

# also intersect ChIP vs. outside as in paper
bedtools intersect -sorted -a $BEDS'/'$NAME'.mA.t190.all.bed' -b $CTCF_ChIP  > $OUT'/'$NAME'.mA.t190.all.overlap.CTCF.bed'
bedtools intersect -sorted -v -a $BEDS'/'$NAME'.mA.t190.all.bed' -b $CTCF_ChIP_5kb_buffer  > $OUT'/'$NAME'.mA.t190.all.overlap.noCTCF.bed'



# sum field 4 for count of A; sum field 5 for mA count
awk -F'\t' '{sum+=$4;}END{print sum;}' $OUT'/'$NAME'.mA.t190.all.overlap.openNoCTCF.bed'   
awk -F'\t' '{sum+=$5;}END{print sum;}' $OUT'/'$NAME'.mA.t190.all.overlap.openNoCTCF.bed'   

awk -F'\t' '{sum+=$4;}END{print sum;}' $OUT'/'$NAME'.mA.t190.all.overlap.openYesCTCF.bed'   
awk -F'\t' '{sum+=$5;}END{print sum;}' $OUT'/'$NAME'.mA.t190.all.overlap.openYesCTCF.bed'   

awk -F'\t' '{sum+=$4;}END{print sum;}' $OUT'/'$NAME'.mA.t190.all.overlap.CTCF.bed'   
awk -F'\t' '{sum+=$5;}END{print sum;}' $OUT'/'$NAME'.mA.t190.all.overlap.CTCF.bed' 

awk -F'\t' '{sum+=$4;}END{print sum;}' $OUT'/'$NAME'.mA.t190.all.overlap.noCTCF.bed'   
awk -F'\t' '{sum+=$5;}END{print sum;}' $OUT'/'$NAME'.mA.t190.all.overlap.noCTCF.bed' 
done