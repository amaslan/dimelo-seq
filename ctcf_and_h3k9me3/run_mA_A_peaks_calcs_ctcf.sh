# script for doing calculations for Supplementary Figure 4a. 

module load bedtools/2.28.0
module load samtools/1.8

BAMS='/bams/megalodon'
BEDS='/beds'
CTCF_ChIP='ENCFF797SDL.chm13.bed'
SCRIPTS='/scripts'

sort -k1,1 -k2,2n $CTCF_ChIP > $CTCF_ChIP'.sorted'

# already have Hia5 file done:

samtools view -h $BAMS/prod_ctcf_mod_mappings_merge.sorted.bam | perl $SCRIPTS/ParseModMappings.hash.mA.pl $BEDS/prod_CTCF.mA.bed
samtools view -h $BAMS/prod_free_Hia5_mod_mappings.sorted.bam | perl $SCRIPTS/ParseModMappings.hash.mA.pl $BEDS/prod_Hia5.mA.bed
samtools view -h $BAMS/prod_IgG_mod_mappings.sorted.bam | perl $SCRIPTS/ParseModMappings.hash.mA.pl $BEDS/prod_IgG.mA.bed

for NAME in CTCF Hia5 IgG;

do 

# bed file with mA probability for every A sequenced
awk '$1=$1' FS="," OFS="\t" $BEDS'/prod_'$NAME'.mA.bed' | awk '{m=$4;for(i=4;i<=NF;i++)if($i>m)m=$i;print $1, $2, $3, m}' | awk -v OFS="\t" '{print $1, $2, $3, $4/255}' > $BEDS'/prod_'$NAME'.mA.max.bed'

# A positions where mA is above threshold
awk -v OFS="\t" '{if ($4 >= 0.9) {print $1, $2, $3, $4}}' $BEDS'/prod_'$NAME'.mA.max.bed' > $BEDS'/'$NAME'_mA_0.9.bed'

# A positions where mA is intersecting cut&run peaks
bedtools intersect -a $CTCF_ChIP -b $BEDS'/'$NAME'_mA_0.9.bed' > $BEDS'/'$NAME'_mA_0.9.overlap.bed'

# need to first sort
sort -k1,1 -k2,2n $BEDS'/prod_'$NAME'.mA.max.bed' > $BEDS'/prod_'$NAME'.mA.max.sorted.bed'

# A positions intersecting cut&run peaks - large to need to sort before intersect to not exceed ram limits
bedtools intersect -sorted -a $CTCF_ChIP'.sorted' -b $BEDS'/prod_'$NAME'.mA.max.sorted.bed' > $BEDS'/'$NAME'_mA.overlap.bed'

echo 'total A genome wide for '$NAME
wc -l $BEDS'/prod_'$NAME'.mA.max.bed'

echo 'total A in ChIP-seq peaks for '$NAME
wc -l $BEDS'/'$NAME'_mA.overlap.bed'

echo 'all mA calls above threshold for '$NAME
wc -l $BEDS'/'$NAME'_mA_0.9.bed'

echo 'mA calls above thresh that fall in CTCF ChIP-seq peaks for '$NAME
wc -l $BEDS'/'$NAME'_mA_0.9.overlap.bed'

done



