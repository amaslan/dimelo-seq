# script for doing calculations for Figure 5a. and 5b.
module load bedtools/2.28.0
module load samtools/1.8

BAMS='/bams/guppywinnow'
BEDS='/beds'
HK9me3_CUTRUN='HG2_K9Me3_ADGCCR09-H_S103_L002_2.chm13_HG002X_HG38Y_peaks.broadPeak'
SCRIPTS='/scripts'
CEN='t2t_cenAnnotation.v2.021921.hg002_X.4.bed'
REF='/ref'
HOR='t2t_cenAnnotation.v2.021921.liveHORonly.HG002.4.sorted.bed'

for NAME in H3K9me3 Hia5_HG002 IgG_HG002;

do 

samtools view -h $BAMS'/prod_'$NAME'_winnowmap_guppy_merge.sorted.q10.bam' | perl $SCRIPTS/ParseModMappings.hash.mA.pl $BEDS'/prod_'$NAME'.mA.bed'

# bed file with mA probability for every A sequenced - take max at each position
awk '$1=$1' FS="," OFS="\t" $BEDS'/prod_'$NAME'.mA.bed' | awk '{m=$4;for(i=4;i<=NF;i++)if($i>m)m=$i;print $1, $2, $3, m}' | awk -v OFS="\t" '{print $1, $2, $3, $4/255}' > $BEDS'/prod_'$NAME'.mA.max.bed'

# A positions where mA is above threshold
awk -v OFS="\t" '{if ($4 >= 0.9) {print $1, $2, $3, $4}}' $BEDS'/prod_'$NAME'.mA.max.bed' > $BEDS'/'$NAME'_mA_0.9.bed'

# A positions where mA is intersecting cut&run peaks
bedtools intersect -a $HK9me3_CUTRUN -b $BEDS'/'$NAME'_mA_0.9.bed' > $BEDS'/'$NAME'_mA_0.9.overlap.bed'

# need to first sort
sort -k1,1 -k2,2n $BEDS'/prod_'$NAME'.mA.max.bed' > $BEDS'/prod_'$NAME'.mA.max.sorted.bed'

# A positions intersecting cut&run peaks - large to need to sort before intersect to not exceed ram limits
bedtools intersect -sorted -a $HK9me3_CUTRUN'.sorted' -b $BEDS'/prod_'$NAME'.mA.max.sorted.bed' > $BEDS'/'$NAME'_mA.overlap.bed'

echo 'total A genome wide for '$NAME
wc -l $BEDS'/prod_'$NAME'.mA.max.bed'

echo 'total A in ChIP-seq peaks for '$NAME
wc -l $BEDS'/'$NAME'_mA.overlap.bed'

echo 'all mA calls above threshold for '$NAME
wc -l $BEDS'/'$NAME'_mA_0.9.bed'

echo 'mA calls above thresh that fall in CTCF ChIP-seq peaks for '$NAME
wc -l $BEDS'/'$NAME'_mA_0.9.overlap.bed'

# repeat analysis with centromere bed file overlap rather than cut&run to show enrichment in centromeres relative to controls

# A positions where mA is intersecting centromeres
bedtools intersect -a $CEN -b $BEDS'/'$NAME'_mA_0.9.bed' > $BEDS'/'$NAME'_mA_0.9.overlap.cen.bed'

# A positions intersecting centromeres
bedtools intersect -sorted -a $CEN'.sorted' -b $BEDS'/prod_'$NAME'.mA.max.sorted.bed' > $BEDS'/'$NAME'_mA.overlap.cen.bed' 

echo 'total A in centromeres for '$NAME
wc -l $BEDS'/'$NAME'_mA.overlap.cen.bed'

echo 'mA calls above thresh that fall in centromeres for '$NAME
wc -l $BEDS'/'$NAME'_mA_0.9.overlap.cen.bed'


# change denominator to instead of being genome wide minus peaks be genome wide minus peaks + 10kb flanking each side of peak

bedtools slop -i $HK9me3_CUTRUN'.sorted' -g $REF'/franken.chromsizes.txt' -b 10000 > $HK9me3_CUTRUN'.sorted.10kbslop'
bedtools merge -i $HK9me3_CUTRUN'.sorted.10kbslop' > $HK9me3_CUTRUN'.sorted.10kbslop.merged'
bedtools intersect -sorted -a $HK9me3_CUTRUN'.sorted.10kbslop.merged' -b $BEDS'/prod_'$NAME'.mA.max.sorted.bed' > $BEDS'/'$NAME'_mA.overlap.10kbslop.bed'
bedtools intersect -a $HK9me3_CUTRUN'.sorted.10kbslop.merged' -b $BEDS'/'$NAME'_mA_0.9.bed' > $BEDS'/'$NAME'_mA_0.9.overlap.10kbslop.bed' # should not have -sorted

echo 'total A in peaks +/- 10kb for '$NAME
wc -l $BEDS'/'$NAME'_mA.overlap.10kbslop.bed'

echo 'mA calls in peaks +/- 10kb for '$NAME
wc -l $BEDS'/'$NAME'_mA_0.9.overlap.10kbslop.bed'


# intersect with HORs
# A positions where mA is intersecting centromeres
bedtools intersect -a $HOR -b $BEDS'/'$NAME'_mA_0.9.bed' > $BEDS'/'$NAME'_mA_0.9.overlap.hor.bed'

# A positions intersecting centromeres
bedtools intersect -sorted -a $HOR -b $BEDS'/prod_'$NAME'.mA.max.sorted.bed' > $BEDS'/'$NAME'_mA.overlap.hor.bed'

echo 'total A in HOR for '$NAME
wc -l $BEDS'/'$NAME'_mA.overlap.hor.bed'

echo 'mA calls above thresh that fall in HOR for '$NAME
wc -l $BEDS'/'$NAME'_mA_0.9.overlap.hor.bed'


done


