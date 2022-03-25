


BEDS='/beds'
HK9me3_CUTRUN='HG2_K9Me3_ADGCCR09-H_S103_L002_2.chm13_HG002X_HG38Y_peaks.broadPeak'
CEN='t2t_cenAnnotation.v2.021921.hg002_X.4.bed'
HOR='t2t_cenAnnotation.v2.021921.liveHORonly.HG002.4.sorted.bed'
chromsizes='franken.chromsizes.txt'
OUT='/out'

bedtools slop -i $HK9me3_CUTRUN'.sorted' -g $chromsizes -b 10000 > $HK9me3_CUTRUN'.sorted.10kbslop'
bedtools merge -i $HK9me3_CUTRUN'.sorted.10kbslop' > $HK9me3_CUTRUN'.sorted.10kbslop.merged'

for NAME in H3K9me3 Hia5_HG002 IgG_HG002;

do 

# bed file with mA probability for every A sequenced 
awk '$1=$1' FS="," OFS="\t" $BEDS'/prod_'$NAME'.mA.bed' | awk '{meth=0;for(i=4;i<=NF;i++)if($i>=153)meth=meth+1;print $1, $2, $3, NF-3, meth}' | awk -v OFS="\t" '{print $1, $2, $3, $4, $5}' > $OUT'/prod_'$NAME'.mA.t153.all.bed'
sort -k1,1 -k2,2n $OUT'/prod_'$NAME'.mA.t153.all.bed' > $OUT'/prod_'$NAME'.mA.t153.all.sorted.bed'

# intersect h3k9me3 cut&run
bedtools intersect -sorted -u -a $OUT'/prod_'$NAME'.mA.t153.all.sorted.bed' -b $HK9me3_CUTRUN'.sorted'  > $OUT'/'$NAME'.mA.t153.all.overlap.cutrun.bed'

# intersect centromeres
bedtools intersect -sorted -u -a $OUT'/prod_'$NAME'.mA.t153.all.sorted.bed' -b $CEN'.sorted' > $OUT'/'$NAME'.mA.t153.all.overlap.cen.bed'

# denominator is genome wide minus peaks + 10kb flanking each side of peak
bedtools intersect -sorted -v -a $OUT'/prod_'$NAME'.mA.t153.all.sorted.bed' -b $HK9me3_CUTRUN'.sorted.10kbslop.merged'  > $OUT'/'$NAME'.mA.t153.all.overlap.nocutrun.bed'

# intersect HORs
bedtools intersect -u -a $OUT'/prod_'$NAME'.mA.t153.all.sorted.bed' -b $HOR > $OUT'/'$NAME'.mA.t153.all.overlap.hor.bed'

# outside centromeres
bedtools intersect -sorted -v -a $OUT'/prod_'$NAME'.mA.t153.all.sorted.bed' -b $CEN'.sorted' > $OUT'/'$NAME'.mA.t153.all.overlap.nocen.bed'



awk -F'\t' '{sum+=$4;}END{print sum;}' $OUT'/'$NAME'.mA.t153.all.overlap.cutrun.bed'   
awk -F'\t' '{sum+=$5;}END{print sum;}' $OUT'/'$NAME'.mA.t153.all.overlap.cutrun.bed'  

awk -F'\t' '{sum+=$4;}END{print sum;}' $OUT'/'$NAME'.mA.t153.all.overlap.cen.bed'   
awk -F'\t' '{sum+=$5;}END{print sum;}' $OUT'/'$NAME'.mA.t153.all.overlap.cen.bed'  

awk -F'\t' '{sum+=$4;}END{print sum;}' $OUT'/'$NAME'.mA.t153.all.overlap.nocutrun.bed'   
awk -F'\t' '{sum+=$5;}END{print sum;}' $OUT'/'$NAME'.mA.t153.all.overlap.nocutrun.bed'  

awk -F'\t' '{sum+=$4;}END{print sum;}' $OUT'/'$NAME'.mA.t153.all.overlap.hor.bed'   
awk -F'\t' '{sum+=$5;}END{print sum;}' $OUT'/'$NAME'.mA.t153.all.overlap.hor.bed'  

awk -F'\t' '{sum+=$4;}END{print sum;}' $OUT'/'$NAME'.mA.t153.all.overlap.nocen.bed'   
awk -F'\t' '{sum+=$5;}END{print sum;}' $OUT'/'$NAME'.mA.t153.all.overlap.nocen.bed' 

done


