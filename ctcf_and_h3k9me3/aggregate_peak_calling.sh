# aggreate peak calling for Extended Data Fig 5f
# unique dimelo peaks detected vs. depth (measurement is independent of CTCF) - 5X, 10X, 15X, 20X, 25X 

# subsample bam file
module load samtools 
BAM=deep_ctcf_mod_mappings_merge.sorted.bam
for FRAC in 0.2 0.4 0.6 0.8;
do 
samtools view -s $FRAC -b $BAM > $BAM'.'$FRAC'.bam' &
done

# create bed file for each bam file
SCRIPTS='/path/to/scripts'
for FRAC in 0.2 0.4 0.6 0.8;
do
samtools view -h $BAM'.'$FRAC'.bam' | perl $SCRIPTS/ParseModMappings.hash.mA.pl $BAM'.'$FRAC'.mA.bed'
done

# call peaks with best peak calling method
# do average for now --> (old: slurm-10339335.out) --> new: slurm-10443392.out (re-run for 0.6: slurm-10446668.out)

module load bedtools
BAM=deep_ctcf_mod_mappings_merge.sorted.bam
windowsS200='windows_s200.sorted.bed' # bed file with 200 bp sliding windows in 20 bp increments genome-wide
windows200='windows_200.sorted.bed' # bed file with 200 bp windows genome-wide 
CTCF_ChIP='ENCFF797SDL.chm13.bed'

bedtools intersect -u -a $windows200 -b $CTCF_ChIP > windows_200.pos.bed
bedtools intersect -v -a $windows200 -b $CTCF_ChIP > windows_200.neg.bed

for FRAC in 0.2 0.4 0.6 0.8;
do
awk '$1=$1' FS="," OFS="\t" $BAM'.'$FRAC'.mA.bed' | awk '{sum=0;num=0;for(i=4;i<=NF;i++){sum+=$i;num++};if(num>0) ret=sum/num; else ret=0; print $1, $2, $3, ret}' | awk -v OFS="\t" '{print $1, $2, $3, $4}' > $BAM'.'$FRAC'.mA.mean.bed'
bedtools map -a $windowsS200 -b $BAM'.'$FRAC'.mA.mean.bed' -c 4 -o mean > $BAM'.'$FRAC'.mA.mean.200.bed'

for CUTOFF in 0 2.5 5 7.5 10 12.5 15 17.5 20 22.5 25 30 40 50;
do
DiMeLo=$BAM'.'$FRAC'.mA.mean.200_'$CUTOFF'.bed'
# drop entries with no overlap in windows (.)
awk '$1=$1' OFS="\t" $BAM'.'$FRAC'.mA.mean.200.bed' | awk -v var="$CUTOFF" '{if($4>=var && $4 != "."){print $1, $2, $3, $4}}' | awk -v OFS="\t" '{print $1, $2, $3, $4}' > $DiMeLo
# TP = in both DiMeLo and ChIP-seq
bedtools intersect -u -sorted -a windows_200.pos.bed -b $DiMeLo > $DiMeLo'.tp.bed'
# FN = in ChIP-seq only
bedtools intersect -v -sorted -a windows_200.pos.bed  -b $DiMeLo > $DiMeLo'.fn.bed'
# FP = in DiMeLo only
bedtools intersect -u -sorted -a windows_200.neg.bed -b $DiMeLo > $DiMeLo'.fp.bed'
# TN
bedtools intersect -v -sorted -a windows_200.neg.bed -b $DiMeLo > $DiMeLo'.tn.bed'
echo 'TP' $CUTOFF
wc -l $DiMeLo'.tp.bed'
echo 'FN' $CUTOFF
wc -l $DiMeLo'.fn.bed'
echo 'FP' $CUTOFF
wc -l $DiMeLo'.fp.bed'
echo 'TN' $CUTOFF
wc -l $DiMeLo'.tn.bed'
done
done

# also do for non-downsampled with same range
FULL='deep_ctcf.mA.bed'
awk '$1=$1' FS="," OFS="\t" $FULL | awk '{sum=0;num=0;for(i=4;i<=NF;i++){sum+=$i;num++};if(num>0) ret=sum/num; else ret=0; print $1, $2, $3, ret}' | awk -v OFS="\t" '{print $1, $2, $3, $4}' > $FULL'.mA.mean.bed'
bedtools map -a $windowsS200 -b $FULL'.mA.mean.bed' -c 4 -o mean > $FULL'.mA.mean.200.bed'

for CUTOFF in 0 2.5 5 7.5 10 12.5 15 17.5 20 22.5 25 30 40 50;
do
DiMeLo=$FULL'.mA.mean.200.bed_'$CUTOFF'.bed'
# drop entries with no overlap in windows (.)
awk '$1=$1' OFS="\t" $FULL'.mA.mean.200.bed' | awk -v var="$CUTOFF" '{if($4>=var && $4 != "."){print $1, $2, $3, $4}}' | awk -v OFS="\t" '{print $1, $2, $3, $4}' > $DiMeLo
# TP = in both DiMeLo and ChIP-seq
bedtools intersect -u -sorted -a windows_200.pos.bed -b $DiMeLo > $DiMeLo'.tp.bed'
# FN = in ChIP-seq only
bedtools intersect -v -sorted -a windows_200.pos.bed  -b $DiMeLo > $DiMeLo'.fn.bed'
# FP = in DiMeLo only
bedtools intersect -u -sorted -a windows_200.neg.bed -b $DiMeLo > $DiMeLo'.fp.bed'
# TN
bedtools intersect -v -sorted -a windows_200.neg.bed -b $DiMeLo > $DiMeLo'.tn.bed'
echo 'TP' $CUTOFF
wc -l $DiMeLo'.tp.bed'
echo 'FN' $CUTOFF
wc -l $DiMeLo'.fn.bed'
echo 'FP' $CUTOFF
wc -l $DiMeLo'.fp.bed'
echo 'TN' $CUTOFF
wc -l $DiMeLo'.tn.bed'
done