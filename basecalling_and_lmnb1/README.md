#basecalling LMNB1 DiMeLo-seq data on AWS

#log in to AWS g4dn.metal instance running Ubuntu 18
ssh -i "keys.pem" ec2-3-85-231-198.compute-1.amazonaws.com -l ubuntu

#initialize local SSD disks
sudo fdisk /dev/nvme1n1
n
p
w
sudo partprobe /dev/nvme1n1
sudo mkfs.ext4 /dev/nvme1n1

sudo mount /dev/nvme1n1 /Data1
sudo chown -R ubuntu /Data1

sudo fdisk /dev/nvme2n1
n
p
w
sudo partprobe /dev/nvme2n1
sudo mkfs.ext4 /dev/nvme2n1

sudo mount /dev/nvme2n1 /Data2
sudo chown -R ubuntu /Data2


#install/update guppy, samtools, bedtools
sudo apt-get update
sudo apt-get install wget lsb-release
export PLATFORM=$(lsb_release -cs)
wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
sudo apt-get update
sudo apt update
sudo apt install samtools
sudo apt-get install bedtools
sudo apt install ont-guppy --no-install-recommends

#obtain the rerio all-contexts model file
git clone https://github.com/nanoporetech/rerio
rerio/download_model.py rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001
sudo cp /Data2/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001 /opt/ont/guppy/data/
sudo cp /Data2/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.cfg /opt/ont/guppy/data/
sudo cp /Data2/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.jsn /opt/ont/guppy/data/
sudo cp /Data2/software/rerio/basecall_models/barcoding/* /opt/ont/guppy/data/barcoding/
 
#install/update megalodon
pip install megalodon 
pip install ont_pyguppy_client_lib 


#install winnowmap
cd /Data2
mkdir software
cd software
wget https://github.com/marbl/Winnowmap/archive/refs/tags/v2.03.tar.gz
tar zxvf v2.03.tar.gz
cd Winnowmap-2.03
make -j8
cd ..
mv Winnowmap-2.03 Winnowmap


#download reference (CHM13v1 can be obtained here: https://www.ncbi.nlm.nih.gov/assembly/GCA_009914755.2, HG002 ChrX v0.7 can be obtained here: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/HG002.chrX_v0.7.fasta.gz)
mkdir /Data2/reference
cd /Data2/reference
mkdir t2t
cd t2t
#download reference fasta here
#also place ciLAD_gold_final_100kb.chm13.sorted2.bed and cLAD_gold_final_100kb.chm13.sorted2.bed chm13.window100kb.sorted.bed chm13.chromsizes.sorted.txt here

#optional: if not yet done, index reference(s) for megalodon and winnowmap
/Data2/software/Winnowmap/bin/meryl count k=15 threads=92 output chm13_merylDB chm13.draft_v1.0.fasta
/Data2/software/Winnowmap/bin/meryl print greater-than distinct=0.9998 chm13_merylDB > chm13.repetitive_k15.txt
/Data2/software/Winnowmap/bin/winnowmap -W /Data2/reference/t2t/chm13.repetitive_k15.txt -d chm13.draft_v1.0.idx chm13.draft_v1.0.fasta
curl -L https://github.com/lh3/minimap2/releases/download/v2.20/minimap2-2.20_x64-linux.tar.bz2 | tar -jxvf -
/Data2/software/minimap2-2.20_x64-linux/minimap2 -d chm13.draft_v1.0.mmi chm13.draft_v1.0.fasta


#download raw fast5 data (do this in a separate screen in case of connection failure; ctrl+a then d to exit screen)
mkdir /Data2/fast5
cd /Data2/fast5
screen -L -R one
aws s3 sync s3://bucket-name/fast5-directory .


#run guppy, parallelize using basecall server
mkdir /Data1/log
#with_demux
#initiate server in one screen
guppy_basecall_server --config res_dna_r941_min_modbases-all-context_v001.cfg --log_path /Data1/log --port 5555 --ipc_threads 16 --device cuda:all --num_callers 16 --gpu_runners_per_device 24 --chunks_per_runner 1024 --chunk_size 2000 --barcode_kits "SQK-NBD110-24" --trim_barcodes
#call supervisor in another screen once server is running
guppy_basecaller_supervisor --num_clients 16 --input_path /Data2/fast5 --save_path /Data1 --fast5_out --bam_out --bam_methylation_threshold 0 --config res_dna_r941_min_modbases-all-context_v001.cfg --port 5555 --barcode_kits "SQK-NBD110-24" --trim_barcodes

#without_demux
#initiate server in one screen
guppy_basecall_server --config res_dna_r941_min_modbases-all-context_v001.cfg --log_path /Data1/log --port 5555 --ipc_threads 16 --device cuda:all --num_callers 16 --gpu_runners_per_device 24 --chunks_per_runner 1024 --chunk_size 2000
#call supervisor in another screen once server is running
guppy_basecaller_supervisor --num_clients 16 --input_path /Data2/fast5 --save_path /Data1 --bam_out --bam_methylation_threshold 0 --fast5_out --config res_dna_r941_min_modbases-all-context_v001.cfg --port 5555 

#high accuracy calling (no mod calling)--need to change chunk and runner settings
guppy_basecall_server --config dna_r9.4.1_450bps_sup.cfg --log_path /Data1/log --port 5555 --ipc_threads 16 --device cuda:all --num_callers 16 --gpu_runners_per_device 24 --chunks_per_runner 256 --chunk_size 2000
guppy_basecaller_supervisor --num_clients 16 --input_path /Data2/fast5 --save_path /Data1/superaccuracy --config dna_r9.4.1_450bps_sup.cfg --port 5555 


#check GPU status while guppy is running
nvidia-smi

#merge sequencing summary files
cd /Data1
awk 'FNR==1 && NR!=1{next;}{print}' sequencing_summary_*.txt >sequencing_summary_ALL.txt


#extract basemod calls (run commands in python)
python

from ont_fast5_api.fast5_interface import get_fast5_file
import numpy as np
from ont_fast5_api.analysis_tools.base_tool import BaseTool
import re
import time
import glob
import csv
from joblib import Parallel, delayed
import multiprocessing
import os

num_cores = multiprocessing.cpu_count()
num_cores

summfile = '/Data1/sequencing_summary_ALL.txt'
directory = r'/Data1/workspace/Data2/fast5/*.fast5'
outfile = '/Data1/output_methylAthresh.txt'


barcodeDict = dict()
with open(summfile, newline='') as csvfile:
	reader = csv.DictReader(csvfile,delimiter="\t")
	for row in reader:
		barcodeDict[row['read_id']]=row['barcode_arrangement']


qthreshInt=10
qthresh = chr(qthreshInt+33)
mAthresh1=252 #255*0.99
mAthresh2=250 #255*0.98
mAthresh3=242 #255*0.95
mAthresh4=230 #255*0.9
mAthresh5=204 #255*0.8
mAthresh6=128 #255*0.5

totA=0
totbases=0

def processFile(fast5_filepath):
	outputstring=""
	with get_fast5_file(fast5_filepath, mode="r") as f5:
		for read_id in f5.get_read_ids():
			read = f5.get_read(read_id)
			latest_basecall = read.get_latest_analysis('Basecall_1D')
			mod_base_table = np.array(read.get_analysis_dataset(latest_basecall, 'BaseCalled_template/ModBaseProbs').transpose()[1,:])
			name, sequence, _, qstring = read.get_analysis_dataset(latest_basecall,'BaseCalled_template/Fastq').strip().split('\n')
			Athreshvec = np.logical_and(np.array(tuple(qstring)) >= qthresh, np.array(tuple(sequence)) == 'A')
			mAthreshvec = np.logical_and(Athreshvec, mod_base_table >= mAthresh4) #decide threshold here				
			printstr = np.array2string(np.flatnonzero(mAthreshvec),separator=',',max_line_width=10000000)
			readlen = len(sequence)
			Acount = np.sum(Athreshvec)
			mAcount1 = np.sum(mAthreshvec)
			mAcount2 = np.sum(np.logical_and(Athreshvec, mod_base_table >= mAthresh2))
			mAcount3 = np.sum(np.logical_and(Athreshvec, mod_base_table >= mAthresh3))
			mAcount4 = np.sum(np.logical_and(Athreshvec, mod_base_table >= mAthresh4))
			mAcount5 = np.sum(np.logical_and(Athreshvec, mod_base_table >= mAthresh5))
			mAcount6 = np.sum(np.logical_and(Athreshvec, mod_base_table >= mAthresh6))
			barcode = barcodeDict[read_id]
			outputstring = outputstring+str(read_id)+'\t'+str(barcode)+'\t'+str(readlen)+'\t'+str(Acount)+'\t'+str(mAcount1)+'\t'+str(mAcount2)+'\t'+str(mAcount3)+'\t'+str(mAcount4)+'\t'+str(mAcount5)+'\t'+str(mAcount6)+'\t'+str(printstr)+'\n'
	return outputstring


files = list(glob.iglob(directory))
start = time.time()
results = Parallel(n_jobs=num_cores)(delayed(processFile)(files[i]) for i in range(len(files)))
stop = time.time()
elapsed=stop-start
print(elapsed)


fhout = open(outfile, 'w')
fhout.write(''.join(results))
fhout.close()


#alternative: extract guppy basemod scores, print *all scores* to text file for mCpG and mA, not just counts above thresholds
python

from ont_fast5_api.fast5_interface import get_fast5_file
import numpy as np
from ont_fast5_api.analysis_tools.base_tool import BaseTool
import re
import time
import glob
import csv
from joblib import Parallel, delayed
import multiprocessing

directory = r'/Data1/workspace/Data2/fast5/*.fast5'
outfile = '/Data1/Output_methylAll.txt'

qthreshInt=10
qthresh = chr(qthreshInt+33)

totA=0
totbases=0
num_cores = multiprocessing.cpu_count()

def processFile(fast5_filepath):
	outputstring=""
	with get_fast5_file(fast5_filepath, mode="r") as f5:
		for read_id in f5.get_read_ids():
			read = f5.get_read(read_id)
			latest_basecall = read.get_latest_analysis('Basecall_1D')
			mod_base_table = np.array(read.get_analysis_dataset(latest_basecall, 'BaseCalled_template/ModBaseProbs').transpose()[1,:])
			mod_base_tableC = np.array(read.get_analysis_dataset(latest_basecall, 'BaseCalled_template/ModBaseProbs').transpose()[3,:])
			name, sequence, _, qstring = read.get_analysis_dataset(latest_basecall,'BaseCalled_template/Fastq').strip().split('\n')
			qarray = np.array(tuple(qstring))
			seqarrayA = np.array(tuple(sequence),dtype=np.object)
			seqarrayC = seqarrayA + np.array(tuple(sequence[1:len(sequence)]+'N'),dtype=np.object)
			Athreshvec = np.logical_and(qarray >= qthresh, seqarrayA == 'A')
			Cthreshvec = np.logical_and(qarray >= qthresh, seqarrayC == 'CG')
			mAthreshvec = np.logical_and(Athreshvec, mod_base_table > 0)
			mCthreshvec = np.logical_and(Cthreshvec, mod_base_table > 0)
			readlen = len(sequence)
			Acount = np.sum(Athreshvec)
			mAcount = np.sum(mAthreshvec)
			Ccount = np.sum(Cthreshvec)
			mCcount = np.sum(mCthreshvec)
			printstr = re.sub('[\[\]\s]','',np.array2string(mod_base_table[mAthreshvec],separator=',',threshold=1000000000))
			printstr2=  re.sub('[\[\]\s]','',np.array2string(np.flatnonzero(mAthreshvec),separator=',',threshold=1000000000))
			printstrC = re.sub('[\[\]\s]','',np.array2string(mod_base_table[mCthreshvec],separator=',',threshold=1000000000))
			printstr2C=  re.sub('[\[\]\s]','',np.array2string(np.flatnonzero(mCthreshvec),separator=',',threshold=1000000000))
			outputstring = outputstring+str(read_id)+'\tA'+'\t'+str(readlen)+'\t'+str(Acount)+'\t'+str(mAcount)+'\t'+str(printstr)+'\t'+str(printstr2)+'\n'
			outputstring = outputstring+str(read_id)+'\tCG'+'\t'+str(readlen)+'\t'+str(Ccount)+'\t'+str(mCcount)+'\t'+str(printstrC)+'\t'+str(printstr2C)+'\n'
	return outputstring


files = list(glob.iglob(directory))
start = time.time()
results=""
results = Parallel(n_jobs=num_cores)(delayed(processFile)(files[i]) for i in range(len(files)))
stop = time.time()
elapsed=stop-start
print(elapsed)

fhout = open(outfile, 'w')
fhout.write(''.join(results))
fhout.close()



##run winnowmap (replace 01 02 03 with list of barcode numbers), submits each barcode to run in parallel. remove '&' from end of line to run serially
cd /Data1
for i in 01 02 03; do
	echo $i
	/Data2/software/Winnowmap/bin/winnowmap -W /Data2/reference/t2t/chm13.repetitive_k15.txt -ax map-ont /Data2/reference/t2t/chm13.draft_v1.0.fasta /Data1/pass/barcode$i/*.fastq | samtools view -b >/Data1/winnowmap.barcode.$i.bam &
done
for i in 01 02 03; do
	samtools sort winnowmap.barcode.$i.bam >winnowmap.barcode.$i.sorted.bam; samtools index -m 500M -@ 7 winnowmap.barcode.$i.sorted.bam; bedtools bamtobed -i winnowmap.barcode.$i.sorted.bam >winnowmap.barcode.$i.sorted.bed &
done


####merge guppy and winnow info
#copy CombineGuppyWinnowBarcode.V2.pl here
perl CombineGuppyWinnowBarcode.V2.pl /Data1 output_methylAthresh.txt 01 02 03


####get c/ci lad proportions
cd /Data1
for i in 01 02 03;do
	echo $i
	cut -f 1,2,3,9,10,11,12,13,14,15 winnowmap.barcode.$i.sorted.methylinfo.q10.bedgraph | bedtools map -F 0.51 -a /Data2/reference/t2t/cLAD_gold_final_100kb.chm13.sorted2.bed -b stdin -c 4,5,6,7,8,9,10 -o sum -g /Data2/reference/t2t/chm13.chromsizes.sorted.txt | awk 'BEGIN {sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;sum6=0;sum7=0} {sum1= sum1+$4;sum2= sum2+$5;sum3= sum3+$6;sum4= sum4+$7;sum5= sum5+$8;sum6= sum6+$9;sum7=sum7+$10} END {print sum1,sum2,sum3,sum4,sum5,sum6,sum7}'
	cut -f 1,2,3,9,10,11,12,13,14,15 winnowmap.barcode.$i.sorted.methylinfo.q10.bedgraph | bedtools map -F 0.51 -a /Data2/reference/t2t/ciLAD_gold_final_100kb.chm13.sorted2.bed -b stdin -c 4,5,6,7,8,9,10 -o sum -g /Data2/reference/t2t/chm13.chromsizes.sorted.txt | awk 'BEGIN {sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;sum6=0;sum7=0} {sum1= sum1+$4;sum2= sum2+$5;sum3= sum3+$6;sum4= sum4+$7;sum5= sum5+$8;sum6= sum6+$9;sum7=sum7+$10} END {print sum1,sum2,sum3,sum4,sum5,sum6,sum7}'
	cut -f 1,2,3,9,10,11,12,13,14,15 winnowmap.barcode.$i.sorted.methylinfo.q10.bedgraph | bedtools map -F 0.51 -a /Data2/reference/t2t/chm13.window100kb.sorted.bed -b stdin -c 4,5,6,7,8,9,10 -o sum -g /Data2/reference/t2t/chm13.chromsizes.sorted.txt | awk 'BEGIN {sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;sum6=0;sum7=0} {sum1= sum1+$4;sum2= sum2+$5;sum3= sum3+$6;sum4= sum4+$7;sum5= sum5+$8;sum6= sum6+$9;sum7=sum7+$10} END {print sum1,sum2,sum3,sum4,sum5,sum6,sum7}'
done


#make bedgraph files for plotting
for i in 01 02 03;do
	echo $i
	cut -f 1,2,3 winnowmap.barcode.$i.sorted.methylinfo.bedgraph | bedtools coverage -a /Data2/reference/t2t/chm13.window100kb.sorted.bed -mean -b stdin -bed -sorted -g /Data2/reference/t2t/chm13.chromsizes.sorted.txt >CoverageWinnow.100kb.barcode.$i.bedgraph
	cut -f 1,2,3 winnowmap.barcode.$i.sorted.methylinfo.q10.bedgraph | bedtools coverage -a /Data2/reference/t2t/chm13.window100kb.sorted.bed -mean -b stdin -bed -sorted -g /Data2/reference/t2t/chm13.chromsizes.sorted.txt >CoverageWinnow.100kb.barcode.$i.q10.bedgraph
done
for i in 01 02 03;do
	echo $i
	cut -f 1,2,3,9,14 winnowmap.barcode.$i.sorted.methylinfo.q10.bedgraph | bedtools map -a /Data2/reference/t2t/chm13.window100kb.sorted.bed -b stdin -c 4,5 -o sum -g /Data2/reference/t2t/chm13.chromsizes.sorted.txt | awk -v OFS="\t" '{print $1,$2,$3,$5/($4+0.00001),$4,$5}' >winnowmap.barcode.$i.sorted.methylinfo.q10.thresh80prop.bedgraph
	cut -f 1,2,3,9,13 winnowmap.barcode.$i.sorted.methylinfo.q10.bedgraph | bedtools map -a /Data2/reference/t2t/chm13.window100kb.sorted.bed -b stdin -c 4,5 -o sum -g /Data2/reference/t2t/chm13.chromsizes.sorted.txt | awk -v OFS="\t" '{print $1,$2,$3,$5/($4+0.00001),$4,$5}' >winnowmap.barcode.$i.sorted.methylinfo.q10.thresh90prop.bedgraph
	cut -f 1,2,3,9,12 winnowmap.barcode.$i.sorted.methylinfo.q10.bedgraph | bedtools map -a /Data2/reference/t2t/chm13.window100kb.sorted.bed -b stdin -c 4,5 -o sum -g /Data2/reference/t2t/chm13.chromsizes.sorted.txt | awk -v OFS="\t" '{print $1,$2,$3,$5/($4+0.00001),$4,$5}' >winnowmap.barcode.$i.sorted.methylinfo.q10.thresh95prop.bedgraph
done
for i in 01 02 03;do
	echo $i
	cut -f 1,2,3,9,10 winnowmap.barcode.$i.sorted.methylinfo.bedgraph | bedtools map -a /Data2/reference/t2t/chm13.window100kb.sorted.bed -b stdin -c 4,5 -o sum -g /Data2/reference/t2t/chm13.chromsizes.sorted.txt | awk -v OFS="\t" '{print $1,$2,$3,$5/($4+0.00001),$4,$5}' >winnowmap.barcode.$i.sorted.methylinfo.thresh50prop.bedgraph
	cut -f 1,2,3,9,14 winnowmap.barcode.$i.sorted.methylinfo.bedgraph | bedtools map -a /Data2/reference/t2t/chm13.window100kb.sorted.bed -b stdin -c 4,5 -o sum -g /Data2/reference/t2t/chm13.chromsizes.sorted.txt | awk -v OFS="\t" '{print $1,$2,$3,$5/($4+0.00001),$4,$5}' >winnowmap.barcode.$i.sorted.methylinfo.thresh80prop.bedgraph
	cut -f 1,2,3,9,13 winnowmap.barcode.$i.sorted.methylinfo.bedgraph | bedtools map -a /Data2/reference/t2t/chm13.window100kb.sorted.bed -b stdin -c 4,5 -o sum -g /Data2/reference/t2t/chm13.chromsizes.sorted.txt | awk -v OFS="\t" '{print $1,$2,$3,$5/($4+0.00001),$4,$5}' >winnowmap.barcode.$i.sorted.methylinfo.thresh90prop.bedgraph
	cut -f 1,2,3,9,12 winnowmap.barcode.$i.sorted.methylinfo.bedgraph | bedtools map -a /Data2/reference/t2t/chm13.window100kb.sorted.bed -b stdin -c 4,5 -o sum -g /Data2/reference/t2t/chm13.chromsizes.sorted.txt | awk -v OFS="\t" '{print $1,$2,$3,$5/($4+0.00001),$4,$5}' >winnowmap.barcode.$i.sorted.methylinfo.thresh95prop.bedgraph
done
for i in 01 02 03;do
	bgzip winnowmap.barcode.$i.sorted.methylinfo.thresh50prop.bedgraph
	tabix -f -p bed winnowmap.barcode.$i.sorted.methylinfo.thresh50prop.bedgraph.gz
	bgzip winnowmap.barcode.$i.sorted.methylinfo.thresh80prop.bedgraph
	tabix -f -p bed winnowmap.barcode.$i.sorted.methylinfo.thresh80prop.bedgraph.gz
	bgzip winnowmap.barcode.$i.sorted.methylinfo.thresh90prop.bedgraph
	tabix -f -p bed winnowmap.barcode.$i.sorted.methylinfo.thresh90prop.bedgraph.gz
	bgzip winnowmap.barcode.$i.sorted.methylinfo.thresh95prop.bedgraph
	tabix -f -p bed winnowmap.barcode.$i.sorted.methylinfo.thresh95prop.bedgraph.gz
done



#optional: split fast5 files by barcode
cd /Data1
for i in 01 02 03; do
	grep "barcode$i" /Data1/sequencing_summary_ALL.txt | cut -f 2 > /Data1/barcode$i.reads.txt
done
for i in 01 02 03; do
	mkdir /Data2/fast5/split/barcode$i
	fast5_subset -i /Data1/workspace/Data2/fast5 -s /Data2/fast5/split/barcode$i -l /Data1/barcode$i.reads.txt -n 4000 -t 1 &
done


#optional: run megalodon
for i in 01 02 03; do
	megalodon /Data2/fast5/split/barcode$i --output-directory /Data1/megalodon/barcode$i --overwrite --guppy-params "-d /opt/ont/guppy/data" --guppy-server-path /usr/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg  --outputs mappings mod_mappings --reference /Data2/reference/t2t/chm13.draft_v1.0.fasta --device cuda:all --processes 92 --mod-min-prob 0
done
#sort output
for i in 01 02 03;do
	samtools sort -m 7000M -@ 4 -o /Data1/megalodon/barcode$i/mod_mappings.$i.sorted.bam /Data1/megalodon/barcode$i/mod_mappings.bam 
done
#make a "pileup" file, which lists all the mod probability scores from all reads overlapping each A or T in the reference
for i in 01 02 03; do
	samtools view -h /Data1/megalodon/barcode$i/mod_mappings.sorted.bam | perl ParseModMappings.hash.pl pileup.$i.bed
	bgzip pileup.$i.bed
done
