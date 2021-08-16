#ParseModMappings.hash.pl
#Nicolas Altemose, 2020
#purpose: extract a list of all mA mod probabilities from the mod_mappings bam file output by megalodon
#note: this is likely to be highly memory intensive
#output: a bed file with one line per A base in the covered genome, with the 4th column containing a comma-separated list of mod probabilities
#sample output line:
#chr1	300000	300001	245,130,0,0,0
#^this A base was overlapped by 5 reads, with the following mod probabilities per read: 245,130,0,0,0
#input: pipe in sorted sam lines with header, see usage below
#it's perl, I'm sorry if it's not very readable

use warnings;
use strict;
my $tic = time; #start runtime clock

my $usage = 'samtools view -h <sorted_mod_mappings_bam_file> | perl ParseModMappings.hash.pl <output_filename>';
my $outfile = "Pileup.bed";

#check for required outfile name
if(defined $ARGV[0]){
	$outfile = $ARGV[0];
	chomp($outfile);
}else{
	die "$usage\n\n";
}

#define subroutines to manipulate sequence
# Complement a sequence
sub comp {
    my ($seq)=@_;
    $seq =~ tr/ACGTUMRWSYKVHDBN/TGCAAKYWSRMBDHVN/;
    return $seq;
}

# Reverse complement a sequence
sub rc {
    my ($seq)=@_;
    $seq = reverse($seq);
    return comp($seq);
}


my %chrsizes;
my %chrarray=();

open(OUT,'>'.$outfile);
my $prevchr="";
my $readnum=0;

#parse sam lines being piped in from STDIN
while(my $line = <STDIN>){
	chomp($line);
	if($line=~/^\@/){ #read in chromosome sizes from header
		if($line=~/^\@SQ\s+SN:(\S+)\s+LN:(\d+)$/){
			$chrsizes{$1}=$2;
		}
		next;
	}
	
	$readnum++;
	my ($qname,$flag,$chr,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual,@aux) = split("\t",$line); #store each sam field into separate variables 

	#once a chromosome has been entirely scanned through, print out all A positions covered by one or more reads
	if($chr ne $prevchr){
		if($prevchr ne ""){
			my @sortedkeys = sort {$a <=> $b} keys %chrarray;
			foreach my $key(@sortedkeys){
				my $end = $key+1;
				$chrarray{$key} =~ s/\,$//; #remove trailing commas
				print OUT "$prevchr\t$key\t$end\t$chrarray{$key}\n";
			}
		}
		%chrarray=();
	}
	
	
	#collect mod names and probabilities 
    my ($mod_str)="@aux"=~m/M[mM]:Z:(\S+)/;
    my ($mod_prob)="@aux"=~m/M[lL]:B:C,(\S+)/;
	
	
	#interpret the strand & rc if needed 
	# nb: All counting is from 5' end irrespective of BAM FLAG 0x10
	my $dir = $flag & 0x10 ? "-" : "+";
	$seq = rc($seq) if ($dir eq "-"); #rc if needed
	my $seqlen = length($seq);
    my @seq = split("", $seq);                 # array of seq bases
    my @seqt = @seq;                           # orientation as shown in SAM
	
   
	if($readnum % 1000 == 0){
		my $toc = time;
		my $elapsed = $toc-$tic;
		printf("\n\n$readnum completed in: %02d:%02d:%02d\n\n (currently processing $chr and printing to $outfile)", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));
	}
	
    $mod_str =~ s/;$//;                        # trim last ; to aid split()
    my @mods = split(";", $mod_str);
    my @probs = split(",", $mod_prob);
    my $pnum = 0;

	#loop through mods and probabilities; ignore mC entirely for now
    foreach (@mods) {
		my ($base, $strand, $types, $poslist) = $_ =~ m/([A-Z])([-+])([^,]+),(.*)/;
	
		if($base eq "A"){
			my $i = 0; # I^{th} bosition in sequence
			foreach my $delta (split(",", $poslist)) {
	 	   		#skip $delta occurences of $base; this is the trickiest part of parsing
	 	   		#it figures out which mod probability goes with which position in the reference
	 	 		do {
					if ($base eq $seq[$i]){
						$delta--;
						if($delta>=0){
							if($dir eq "-"){
								$chrarray{$pos+$seqlen-$i-2}.="0,"; #assign probability of 0 to skipped bases
							}else{
								$chrarray{$pos+$i-1}.="0,";
							}
						}
					}
					$i++;
	   	 		}while ($delta >= 0);
	   		 	$i--;

				my $pval = $probs[$pnum++]; #get the probability for this position
	   
	  	 		# Add to top or bottom seq
			   	if ($dir eq "-") {
					$chrarray{$pos+$seqlen-$i-2}.="$pval,"; #add this probability to the string in memory for that position (or initialize this string) 
	   	 		}else {
					$chrarray{$pos+$i-1}.="$pval,";
	   	 		}
	   		 	$i++;
			}
		}else{
			$pnum += scalar(split(",",$poslist));
		}
    }
	$prevchr = $chr;
}

			#print last chromosome
			my @sortedkeys = sort {$a <=> $b} keys %chrarray;
			foreach my $key(@sortedkeys){
				my $end = $key+1;
				$chrarray{$key} =~ s/\,$//;
				print OUT "$prevchr\t$key\t$end\t$chrarray{$key}\n";
			}

close OUT;


	#print runtime info to STDOUT:
	my $toc = time;
	my $elapsed = $toc-$tic;
	printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));

