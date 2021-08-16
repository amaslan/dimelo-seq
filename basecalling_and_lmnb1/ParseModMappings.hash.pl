use warnings;
use strict;
my $tic = time; 

my $outfile = "Pileup.bed";

if(defined $ARGV[0]){
	$outfile = $ARGV[0];
	chomp($outfile);
}


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
while(my $line = <STDIN>){
	chomp($line);
	if($line=~/^\@/){
		if($line=~/^\@SQ\s+SN:(\S+)\s+LN:(\d+)$/){
			$chrsizes{$1}=$2;
		}
		next;
	}
	
	$readnum++;
	my ($qname,$flag,$chr,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual,@aux) = split("\t",$line);

	if($chr ne $prevchr){
		if($prevchr ne ""){
			my @sortedkeys = sort {$a <=> $b} keys %chrarray;
			foreach my $key(@sortedkeys){
				my $end = $key+1;
				$chrarray{$key} =~ s/\,$//;
				print OUT "$prevchr\t$key\t$end\t$chrarray{$key}\n";
			}
		}
		%chrarray=();
	}
	
    my $dir = $flag & 0x10 ? "-" : "+";
    my ($mod_str)="@aux"=~m/M[mM]:Z:(\S+)/;
    my ($mod_prob)="@aux"=~m/M[lL]:B:C,(\S+)/;

    # All counting is from 5' end irrespective of BAM FLAG 0x10
    $seq = rc($seq) if ($dir eq "-");
	my $seqlen = length($seq);
    my @seq = split("", $seq);                 # array of seq bases
    my @seqt = @seq;                           # orientation as shown in SAM
#    my @seqb = split("", comp($seq));          # plus a complemented copy
	if($readnum % 1000 == 0){
		my $toc = time;
		my $elapsed = $toc-$tic;
		printf("\n\n$readnum completed in: %02d:%02d:%02d\n\n (currently processing $chr and printing to $outfile)", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));
	}
    $mod_str =~ s/;$//;                        # trim last ; to aid split()
    my @mods = split(";", $mod_str);
    my @probs = split(",", $mod_prob);
    my $pnum = 0;
#    my $psize = scalar(@probs);

    foreach (@mods) {
		my ($base, $strand, $types, $poslist) = $_ =~ m/([A-Z])([-+])([^,]+),(.*)/;
	
		if($base eq "A"){
			my $i = 0; # I^{th} bosition in sequence
			foreach my $delta (split(",", $poslist)) {
	 	   		# Skip $delta occurences of $base
	 	 		do {
					if ($base eq $seq[$i]){
						$delta--;
						if($delta>=0){
							if($dir eq "-"){
								$chrarray{$pos+$seqlen-$i-2}.="0,";
							}else{
								$chrarray{$pos+$i-1}.="0,";
							}
						}
					}
					$i++;
	   	 		}while ($delta >= 0);
	   		 	$i--;

				my $pval = $probs[$pnum++];
	   
	  	 		# Add to top or bottom seq
			   	if ($dir eq "-") {
					$chrarray{$pos+$seqlen-$i-2}.="$pval,";
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
	#print "$pnum of $psize recorded\n";
}
#print "here!\n";

			my @sortedkeys = sort {$a <=> $b} keys %chrarray;
			foreach my $key(@sortedkeys){
				my $end = $key+1;
				$chrarray{$key} =~ s/\,$//;
				print OUT "$prevchr\t$key\t$end\t$chrarray{$key}\n";
			}

close OUT;


	#put this at the end of your code:
	my $toc = time;
	my $elapsed = $toc-$tic;
	printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));

