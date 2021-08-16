use warnings;
use strict;

my $date = "20200320";
my @barcodelist;
if(defined $ARGV[0]){
	($date,@barcodelist) = @ARGV;
}else{
	die "provide list of barcodes\n";
}
my $basedir = "/Data1";


my $infile1 = "$basedir/$date"."_methylA99.txt";
print "opening $infile1\n";


my %methylmemory;
open(IN, $infile1);
while(my $line = <IN>){
	chomp($line);
	if($line=~/^(\S+)\s+(.+)$/){
		$methylmemory{$1}=$2;
	}
}
close IN;


my %barcodes;
my %barcodememory;
foreach my $barcode(@barcodelist){
	chomp($barcode);
	$barcodes{$barcode}=0;
	my $infile2 = "$basedir/winnowmap.barcode.$barcode.sorted.bed";
	print "opening $infile2\n";
	my $outfilehere = "$basedir/winnowmap.barcode.$barcode.sorted.methylinfo.bedgraph";
	my $outfilehere10 = "$basedir/winnowmap.barcode.$barcode.sorted.methylinfo.q10.bedgraph";
	open(OUT,'>'.$outfilehere);
	open(OUT10,'>'.$outfilehere10);
	open(IN, $infile2);
	while(my $line = <IN>){
		chomp($line);
		if($line=~/^\S+\s+\S+\s+\S+\s+(\S+)\s+(\d+)\s+/){
			if(exists $methylmemory{$1}){
				print OUT "$line\t$methylmemory{$1}\n";
				if($2>=10){
					print OUT10 "$line\t$methylmemory{$1}\n";
				}
			}
		}
	}
	close IN;
	close OUT;
	close OUT10;
}



