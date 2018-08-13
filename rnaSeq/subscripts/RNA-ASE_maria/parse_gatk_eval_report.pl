#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;

# usage: perl parse_gatk_eval_report.pl outFile allEvalFilesToParse
# i.e.: perl parse_gatk_eval_report.pl summary_eval_1KGP.txt /medpop/srlab/cd4timelinepilot/rnaseq/141125_PR1504/variantCalling/eval.*_1KGP

my @files = @ARGV;
my $outfile = $files[0];

open (OUT, ">$outfile") || die ("Cannot open new file $outfile\n");

print OUT "File\tEvalVariants\tnovelSites\tnVariantsAtComp\tcompRate\tnConcordant\tconcordantRate\n";

for(my $i = 1; $i < @files; $i++){
	chomp($files[$i]);
	my @path = split /\//, $files[$i];
	my $last = scalar(@path);
	my $fileName = $path[$last - 1] ;
	open(FH, "$files[$i]") || die ("Cannot open file $files[$i]\n");
	while(my $line = <FH>){
		chomp $line;
		if($line =~ /CompOverlap\s+dbsnp/){
			my @ts = split /\s+/, $line;
			print OUT "$fileName\t$ts[5]\t$ts[6]\t$ts[7]\t$ts[8]\t$ts[9]\t$ts[10]\n";
			last;
		}
	}
}
