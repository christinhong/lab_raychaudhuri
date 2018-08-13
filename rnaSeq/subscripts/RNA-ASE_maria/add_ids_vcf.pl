#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;

################################################################################
################################################################################

my $vcf = $ARGV[0];

open(FH, $vcf) || die("Cannot open input vcf file\n");

while (my $line = <FH>) {
	chomp $line;
	if ($line =~ /^#/) {
		print "$line\n";
	}
	else{
		my @ts = split /\t/, $line;
		my $id = "snp_" . $ts[0] . "_" . $ts[1];
		print "$ts[0]\t$ts[1]\t$id\t$ts[3]\t$ts[4]\t$ts[5]\t$ts[6]\t$ts[7]\t$ts[8]\t$ts[9]\n";
    }
}      
close FH;
