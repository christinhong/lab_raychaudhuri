#!/usr/bin/perl

use strict;
#use warnings;
#use diagnostics;

################################################################################
################################################################################


my $vcf = $ARGV[0];
my $het = "0/1";
my $flag = 0;

open(FH, $vcf) || die("Cannot open input vcf file\n");

while (my $line = <FH>) {
	chomp $line;
	if ($line =~ /^#/) {
		print "$line\n";
		if($flag == 0 and $line =~ /^##FILTER=<ID=LowQual/) {
			print "##FILTER=<ID=Hom,Description=\"SNP is not heterozygous\">\n";
			$flag = 1;
		}
	}
	else{
		my @ts = split /\t/, $line;
		my @data = split /:/, $ts[9];
		my $G = $data[0];
		if($G eq $het){
			print "$ts[0]\t$ts[1]\t$ts[2]\t$ts[3]\t$ts[4]\t$ts[5]\t$ts[6]\t$ts[7]\t$ts[8]\t$ts[9]\n";
		}
		else {
			print "$ts[0]\t$ts[1]\t$ts[2]\t$ts[3]\t$ts[4]\t$ts[5]\tHom\t$ts[7]\t$ts[8]\t$ts[9]\n";
		}
    }
}      
close FH;
