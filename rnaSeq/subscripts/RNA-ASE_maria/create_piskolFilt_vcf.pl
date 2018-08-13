#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;

################################################################################
################################################################################

my $listPassSNPs = $ARGV[0];
my $vcf = $ARGV[1];

my %SNPs;

open(FH, $listPassSNPs) || die("Cannot open list of SNPs with chr & position as first and second columns\n");
while (my $line = <FH>) {
	chomp $line;
	my @ts = split /\t/, $line;
	my $snp = $ts[0] . "-" . $ts[1];
	$SNPs{$snp} = 1;
}
close FH;


open(FH, $vcf) || die("Cannot open input vcf file\n");

my $flag = 0;
while (my $line = <FH>) {
	chomp $line;
	if ($line =~ /^#/) {
		print "$line\n";
		if($flag == 0 and $line =~ /^##FILTER=<ID=SnpCluster/) {
			print "##FILTER=<ID=SNPiR,Description=\"SNP discarded by Piskol et al filters (rmIndel, rmintron, rmsk, rmhom, rmedt)\">\n";
			$flag = 1;
		}
	}
	else{
		my @ts = split /\t/, $line;
		my $snp = $ts[0] . "-" . $ts[1];
		if(exists($SNPs{$snp})){
			print "$line\n";
		}
		else {#if SNP is not in list then it was filtered out by piskol et al filters
			print "$ts[0]\t$ts[1]\t$ts[2]\t$ts[3]\t$ts[4]\t$ts[5]\tSNPiR\t$ts[7]\t$ts[8]\t$ts[9]\n";
		}
    }
}      
close FH;
