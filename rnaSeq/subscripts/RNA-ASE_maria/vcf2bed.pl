#!/usr/bin/perl
use strict;

open(FH, $ARGV[0]) || die("Cannot open input file\n");
while (my $line = <FH>) {
    chomp $line;
    my @ts = split /\t/, $line;
    my $snpID = $ts[2];
    my $chr = $ts[0];
    my $end = $ts[1];
    my $start = $end - 1;
    
    print "$chr\t$start\t$end\t$snpID\n";
}
close(FH);

print STDERR "Done!\n";
