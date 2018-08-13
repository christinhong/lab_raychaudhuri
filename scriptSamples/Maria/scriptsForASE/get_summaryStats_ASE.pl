#!/usr/bin/perl
use strict;

#usage: perl get_allelesSeen_access.pl hets_vcf_file_with_ids all_access_ase_files

my $out_bas = "summary_both_alleles_seen.txt";
my $out_cov = "summary_total_count.txt";
my $out_refRatio = "summary_refRatio.txt";

my @files = @ARGV;#the @ARGV vector are all the files listed, for example /medpop/rnaseq/tcells/G72282/SC*/v1/subread/log.txt
my $vcf = $files[0];

my %VARS;
my %BAS;
my %COV;
my %REFRATIO;
my @samples;

open(VCF, "$vcf") || die ("Cannot open file with het variants\n");
while(my $line = <VCF>){
	chomp $line;
	unless($line =~ /^#/) {
		my @ts = split /\t/, $line;
		push(@{$VARS{$ts[2]}}, 1);
	}
}      
close VCF;
push(@samples, 1);

for(my $i = 1; $i < @files; $i++){
	my %TEMP;
	chomp($files[$i]);
	my @path = split /\//, $files[$i];
	my $last = scalar(@path);
	my $sample_id = $path[$last - 1];
	push(@samples, $sample_id);	
	open(FH, "$files[$i]") || die ("Cannot open file $files[$i]\n");
	my $header = <FH>;
	while(my $line = <FH>){
		chomp $line;		
		my @ts = split /\t/, $line;
		my $bas = $ts[5];
		my $totalCount = $ts[8];
		my $refCount = $ts[6];
		my $refRatio = $refCount/$totalCount;
		my $varID = $ts[1];
		push(@{$TEMP{$varID}}, $bas);
		push(@{$TEMP{$varID}}, $totalCount);
		push(@{$TEMP{$varID}}, $refRatio);
	}
	close FH;
	foreach my $var (keys %VARS){
		if(exists($TEMP{$var})){
			push(@{$BAS{$var}}, $TEMP{$var}[0]);
			push(@{$COV{$var}}, $TEMP{$var}[1]);
			push(@{$REFRATIO{$var}}, $TEMP{$var}[2]);
		}else{
			push(@{$BAS{$var}}, "NA");
			push(@{$COV{$var}}, "NA");
			push(@{$REFRATIO{$var}}, "NA");
		}
	}
}

#Print hashes

open(BAS, ">$out_bas") || die ("Cannot open new file $out_bas\n");
open(COV, ">$out_cov") || die ("Cannot open new file $out_cov\n");
open(REF, ">$out_refRatio") || die ("Cannot open new file $out_refRatio\n");

print BAS "SNP_ID";
print COV "SNP_ID";
print REF "SNP_ID";

for(my $i = 1; $i < @samples; $i++){
	print BAS "\t$samples[$i]";
	print COV "\t$samples[$i]";
	print REF "\t$samples[$i]";
}
print BAS "\n";
print COV "\n";
print REF "\n";

foreach my $var (sort keys %VARS){
	print BAS "$var";
	print COV "$var";
	print REF "$var";
	for(my $i = 0; $i < @{$BAS{$var}}; $i++){
		print BAS "\t$BAS{$var}[$i]";
		print COV "\t$COV{$var}[$i]";
		print REF "\t$REFRATIO{$var}[$i]";
	}
	print BAS "\n";
	print COV "\n";
	print REF "\n";
}

close BAS;
close COV;
close REF;
