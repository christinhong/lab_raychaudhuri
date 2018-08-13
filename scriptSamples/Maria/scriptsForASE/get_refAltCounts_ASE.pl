#!/usr/bin/perl
use strict;

#usage: perl get_allelesSeen_access.pl hets_vcf_file_with_ids all_access_ase_files

my $out_ref = "summary_ref_counts.txt";
my $out_alt = "summary_alt_counts.txt";

my @files = @ARGV;#the @ARGV vector are all the files listed, for example /medpop/rnaseq/tcells/G72282/SC*/v1/subread/log.txt
my $vcf = $files[0];

my %VARS;
my %REF;
my %ALT;
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
		#my $bas = $ts[5];
		#my $totalCount = $ts[8];
		my $refCount = $ts[6];
		my $altCount = $ts[7];
		#my $refRatio = $refCount/$totalCount;
		my $varID = $ts[1];
		push(@{$TEMP{$varID}}, $refCount);
		push(@{$TEMP{$varID}}, $altCount);
	}
	close FH;
	foreach my $var (keys %VARS){
		if(exists($TEMP{$var})){
			push(@{$REF{$var}}, $TEMP{$var}[0]);
			push(@{$ALT{$var}}, $TEMP{$var}[1]);
		}else{
			push(@{$REF{$var}}, "NA");
			push(@{$ALT{$var}}, "NA");
		}
	}
}

#Print hashes

open(REF, ">$out_ref") || die ("Cannot open new file $out_ref\n");
open(ALT, ">$out_alt") || die ("Cannot open new file $out_alt\n");


print REF "SNP_ID";
print ALT "SNP_ID";


for(my $i = 1; $i < @samples; $i++){
	print REF "\t$samples[$i]";
	print ALT "\t$samples[$i]";
}
print REF "\n";
print ALT "\n";


foreach my $var (sort keys %VARS){
	print REF "$var";
	print ALT "$var";
	for(my $i = 0; $i < @{$REF{$var}}; $i++){
		print REF "\t$REF{$var}[$i]";
		print ALT "\t$ALT{$var}[$i]";
	}
	print REF "\n";
	print ALT "\n";
}

close REF;
close ALT;
