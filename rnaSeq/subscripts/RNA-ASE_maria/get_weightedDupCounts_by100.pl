#!/usr/bin/perl
use strict;

#usage: perl get_weightedDupCounts.pl ase_output_file bam_file > weighted_dup_counts_output_file
#needs: samtools view
#you may need to change these values:
my $minMapQ = 20; #minimum mapping quality
my $base_quality_threshold = 10; #minimum base quality

my $datestring = localtime();
print STDERR "Local date and time $datestring\n";

my %noALT;
my %notCALLED;

#get header of bam file
my $bam = $ARGV[1];
my @parts = split /\//, $bam;
my $size = scalar(@parts);
my @subparts = split /\./, $parts[$size -1];
print STDERR "The suffix for temp files is $subparts[0].\n";

system("samtools view -H -o header.$subparts[0].temp.sam $bam");

#open ase output file
open(ASE, $ARGV[0]) || die("Cannot open input file\n");
my $discarded_header = <ASE>;

print "INDIVIDUAL\tRSID\tBOTH_ALLELES_SEEN\tREF_ALLELE\tALT_ALLELE\tDUP_WEIGHTED_REF_COUNT\tDUP_WEIGHTED_ALT_COUNT\tDUP_WEIGHTED_REF_RATIO\n";

#for each het SNP
while (my $line = <ASE>) {
    my %FWD; #hash for storing reads in forward strand, by 5' position (all dups together in a vector)
	my %REV; #hash for storing reads in reverse strand, by 5' position (all dups together in a vector)

    chomp $line;
    my @ts = split /\t/, $line;

    my $indiv = $ts[0];
    my $snpID = $ts[1];
    my $chr = "chr" . $ts[2];
    my $pos = $ts[3];
    my $bas = $ts[5];
    
    my @calledAlleles = split /\//, $ts[4]; 
    
    if($bas == 1 or $bas == 0){
	    
	    #define $refAllele and $altAllele and make sure they are OK, if not, report in hash and print at the end
	    my $refAllele = $ts[13];
	    my $altAllele = $ts[14];
	    
	    if(!($refAllele =~ /[ATCG]/)){
	    	die ("Did not find ref allele for variant $snpID !\n");
	    }elsif(!($altAllele =~ /[ATCG]/)){
	    	$noALT{$snpID} = 1;
	    	$altAllele = $calledAlleles[1]; #take info of alternative allele from called variant
	    }
	    if($refAllele ne $calledAlleles[0] or $altAllele ne $calledAlleles[1]){
	    	$notCALLED{$snpID} = 1;
	    }
	    
	    #print STDERR "\nSTARTING WITH VARIANT $snpID\nMy ref allele is $refAllele, my alt allele is $altAllele\n";
	    
		#do samtools view on the position (like: samtools view chr:pos-pos -o tmp.sam)
		#requiring here already minimum mapping quality
		system("samtools view -q $minMapQ -o $subparts[0].temp.sam $bam $chr:$pos-$pos");
		
		#open tmp.sam
		open(SAM, "$subparts[0].temp.sam") || die("Cannot open temp sam file\n");
		my $cntReadOK = 0;
		my $cntAlleleOK = 0;
		my $cntBadNt = 0;
		my $cntA = 0;
		my $cntC = 0;
		my $cntG = 0;
		my $cntT = 0;
		
		#FOR EACH READ
		while (my $read_line = <SAM>) {
			
			#print STDERR "I am in this read line from temp.sam:\n$read_line";
	    	
	    	chomp $read_line;
	    	my @ts = split /\t/, $read_line;
	    	
	    	my $start_pos = $ts[3];
	    	my $cigar = $ts[5];
	    	my $read_seq = $ts[9];
	    	my $base_quals = $ts[10];
	    	
	    	my @read = split //, $read_seq;
	    	my @bq = split //, $base_quals;    	
	        my @cigar_numbers = split /[A-Z]/, $cigar;#split by letters to get numbers
			my @cigar_letters = split /\d+/, $cigar; #split by numbers to get letters
			shift @cigar_letters;#shifted because first element of array will be empty because the first character of cigar is number and the split function creates an empty element		
	
			### EXTRACT BASE AND BASE QUALITY FROM DESIRED VARIANT POSITION
			
			my $read_idx = 0;
			my $ref_idx = 0;
			my $flag = 0;
			my $final_read_pos;
			
			my $ref_objective = $pos - $start_pos; # SNP position - read start position in sam file
			#print STDERR "My final read position BEFORE reading cigar is $final_read_pos\n";
			for(my $cigar_idx = 0; $cigar_idx < @cigar_numbers; $cigar_idx++){
			#	print STDERR "go\t";
				if($ref_idx == $ref_objective){
					$final_read_pos = $read_idx;
					last;
				}
				if($cigar_letters[$cigar_idx] eq "S"){#soft clipped bases
					$read_idx = $read_idx + $cigar_numbers[$cigar_idx];
					#ref index does not change
				}
				elsif($cigar_letters[$cigar_idx] eq "M"){#match/mismatch
					
					#for 5'coord
					#$end_pos = $end_pos + $cigar_numbers[$i];#advance end_pos by number of matches/mismatches
					
					#for extracting exact base
					for(my $m_idx = 0; $m_idx < $cigar_numbers[$cigar_idx]; $m_idx++){				
						$read_idx++;
						$ref_idx++;
						if($ref_idx == $ref_objective){
							$final_read_pos = $read_idx;
							$flag = 1;
							last;#exit from match/mismatch loop
						}
					}
				}
				elsif($cigar_letters[$cigar_idx] eq "I"){#insertion
					#no change in end_pos
					$read_idx = $read_idx + $cigar_numbers[$cigar_idx];
					#ref index does not change
				}
				elsif($cigar_letters[$cigar_idx] eq "D"){#deletion
					
					#for 5'coord
					#$end_pos = $end_pos + $cigar_numbers[$i];#advance end_pos by number of bases deleted with respect to ref
					
					$ref_idx = $ref_idx + $cigar_numbers[$cigar_idx];
					if($ref_idx == $ref_objective){
						$final_read_pos = $read_idx;
						$flag = 1;
					}
					elsif($ref_idx > $ref_objective){#SNP position of interest in deleted from read, so no value assigned to final_read_pos
						$flag = 1;
					}
				}
				else{
					print STDERR "For this variant: $snpID\n";
					print STDERR "For this sam file line: $read_line\n";
					die ("Not recognized letter in cigar string! Letter found: $cigar_letters[$cigar_idx]\n");
				}
				if($flag == 1){
					last;#exit from cigar loop
				}
			}
			#print STDERR "\nMy final read position is $final_read_pos\n\n\n";
			unless(!defined($final_read_pos)){
				
				$cntAlleleOK++;
				my $nucleotide = uc($read[$final_read_pos]);
				my $base_quality = ord($bq[$final_read_pos]) - 33;
				
				if($nucleotide ne $refAllele and $nucleotide ne $altAllele){
	        	    #print STDERR "\tSkipping read because not ref nor alt allele: $nucleotide...\n";
	        	    #print STDERR "\t$read_seq\n";
	        	    #print STDERR "\tfinal read pos is: $final_read_pos\n";
	        	    #print STDERR "\tone before $final_read_pos is $read[$final_read_pos - 1]\n";
	        	    #print STDERR "\tone after $final_read_pos is $read[$final_read_pos + 1]\n";
	        	    $cntBadNt++;
	        	    next; #skip this read because it has another allele
	        	}
	        	 #if base qual OK
	        	if ($base_quality >= $base_quality_threshold) {
	        		#print STDERR "**Read allele OK and base quality OK**\n";
	        		$cntReadOK++;
					#1# Infer strand from flag
					my $strand;
					if($ts[1] & 0x0010){
						$strand = "-";
					}else{
						$strand = "+";
					}
		
		        	#2# If rev strand, infer 5' coord from cigar string
					
					#if strand is fwd, 
						#push to has %REV with key end and value base
					my $end_pos = $start_pos;		
					for ( my $i = 0; $i<@cigar_numbers; $i++){
						if($cigar_letters[$i] eq "S"){#sof clipped bases
							#no change in end_pos
						}
						elsif($cigar_letters[$i] eq "M"){#match/mismatch
							$end_pos = $end_pos + $cigar_numbers[$i];#advance end_pos by number of matches/mismatches
						}
						elsif($cigar_letters[$i] eq "I"){#insertion
							#no change in end_pos
						}
						elsif($cigar_letters[$i] eq "D"){#deletion
							$end_pos = $end_pos + $cigar_numbers[$i];#advance end_pos by number of bases deleted with respect to ref
						}
						else{
							die ("Not recognized letter in cigar string! Letter found: $cigar_letters[$i]\n");
						}
					}
					my $coords = $start_pos . "-" . $end_pos;
					#print STDERR "Coordinates of read with cigar $cigar are: $coords\n";
					if($strand eq "+") {
						#push to hash %FWD with key start and value base
						push(@{$FWD{$coords}}, $nucleotide);
			#			print STDERR "\nThe read is in FWD strand\n";
					}else{#if strand is rev, 
						push(@{$REV{$coords}}, $nucleotide);
				#		print STDERR "\nThe read is in REV strand\n";					
					}			
				}
			}
			
		}
		
		close SAM;
		
		#for each my $dupGroup keys %FWD
		# define $weighted_ref $weighted_alt
		my $weighted_ref = 0;
		my $weighted_alt = 0;
		
		my @keys_fwd = keys %FWD;
		my $num_keys_fwd = scalar(@keys_fwd);
		
		#print STDERR "\nFinished going through each read...\n";
		unless($num_keys_fwd == 0){
	
			#print STDERR "\nNow going through hash with dup groups FWD strand\n";
			foreach my $dupGroup (keys %FWD) {
				#print STDERR "MY KEY IS $dupGroup\n";
				my $refCnt = 0;
				my $altCnt = 0;
				my $value_ref;
				my $value_alt;
				for(my $i = 0; $i < @{$FWD{$dupGroup}}; $i++){
					if($FWD{$dupGroup}[$i] eq $refAllele){
						$refCnt++;
						#print STDERR "\tmy index is $i, found ref allele: $FWD{$dupGroup}[$i] \n";
					}
					elsif( $FWD{$dupGroup}[$i] eq $altAllele ){
						$altCnt++;
						#print STDERR "\tmy index is $i, found alt allele: $FWD{$dupGroup}[$i] \n";
					}
				}
				#print STDERR "Finished with this dupGroup. My ref count is $refCnt, my alt count is $altCnt\n";
				
				#multiplying each ref/alt proportion by 100 (making sure they are rounded up and together sum up 100 counts)
				if($refCnt > $altCnt){
					$value_ref = int(( ( $refCnt / ($refCnt + $altCnt) ) * 100 ) + 0.5);
					$value_alt = 100 - $value_ref;
				}else{
					$value_alt = int(( ( $altCnt / ($refCnt + $altCnt) ) * 100 ) + 0.5);
					$value_ref = 100 - $value_alt;
				}
				my $sum = $value_ref + $value_alt;
				#print STDERR "My value alt is $value_alt, value ref is $value_ref, my sum is $sum\n";
				if( $sum != 100){
					print STDERR "The value of ref is $value_ref, of alt is $value_alt\n";
					print STDERR "Ref Count is $refCnt, alt count is $altCnt\n";
					die("Ref plus alt counts/proportions did not sum up to 100\n");
				}
				$weighted_ref = $weighted_ref + $value_ref;
				$weighted_alt = $weighted_alt + $value_alt;
			#	print STDERR "After summing, my weighted_ref is $weighted_ref, and my weighted_alt is $weighted_alt\n";
			}
		}	
		
		##SAME FOR FOR %REV
		my @keys_rev = keys %REV;
		my $num_keys_rev = scalar(@keys_rev);
		unless($num_keys_rev == 0){
			#print STDERR "\nNow going through hash with dup groups REV strand\n";
			foreach my $dupGroup (keys %REV) {
				#print STDERR "MY KEY IS $dupGroup\n";
				my $refCnt = 0;
				my $altCnt = 0;
				my $value_ref;
				my $value_alt;
				for(my $i = 0; $i < @{$REV{$dupGroup}}; $i++){
					if($REV{$dupGroup}[$i] eq $refAllele){
						$refCnt++;
						#print STDERR "\tmy index is $i, found ref allele: $REV{$dupGroup}[$i] \n";
					}
					elsif( $REV{$dupGroup}[$i] eq $altAllele ){
						$altCnt++;
						#print STDERR "\tmy index is $i, found alt allele: $REV{$dupGroup}[$i] \n";
					}
				}
				#print STDERR "Finished with this dupGroup. My ref count is $refCnt, my alt count is $altCnt\n";
				#multiplying each ref/alt proportion by 100 (making sure they are rounded up and together sum up 100 counts)
				if($refCnt > $altCnt){
					$value_ref = int(( ( $refCnt / ($refCnt + $altCnt) ) * 100 ) + 0.5);
					$value_alt = 100 - $value_ref;
				}else{
					$value_alt = int(( ( $altCnt / ($refCnt + $altCnt) ) * 100 ) + 0.5);
					$value_ref = 100 - $value_alt;
				}
				if(($value_ref + $value_alt) != 100){
					print STDERR "The value of ref is $value_ref, of alt is $value_alt\n";
					print STDERR "Ref Count is $refCnt, alt count is $altCnt\n";
					die("Ref plus alt counts/proportions did not sum up to 101\n");
				}
				$weighted_ref = $weighted_ref + $value_ref;
				$weighted_alt = $weighted_alt + $value_alt;
			#	print STDERR "After summing, my weighted_ref is $weighted_ref, and my weighted_alt is $weighted_alt\n";
			}
		}
		
		my $weighted_ref_ratio = $weighted_ref / ($weighted_ref + $weighted_alt);
		#rm temp files?
		#system("rm *temp*");
		#print "INDIVIDUAL\tRSID\tBOTH_ALLELES_SEEN\tREF_ALLELE\tALT_ALLELE\tDUP_WEIGHTED_REF_COUNT\tDUP_WEIGHTED_ALT_COUNT\tDUP_WEIGHTED_REF_RATIO\n";
		print "$indiv\t$snpID\t$bas\t$refAllele\t$altAllele\t$weighted_ref\t$weighted_alt\t$weighted_ref_ratio\n";
		
		#print STDERR "\nPRINTED DATA FOR THIS VARIANT!!\n";
		#print STDERR "Number of reads with final pos defined: $cntAlleleOK\n";
		#print STDERR "Number of reads with allele not in ref or alt: $cntBadNt\n";
		#print STDERR "Number of reads with allele and base qual OK is: $cntReadOK\n";
		#print STDERR "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";	
	}
}
close(ASE);

print STDERR "Done!\n\n";

my @snpsNoAlt = keys %noALT;
my @snpsNotCalled = keys %notCALLED;
my $numAltNotFound = scalar(@snpsNoAlt);
my $numAlleleNotCalled = scalar ( @snpsNotCalled );

print STDERR "Found $numAltNotFound cases for which the alternative allele was not reported.\nHere are some examples:\n";
print STDERR "$snpsNoAlt[0]\t$snpsNoAlt[1]\t$snpsNoAlt[2]\t$snpsNoAlt[3]\t$snpsNoAlt[4]\n";

print STDERR "Found $numAlleleNotCalled cases for which at least one of the alleles was not in the original variant call.\nHere are some examples:\n";
print STDERR "$snpsNotCalled[0]\t$snpsNotCalled[1]\t$snpsNotCalled[2]\t$snpsNotCalled[3]\t$snpsNotCalled[4]\n";

$datestring = localtime();
print STDERR "Local date and time $datestring\n";
