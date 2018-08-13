#!/bin/bash

# Christin M. Hong
# Last modified: 2015-10
# Lab of Soumya Raychaudhuri, Harvard Medical School

# Assuming that some PCR "duplicates" are actually true reads.  Normalizing duplicates at heterozygous SNP sites to 100 (using 100 for convenience) in order to account for the proportion of duplicate reads at a site that come from the reference allele vs. the proportion from the alternate allele.


#####

# VCF file with filtered variants with IDs
varsFilt=ASE02f10_${pathInd}_hets_allFilts.ids.vcf


# Input and output identifiers
regexASE="ASE02f14_5SamASE_${pathInd}_(SR.*).ASE.txt"
[[ ${1} =~ ${regexASE} ]]

input2=ASE02f14_1SamSort_${pathInd}_"${BASH_REMATCH[1]}".sorted.bam
out1=ASE02f15_dupGroupsScaledTo100_${pathInd}_"${BASH_REMATCH[1]}".txt


#####

# Scale reads overlapping a het SNP site with a PCR "duplicate" block to 100 counts per event, with the entire PCR block counting as 1 event, and every other unique read overlapping that site counting as 1 event.
perl ${pathSubscriptsMaria}/get_weightedDupCounts_by100.pl \
     ${1} \
     ${input2} \
     > ${out1}


