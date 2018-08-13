#!/bin/bash -l

# Christin M. Hong
# Lab PI: Soumya Raychaudhuri, Harvard Medical School
# Last updated: 2017-05

# Note: Parallel on Broad doesn't seem to like line breaks, even with \ preceeding them.

#####

source /broad/software/scripts/useuse

reuse Bcftools

bcftools view -O v ${G1K}/ALL.chr${1}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -R ${proj}/analyses/SNPs-for-subsetting.txt > ${proj}/analyses/temp-shared-G1K-chr${1}.vcf

