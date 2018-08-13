!/bin/bash -l

# Christin M. Hong
# Lab PI: Soumya Raychaudhuri, Harvard Medical School
# Last updated: 2017-05


####

source /broad/software/scripts/useuse

reuse Bcftools

bcftools index -f ${G1K}/ALL.chr${1}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
