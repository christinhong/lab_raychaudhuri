#!/bin/bash

# Christin M. Hong
# Last modified: 2015-09
# Lab of Soumya Raychaudhuri, Harvard Medical School


#####

# Merge all files from the individual for the input chromosome

cat ASE02_aa_bqsrBams.txt | {
    read firstbam
    samtools view -h "$firstbam" ${1}
    while IFS= read -r bam; do
        samtools view "$bam" ${1}
    done
} | samtools view -ubS - | samtools sort - ASE02f04_${pathInd}_${1}_merged


# Notes
     # Using this instead of samtools merge because samtools merge kept giving me segmentation fault errors.  (Maybe due to header incompatibilities?)
     # Checked for success using samtools view -c (count) and confirming that the merged file had the correct number of sequences (sum of individual files).  (Thanks to Maria for the tip!)
