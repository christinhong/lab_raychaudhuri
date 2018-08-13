#!/bin/bash

# Christin M. Hong
# Last modified: 2015-09
# Lab of Soumya Raychaudhuri, Harvard Medical School

# 1st pass of STAR mapping of reads to reference genome


#####

# Tagging output with sample identifier.
regexFastq="(SR.*).sra.fastq"
[[ ${1} =~ ${regexFastq} ]]
outBAM=ASE02f01_starP1_"${BASH_REMATCH[1]}"_

echo ${outBAM}

${STAR} \
     --readFilesIn ${1} \
     --genomeDir ${pathStarRef} \
     --genomeLoad LoadAndKeep \
     --outFileNamePrefix ${outBAM} \
     --outSAMtype BAM Unsorted


# Notes
     # Decided against specifying number of threads for multithreading since already running with GNU parallel.
     # The Partners cluster only allows creation of new directories if the full path is specified.  Trying "mkdir" without coding the full path from within a scriptS results in "mkdir: cannot create directory '': Permission denied".
