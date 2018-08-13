#!/bin/bash

# Modified by Christin Hong on 2015-10-19 from:
#@!/usr/bin/env bash
# run_subread-align.sh
# Kamil Slowikowski / Maria Gutierrez-Arcelus 
# June  9, 2014
#
# Run the subread aligner to map reads to a personalized genome (called variants are masked to reduce mapping biases due to individual SNPs).


#####

# Input and output identifiers
regex="(SR.*).fastq"
[[ $1 =~ $regex ]]
out=ASE02f13_3subreadAligned_"${BASH_REMATCH[1]}".bam
log=${out}_log.txt


# Map libraries with subread align to personalized masked genome
${subread}/subread-align \
     -i hg19-masked-subread \
     -r ${1} \
     -o ${out} \
     -u \
     -Q \
     --BAMoutput \
     -d 0 \
     -D 1000000 \
     -T ${intCPUs} \
     &> ${log}



