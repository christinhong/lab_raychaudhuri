#!/bin/bash

# Christin M. Hong
# Last modified: 2015-09
# Lab of Soumya Raychaudhuri, Harvard Medical School

# Make a subdirectory for each individual and transfer in their FASTQ files


#####

sub="${pathWd}/all${1}"

mkdir $sub

while IFS= read -r line ; do mv "${line}" $sub ; done < ${1}
