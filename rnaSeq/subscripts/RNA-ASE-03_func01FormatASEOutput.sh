#!/bin/bash

# Christin M. Hong
# Last modified: 2015-10
# Lab of Soumya Raychaudhuri, Harvard Medical School

# Cleaning up allelic-specific expression output for downstream analysis in R.


#####

# Input and output identifiers
regex="ASE02f14_5SamASE_(.*).sra.ASE.txt"
[[ ${1} =~ ${regex} ]]

out1=ASE03f01_trimmed_"${BASH_REMATCH[1]}".txt


#####

cut -f1-9 < ${1}  > ${out1}


