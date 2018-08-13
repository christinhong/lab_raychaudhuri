#!/bin/bash

# Christin M. Hong
# Last modified: 2015-09
# Lab of Soumya Raychaudhuri, Harvard Medical School


#####

# Input and output identifiers
regex="(.*)output.all_chrs.(.*).vcf"
[[ ${1} =~ $regex ]]
base="${BASH_REMATCH[2]}"


nice -n20 java -jar ${GATK} -T VariantEval \
     -R ${fileRef} \
     -o ASE02f11_${pathInd}_eval.${base}.gatkreport_dbsnp \
     --dbsnp ${fileDbsnp} \
     --eval:set1 ${1}

nice -n20 java -jar ${GATK} -T VariantEval \
     -R ${fileRef} \
     -o ASE02f11_${pathInd}_eval.$base.gatkreport_1KGP \
     --dbsnp ${file1000G} \
     --eval:set1 ${1}
    

