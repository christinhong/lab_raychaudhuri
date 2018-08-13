#!/bin/bash

# Christin M. Hong
# Last modified: 2015-09
# Soumya Raychaudhuri Lab, Harvard Medical School

# 2nd pass of STAR mapping of reads to reference genome with splice junctions collected from all samples during 1st pass


#####

# Tagging output with sample identifier.
regexFastq="(SR.*).sra.fastq"
[[ ${1} =~ ${regexFastq} ]]
outBAM=ASE02f02_starP2_"${BASH_REMATCH[1]}"_

echo ${outBAM}

# Store SJ files in a bash array.
readarray -t sj < ASE02_aa_sj.txt

# Run STAR mapping with SJ files.
${STAR} \
     --readFilesIn ${1} \
     --genomeDir ${pathStarRef} \
     --outFileNamePrefix ${outBAM} \
     --outSAMtype BAM Unsorted \
     --genomeLoad NoSharedMemory \
     --sjdbFileChrStartEnd ${sj[*]}


# NOTES
     # --genomeLoad LoadAndKeep is incompatible with on-the-fly junction insertion (sjdbFileChrStartEnd) - need to run with --genomeLoad NoSharedMemory.

     # My understanding is that:
          # "${myarray[@]}" leads to each element of the array being treated as if they were separated by newlines, and 
          # "${myarray[*]}" results in the elements of the array being treated as if they were separated by spaces (or whatever the first character of IFS is).  See http://stackoverflow.com/a/3355375
     # Like subjunc "" within command aren't supported.  When I tried for ${sj[*]}, I think it made the command try to find the full expanded array as one filename instead of files separated by whitespace.

