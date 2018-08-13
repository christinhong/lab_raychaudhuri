#!/usr/bin/env bash

# Kamil Slowikowski
# June  9, 2014
#
#Modif Maria Nov, 2014

fq="(\.fq|\.fastq)$"
if [[ $# != 2 || ! "$1" =~ $fq ]]
then
    echo "Received: $0 $@"
    echo "Usage: $0 SC10054_A01.end1.fq regExprFile"
    echo
    exit 1
fi
#[[ ! -e "$1" ]] && echo "FASTQ file not found $1" && exit 1

regex=$2
#/medpop/srlab/cd4timelinepilot/rnaseq/141125_PR1504/input/20141125-0hr-NT-20141014-IDX2-PR1504_S1_R1.fastq
[[ $1 =~ $regex ]]
basedir="${BASH_REMATCH[1]}"
library="${BASH_REMATCH[2]}"

echo $basedir
echo $library

out=$basedir/mapping/$library/starPass1
[[ ! -d $out ]] && mkdir -p $out

# STAR generates files in the current directory.
cd $out

opt=(
--genomeDir /medpop/rnaseq/access/star_genome3
--readFilesIn $1
--runThreadN 4
--genomeLoad LoadAndKeep
)

if [[ "$DRY_RUN" == "" ]]
then
    STAR ${opt[*]} &> log.txt
    # Convert SAM to BAM.
    sam=Aligned.out.sam
    bam=Aligned.out.bam
	[[ -e $sam ]] && samtools view -Sb -o $bam $sam && rm -f $sam
else
    echo "STAR ${opt[*]} &> log.txt"
    echo "samtools view -Sb -o Aligned.out.bam Aligned.out.sam"
    echo "rm -f Aligned.out.sam"
fi


