#!/usr/bin/env bash
# run_subread-align.sh
# Kamil Slowikowski
# October 6, 2014
#
# Run the subread aligner.
#
# Usage:
#For Broad
# parallel --env _ -S 10/piglet,8/tigger --colsep '\t' /home/unix/slowikow/work/rnaseq/scripts/run_subread-align.sh {1} {2} :::: <(head -n1 /medpop/rnaseq/tcells/fastqs.txt)
#
# parallel --env _ -S 10/piglet,10/tigger --colsep '\t' /home/unix/slowikow/work/rnaseq/scripts/run_subread-align.sh {1} {2} :::: /medpop/rnaseq/tcells/fastqs.txt

# $RAM
# Human genome ref file and copying to memory - less reads from disk = faster
#source /home/unix/slowikow/work/rnaseq/scripts/copy_genome_to_mem.sh

#prog='/home/unix/slowikow/src/subread-1.4.5-p1-source/bin/subread-align'

fq="(\.fq|\.fastq|\.fq.gz|\.fastq.gz)$"
#two arguments (forward and reverse), both ending in same extension
if [[ $# != 2 || ! "$1" =~ $fq || ! "$2" =~ $fq ]]
then
    echo "Received: $0 $@"
    echo "Usage: $0 file.end1.fq[.gz] file.end2.fq[.gz]"
    echo
    exit 1
fi
#do the files exist
[[ ! -e "$1" ]] && echo "FASTQ file not found $1" && exit 1
[[ ! -e "$2" ]] && echo "FASTQ file not found $2" && exit 1

out=/data/srlab/pfizer/Subread
#create output files in same directory as fasq
[[ ! -d $out ]] && mkdir -p $out
log=$out/subread-align-log-$(basename ${1%_*}).txt
#text file containing status
opt=(
    --index /data/srlab/external-data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/SubreadIndex/hg19-subread
    --read $1
    --read2 $2
    --gzFASTQinput
    --subreads 20
    --BAMoutput
    --output $out/$(basename ${1%_*}).bam
    --order fr
    --mindist 0
    --maxdist 1000000
    --multi 2
    --quality
    --reportFusions
    --threads 4
)

(time subread-align ${opt[*]}) &> $log
