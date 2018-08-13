#!/usr/bin/env bash
# run_subread-align.sh
# Kamil Slowikowski
# June  9, 2014
#
# Run the subread aligner.

# $RAM
#source /home/unix/mgutierr/work/rnaseq/SC/scripts/copy_genome_to_mem.sh

fq="(\.fastq)$"
if [[ $# != 1 || ! "$1" =~ $fq ]]
then
    echo "Received: $0 $@"
    echo "Usage: $0 SC10054_A01.end1.fq"
    echo
    exit 1
fi
[[ ! -e "$1" ]] && echo "FASTQ file not found $1" && exit 1

regex="(/medpop/.*)/(SRR.*)\.fastq"
#/medpop/srlab/cd4timelinepilot/rnaseq/141125_PR1504/input/20141125-0hr-NT-20141014-IDX2-PR1504_S1_R1.fastq

[[ $1 =~ $regex ]]
basedir="${BASH_REMATCH[1]}"
library="${BASH_REMATCH[2]}"

echo $basedir
echo $library

out=$basedir/mapping/$library/subread
[[ ! -d $out ]] && mkdir -p $out
log=$out/log.txt

opt=(
    -i /medpop/rnaseq/genome/subread/hg19-subread
    -u
    -Q
    -D 100000
    -r $1
    --BAMoutput
    -o $out/${library}.bam
    -T 4
)

if [[ "$DRY_RUN" == "" ]]
then
    nice -n20 /home/unix/slowikow/src/subread-1.4.5-p1-source/bin/subread-align ${opt[*]} &> $log
else
    echo "/home/unix/slowikow/src/subread-1.4.5-p1-source/bin/subread-align ${opt[*]} &> $log"
fi

#/home/unix/slowikow/work/rnaseq/scripts/run_subread-align.sh /medpop/rnaseq/tcells/G71162/SCC10055_A01/v1/SCC10055_A01.end1.fq /medpop/rnaseq/tcells/G71162/SCC10055_A01/v1/SCC10055_A01.end2.fq

# subread-align
# Version 1.4.5-p1
# 
# Usage:
# 
#  ./subread-align [options] -i <index_name> -r <input> -o <output>
# 
# Required arguments:
# 
#     -i --index     <index>  base name of the index.
# 
#     -r --read      <input>  name of the input file(FASTQ/FASTA format by default
#                             . See below for more supported formats). Both base-
#                             space and color-space read data are supported. For
#                             paired-end reads, this gives the first read file
#                             and the other read file should be specified using
#                             the -R option.
# 
# Optional general arguments:
# 
#     -o --output    <output> name of the output file(SAM format by default). If
#                             not provided, mapping results will be output to the
#                             standard output (stdout).
# 
#     -n --subreads  <int>    number of selected subreads, 10 by default.
# 
#     -m --minmatch  <int>    consensus threshold (minimal number of consensus
#                             subreads required) for reporting a hit. If paired-
#                             end read data are provided, this gives the consensus
#                             threshold for the read which receives more votes
#                             than the other read from the same pair. 3 by default
# 
#     -T --threads   <int>    number of threads, 1 by default.
# 
#     -I --indel     <int>    number of indels allowed, 5 by default. Indels of up
#                             to 200bp long can be detected.
# 
#     -B --multi     <int>    Specify the maximal number of equally-best mapping
#                             locations allowed to be reported for each read. 1
#                             by default. Allowed values are between 1 to 16
#                             (inclusive). 'NH' tag is used to indicate how many
#                             alignments are reported for the read and 'HI' tag
#                             is used for numbering the alignments reported for
#                             the same read, in the output. Note that -u option
#                             takes precedence over -B.
# 
#     -P --phred     <3:6>    the format of Phred scores in input files, '3' for
#                             phred+33 and '6' for phred+64. '3' by default.
# 
#     -u --unique             only uniquely mapped reads will be reported (reads
#                             mapped to multiple locations in the reference genome
#                             will not be reported). This option can be used
#                             together with option '-H' or '-Q'.
# 
#     -Q --quality            using mapping quality scores to break ties when more
#                             than one best mapping location is found.
# 
#     -H --hamming            using Hamming distance to break ties when more than
#                             one best mapping location is found.
# 
#     -b --color-convert      convert color-space read bases to base-space read
#                             bases in the mapping output. Note that the mapping
#                             itself will still be performed at color-space.
# 
#     -M --maxMismatches <int> Specify the maximum number of mis-matched bases
#                             allowed in the alignment. 10 by default. Mis-matches
#                             found in soft-clipped bases are not counted.
# 
#        --reportFusions      report discovered genomic fusion events such as
#                             chimeras. Discovered fusions will be saved to a file
#                             (*.fusions.txt). Detailed mapping results for fusion
#                             reads will be saved to the SAM/BAM output file as
#                             well. Secondary alignments of fusion reads will be
#                             saved to the following optional fields: CC(Chr),
#                             CP(Position), CG(CIGAR) and CT(strand). Note that
#                             each fusion read occupies only one row in the
#                             SAM/BAM output file.
# 
#        --trim5     <int>    trim off <int> number of bases from 5' end of each
#                             read. 0 by default.
# 
#        --trim3     <int>    trim off <int> number of bases from 3' end of each
#                             read. 0 by default.
# 
#        --rg-id     <string> specify the read group ID. If specified,the read
#                             group ID will be added to the read group header
#                             field and also to each read in the mapping output.
# 
#        --rg        <string> add a <tag:value> to the read group (RG) header in
#                             in the mapping output.
# 
#        --gzFASTQinput       specify that the input read data is in gzipped
#                             FASTQ/FASTA format.
# 
#        --SAMinput           specify that the input read data is in SAM format.
# 
#        --BAMinput           specify that the input read data is in BAM format.
# 
#        --BAMoutput          specify that mapping results are saved into a BAM
#                             format file.
# 
#        --DPGapOpen  <int>   a numeric value giving the penalty for opening a
#                             gap when using the Smith-Waterman dynamic
#                             programming algorithm to detect insertions and
#                             deletions. The Smith-Waterman algorithm is only
#                             applied for those reads which are found to contain
#                             insertions or deletions. -1 by default.
# 
#        --DPGapExt   <int>   a numeric value giving the penalty for extending the
#                             gap, used by the Smith-Waterman algorithm. 0 by
#                             default.
# 
#        --DPMismatch <int>   a numeric value giving the penalty for mismatches,
#                             used by the Smith-Waterman algorithm. 0 by default.
# 
#        --DPMatch    <int>   a numeric value giving the score for matches used by
#                             the Smith-Waterman algorithm. 2 by default.
# 
#     -v                      output version of the program.
# 
# 
# Optional arguments for paired-end reads:
# 
#     -R --read2     <input>  name of the second input file. The program will then
#                             be switched to the paired-end read mapping mode.
# 
#     -p --minmatch2 <int>    consensus threshold for the read which receives less
#                             votes than the other read from the same pair, 1 by
#                             default.
# 
#     -d --mindist   <int>    minimum fragment/template length, 50bp by default.
# 
#     -D --maxdist   <int>    maximum fragment/template length, 600bp by default.
# 
#     -S --order     <ff:fr:rf> orientation of the two read from the same pair,
#                             'fr' by default.
# 
# 
# For more information about these arguments, please refer to the User Manual.
