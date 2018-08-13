#!/usr/bin/env bash
# run_rnaseq.sh
# Kamil Slowikowski <kslowikowski@fas.harvard.edu>
#
# Execute the following steps:
#
# 1. Map RNA-seq reads from 1 (single-end) or 2 (paired-end) FASTQ files.
# 2. Count genes using the GENCODE v19 reference.
# 3. Compute metrics on the BAM file using Picard tools.
#
# Usage:
#
#   /PHShome/ks38/work/ngs/scripts/run_rnaseq.sh DIR FQ1 FQ2
#
#   DIR     A folder called "subjunc" will be written here.
#   FQ1     Read 1 of a paired-end RNA-seq experiment.
#   FQ2     Read 2 of a paired-end RNA-seq experiment.
#
# A BAM file will be written to "DIR/subjunc/FILE.bam"
#
# The name of the BAM file depends on the name of the FQ1 file:
#
#   Input FASTQ:    /path/to/C123_A1-1_B1_ASDF_R1.fastq.gz
#   Output BAM:     C123_A1-1_B1_ASDF_R1.bam
#

# Global variables -----------------------------------------------------------

THREADS=4

##############################################################################
# Ensembl 75 + Gencode 19 ----------------------------------------------------

SUBREAD_INDEX="/data/srlab/external-data/Ensembl/pub/release-75/fasta/homo_sapiens/dna/chr/repbase/subread/Homo_sapiens.GRCh37.75.dna.primary_assembly.repbase"
GTF="/data/srlab/external-data/GENCODE/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf"

##############################################################################
# UCSC hg19 + Gencode 19 -----------------------------------------------------

# SUBREAD_INDEX="/data/srlab/external-data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/SubreadIndex/hg19-subread"
# GTF="/data/srlab/external-data/GENCODE/Gencode_human/release_19/gencode.v19.annotation.gtf"
# PICARDMETRICS_CONF="/PHShome/ks38/picardmetrics.conf"

# GTF_INTRONS="/medpop/rnaseq/genome/gencode.v19.annotation.FILT2.0.introns.gtf"
# [[ ! -s "$GTF_INTRONS" ]] && echo "GTF not found: $GTF_INTRONS" && exit 1

##############################################################################
# Gencode GRCh37.p13 with alternative haplotypes and unplaced contigs --------

# SUBREAD_INDEX="/data/srlab/external-data/GENCODE/Gencode_human/release_19/subread/GRCh37.p13"
# GTF="/data/srlab/external-data/GENCODE/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf"
# PICARDMETRICS_CONF="/PHShome/ks38/picardmetrics-GRCh37.p13.conf"

# Main function --------------------------------------------------------------

main() {

  out="$1"
  fq1="$2"
  fq2="$3"

  if ! mkdir -p "$out"; then
    echo "ERROR: Failed to create '$out'" 1>&2
    exit 1
  fi

  if ! is_fastq "$2"; then
    echo "ERROR: Not a FASTQ file: '$2'" 1>&2
    exit 1
  fi

  if ! is_fastq "$3"; then
    echo "ERROR: Not a FASTQ file: '$3'" 1>&2
    exit 1
  fi

  echo "$(date) -- START -- $fq1"

  # Map reads ----------------------------------------------------------------

  sample=$(basename ${fq1%%.*})

  # Create and then enter the subread/ directory.
  progs=(
    subjunc_ensembl75_repbase
  )
    #star
    #subjunc
    #subread

  for prog in ${progs[*]}; do

    if ! mkdir -p "$out/$prog"; then
      echo "ERROR: Failed to create '$out'" 1>&2
      exit 1
    fi

    if ! cd "$out/$prog"; then
      echo "ERROR: cd '$out/$prog'" 1>&2
      exit 1
    fi

    local bam_mapped="${sample}.bam"

    # Check if the BAM exists and also check if it is incomplete (corrupt).
    if [[ "$fq1" -nt "$bam_mapped" ]] || bam_incomplete "$bam_mapped"; then
      subjunc_align "$fq1" "$fq2" "$bam_mapped"
    fi

    # opt=(
    #   -f "$PICARDMETRICS_CONF"
    #   -o picardmetrics
    #   -r
    #   "$bam_mapped"
    # )
    # log=${sample}.picardmetrics.log
    # echo "$(date) -- Running picardmetrics"
    # if ! picardmetrics run ${opt[*]} &> $log; then
    #   echo "ERROR: picardmetrics failed. See '$log'"
    #   cat $log >&2
    #   exit 1
    # fi

    # Count reads mapped to genes with featureCounts ---------------------------

    # local counts=${sample}.featureCounts-gene.tsv
    # if [[ "$bam_mapped" -nt "$counts" || ! -s "$counts" ]]; then
    #   count_genes_by_gene "$bam_mapped" "$counts"
    # fi

    local counts=${sample}.featureCounts-exon.tsv
    if [[ "$bam_mapped" -nt "$counts" || ! -s "$counts" ]]; then
      count_genes_by_exon "$bam_mapped" "$counts"
    fi
  done

  echo "$(date) -- DONE  -- $fq1"
}

# Extra functions ------------------------------------------------------------

is_bam() {
  re="\.bam$"
  if [[ ! -s "$1" ]]; then
    echo "ERROR: File not found: '$1'" >&2
    return 1
  fi
  [[ "$1" =~ $re ]]
}

is_fastq() {
  re="\.(fq|fq.gz|fastq|fastq.gz)$"
  if [[ ! -s "$1" ]]; then
    echo "ERROR: File not found: '$1'" >&2
    return 1
  fi
  [[ "$1" =~ $re ]]
}

bam_incomplete() {
  if [[ ! -s "$1" ]]; then
    return 0
  fi
  samtools view -H "$1" 2>&1 | grep -q 'EOF marker is absent'
  return "$?"
}

function fastq_readlength {
  zcat $1 | head -n4 | awk 'NR % 4 == 2 {print length($1)}'
}

subjunc_align() {
  local fq1="$1" fq2="$2" bam="$3"
  # I think 20 subreads might be causing the CIGAR error ... not sure.
  local common_opt=(
      --subreads 10
      --index $SUBREAD_INDEX
      --BAMoutput
      --output "$bam"
      --multi 16
      --quality
      --threads $THREADS
      --maxMismatches 10
      --allJunctions
  )
  local opt
  # If a second FASTQ file is present, assume paired-end sequencing.
  if is_fastq "$fq1" && is_fastq "$fq2"; then
    opt=(
      "${common_opt[*]}"
      --gzFASTQinput
      --read "$fq1"
      --read2 "$fq2"
      --order fr
      --mindist 0
      --maxdist 1000000
    )
  # If only the first FASTQ file is present, assume single-end sequencing.
  elif is_fastq "$fq1"; then
    opt=(
      "${common_opt[*]}"
      --gzFASTQinput
      --read "$fq1"
    )
  # If a BAM file is present, assume paired-end sequencing.
  elif is_bam "$fq1"; then
    opt=(
      "${common_opt[*]}"
      --BAMinput
      --read "$fq1"
      --order fr
      --mindist 0
      --maxdist 1000000
    )
  else
    echo "ERROR: Bad arguments: fq1='$fq1' fq2='$fq2' bam='$bam'" 1>&2
    exit 1
  fi
  # Run subjunc.
  echo "$(date) -- Running subjunc"
  log=$(basename "$bam" .bam).subjunc-align.log
  if ! ( subjunc ${opt[*]} ) &> $log; then
    echo "ERROR: subjunc failed. See '$log'" >&2
    cat $log >&2
    exit 1
  fi
}

# featureCounts --------------------------------------------------------------

count_genes_by_exon() {
  local bam="$1" counts="$2"
  # See explanations of each option below.
  local opt=(
    -T $THREADS
    -s 0
    -Q 40
    -p
    -F GTF
    -t exon
    -g gene_id
    -a "$GTF"
    -o "$counts"
    "$bam"
  )
  echo "$(date) -- Running featureCounts on exons"
  # featureCounts creates temporary files in the current directory.
  log=${counts%.*}.log
  if ! ( featureCounts ${opt[*]} ) &> $log; then
    echo "ERROR: count_exons failed. See '$log'" >&2
    cat $log >&2
    exit 1
  fi
}

main "$@"
