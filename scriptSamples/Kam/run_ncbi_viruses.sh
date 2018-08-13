#!/usr/bin/env bash
# run_rnaseq.sh
# Kamil Slowikowski <kslowikowski@fas.harvard.edu>

# Global variables -----------------------------------------------------------

THREADS=4

##############################################################################
# NCBI Viruses ---------------------------------------------------------------

#NCBI_VIRUSES="/data/srlab/external-data/NCBI/viruses/subread/ncbi-viruses-refseq"

##############################################################################
# Ensembl 75 + NCBI Viruses --------------------------------------------------

NCBI_VIRUSES="/data/srlab/external-data/Ensembl/pub/release-75/fasta/homo_sapiens/dna/chr/viruses/subread/Homo_sapiens.GRCh37.75.dna.primary_assembly.viruses"

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
  prog=subjunc_ensembl75_ncbi_viruses

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
  subjunc_align "$fq1" "$fq2" "$bam_mapped"

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
  local common_opt=(
      --subreads 10
      --index $NCBI_VIRUSES
      --BAMoutput
      --output "$bam"
      --multi 8
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

main "$@"
