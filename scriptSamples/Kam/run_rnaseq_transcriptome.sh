#!/usr/bin/env bash
# run_rnaseq_transcriptome.sh
# Kamil Slowikowski
# June 18, 2015
#
# Usage:
#
#   /PHShome/ks38/work/ngs/scripts/run_rnaseq.sh DIR FQ1 FQ2
#
#   DIR     A folder called "subjunc" will be written here.
#   FQ1     Read 1 of a paired-end RNA-seq experiment.
#   FQ2     Read 2 of a paired-end RNA-seq experiment.
#

##############################################################################
# Ensembl 75 Transcriptome ---------------------------------------------------

EXPRESS_TARGETS="/data/srlab/external-data/Ensembl/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa"
SUBREAD_INDEX="/data/srlab/external-data/Ensembl/pub/release-75/fasta/homo_sapiens/cdna/subread/Homo_sapiens.GRCh37.75.cdna.all"
BOWTIE2_INDEX="/data/srlab/external-data/Ensembl/pub/release-75/fasta/homo_sapiens/cdna/bowtie2/Homo_sapiens.GRCh37.75.cdna.all"

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

  sample=$(basename ${fq1%%.*})
  out="$out/$sample"
  log="$out/express.log"

  if ! mkdir -p "$out"; then
    echo "ERROR: Failed to create '$out'" 1>&2
    exit 1
  fi

  (echo "$(date) -- START -- $sample"

  # Map reads ----------------------------------------------------------------

    #--upto 100
  bowtie_opt=(
    --threads 4
    --sensitive-local
    -a
    -X 800
    -x $BOWTIE2_INDEX
    -1 $fq1
    -2 $fq2
  )

  express_opt=(
    --no-update-check
    --logtostderr
    --output-dir $out
    $EXPRESS_TARGETS
  )

  bowtie2 ${bowtie_opt[*]} | express ${express_opt[*]}

  echo "$(date) -- DONE  -- $sample") &> $log
}

# Functions ------------------------------------------------------------------

is_bam() {
  re="\.bam$"
  [[ -s "$1" && "$1" =~ $re ]]
}

is_fastq() {
  re="\.(fq|fq.gz|fastq|fastq.gz)$"
  [[ -s "$1" && "$1" =~ $re ]]
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
      --subreads 20
      --index $SUBREAD_INDEX
      --BAMoutput
      --output "$bam"
      --multi 16
      --quality
      --threads 4
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
    echo "ERROR: subjunc failed" 1>&2
    exit 1
  fi
}

main "$@"
