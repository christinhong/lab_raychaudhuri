#!/usr/bin/env bash
# drop_long_skipped_reads.sh
# Kamil Slowikowski
# July 19, 2015
#
# Drop reads from a BAM file if they have long skipped regions.
#
# Example:
#
#   ./drop_long_skipped_reads.sh FILE.bam 1000000
#
# This will produce two files:
#
#   FILE.no_long_skipped.bam
#
#       The clean BAM file, omitting reads with long skipped regions.
#
#   FILE.long_skipped.bam
#
#       The subset of reads that contain long skipped regions.

if [[ $# -ne 2 ]]; then
  echo "usage: $0 BAM N"
  exit 1
fi

IN="$1"
N="$2"

if [[ ! -f "$IN" ]]; then
  echo "'$IN' must be a BAM file."
  exit 1
fi

if [[ ! "$N" =~ ^[0-9]+$ ]]; then
  echo "N must be an integer."
  exit 1
fi

OUT1="${IN%.bam}.no_long_skipped.bam"
OUT2="${IN%.bam}.long_skipped.bam"

samtools view -h $IN \
  | perl -lane '
  if (/^\@/) {
    print STDOUT;
    print STDERR;
  } elsif ($F[5] =~ /([0-9]{'${#N}',})N/) {
    if ($1 <= '$N') {
      print STDOUT;
    } else {
      print STDERR;
    }
  } else {
    print STDOUT;
  }
  ' \
  1> >(samtools view -bS - > "$OUT1") \
  2> >(samtools view -bS - > "$OUT2")

