#!/usr/bin/env bash
# collate_express.sh
# Kamil Slowikowski
# July 21, 2015
#
# Collate multiple expression results files produced by kallisto.

files=($@)

(
  # The first file provides the header.
  first=${files[0]}
  sample=$(basename $(dirname $first))
  perl -pe 'if ($.==1) { s/^/SAMPLE\t/ } else { s/^/'$sample'\t/ }' < $first

  # Skip the header in the remaining files.
  for f in ${files[*]:1}; do
	sample=$(basename $(dirname $f))
	tail -n+2 $f | perl -pe 's/^/'$sample'\t/'
  done
) | pigz > abundance.txt.gz
