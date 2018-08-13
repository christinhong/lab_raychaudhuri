#!/usr/bin/env bash
# collate_express.sh
# Kamil Slowikowski
# July 21, 2015
#
# Collate multiple expression results files produced by eXpress.

express_files="$@"

for f in ${express_files[*]}; do
  sample=$(basename $(dirname $f))
  perl -pe 'if ($.==1) { s/^/SAMPLE\t/ } else { s/^/'$sample'\t/ }' < $f
done | pigz > results.xprs.gz
