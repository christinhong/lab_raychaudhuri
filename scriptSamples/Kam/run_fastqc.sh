#!/usr/bin/env bash
# run_fastqc.sh
# Kamil Slowikowski
# April 27, 2015
#
# Submit a job to LSF to run FASTQC on a set of files.

if [[ "$1" == "" || "$2" == "" ]]; then
	echo "usage: $0 outdir file.fq.gz [file.fq.gz...]"
	exit 1
fi

OUT="$1"; shift
mkdir -p "$OUT"

QUEUE="short"
# QUEUE="medium"
JOB="fastqc_$(date +%Y-%m-%d_%H-%M)"

DIR=/PHShome/ks38/src/fastqc-0.11.3

opt=(
--extract
--nogroup
--contaminants $DIR/Configuration/contaminant_list.txt
--adapters $DIR/Configuration/adapter_list.txt
--threads 4
--out $OUT
)

for f in $@; do
	if [[ ! -f $f ]]; then
		echo "Not a file: '$f'"
		exit 1
	fi

	# Skip if we already have output.
	html="$OUT/$(basename $f .gz)_fastqc.html"
	if [[ ! -e $html ]]; then
		echo "fastqc ${opt[@]} $f"
	fi
done | \
asub -M 8 \
	-R "select[type==X86_64&&mem>800] rusage[mem=800]" \
	-n 4 -q $QUEUE -j $JOB

