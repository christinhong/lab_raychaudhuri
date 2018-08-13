#!/usr/bin/env bash
# run_kallisto.sh
# Kamil Slowikowski
# July 21, 2015

set -eou pipefail

KALLISTO_INDEX=/data/srlab/external-data/Ensembl/pub/release-81/fasta/homo_sapiens/cdna/kallisto-repbase/Homo_sapiens.GRCh38.cdna.all

paths=(
	/data/srlab/DeepakRao/Rao_SSF-700/H7F3VBGXX
	/data/srlab/DeepakRao/Rao_SSF-876/HCHT7BGXX
	/data/srlab/MichaelBrenner/Brenner_SSF-686/H5WYVBGXX
)
	# /data/srlab/AMP/SN0060355/fastq
script="/PHShome/ks38/work/ngs/scripts/run_rnaseq.sh"

for path in ${paths[*]}; do

  sleep 2

  mkdir -p $path/jobs
  cd $path/jobs

  job="kallisto_ensembl81_repbase_$(date +%Y-%m-%d_%H-%M)"

  ls $path/*/*.fastq.gz | sort | while read fq1; read fq2; do
	  # Ensure we have forward and reverse reads from the same sample.
	  sample1=$(basename ${fq1%%.*})
	  sample2=$(basename ${fq2%%.*})

	  if [[ "$sample1" != "$sample2" ]]; then
		  echo "ERROR: Mismatched samples '$sample1' and '$sample2'"
		  echo "       '$fq1'"
		  echo "       '$fq2'"
		  exit 1
	  fi

	  out=$path/kallisto_ensembl81_repbase/$sample1
	  mkdir -p $out

	  log=$out/kallisto-quant.log

	  opt=(
		quant
		-i $KALLISTO_INDEX
		-o $out
		$fq1
		$fq2
	  )

	  echo "(time kallisto ${opt[*]}) &> $log"
  done | \
  asub \
	  -R "rusage[mem=4000]" \
	  -q normal -j $job
  
done
