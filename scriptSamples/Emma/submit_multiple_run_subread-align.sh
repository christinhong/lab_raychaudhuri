#!/usr/bin/env bash

while read f1
read f2
do
	if [ $(basename "${f1%_*}") == $(basename "${f2%_*}") ]
	then	
	bsub -q normal -J "run_subread_test" -o "/PHShome/ed686/cluster_output/run_subread_test-%J.out" -e "/PHShome/ed686/cluster_output/run_subread_test-%J.err" "sh /PHShome/ed686/bin/run_subread-align.sh $f1 $f2"
	fi
done < /PHShome/ed686/bin/fastq-files_repeat.txt
