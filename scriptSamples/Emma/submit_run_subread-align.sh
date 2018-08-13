#!/bin/bash
#BSUB -J run_subread_test
#BSUB -o cluster_output/run_subread_test-%J.out
#BSUB -e cluster_output/run_subread_test-%J.err

sh ~/bin/run_subread-align.sh US-1425956_1.fq.gz US-1425956_2.fq.gz
