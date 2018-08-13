#!/bin/bash

#Variant Caller
#runs star pass 2 mapping
#runs pre-processing for GATK
#runs hapCaller by merging files by chromosome


#parameters to change
workingDir=/home/unix/mgutierr/work/rnaseq/cd4timelinepilot/tofiPub
mappingDir=/medpop/srlab/mgutierr/tofiPub/mapping
suffixSamples=SRR
listFastq=list_fastq.txt#assumes it is in working dir
numJobs=3
regExpr=regExpr.txt
#check ref star genome is still OK for run_star_pass1.sh

#shortcuts
parallel=/home/unix/slowikow/src/parallel-20140422/src/parallel

#STAR
#$parallel --xapply echo :::: list_fastq.txt $regExpr

#Align reads with STAR for each sample
$parallel -j$numJobs --xapply --eta 'sh run_star_pass1.sh' :::: $listFastq $regExpr

#Merge the exon-exon junctions found for each sample into one file
nice -n20 perl ~/work/rnaseq/cd4timelinepilot/merge_SJ_pass1.pl $mappingDir/$suffixSamples*/starPass1/SJ.out.tab > $mappingDir/SJ.out.merged.tab

#Prepare new reference genome for star using exon junctions
sh run_star-prepare-reference-mergedSJs.sh $mappingDir $mappingDir/SJ.out.merged.tab

#Run STAR again for each sample using new ref genome	
$parallel -j$numJobs --xapply --eta 'sh run_star_pass2.sh' :::: $listFastq $regExpr

#Run GATK data pre-processing
ls -1 $mappingDir/$suffixSamples*/starPass2/*.bam > list_bams_starPass2.txt 

parallel -j3 --eta 'sh run_gatk_markDups.sh' :::: list_bams_starPass2.txt
#/home/unix/slowikow/src/parallel-20140422/src/parallel echo ::: A B C > abcFile

