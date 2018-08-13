#!/bin/bash

# Christin M. Hong
# Last modified: 2015-09
# Lab of Soumya Raychaudhuri, Harvard Medical School

# Picard and GATK data processing pipeline to prep for Haplotype Caller


#####

# Input and output identifiers
regexBam="ASE02f02_starP2_(SR.*)_Aligned.out.bam"
[[ ${1} =~ ${regexBam} ]]

out1=ASE02f03_1PicRG_"${BASH_REMATCH[1]}"
out2=ASE02f03_2PicMD_"${BASH_REMATCH[1]}"
out3=ASE02f03_3PicReO_"${BASH_REMATCH[1]}"
out4=ASE02f03_4GatkSplit_"${BASH_REMATCH[1]}"
out5=ASE02f03_5GatkRecal_"${BASH_REMATCH[1]}"
out6=ASE02f03_6GatkBQSR_"${BASH_REMATCH[1]}"

echo "Checking BASH_REMATCH in first and last output variables: ${out1} ${out6}"


#####

# Picard: Add read groups, sort, mark duplicates, create index, and reorder to match hg19 reference DICT file.
java -jar ${picard} AddOrReplaceReadGroups \
     I=${1} \
     O=${out1}.bam \
     SORT_ORDER=coordinate \
     RGID=readGroupID \
     RGLB=libraryID \
     RGPL=illumina \
     RGPU=HiSeq \
     RGSM=sampleID \
     CREATE_INDEX=true

java -jar ${picard} MarkDuplicates \
     I=${out1}.bam \
     O=${out2}.bam \
     CREATE_INDEX=true \
     VALIDATION_STRINGENCY=SILENT \
     M=ASE02f03_2PicMD_outMetrics.txt

java -jar ${picard} ReorderSam \
     I=${out2}.bam \
     R=${fileRef} \
     O=${out3}.bam \
     CREATE_INDEX=true


# GATK: Split'N'Trim and reassign mapping qualities (to keep unknown quality="255" reads from STAR mapping)
java -jar ${GATK} -T SplitNCigarReads \
     -R ${fileRef} \
     -I ${out3}.bam \
     -o ${out4}.bam \
     --read_filter ReassignOneMappingQuality \
     -RMQF 255 \
     -RMQT 60 \
     --unsafe ALLOW_N_CIGAR_READS
     # --unsafe ALLOW_N_CIGAR_READS is currently necessary for working with RNA-seq data


# GATK: Base recalibration for adjusting base quality scores (based largely on sequencing platform, AT/GC content, and expected haplotypes.  Makes allowances for known sites of variance).
java -jar ${GATK} -T BaseRecalibrator \
     -I ${out4}.bam \
     -R ${fileRef} \
     -knownSites ${fileDbsnp} \
     -knownSites ${file1000G} \
     -knownSites ${pathRef}/Mills_and_1000G_gold_standard.indels.hg19.vcf \
     -o ${out5}.table


# GATK: Reassign base quality scores and write processed read data
java -jar ${GATK} -T PrintReads \
     -R ${fileRef} \
     -I ${out4}.bam \
     -BQSR ${out5}.table \
     -o ${out6}.bam


# GATK: Create plots to visualize results of base recalibration after building a second model of base recalibration.  Optional but highly recommended for QC.
#     -T AnalyzeCovariates \
# java -jar GenomeAnalysisTK.jar \
#     -R ${fileRef} \
#     -before recal2.table \
#     -after recal3.table \
#     -plots recalQC.pdf
