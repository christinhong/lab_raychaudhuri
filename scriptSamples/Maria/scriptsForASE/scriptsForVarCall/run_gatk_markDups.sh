
gatk=/home/unix/mgutierr/src/GATK
picard=/home/unix/mgutierr/src/picard-tools-1.119
genome=/humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta

infile=$1
#/medpop/srlab/cd4timelinepilot/rnaseq/141125_PR1504/mapping/20141125-0hr-NT-20141014-IDX2-PR1504_S1_R1/starPass2/Aligned.out.bam
#/medpop/rnaseq/access/mapping/140509_SN7001282_0537_BH0RGBADXX_Analysis/CD4_S5_L001/starPass2/Aligned.out.bam

out=$(dirname $1)/gatk
[[ ! -d $out ]] && mkdir -p $out
cd $out

#2. Add read groups, sort, mark duplicates, and create index
java -jar $picard/AddOrReplaceReadGroups.jar I=$infile O=rg_added_sorted.bam SO=coordinate RGID=readGroupID RGLB=libraryID RGPL=illumina RGPU=HiSeq RGSM=sampleID 

java -jar $picard/MarkDuplicates.jar I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

#need to index since MarkDups is not being ran in this case
#samtools index rg_added_sorted.bam

#3. Split'N'Trim and reassign mapping qualities

java -jar $gatk/GenomeAnalysisTK.jar -T SplitNCigarReads -R $genome -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

#5. Base Recalibration
java -Xmx4g -jar $gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -I split.bam -R $genome -knownSites /humgen/gsa-hpprojects/GATK/bundle/current/hg19/dbsnp_138.hg19.vcf -knownSites /humgen/gsa-hpprojects/GATK/bundle/current/hg19/1000G_phase1.snps.high_confidence.hg19.vcf -knownSites /humgen/gsa-hpprojects/GATK/bundle/current/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -o recal_data.table -nct 4

java -Xmx60g -jar $gatk/GenomeAnalysisTK.jar -T PrintReads -R $genome -I split.bam -BQSR recal_data.table -o split-BQSR.bam -nct 4

