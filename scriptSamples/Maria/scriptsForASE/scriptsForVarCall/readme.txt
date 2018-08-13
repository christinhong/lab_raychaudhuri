
Software you will need:

- GATK
- Picard tools
- vcf tools
- samtools
- SNPiR scripts (from Piskol et al)

Things to edit or keep in mind in scripts:

- paths to software, input and output files, etc.
- I often use "gnu parallel" to parallelize jobs...not sure if you can use that in partners, may be yes as long as you include in the parallel cmd the bjob? (you can ask Kam)
- the scripts often depend on the directory structure in which I work. For example:
path_of_project/mapping/sample_name/subread(or mapper_name)/sample_name.bam 

##############################################
######          CALLING VARIANTS        ######
##############################################

#### Commands to call variants 76bp single-end bulk RNA-seq ####

## Working directory
cd /medpop/srlab/mgutierr/tofiPub/varCall

## Output directory
/medpop/srlab/mgutierr/tofiPub/varCall

## Copying and slightly modifying scripts from
/home/unix/mgutierr/work/rnaseq/cd4timelinepilot/variantCalling
#/medpop/srlab/cd4timelinepilot/rnaseq/141125_PR1504/variantCalling
#/home/unix/mgutierr/work/rnaseq/testCallVar

## Taking fastq files from
/medpop/srlab/mgutierr/tofiPub/*.fastq

## 1. Star 2-pass mapping
#		a. Using STAR generated genome in: /medpop/rnaseq/access/star_genome3
				(was generataed using sh run_star-prepare-reference.sh)		
#	 	b. First pass: Map reads to reference genome (output will be saved in star_pass1 folder, within each library folder)
#		c. Create reference genome with merged junctions from first pass of all samples
#		d. Second pass: map reads to new reference genome that has all observed splice junctions
##
## 2. Data pre-processing with Picard tools and GATK
# 		a. Add read groups, sort, and create index
#		b. Split'N'Trim and reassign mapping qualities
#		c. Base recalibration

sh try_pipeline_star2pass.sh


## 3. Merging files and calling variants with Haplotype Caller
# (in the files that look like list_chrs.txt there are the chromosome numbers or letters)
# (there are many cmds because I tried different approaches, since it takes a long time for each chr)
		a. Merge per chromosome, sort, index and call variants
		#tried all first for chr22
		#samtools merge -R chr22 /medpop/srlab/cd4timelinepilot/rnaseq/141125_PR1504/mapping/gatkMergedChrs/chr22.bam /medpop/srlab/cd4timelinepilot/rnaseq/141125_PR1504/mapping/201*/starPass2/gatk/split-BQSR.bam
		
		#lorax
		screen
		#parallel -j2 --eta 'sh run_merge_chr_and_HapCall.sh' :::: list_chrs_1.txt
	
		parallel -S 3/tigger,2/lorax,3/piglet --eta 'sh run_merge_chr_and_HapCall.sh' :::: list_chrs.txt
		
		parallel -j2 --eta 'sh run_merge_chr_and_HapCall.sh' :::: list_chrs.txt
		
		#lorax
		screen
		parallel -j2 --eta 'sh run_merge_chr_and_HapCall.sh' :::: list_chrs_1.txt
		
		#piglet
		screen
		parallel -j3 --eta 'sh run_merge_chr_and_HapCall.sh' :::: list_chrs_2.txt
		
		#tigger
		screen
		parallel -j3 --eta 'sh run_merge_chr_and_HapCall.sh' :::: list_chrs_3.txt
		
		b. Move output to common folder
		
		cd /medpop/srlab/mgutierr/tofiPub/mapping/gatkMergedChrs
		mv output.chr* ../../varCall/rawPerChr/
		mv hapCall.chr* ../../varCall/rawPerChr/

## Working directory
cd /medpop/srlab/mgutierr/tofiPub/varCall
		
#		c. Concatenate vcf files of all chrs into one

		vcf-concat rawPerChr/output.chr1.vcf rawPerChr/output.chr2.vcf rawPerChr/output.chr3.vcf rawPerChr/output.chr4.vcf rawPerChr/output.chr5.vcf rawPerChr/output.chr6.vcf rawPerChr/output.chr7.vcf rawPerChr/output.chr8.vcf rawPerChr/output.chr9.vcf rawPerChr/output.chr1{0,1,2,3,4,5,6,7,8,9}.vcf rawPerChr/output.chr20.vcf rawPerChr/output.chr21.vcf rawPerChr/output.chr22.vcf rawPerChr/output.chrX.vcf rawPerChr/output.chrY.vcf > output.all_chrs.raw.sorted.vcf

		#vcf-concat rawPerChr/output.chr*.vcf > output.all_chrs.raw.vcf

## 4. Filtering

	sh try_pipeline_filtVars.sh


#######

