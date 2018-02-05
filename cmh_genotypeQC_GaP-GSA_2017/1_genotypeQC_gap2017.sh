 #!/usr/local/bin/bash 

# Christin M. Hong
# Lab PI: Soumya Raychaudhuri, Harvard Medical School
# Last updated: 2017-07


#### Project Notes ####

# Source: GaP Long Island registry
# Platform: QC Illumina Global Screening Array (genotype data)
# 
# Script purpose(s): 
# 	1. Process and QC GaP genotype data, do PCA.
# 	2. Merge with cleaned genotypes from 1000 Genomes, do PCA.
# 
# Script started off Anderson et al. Nat Protoc. 2010, per Emma and Harm-Jan's advice
# Merged with Yang Luo's genotyping QC script for LIIMA
# Coded with suggestions from Yang, Maria, and Harm-Jan


#### Infrastructure ####

## Log onto interactive shell for compute-intensive work at Broad with Broad Mac
# from https://intranet.broadinstitute.org/bits/service-catalog/scientific-computing/login-servers/login-servers

ssh silver.broadinstitute.org # login node, as are gold and platinum
ssh tigger 
#or "ssh piglet" for lab nodes with 64 CPU + 250 GB RAM each.
# On tigger/lab servers, need to prepend ALL commands with "nice -n20" to avoid locking up login node!!


## Set up interactive session
use UGER
ish # expires at 36h, default is 1G memory
#ish -l h_vmem=4G # requesting 4G memory

## Load useful Broad dotkit tools
reuse UGER
#use -l # list of available tools
use PLINK2
use R-3.3
use Perl-5.20
use Bcftools
use VCFtools



## Set project-specific variables (remember that target names with underscores need to be within quotes)
export proj="/medpop/srlab/chong/immunomics/genotypeQC_GaP-GSA_2017"

export ORIGDATA=${proj}"/data/orig/GAP_July17_2017"

export OUTPUT="gap201707"

export DATA=${proj}"/data/"${OUTPUT}

export OUT1KG=${OUTPUT}"-1kg"


## Set general variables
export resources="/medpop/srlab/chong/aaa_resources"

export parallel=${resources}"/toolkits/parallel/parallel-20170522/src/parallel"

export intCPUs=22

export fileChrList=${resources}"/ref/fileChrList.txt"

export G1K="/medpop/srlab/external-data/1000genomes/ftp/release/20130502"



#### Set up project directory ####
mkdir ${proj}
cd ${proj}

mkdir scripts
mkdir analysis
mkdir -p data/orig



#### Process original data ####
# Copy original data (including hidden files/folders) to project data folder
cp -R /medpop/srlab/Gap2017/. ${proj}/data/orig/

# Convert PED and MAP files to Plink BED, BIM, and FAM.  Only getting SNPs (no indels).
cd ${proj}/data

nice -n20 plink --file ${ORIGDATA} --make-bed --snps-only --out ${proj}/data/${OUTPUT}

# Output: 698411 variants loaded from .bim file.  1835 people (0 males, 0 females, 1835 ambiguous) loaded from .fam.  Total genotyping rate is 0.991331.
	# 1835 people here vs. 2504 in 1000 Genomes = should be okay to do PCA without EIGENSTRAT/weighting variants



#### Check build ####
## Using HJ's CD28/CTLA4 SNPs: rs117701653 (not in data, probably imputed), rs3087243
grep rs3087243 ${OUTPUT}.bim

# Output: 2	rs3087243	0	204738919	A	G
# From https://www.immunobase.org/marker/rs3087243/, can see this is on genome build GRCh37 = hg19, so no need to liftover.


## If need to lift over (use appropriate chain file - this is for hg18 to hg19):
# ${resources}/tools/liftOverPlink.py -m ${ORIGDATA}.map -p ${ORIGDATA}.ped -o ${proj}/data/${OUTPUT} -e ${resources}/tools/liftOver -c ${resources}/ref/hg18ToHg19.over.chain.gz 

# plink --file ${proj}/data/${OUTPUT} --make-bed --snps-only --out ${proj}/data/${OUTPUT}


#### 0 vs. 1-based position: If necessary, add +1 to genomic positions in BIM. ####
# cp ${OUTPUT}.bim ${OUTPUT}-before-adding-1-to-position.bim
# 
# head ${OUTPUT}-before-adding-1-to-position.bim
# 
# awk -F $'\t' 'BEGIN { OFS = FS } { $4 = $4 + 1; print }' ${OUTPUT}-before-adding-1-to-position.bim > ${OUTPUT}.bim



#### Check for duplicate IIDs ####
cat ${DATA}.fam | awk '{ print $1, $2 }' | sort | uniq -d > dup-iids.txt

wc dup-iids.txt

# No duplicate individuals in this data.


# If there are duplicate IIDs - which break Plink - identify and make unique by appending row number
# cat ${proj}/data/${OUTPUT}.fam  | awk '{ print $2, NR }' > iids-rownumbers.txt
# 
# cat iids-rownumbers.txt | awk '{ print $1"-"$2 }' > iids-updated.txt
#  
# cat iids-updated.txt | awk '{ print $1 }' | less


# Manually replace old IID column with the new IID column (again, because Plink breaks when trying to handle dup IID)
# cp ${OUTPUT}.fam ${OUTPUT}-withDupIIDs.fam

# awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' iids-updated.txt ${OUTPUT}-withDupIIDs.fam > ${OUTPUT}.fam

# head ${OUTPUT}.fam # OK, nice!


# I think dup IIDs are more likely to be mislabeled.  
# Since we want to select samples for callback based on the results of this analysis, extract duplicate IIDs (now with append row number) to store in a fail filter for later removal from dataset
# cat dup-iids.txt | awk '{ print $2 }' | sort > dup-iids-only.txt

# grep -f ${proj}/data/dup-iids-only.txt ${proj}/data/${OUTPUT}.fam | sort |  awk '{ print $1, $2 }' > ${proj}/analysis/fail-dup-iids.txt




#### Identify individuals with discordant sex information ####
cd ${proj}/analysis/

# Calculate mean homozygosity rate across X chromosome markers for each individual
nice -n20 plink --bfile ${DATA} --check-sex --out ${OUTPUT}

# All ambiguous -> won't do sexcheck QC.


# If doing sexcheck QC:
# Generate list of individuals with discordant sex data
# grep PROBLEM ${OUTPUT}.sexcheck > ${OUTPUT}.sexprobs

# Open list (see explanation of list format in Anderson 2010 tutorial)
# nano ${OUTPUT}.sexcheck

# Add family ID (FID, column 1) and individual ID (IID, column 2) to file 'fail-sexcheck-qc.txt' (one individual per line, tab delimited).
# awk -v OFS='\t' '{print $1,$2}' ${OUTPUT}.sexprobs > fail-sexcheck-qc.txt



#### Identify individuals with elevated missing data rates / outlying heterozygosity rate ####

# Calculate missing genotypes (can't resolve as aa, Aa, or AA from raw intensity sequencing data) and store in IMISS (by individual) and LMISS (by locus) files.
nice -n20 plink --bfile ${DATA} --missing --out ${OUTPUT} --allow-no-sex

# Create a HET file that records number of homozygous genotypes [O(Hom)] in 3rd column and number of non-missing genotypes [N(NM)] per individual in 5th column.
nice -n20 plink --bfile ${DATA} --het --out ${OUTPUT} --allow-no-sex



#### 3. Analyze and remove samples with elevated missingness or outlying heterozygosity
## Calculate the observed homozygosity rate per individual (of non-missing genotypes) using formula
# (N(NM) - O(Hom))/N(NM) (but Yang uses E.HOM instead of N.NM.  Not sure why, but the results look the same when plotted.)

# Plot the heterozygosity rate (X) vs. missingness (Y) per individual using cmh_imiss-vs-het.R, which generates "fail-imisshet-qc.txt" and "${OUTPUT}.imiss-vs-het.pdf"
	# a is the x standard deviation threshold for heterozygosity
	# b is the missingness threshold (lower value = less missingness accepted)

nice -n20 R CMD BATCH '--args filename="gap201707" a=3 b=0.05' ${resources}/tools/cmh-imiss-vs-het.R

# REQUIRES MANUAL UPDATE PER PROJECT DUE TO ${OUTPUT} FILENAME


# Check PDFs by rsync from home (not ssh) terminal.  See end of script.
	# Handful of samples removed for missingness >0.05, none removed for excess heterozygosity



#### Reduce computational complexity by accounting for LD and duplicated/related individuals ####

## Generate list of SNPs to be kept for analysis by removing SNPs with high LD or in regions with high LD.
# High LD: Pair of SNPs within bp window (size?) with r^2 > 0.2.
# Copied the high-LD-regions.txt from Yang (who I presume found it from the Anderson 2010 tutorial)

nice -n20 plink --bfile ${DATA} --exclude range ${resources}/ref/high-LD-regions.txt --indep-pairwise 50 5 0.2 --out ${OUTPUT} --allow-no-sex

# 376437 of 670452 variants removed (56%)


## Generate pairwise IBS (identity by state) matrix for all individuals based on retained SNPs
# Only took ~5 minutes for GaP samples, but can be a long process -> prefix with "nohup" to allow command to continue after logout.  Adding & to end allows continuing use of terminal while command runs in background.
	# Generates files nohup.out (log of usual terminal output) and OUTPUT.genome (file of interest)

nice -n20 nohup plink --bfile ${DATA} --allow-no-sex --extract ${OUTPUT}.prune.in --genome --out ${OUTPUT} &

pid=$!

ls

# Check for when done
#wait $pid

# Check and clear log
cat nohup.out
rm nohup.out


## Identify all pairs of individuals with an IBD >0.185 (indicates relatedness closer than third-degree relatives)
# Uses .genome and .imiss files.  Outputs "fail-IBD-QC.txt"

nice -n20 perl ${resources}/tools/run-IBD-QC.pl ${OUTPUT}


# Consolidate remaining SNPs
nice -n20 plink --bfile ${DATA} --extract ${OUTPUT}.prune.in --make-bed --out ${OUTPUT}-pruned

# Total genotyping rate is 0.994154.  294015 variants and 1835 people pass filters and QC.



#### Remove SNPs that fall below missingness and MAF thresholds, and are out of Hardy-Weinberg equilibrium #### 
nice -n20 plink --bfile ${OUTPUT}-pruned --mind 0.1 --maf 0.05 --geno 0.05 --hwe 0.00001 --make-bed --out ${OUTPUT}-hwe

	#12 people removed due to missing genotype data (--mind).
#IDs written to gap201707-hwe.irem .
#Before main variant filters, 1823 founders and 0 nonfounders present.
#Total genotyping rate in remaining samples is 0.995868.
#4690 variants removed due to missing genotype data (--geno).
#--hwe: 13512 variants removed due to Hardy-Weinberg exact test.
#169845 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#105968 variants and 1823 people pass filters and QC.



#### Remove duplicate variants and individuals that failed QC ####
# Data processing: Rename variants as CHR:BP ("old name" "new name")
awk '{print $2,$1":"$4}' ${OUTPUT}-hwe.bim > ${OUTPUT}-newname-map.txt

# Remake bim with new names
nice -n20 plink --bfile ${OUTPUT}-hwe --update-name ${OUTPUT}-newname-map.txt --make-bed --out ${OUTPUT}-hwe-renamed

# Find duplicate SNPs
awk '{print $2}' ${OUTPUT}-hwe-renamed.bim | uniq -d > dup-SNPs.txt

# Remake files without duplicate SNPs
nice -n20 plink --bfile ${OUTPUT}-hwe-renamed --exclude dup-SNPs.txt --make-bed --out ${OUTPUT}-dedup

# Total genotyping rate is 0.997212.  105968 variants and 1823 people pass filters and QC.


# Find individuals that failed QC
cat fail-* | sort -k1 | uniq > inds-fail-qc.txt

# Remake files without problematic individuals
nice -n20 plink --bfile ${OUTPUT}-dedup --remove inds-fail-qc.txt --make-bed --out ${OUTPUT}-clean

# Total genotyping rate in remaining samples is 0.998371.  105968 variants and 1752 people pass filters and QC.



#### Remove variants under minor allele frequency threshold and calculate first 20 PCs for ancestry analysis ####
## Pruning and extracting variants that pass MAF and indep-pairwise thresholds

nice -n20 plink --bfile ${OUTPUT}-clean --maf 0.05 --indep-pairwise 50 10 0.1 --out ${OUTPUT}-clean-maf 

nice -n20 plink --bfile ${OUTPUT}-clean --extract ${OUTPUT}-clean-maf.prune.in --make-bed --out ${OUTPUT}-clean-maf

# Total genotyping rate is 0.998336.  67215 variants and 1752 people pass filters and QC.


# Run PCA
nice -n20 plink --bfile ${OUTPUT}-clean-maf --pca --out pca-${OUTPUT}



#### Identify ancestry outliers by comparing with 1000 Genomes ####

## Subset 1000 Genome data with Bcftools for overlapping variants
# Subset QC'd BIM such that column 1 = chr and 2 = position, tab-delimited, since that's the subsetting format accepted by Bcftools
awk -v OFS='\t' '{print $1, $4}' ${OUTPUT}-clean-maf.bim > SNPs-for-subsetting.txt


# If necessary, generate CSI index (instead of TBI) for the genotype data from 1000 Genomes.  (Bcftools has issues with TBI, but this only needs to be done once.  Took ~1 hour.)



#### Subset 1KG VCFs (which are available by chr) based on SNPs common with project data ####
# Time: ~10 minutes)
nice -n20 nohup ${parallel} -j ${intCPUs} -k 'sh ${resources}/tools/bcftools-subset1kG.sh' :::: ${fileChrList} &

pid=$!

ls

wc temp-shared*

# Wait until done
#wait $pid

# Check and clean up
cat nohup.out # Expect that subsetting sex chromosomes failed
rm nohup.out

# Note: Above script excludes sex chromosomes, per Maria's suggestion.  Also, Y chromosomes trigger "different number of samples" error with Bcftools (since not all samples have Y chr).



## Create sorted list of VCF files, concatenate VCFs with Bcftools, and sort merged VCF
nice -n20 find ${proj}/analysis/ -name "temp-shared-G1K-chr[0-9]*.vcf" -type f -printf '%f\n' | sort > ${proj}/analysis/temp-shared-G1K.txt

nice -n20 bcftools concat -o ${proj}/analysis/shared-G1K.vcf -O v -f ${proj}/analysis/temp-shared-G1K.txt

nice -n20 cat ${proj}/analysis/shared-G1K.vcf | vcf-sort > shared-G1K-sorted.vcf



#### Converting VCF of 1000 Genomes subset to plink format ####

nice -n20 vcftools --vcf shared-G1K-sorted.vcf --plink-tped --out ${OUT1KG}

# For above command, I originally used --plink flag instead of --plink-tped flag for above command.  That works on silver.broadinstitute.org, but I'm not allowed to open enough temp files on tigger ("ulimit -n 10000" -> operation not permitted).
	# I've also tried using Plink 1.9b version for converting the VCF to Plink format, but that runs into issues when renaming the variants.  I think VCF names some things '.', which Plink dislikes. 



nice -n20 plink --tfile ${OUT1KG} --make-bed --snps-only --out ${OUT1KG}

# 67027 variants loaded from .bim file.
# 67027 variants and 2504 people pass filters and QC.


# Remove temp VCFs
rm temp-shared-G1K-chr*



#### Filtering 1kG for missingness and outlying heterozygosity ####
nice -n20 plink --bfile ${OUT1KG} --missing --out ${OUT1KG} --allow-no-sex

nice -n20 plink --bfile ${OUT1KG} --het --out ${OUT1KG} --allow-no-sex

nice -n20 R CMD BATCH '--args filename="gap201707-1kg" a=3 b=0.05' ${resources}/tools/cmh-imiss-vs-het.R

# REQUIRES MANUAL UPDATE PER PROJECT DUE TO ${OUT1KG} FILENAME

# Check PDFs by rsync from home (not ssh) terminal.  See end of script.
	# No missingness at all -> breaks R code for plotting missingness vs. outlier het.
	# However, outlying hets are still captured in the "fail-imisshet-qc.txt" generated by the R script, so that's OK.


## Remove samples that failed G1K het QC
nice -n20 plink --bfile ${OUT1KG} --remove fail-imisshet-qc.txt --make-bed --out ${OUT1KG}-qc

# 67027 variants and 2478 people pass filters and QC.



#### Filter out duplicate variants and individuals that failed QC ####
## Note: Need to dedup before renaming is possible with 1kG data...not sure why, but OK.

# Find duplicate SNPs
awk '{print $2}' ${OUT1KG}-qc.bim | uniq -d > dup-SNPs2.txt

# Remake files without duplicate SNPs
nice -n20 plink --bfile ${OUT1KG}-qc --exclude dup-SNPs2.txt --make-bed --out ${OUT1KG}-dedup

# 66981 variants and 2478 people pass filters and QC.


## Rename variants as CHR:BP ("old name" "new name")
awk '{print $2,$1":"$4}' ${OUT1KG}-dedup.bim > ${OUT1KG}-newname-map.txt

# Remake bim with new names
nice -n20 plink --bfile ${OUT1KG}-dedup --update-name ${OUT1KG}-newname-map.txt --make-bed --out ${OUT1KG}-renamed



#### Filter for MAF in 1KG ####
nice -n20 plink --bfile ${OUT1KG}-renamed --maf 0.05 --indep-pairwise 50 10 0.1 --out ${OUT1KG}-renamed-maf 

nice -n20 plink --bfile ${OUT1KG}-renamed --extract ${OUT1KG}-renamed-maf.prune.in --make-bed --out ${OUT1KG}-renamed-maf

# 55796 variants and 2478 people pass filters and QC.



#### Subset project data for common SNPs with G1K ####
awk '{print $2}' ${OUT1KG}-renamed-maf.bim > SNPs-for-subsetting2.txt

nice -n20 plink --bfile ${OUTPUT}-clean-maf --extract SNPs-for-subsetting2.txt --make-bed --out ${OUTPUT}-shared

# Total genotyping rate is 0.998342.  55796 variants and 1752 people pass filters and QC. 


# Check that number of lines = variants is the same
wc ${OUT1KG}-renamed-maf.bim
wc ${OUTPUT}-shared.bim

cat ${OUT1KG}-renamed-maf.bim ${OUTPUT}-shared.bim | awk '{print $2}' | sort | uniq -u
# No unique names, good.



#### Merge common variants from data with 1000G and re-QC ####
nice -n20 plink --bfile ${OUT1KG}-renamed-maf --bmerge ${OUTPUT}-shared --make-bed --out ${OUTPUT}-merged

# Error: 28103 variants with 3+ alleles present. 


# Flip strand due to error log.  See https://www.cog-genomics.org/plink/1.9/data#merge
nice -n20 plink --bfile ${OUTPUT}-shared --flip ${OUTPUT}-merged-merge.missnp --make-bed --out ${OUTPUT}-flip


nice -n20 plink --bfile ${OUT1KG}-renamed-maf --bmerge ${OUTPUT}-flip --make-bed --out merged_trial

# Error: 1 variant with 3+ alleles present.


# Exclude variants that remain problematic after flipping
nice -n20 plink --bfile ${OUT1KG}-renamed-maf --exclude merged_trial-merge.missnp --make-bed --out tmp-${OUT1KG}-renamed

nice -n20 plink --bfile ${OUTPUT}-flip --exclude merged_trial-merge.missnp --make-bed --out tmp-${OUTPUT}-flip


# Final merge
nice -n20 plink --bfile tmp-${OUT1KG}-renamed --bmerge tmp-${OUTPUT}-flip --make-bed --out ${OUTPUT}-merged

# Total genotyping rate is 0.999313.  55795 variants and 4230 people pass filters and QC.


# Clean
rm tmp-${OUT1KG}-renamed.*

rm tmp-${OUTPUT}-flip.*



## Prune and extract variants that pass MAF threshold
nice -n20 plink --bfile ${OUTPUT}-merged --maf 0.05 --indep-pairwise 50 10 0.1 --out ${OUTPUT}-merged-maf 

nice -n20 plink --bfile ${OUTPUT}-merged --extract ${OUTPUT}-merged-maf.prune.in --make-bed --out ${OUTPUT}-merged-maf

# Total genotyping rate is 0.999315.  55608 variants and 4230 people pass filters and QC.



## Identify markers with an excessive missing data rate...fails because G1K data has 0 missingness, so I didn't change this much from Yang's script.
#plink --bfile ${OUTPUT}-merged-maf --missing --out ${OUTPUT}-merged-maf

#R CMD BATCH ${resources}/tools/lmiss-hist.Rscript

# Test markers for different genotype call rates between cases and controls
#plink --bfile ${OUTPUT}-merged-maf --test-missing --out ${OUTPUT}-merged-maf

#awk '$NF<1e-5 {print $2}' clean-ind-LIMAA.missing >fail-diffmiss-qc.txt



#### Final clean! ####
nice -n20 plink --bfile ${OUTPUT}-merged-maf --mind 0.1 --maf 0.05 --geno 0.05 --hwe 0.00001 --make-bed --out ${OUTPUT}-merged-clean

# 0 people removed due to missing genotype data (--mind).
#Before main variant filters, 4230 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.999315.
#0 variants removed due to missing genotype data (--geno).
#--hwe: 17392 variants removed due to Hardy-Weinberg exact test.
#0 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#38216 variants and 4230 people pass filters and QC.



#### Run PCA on merged data ####
nice -n20 plink --bfile ${OUTPUT}-merged-clean --pca --out pca-${OUTPUT}-merged



## Currently plotting and analyzing PCA on R on home computer, but theoretically could run on Broad server.
# Soft-linking to panel for R script to reference while running on remote server:
#ln -s /medpop/srlab/external-data/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel g1k.sample

#R CMD BATCH ${proj}/analysis/pca-ancestry-plots.R



#### Notes ####
## Commands for file transfer (run on home terminal, not Broad server)
# Upload to Broad
rsync -av -e "ssh -l chong" /Users/chong/00_SR-lab-data/aaa_resources/tools/ silver.broadinstitute.org:/medpop/srlab/chong/aaa_resources/tools

# Download to computer (e.g. to run R or look at PDFs)
rsync -av -e "ssh -l chong" silver.broadinstitute.org:/medpop/srlab/chong/immunomics/genotypeQC_GaP-GSA_2017/analysis/ /Users/chong/00_SR-lab-data/immunomics/genotypeQC_GaP-GSA_2017/analysis



## Can check MAF with:
# /medpop/srlab/yang/software/plink --bfile clean-GAP_HG19_rmLDXY_OV1KG_1KGID --freq --out gap
# /medpop/srlab/yang/software/plink --bfile 1KG_OVGAP --freq --out 1kg


# From Samira: Use GCTA to generate PCA for larger number of individuals (faster than EIGENSTRAT)
# only varints with maf > 0.05 and missingness < 0.01 are included in PCA analysis
# gcta=/medpop/srlab/sasgari/Software/gcta/gcta64
# in1=clean-TB-Reich-1kg-clean
# $gcta --bfile $in1 --autosome  --make-grm  --out $in1 --thread-num  30
# $gcta --grm $in1 --thread-num  30 --pca 20  --out $in1
