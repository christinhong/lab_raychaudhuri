#!/bin/bash

## Working directory
workingDir=/medpop/srlab/mgutierr/tofiPub/varCall
genome=/humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta
gatk=/home/unix/mgutierr/src/GATK
mappingDir=/medpop/srlab/mgutierr/tofiPub/mapping
suffixSamples=SRR
scripts=/home/unix/mgutierr/work/rnaseq/cd4timelinepilot/scripts
parallel=/home/unix/slowikow/src/parallel-20140422/src/parallel

cd $workingDir
		
#		c. Concatenate vcf files of all chrs into one

vcf-concat rawPerChr/output.chr1.vcf rawPerChr/output.chr2.vcf rawPerChr/output.chr3.vcf rawPerChr/output.chr4.vcf rawPerChr/output.chr5.vcf rawPerChr/output.chr6.vcf rawPerChr/output.chr7.vcf rawPerChr/output.chr8.vcf rawPerChr/output.chr9.vcf rawPerChr/output.chr1{0,1,2,3,4,5,6,7,8,9}.vcf rawPerChr/output.chr20.vcf rawPerChr/output.chr21.vcf rawPerChr/output.chr22.vcf rawPerChr/output.chrX.vcf rawPerChr/output.chrY.vcf > output.all_chrs.raw.sorted.vcf

## 4. Filtering
		#a. SNPcluster filter
java -jar $gatk/GenomeAnalysisTK.jar -T VariantFiltration -R $genome  -V output.all_chrs.raw.sorted.vcf -window 35 -cluster 3 -o output.all_chrs.filtSNPclstr.vcf &> snpClstr.log
		
#b. Piskol et al filters

if [[ 1 =~ 2 ]]
then

	
cp -r ~/work/rnaseq/access/variantCalling/SNPiR .	
cd SNPiR
		
### Step 2. Convert VCF format to theirs and filter low qual variants

	#rm indels from vcf file

vcftools --vcf ../output.all_chrs.raw.sorted.vcf --remove-indels --recode --recode-INFO-all --out ../output.all_chrs.raw.sorted.rmIndel &> rmIndels.log

sh convertVCF.sh ../output.all_chrs.raw.sorted.rmIndel.recode.vcf ../output.all_chrs.rmIndel.txt 20

#This one did not work las time, so I will skip to next while merging is done (in tigger)
#### Step 3. Remove mismatches in first 6 bp of reads: 
#screen
#samtools merge $mappingDir/gatkMergedChrs/rg_added_sorted.mergedAllChrs.bam $mappingDir/$suffixSamples*/starPass2/gatk/rg_added_sorted.bam
#use bam file in which reads haven't passed by the  split thingy and probably haven't messed up cigar strngs
#perl filter_mismatch_first6bp.pl -infile /home/unix/mgutierr/work/rnaseq/access/variantCalling/piskolFilts/output.gatk.rmIndel.SNPiR.txt -outfile /home/unix/mgutierr/work/rnaseq/access/variantCalling/piskolFilts/output.gatk.rmIndel.SNPiR.rmhex.txt -bamfile /medpop/rnaseq/access/mapping/mergedGATKbams/rg_added_sorted.mergedAllChrs.bam


#Step 4. Use bedtools to remove sites in repetitive regions based on RepeatMasker annotation
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../output.all_chrs.rmIndel.txt > ../output.all_chrs.rmIndel.bed
intersectBed -a ../output.all_chrs.rmIndel.bed -b /medpop/srlab/mgutierr/GFs/RepeatMasker_UCSCtab_hg19_Nov2014.bed -v > ../output.all_chrs.rmIndel.rmsk.bed
awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' ../output.all_chrs.rmIndel.rmsk.bed > ../output.all_chrs.rmIndel.rmsk.txt

echo "Done with removing snps in repeats..."
#Further Filtering
#Step 1. Filter intronic candidates that are within 4 bp of splicing junctions: perl filter_intron_near_splicejuncts.pl
perl filter_intron_near_splicejuncts.pl -infile ../output.all_chrs.rmIndel.rmsk.txt -outfile ../output.all_chrs.rmIndel.rmsk.rmintron.txt -genefile UCSCgenes_knowGene_hg19_Nov2014.SNPiRformat.txt

echo "Done with filtering out snps in intronic near splice junctions..."

#Step 2. Filter candidates in homopolymer runs: perl filter_homopolymer_nucleotides.pl
perl filter_homopolymer_nucleotides.pl -infile ../output.all_chrs.rmIndel.rmsk.rmintron.txt -outfile ../output.all_chrs.rmIndel.rmsk.rmintron.rmhom.txt -refgenome $genome

echo "Done filtering out homopolymer runs..."
	
			##skipping BLAT step

#Step 4. Use bedtools to separate out candidates that are known RNA editing sites
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../output.all_chrs.rmIndel.rmsk.rmintron.rmhom.txt > ../output.all_chrs.rmIndel.rmsk.rmintron.rmhom.bed
intersectBed -a ../output.all_chrs.rmIndel.rmsk.rmintron.rmhom.bed -b /medpop/srlab/mgutierr/GFs/darned_UCSCtab_Nov2014.bed -v > ../output.all_chrs.rmIndel.rmsk.rmintron.rmhom.rmedt.bed
awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' ../output.all_chrs.rmIndel.rmsk.rmintron.rmhom.rmedt.bed > ../output.all_chrs.rmIndel.rmsk.rmintron.rmhom.rmedt.txt

echo "Done filtering out RNA editing sites..."


		
	#done! (except step 3)
			#Apply Piskol filters to vcf file
			
perl $scripts/create_piskolFilt_vcf.pl ../output.all_chrs.rmIndel.rmsk.rmintron.rmhom.rmedt.txt ../output.all_chrs.filtSNPclstr.vcf > ../output.all_chrs.filtSNPclstr.filtPiskol.vcf

echo "Done creating vcf with all used piskol Filts..."
	
		#c. GQ filter? This one I every time I try it doesn't really help, so not necessary
cd ..
perl $scripts/filter_vcf_byGQ.pl output.all_chrs.filtSNPclstr.filtPiskol.vcf > output.all_chrs.filtSNPclstr.filtPiskol.filtGQ.vcf

echo "Done filtering by genotype qual..."
		
#d. Each allele seen twice filter 
		
#need to run ASE pipeline on all biallelic variants
		
#select only bi allelic

java -Xmx2g -jar /home/unix/mgutierr/src/GATK/GenomeAnalysisTK.jar \
   -R /humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta \
   -T SelectVariants \
   --variant output.all_chrs.filtSNPclstr.filtPiskol.filtGQ.vcf \
   -o output.all_chrs.filtSNPclstr.filtPiskol.filtGQ.biallel.vcf \
   -restrictAllelesTo BIALLELIC

echo "Done selecting only biallelic with gatk..."

#keeping only variants that pass filters and are biallelic
java -Xmx2g -jar /home/unix/mgutierr/src/GATK/GenomeAnalysisTK.jar \
   -R /humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta \
   -T SelectVariants \
   --variant output.all_chrs.filtSNPclstr.filtPiskol.filtGQ.vcf \
   -o all_chrs.allFilts.vcf \
   -ef \
   -restrictAllelesTo BIALLELIC

echo "Done ballelic cmd 2..."

		
perl $scripts/filter_vcf_byHets.pl output.all_chrs.filtSNPclstr.filtPiskol.filtGQ.biallel.vcf > output.all_chrs.filtSNPclstr.filtPiskol.filtGQ.biallel.filtHom.vcf
			
echo "Done filtering  vcf by hets..."

fi

java -Xmx2g -jar /home/unix/mgutierr/src/GATK/GenomeAnalysisTK.jar \
   -R /humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta \
   -T SelectVariants \
   --variant output.all_chrs.filtSNPclstr.filtPiskol.vcf \
   -o output.all_chrs.filtSNPclstr.filtPiskol.biallel.vcf \
   -restrictAllelesTo BIALLELIC

#apply hom filter to  file below:
nice -n20 perl ~/work/rnaseq/access/scripts/filter_vcf_byHets.pl output.all_chrs.filtSNPclstr.filtPiskol.biallel.vcf > output.all_chrs.filtSNPclstr.filtPiskol.biallel.filtHom.vcf

#Get Hets with alleles seen in at least 2 samples
#nice -n 20 perl ~/work/rnaseq/access/scripts/create_easFilt_vcf.pl /home/unix/mgutierr/work/rnaseq/cd4timelinepilot/ASE/list_snps_as_min2.txt /medpop/srlab/cd4timelinepilot/rnaseq/141125_PR1504/variantCalling/output.all_chrs.filtSNPclstr.filtPiskol.biallel.filtHom.vcf > /medpop/srlab/cd4timelinepilot/rnaseq/141125_PR1504/variantCalling/output.all_chrs.filtSNPclstr.filtPiskol.biallel.filtHom.filtEAS2.vcf

#Get hets concordant with dbSNP
 java -Xmx2g -jar /home/unix/mgutierr/src/GATK/GenomeAnalysisTK.jar \
 	-R /humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta \
   -T SelectVariants \
   --variant output.all_chrs.filtSNPclstr.filtPiskol.biallel.filtHom.vcf \
   --concordance /humgen/gsa-hpprojects/GATK/bundle/current/hg19/dbsnp_138.hg19.vcf \
   -o output.all_chrs.filtSNPclstr.filtPiskol.biallel.filtHom.dbsnp.vcf


grep "0/1" output.all_chrs.filtSNPclstr.filtPiskol.biallel.filtHom.dbsnp.vcf | grep PASS > hets_allFilts.vcf

perl $scripts/add_ids_vcf.pl hets_allFilts.vcf > hets_allFilts.ids.vcf

echo "Done getting hets with ids."

ls -1 output*.vcf > list_vcfs2eval.txt

$parallel -j7 --eta 'sh /home/unix/mgutierr/work/rnaseq/cd4timelinepilot/scripts/run_gatk_eval.sh' :::: list_vcfs2eval.txt

perl ~/work/rnaseq/access/scripts/parse_gatk_eval_report.pl summary_eval_1KGP.txt eval.*_1KGP
perl ~/work/rnaseq/access/scripts/parse_gatk_eval_report.pl summary_eval_dbsnp.txt eval.*_dbsnp
			

