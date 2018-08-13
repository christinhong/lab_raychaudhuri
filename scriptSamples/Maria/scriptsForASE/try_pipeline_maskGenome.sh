#!/bin/bash

# Mask genome on heterozygous sites and map to masked genome using subread. 

## Working directory
workingDir=/medpop/srlab/mgutierr/tofiPub
genome=/medpop/srlab/mgutierr/ucsc_hg19/genome.fa
#genome=/humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta
hetSites=/medpop/srlab/mgutierr/tofiPub/varCall/hets_allFilts.ids.vcf
outDir=/medpop/srlab/mgutierr/tofiPub
scripts=/home/unix/mgutierr/work/rnaseq/cd4timelinepilot/scripts

#njobs=3
#listFastq=
#suffixSamples=SRR
#parallel=/home/unix/slowikow/src/parallel-20140422/src/parallel

cd $workingDir

regex="(/.*)\.vcf"
[[ $hetSites =~ $regex ]]
baseHets="${BASH_REMATCH[1]}"

echo $baseHets

outMasked=$outDir/maskedGenome
[[ ! -d $out ]] && mkdir -p $out

#convert vcf to bed
echo "Converting vcf to bed..."
perl $scripts/vcf2bed.pl $hetSites > $baseHets.bed

#mask genome in heterozygous sites
echo "Masking genome with bed tools..."
maskFastaFromBed -fi $genome -bed $baseHets.bed -fo $outMasked/ucsc_hg19_genome_maskedHets.fa

#create new ref genome with subread
echo "Creating new subread masked indexed genome..."
cd $outMasked
nice -n20 /home/unix/slowikow/src/subread-1.4.5-p1-source/bin/subread-buildindex -o hg19-masked-subread $outMasked/ucsc_hg19_genome_maskedHets.fa


###### map to masked genome

###NEED TO MODIFY SCRIPT TO ACCEPT NEW MASKED GENOME AS PARAMETER?
# Working Dir
#cd $outDir/mapping

#parallel -j$njobs --eta 'sh run_subread-allign_cd4timelinepilot_u_Q_FL100K_maskedGenome.sh' :::: list_fastq_gz.txt
