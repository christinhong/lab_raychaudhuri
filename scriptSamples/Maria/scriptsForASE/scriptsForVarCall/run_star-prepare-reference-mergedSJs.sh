# run_star-prepare-reference.sh
# Kamil Slowikowski
# July 18, 2014
#
# Modified by Maria, Nov 2014
# Prepare reference files needed to run the STAR read mapper, using observed splice-junctions from pass 1.

# sh run_star-prepare-reference-mergedSJs.sh directoryForOutputRefGenome fileWithMergedSpliceJunctions

#genome=/medpop/rnaseq/access/kam_ucsc_genome/genome.fa

#screen
#cd /medpop/rnaseq/access/mapping/140506_SN7001282_0527_AH0H0GADXX/CD5_S4_L002/hg19_2pass
#STAR --runMode genomeGenerate --genomeDir /medpop/rnaseq/access/mapping/140506_SN7001282_0527_AH0H0GADXX/CD5_S4_L002/hg19_2pass --genomeFastaFiles /humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta --sjdbFileChrStartEnd /medpop/rnaseq/access/mapping/140506_SN7001282_0527_AH0H0GADXX/CD5_S4_L002/SJ.out.tab --sjdbOverhang 75 --runThreadN 4

#regex="(14[a-zA-Z0-9_]+)/([a-zA-Z0-9]+_S[0-9]+_L[0-9]+)"
#/medpop/rnaseq/access/mapping/140506_SN7001282_0527_AH0H0GADXX/CD4_S2_L002/CD4_S2_L002_R1_001.fastq

#[[ $1 =~ $regex ]]
#project="${BASH_REMATCH[1]}"
#library="${BASH_REMATCH[2]}"

genome=/humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta

out=$1/GenomeForStarPass2
[[ ! -d $out ]] && mkdir -p $out

# STAR generates files in the current directory.
cd $out

log=$out/log.txt

opt=(
--runMode genomeGenerate
--genomeDir $out
--genomeFastaFiles $genome
--runThreadN 4
--sjdbFileChrStartEnd $2
--sjdbOverhang 75
)

STAR ${opt[*]} &> $log
