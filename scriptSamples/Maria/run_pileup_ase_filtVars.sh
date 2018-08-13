####WARNING: THIS IS FOR SUBREAD MAPPED READS
outdir=/medpop/srlab/cd4timelinepilot/rnaseq/141125_PR1504/ASE/FiltScPkBaCnc_subread_masked
hetSites=/medpop/srlab/cd4timelinepilot/rnaseq/141125_PR1504/variantCalling/vars_passed_ScPkBaCnc.vcf

ASEscript=/home/unix/mgutierr/work/rnaseq/access/ASE/samase.pl
genome=/humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta
refRatio=0.5
#geneSNPFile=/home/unix/mgutierr/work/rnaseq/access/ASE/gene_snp_mapping_file.txt

mrs=10 #min reads per site
mq=20 #min mapping quality
bq=10 #min base quality

#/medpop/srlab/cd4timelinepilot/rnaseq/141125_PR1504/mapping/20141125-0hr-NT-20141014-IDX2-PR1504_S1_R1/starPass2/Aligned.out.bam

regex="mapping/(2014.*)/subread_masked"
[[ $1 =~ $regex ]]
sampleID="${BASH_REMATCH[1]}"
indir=$(dirname $1)

echo $sampleID
echo $indir

#echo "Creating gene_snp_mapping file..."
#nice -n20 perl $ASEscript --get_gene_snp_mapping_file --use_gencode -is $hetSites > $outdir/gene_snp_mapping_file.txt

echo "Sorting bam file..."
nice -n20 samtools sort $indir/$sampleID.bam $indir/$sampleID.sorted

echo "Doing pileup..."
nice -n20 samtools pileup -s -B -f $genome $indir/$sampleID.sorted.bam -l $hetSites > $outdir/$sampleID.pileup.txt

echo "Parsing pileup..."
nice -n20 perl $ASEscript -parse_pileup -sp $outdir/$sampleID.pileup.txt -vcf $hetSites > $outdir/$sampleID.parsed_pileup.txt

echo "Calculating ASE..."
nice -n20 perl $ASEscript --calculate_ase_normal -pp $outdir/$sampleID.parsed_pileup.txt -is $hetSites -r $refRatio -gf $outdir/gene_snp_mapping_file.txt -mrs $mrs -mq $mq -bq $bq > $outdir/$sampleID.ase.txt
