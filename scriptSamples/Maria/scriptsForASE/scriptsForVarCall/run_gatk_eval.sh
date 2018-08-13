gatk=/home/unix/mgutierr/src/GATK
genome=/humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta
dbsnp=/humgen/gsa-hpprojects/GATK/bundle/current/hg19/dbsnp_138.hg19.vcf
kg=/humgen/gsa-hpprojects/GATK/bundle/current/hg19/1000G_phase1.snps.high_confidence.hg19.vcf

infile=$1
#/home/unix/mgutierr/work/rnaseq/access/variantCalling/output.all_chrs.filtSNPclstr.filtSNPiR.GQ30.biallel.vcf

regex="(.*)/output\.all_chrs\.(.*)\.vcf"
[[ $1 =~ $regex ]]
base="${BASH_REMATCH[2]}"

outdir=$(dirname $1)
cd $outdir

nice -n20 java -Xmx2g -jar $gatk/GenomeAnalysisTK.jar -T VariantEval -R $genome -o eval.$base.gatkreport_dbsnp --dbsnp $dbsnp --eval:set1 $infile

nice -n20 java -Xmx2g -jar $gatk/GenomeAnalysisTK.jar -T VariantEval -R $genome -o eval.$base.gatkreport_1KGP --dbsnp $kg --eval:set1 $infile
