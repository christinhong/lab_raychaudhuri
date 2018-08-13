outdir=/medpop/srlab/mgutierr/tofiPub/mapping/gatkMergedChrs/

#chromosome to process
chr=chr$1

#merge all files from the individual in turn for the input chromosome
samtools merge -R $chr $outdir/$chr.bam /medpop/srlab/mgutierr/tofiPub/mapping/SR*/starPass2/gatk/split-BQSR.bam

#sort and index the chromosome file
samtools sort $outdir/$chr.bam $outdir/$chr.sorted
samtools index $outdir/$chr.sorted.bam

#call variants with HaplotypeCaller
java -Xmx50g -jar  /home/unix/mgutierr/src/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /humgen/gsa-hpprojects/GATK/bundle/current/hg19/ucsc.hg19.fasta -I $outdir/$chr.sorted.bam -recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $outdir/output.$chr.vcf -nct 40 &> $outdir/hapCall.$chr.log
