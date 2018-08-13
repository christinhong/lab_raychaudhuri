mgutierr 2:44 PM To summarize:
1) Mask genome and map to masked genome
2) run ASE pipeline with new bam files
3) weigh duplicates and get new allelic counts
4) identify problematic SNPs from your called het variants (SNPs in low mappability regions or that presented mapping bias in simulations (Panousis et al, Genome Biology)
5) run logistic regression

mgutierr 3:10 PM you can find the files for the problematic SNPs here: /PHShome/mg154 
