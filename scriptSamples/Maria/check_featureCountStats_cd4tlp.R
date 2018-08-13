# ANALYZE GENE EXPRESSION 
# ANALYZE FEATURE COUNT STATS!
wd <- "/Users/mgutierr/Documents/work/rnaseq/cd4timelinepilot/geneExpr"
setwd(wd)
m <- read.table("summaryStats_featureCounts.txt", stringsAsFactors = F, header = T, row.names = 1)

m

tpts <- sapply(rownames(m), function(id){
  strsplit(id, "-")[[1]][2]
})
tpts
rownames(m) <- tpts
###comparing all fragments lengths vs FL > 25K
m2 <- m[c(1,4,5,7,2,3,6), c(5, 6, 4, 2, 1)]

m2
x11(width = 6, height = 6)
mosaicplot(m2, color= c("black", "darkorange2", "gray", "goldenrod1", "forestgreen"), las = 1, cex.axis = 1.2, main = "")
dev.print( "mosaciplot_fc_stats_cd4timelinepilot_mapq20.pdf", dev = pdf)

m2
barplot(rowSums(m2), names = rownames(m2), ylab = "Number of sequenced reads", cex.axis = 1.5, cex.lab = 1.5, cex.name = 1.3)
dev.print("numberOfSequencedReadsPerLib_cd4tlp.pdf", dev = pdf)

barplot(m2[,3] + m2[,4] + m2[,5], names = rownames(m2), ylab = "Number mapped reads with mapq > 20", cex.axis = 1.5, cex.lab = 1.5, cex.name = 1.3)
dev.print("numberOfMappedReadsWithMinMapq20PerLib_cd4tlp.pdf", dev = pdf)



?barplot
