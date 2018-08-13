# Christin M. Hong
# Last modified: 2015-08-30
# Soumya Raychaudhuri Lab, Harvard Medical School


# R script for analyzing human RNA-seq data for differential expression.
    # Input: Output from script RNA-01_sraToFcounts-genes_cmh*.lsf
    # Processes run: TMM {edgeR}, voom {limma}, and graphing {ggplot2}.

#####

# LIBRARIES

library(edgeR)
library(limma)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(gtable)
library(tidyr) # not yet used
library(dplyr) # not yet used


# GLOBAL VARIABLES

vWorkingDir <- "/home/christin/001_Raychaudhuri-Lab/data/2015-08_RNA_Hirahara-OShea_stat1"
setwd(vWorkingDir)

vIntCPUS <- 4
vDate <- format(Sys.Date(), format="%Y%m%d")

vFileSessionInfo <- paste("aa_RNA-03-DE_SessionAndOptions_", vDate, ".txt", sep = "")

# These variables hold values from upstream script (RNA-01_sraToFcounts-genes_cmh*.lsf). WARNING: Changing will break pipe!
vFileTargets <- "aa_RNA-01_Targets.txt"
vFileFcCounts <- "aa_RNA-01_fcCountsNames.txt"
vFileFcStats <- "aa_RNA-01_fcCounts.txt.summary"

# prefix <- "aa_RNA-03-DE_"

#####

# OPTIONS.  See https://stat.ethz.ch/R-manual/R-devel/library/base/html/options.html
options(digits = 3,
        stringsAsFactors = T) # Default value.  Leaving this alone in case packages are built to expect this behavior, but it has enough ramifications to make the option explicit.


#####

# MOSAIC PLOTS OF FEATURECOUNT STATS

# Store featureCounts stats/summary output in data frame
vFcStats <- read.delim(vFileFcStats, stringsAsFactors = T)
# Check
names(vFcStats)
head(vFcStats, n = 3)


# RELATIVE MOSAIC PLOT

# Invert rows and columns
vFcStats.inv <- as.data.frame(t(vFcStats[,-1]))
colnames(vFcStats.inv) <- vFcStats$Status
names(vFcStats.inv)
head(vFcStats.inv)

# Remove columns with all 0 (or NA, when applicable) values
vFcStats.inv.filt <- vFcStats.inv[, colSums(vFcStats.inv != 0, na.rm = TRUE) > 0]
names(vFcStats.inv.filt)
head(vFcStats.inv.filt)

# OPTIONAL: Reverse the order of columns
vFcStats.inv.filt.cRev <- vFcStats.inv.filt[, rev(seq_len(ncol(vFcStats.inv.filt)))]
head(vFcStats.inv.filt.cRev)

x11()
mosaicplot(vFcStats.inv.filt.cRev, shade = FALSE, color=c("#C0C0C0", "#C0C0C0", "#C0C0C0", "#C0C0C0", "#00ccff"), las = 2, cex.axis = .8, main = "featureCount Stats", sub = "From human CD4+CD45RA(hi)CD45RO(lo) cells with and without STAT1 mutations.\nTreated 3 days with anti-CD3/CD28.")

# Save to PDF
pdf("aa_RNA-02-DE_NUM_fcStat rel mosaic.pdf",width = 10, height = 8, useDingbats = F);
mosaicplot(vFcStats.inv.filt.cRev, shade = FALSE, color=c("#C0C0C0", "#C0C0C0", "#C0C0C0", "#C0C0C0", "#00ccff"), las = 2, cex.axis = .8, main = "featureCount mosaic plot for Hirahara-O'Shea STAT1 data", sub = "From human CD4+CD45RA(hi)CD45RO(lo) cells with and without STAT1 mutations.\nTreated 3 days with anti-CD3/CD28.")

dev.off ();



# QUANTITATIVE MOSAIC/STACKED BAR PLOTS

# Filter out rows with all 0 values
vFcStats.filt <- vFcStats[rowSums(vFcStats[, -(1)]) > 0, ]
names(vFcStats.filt)
head(vFcStats.filt)


# Transform stat data from wide to long format
vFcStats.filt.long = melt(vFcStats.filt, id = c("Status"), variable.name = "sampleID")
summary(vFcStats.filt.long)
head(vFcStats.filt.long)

# Restructure as data frame
vFcStats.filt.long.df <- data.frame(vFcStats.filt.long)
names(vFcStats.filt.long.df)
head(vFcStats.filt.long.df)


# Making quantitative bar graphs of featureCount Stats
x11()
ggplot(vFcStats.filt.long.df, aes(x = sampleID, y = value, fill = Status)) + geom_bar(stat="identity",colour="black", position = "stack") + xlab("") + ylab(expression("Number of reads")) + guides(fill=FALSE) + ggtitle("featureCount Stats") + facet_wrap(~ Status) +  theme(axis.text.x = element_text(angle=90, vjust=1), text = element_text(size=12), axis.text=element_text(colour="#525252")) 
#vjust adjust the vertical justification of the labels, which is often useful


x11()
ggplot(vFcStats.filt.long.df, aes(x = sampleID, y = value, fill = Status)) + geom_bar(stat="identity", colour="black", position = "stack") + xlab("") + ylab(expression("Number of reads")) + scale_fill_manual(values = c("#00ccff","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0")) + guides(fill=FALSE) + ggtitle(expression(atop("Quant featureCount stats for Hirahara-O'Shea STAT1 data", atop("Assigned, Ambiguity, No Features, and Unmapped"), ""))) +  theme(text = element_text(size=18), axis.text=element_text(colour="#525252", size=16), axis.title.y=element_text(angle=90), axis.text.x=element_text(angle=90)) 
#+ coord_flip()




# Save to PDF
pdf("aa_RNA-02-DE_NUM_fcStat quant stacked bar graphs.pdf",width = 14, height = 10, useDingbats = F);

ggplot(vFcStats.filt.long.df, aes(x = sampleID, y = value, fill = Status)) + geom_bar(stat="identity",colour="black", position = "stack") + xlab("") + ylab(expression("Number of reads")) + guides(fill=FALSE) + ggtitle("featureCount Stats") + facet_wrap(~ Status) +  theme(axis.text.x = element_text(angle=90, vjust=1), text = element_text(size=14), axis.text=element_text(colour="#525252")) 


ggplot(vFcStats.filt.long.df, aes(x = sampleID, y = value, fill = Status)) + geom_bar(stat="identity", colour="black", position = "stack") + xlab("") + ylab(expression("Number of reads")) + scale_fill_manual(values = c("#00ccff","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0")) + guides(fill=FALSE) + ggtitle(expression(atop("Quant featureCount stats for Hirahara-O'Shea STAT1 data", atop("Assigned, Ambiguity, No Features, and Unmapped"), ""))) +  theme(text = element_text(size=18), axis.text=element_text(colour="#525252", size=16), axis.title.y=element_text(angle=90), axis.text.x=element_text(angle=90)) 
#+ coord_flip()


dev.off ();


#####

# DATA ANALYSIS

# Import Targets file
vTargets <- readTargets(vFileTargets, sep = "\t", path = NULL, row.names = NULL)
# Check and confirm data import with head
head(vTargets, n = 3)
# Check and confirm subset structure (originally columns), which can be specifically called with the $ operator
names(vTargets)

#####

# Add label of interest to Targets file as columns - will work on later

# First label: Cytokine treatment
# vTargets$cytokine <- 0
# head(vTargets)

# Finding rows with specified treatment
# vTargets[grep("nc", vTargets$sampleID, value = F), ]

# Building loop to autofill rows
# if (grepl("_nc", vTargets$sampleID)) {
#    vTargets$cytokine <- as.character(vTargets$cytokine)
#    vTargets$cytokine[vTargets$cytokine == "0"] <- "nc"
#    vTargets$cytokine <- as.factor(vTargets$cytokine)
# }

#which("nc" %in% vTargets$sampleID)
#match("nc", vTargets$sampleID)

#vTargets[
#    grepl("_nc", vTargets$sampleID)


# FEATURECOUNTS COUNT DATA
# featureCounts count output in data frame
# Note that hyphens are converted to periods because R interprets hyphens as minus signs.
vFc <- read.delim(vFileFcCounts, stringsAsFactors = F)
# Check
names(vFc)
head(vFc, n = 3)

# Create a subset with only counts by setting rownames to the geneID and removing Geneid and Length columns
vFc.Counts <- vFc
rownames(vFc.Counts) <- vFc.Counts$Geneid
vFc.Counts$Geneid <- NULL
vFc.Counts$Length <- NULL
# Check
head(vFc.Counts, n = 3)


#####

# TMM and edgeR classic {edgeR}

## Calculating TMM normalization factors

# Filter for genes with > 0 mapped reads for all samples
vGenesToKeep <- rowSums(vFc.Counts) > 0
vFc.Counts.filt = vFc.Counts[vGenesToKeep, ]
vFc.filt = vFc[vGenesToKeep, ]

# Checks
dim(vFc.Counts)
dim(vFc.Counts.filt)
dim(vFc)
dim(vFc.filt)

# Create DGEList object from counts based on cytokine treatment
vGroupCyto <- vTargets$cytokine

vDGE.cyto <- DGEList(counts=vFc.Counts.filt, group = vGroupCyto, genes = vFc.filt[,c("Geneid","Length")])


# Calculate TMM normalization factors
vFcCyto.tmm <- calcNormFactors(vDGE.cyto, method = c("TMM"), refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)

# From http://davetang.org/muse/2012/01/19/the-dgelist-object-in-r/:

# Typical values for the common BCV (square-root-dispersion) for datasets arising from well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates.

vFcCyto.tmm <- estimateCommonDisp(vFcCyto.tmm, verbose = T) 
# Disp = 0.12098 , BCV = 0.3478

vFcCyto.tmm <- estimateTagwiseDisp(vFcCyto.tmm)
vFcCyto.tgw <- exactTest(vFcCyto.tmm)


# https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html and edgeR tutorial:

#   The function topTags() takes the output from exactTest(), adjusts the raw p-values using the False Discovery Rate (FDR) correction, and returns the top differentially expressed genes.  The output is similar to that of exactTest() but with a column of adjusted p-values and sorted by increasing p-value.  The sort.by argument allows you to sort the table by p-value, concentration or fold-change if desired.  We can also use topTags() to return the original counts of the top differentially expressed genes.  By setting the n parameter to the total number of genes, we can save the entire topTags() results table

topTags(vFcCyto.tgw, n = 20, sort.by = "p.value")

# Store full topTags results table
vFcCyto.tgw.TTTable <- topTags(vFcCyto.tgw, n = nrow (vFcCyto.tgw$table), sort.by = "p.value")$table

# Store topTags results in csv:
write.table(vFcCyto.tgw.TTTable, file="aa_RNA-03-DE_TopTagsTable.csv", append = FALSE, quote = TRUE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = NA, qmethod = c("escape"), fileEncoding = "")



# Default: decideTestsDGE(object, adjust.method="BH", p.value=0.05)
vFcCyto.tgw.bh <- decideTestsDGE(vFcCyto.tgw, adjust.method="BH", p.value=0.05)

summary(vFcCyto.tgw.bh)
#   [,1] 
# -1     3 = -1 for down-regulated genes
# 0  40039 = 0 for non-DE genes
# 1     15 = 1 for up-regulated

# differentially expressed tags from the naive method in d1
vFcCyto.tmm.de <- rownames(vFcCyto.tmm)[as.logical(vFcCyto.tgw.bh)] 

x11()
plotSmear(vFcCyto.tgw, de.tags=vFcCyto.tmm.de) 
# plotSmear is a more sophisticated and superior way to produce an 'MA plot'. plotSmear resolves the problem of plotting tags that have a total count of zero for one of the groups by adding the 'smear' of points at low A value. The points to be smeared are identified as being equal to the minimum estimated concentration in one of the two groups. The smear is created by using random uniform numbers of width smearWidth to the left of the minimum A. plotSmear also allows easy highlighting of differentially expressed (DE) tags. 
abline(h = c(-2, 2), col = "blue")
dev.print(png,filename="Data/2015-Hirahara_OShea-stat1/Images/14-TMM.png",width=800);
dev.off ();


pdf("aa_RNA-03-DE_NUM_vFcCyto.tgw.bh.pdf",width = 8, height = 6, useDingbats = F);
plotSmear(vFcCyto.tgw, de.tags=vFcCyto.tmm.de) 
abline(h = c(-2, 2), col = "blue")
dev.off();



# Plotting the tagwise dispersion calculated from this
names(vFcCyto.tmm)
vFc.Counts.filt$twd <- vFcCyto.tmm$tagwise.dispersion
head(vFc.Counts.filt)

pdf("aa_RNA-03-DE_NUM_hist-vFcCyto.tmm.tgw.pdf",width = 8, height = 6, useDingbats = F);
hist(vFc.Counts.filt$twd, breaks=20, xlim=c(0,3))
dev.off()

#now if you wanted to save the fold changes
names(vFcCyto.tgw)
head(vFcCyto.tgw$table)
head(vFcCyto.tgw$genes)
vFc.Counts.filt.changes <- cbind(vFc.Counts.filt, vFcCyto.tgw$table)
head(vFc.Counts.filt.changes)

#I want the fdr
vFc.Counts.filt.changes$PValue_fdr <- p.adjust(method="fdr",p=vFc.Counts.filt.changes$PValue)
head(vFc.Counts.filt.changes)

#sanity check with the decideTestsDGE() call
table(vFc.Counts.filt.changes$PValue_fdr<0.05) # TRUE should be sum of -1 and 1 from summary(de).  Works here: 3 + 15 = 18

table(vFc.Counts.filt.changes$PValue_fdr<0.01)


write.table(vFc.Counts.filt.changes, file="aa_RNA-03-DE_NUM_edgeR-TMM-DE.csv", append = FALSE, quote = TRUE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = NA, qmethod = c("escape"), fileEncoding = "")


# Trying to make plot with top 20 hits labeled with gene ID, from https://support.bioconductor.org/p/56590/:
de.name <- rownames(topTags(vFcCyto.tgw, n = 20)$table)

x11()
plotSmear(vFcCyto.tgw, de.tags=de.name, cex = 0.5)
vFcCyto.tgw$genes <- rownames(vFcCyto.tmm$counts)

ids <- c("ENSG00000110324.5_IL10RA", "ENSG00000026751.12_SLAMF7"
, "ENSG00000181847.7_TIGIT")

gene.labels <- vFcCyto.tgw$table[vFcCyto.tgw$genes %in% ids,]
abline(h = c(-2, 2), col = "blue")
points(x=gene.labels$logCPM, y=gene.labels$logFC, cex=2, col="red")
text(x = gene.labels$logCPM, y = gene.labels$logFC, labels = rownames(gene.labels), cex = 1.5, pos = 2)




pdf("aa_RNA-03-DE_NUM_MA-vFcCyto.tmm.tgw.pdf",width = 8, height = 6, useDingbats = F);
plotSmear(vFcCyto.tgw, de.tags=de.name, cex = 0.5)
vFcCyto.tgw$genes <- rownames(vFcCyto.tmm$counts)

ids <- c("ENSG00000110324.5_IL10RA", "ENSG00000026751.12_SLAMF7"
         , "ENSG00000181847.7_TIGIT")

gene.labels <- vFcCyto.tgw$table[vFcCyto.tgw$genes %in% ids,]
abline(h = c(-2, 2), col = "blue")
points(x=gene.labels$logCPM, y=gene.labels$logFC, cex=2, col="red")
text(x = gene.labels$logCPM, y = gene.labels$logFC, labels = rownames(gene.labels), cex = 1.5, pos = 2)

dev.off()
# Yay!


# Looking up genes on ensembl.org and saving in spreadsheet.


# Building volcano plot!
head(fc.tgw$table)

tmm.table = data.frame(logFC = fc.tgw$table[, 1], negLogPval = -log10(fc.tgw$table[, 3]))

head(tmm.table)

x11()
par(mar = c(5, 4, 4, 4))
plot(tmm.table, pch = 16, cex  = 0.6, xlab = expression (log[2]~fold~change), ylab = expression(~log[10]~pvalue))

## Set log-fold change and p-value cutoffs
lfc = 2
pval = 0.01

## Selecting interest genes
sigGenes = (abs(tmm.table$logFC) > lfc & tmm.table$negLogPval > -log10(pval))

## Identifying selected genes
points(tmm.table[sigGenes, ], pch = 16, cex = 0.8, col = "red")
abline(h = -log10(pval), col = "blue", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)

text(x = gene.labels$logCPM, y = gene.labels$logFC, labels = rownames(gene.labels), cex = 1, pos = 1)

dev.print(png,filename="Data/2015-Hirahara_OShea-stat1/Images/20-TMM volcano.png",width=800);
dev.off ();

# Can learn how to build volcano plot with ggplot here: http://bioinformatics.knowledgeblog.org/2011/06/21/volcano-plots-of-microarray-data/


#####

# VOOM {limma}

# Store count data by sample ID
vSampleID <- factor(vTargets$sampleID)

# Build matrix of experiment design by group
vDesign <- model.matrix(~vSampleID)


# Transfer counts and annotations to DGE list {edgeR}.  See http://davetang.org/muse/2012/01/19/the-dgelist-object-in-r/
vDCounts <- DGEList(counts=vFc.Counts, genes=vFc[,c("Geneid","Length")])

# Subset for genes with > 10 mapped reads/million in at least 2 libraries - NEED TO FIX THIS IS AFTER RPKM WHICH IM NOT DOING HERE
vGenesToKeep <- rowSums(cpm(vDCounts) > 10) >= 2
vDCounts.f <- vDCounts[vGenesToKeep,]

# Normalizing with TMM {edgeR}
vDCounts.f.tmm <- calcNormFactors(vDCounts.f, method = c("TMM"), refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)

# Re-normalizing and accounting for mean-variance relationship with voom {limma}
    # Note that voom also has an internal calculation for normalizing to library size, so I'm not sure why it's standard practice to also run TMM?
    # Saving results in external PDF (which is a vector file format)
pdf("bb_RNA-02-DE_NUM-filt-TMM-voom.pdf", width=8, height=6, useDingbats=FALSE)
vVoom <- voom(vDCounts.f.tmm, vDesign, plot = T)
dev.off()

# Sample clustering by MDS
pdf("bb_RNA-02-DE_NUM-MDS-voom.pdf", width=8, height=6, useDingbats=FALSE)
plotMDS(vVoom,xlim=c(-2.5,2.5))
dev.off()

# Fitting to linear model and assessing differential expression using eBayes moderated t statistic.
vVoom.eBayes <- eBayes(lmFit(vVoom, vDesign))
vVoom.eBayes.top <- topTable(vVoom.eBayes, coef=2)


write.table(vVoom.eBayes.top,
            file = "bb_RNA-02-DE_NUM-topTable-voom-eBayes.txt",
            sep = "\t", 
            quote=FALSE,
            row.names = FALSE)



# Open tabs
    # http://bioinf.wehi.edu.au/RNAseqCaseStudy/
    # http://davetang.org/muse/2012/01/19/the-dgelist-object-in-r/
    # https://www.biostars.org/p/71983/
    # https://rstudio-pubs-static.s3.amazonaws.com/85101_ad12a157e52d4fbb8361219c0c8b50c1.html

#####

sink(vFileSessionInfo, append = F, type = "output", split = F)
sessionInfo()
.Options
sink()