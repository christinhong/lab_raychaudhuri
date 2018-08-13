library(limma)
library(edgeR)
library(stringr)

wd <- "/Users/mgutierr/Documents/work/rnaseq/cd4timelinepilot/tofiPub/geneExpression/diffExpr"
suffix <- "_tofiPubAllDonors"
setwd(wd)

counts <- read.table("../gene_counts_tofiPubAllDonors_oneChrStrandMinStartMaxEnd.txt", head = T, stringsAsFactors = FALSE)
cov <- read.csv("master_covariate.csv", head = T, stringsAsFactors = FALSE)
head(cov)

des <- cov
head(des)

head(counts)
summary(rowMeans(counts[,7:ncol(counts)]))

shortN <- sapply(colnames(counts)[7:ncol(counts)], function(sample){
  cov$ShortName[which(cov$Library %in% sample)]
})
shortN

colnames(counts)[7:ncol(counts)] <- shortN
head(counts)
#read gene id-name match
match <- read.table("/Users/mgutierr/Documents/work/rnaseq/access/diffExpr/gene_ID_name_match.txt", head = T, stringsAsFactors = F)
head(match)
dim(match)

# generate RPKM values if you need them
x_rpkm <- rpkm(counts[,7:ncol(counts)],counts$Length)
rownames(x_rpkm) <- counts$Geneid
head(x_rpkm)
hist(x_rpkm[which(x_rpkm[,1] > 1),1], breaks = 100)
ln <- log1p(x_rpkm)# ln( 1 + rpkm)
hist(ln[which(ln[,1] > 1),1], breaks = 100)
abline(v = ln[which(rownames(ln) %in% "ENSG00000049768.10"),1], col = "red", lwd = 4)
boxplot()
density()
plot(runif(100), runif(100), pch = 20, cex = 2, cex.axis = 2)
?par
barplot(ln[which(rownames(ln) %in% "ENSG00000049768.10"),])

x11(width = 9, height = 7)
barplot(ln[which(rownames(ln) %in% "ENSG00000049768.10"),], las = 2)


#write.table(x_rpkm, file = "rpkm_effectiveLength_access_sampleNames.txt", sep = "\t", quote = FALSE)
#summary(counts$Length)

# generate DGElist with counts
x <- DGEList(counts=counts[,7:ncol(counts)], genes=counts[,c("Geneid","Length")])
isexpr <- rowSums(cpm(x) > 10) >= 2
dim(x)
x <- x[isexpr,]
dim(x)

#perform tmm norm (I think norm factors coded separately, and eventually are multiplied by expr values for each sample)
x.tmm <- calcNormFactors(x)
head(x.tmm)
attributes(x.tmm)
x.tmm$samples
summary(x.tmm$samples$norm.factors)

# generate FPKM values with effective lib size and gene length
x_fpkm <- rpkm(x.tmm)
head(x_fpkm)
class(x_fpkm)
dim(x_fpkm)
x_fpkm2 <- cbind(x.tmm$genes, x_fpkm)
head(x_fpkm2)
dim(x_fpkm2)
mat <- merge(match, x_fpkm2, by = 1)
dim(mat)
head(mat)
colnames(mat)[4] <- "Effective_Length"
head(mat)

write.table(mat, file = "gene_FPKMs_TMMnorm_effectiveLength_exprCMP10in2_tofiAllDonors_valuesUsedForLimma.txt", row.names = F, quote = F, sep = "\t")
#head(x_rpkm)
#write.table(x_rpkm, file = "rpkm_effectiveLength_access_sampleNames.txt", sep = "\t", quote = FALSE)
#summary(counts$Length)

#create a desing matrix
head(des)
pos <- sapply(colnames(x.tmm), function(s){
  which(des[,4] %in% s)
})
pos
des[12,]
des.ord <- des[pos,]
head(des.ord)

celltype <- factor(des.ord[,"Tofi"])
design <- model.matrix(~0+celltype)
design

#x11()
y <- voom(x.tmm,design,plot=TRUE)

dev.print("plot_voom_meanVarianceTrend_tmmNorm_voom_desTofi.pdf", dev = pdf)
#dev.off()
# cluster libraries
#x11()
plotMDS(y)
dev.print("mds_default_tofiAllDonors_tmmVoom_desTofi.pdf", dev = pdf)
#head(des[,"CellType"])

#x11()
#plotMDS(y,labels =  des[pos,2] ,top = 50, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 50 genes") #, prior.count = 5)
#plotMDS(y,labels =  des[pos,2] ,top = 100, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 100 genes")
#plotMDS(y,labels =  des[pos,2] ,top = 200, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 200 genes")
#plotMDS(y,labels =  des[pos,2] ,top = 300, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 300 genes")
#plotMDS(y,labels =  des[pos,2] ,top = 500, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 500 genes") #, prior.count = 5)
#plotMDS(y,labels =  des[pos,2] ,top = 700, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 700 genes") #, prior.count = 5)

#colors <- rep("black", length(pos))
#l <- levels(factor(des[pos,"RunID"]))
#colors[which(des[pos,"RunID"] %in% l[1])] <- "red"
#colors[which(des[pos,"RunID"] %in% l[2])] <- "blue"
#colors[which(des[pos,"RunID"] %in% l[3])] <- "forestgreen"
#pdf("mdsPlots_tmm_voom_norm_counts_colByRun_CORRECTEDnumTop.pdf")
#plotMDS(y,labels =  rownames(des)[pos], col = colors ,top = 700, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 700 genes")
#plotMDS(y,labels =  des[pos,"CellType"], col = colors ,top = 700, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 700 genes")
#plotMDS(y,labels =  rownames(des)[pos], col = colors ,top = 500, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 500 genes")
#plotMDS(y,labels =  des[pos,"CellType"], col = colors ,top = 500, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 500 genes")
#plotMDS(y,labels =  rownames(des)[pos], col = colors ,top = 300, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 300 genes")
#plotMDS(y,labels =  des[pos,"CellType"], col = colors ,top = 300, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 300 genes")
#plotMDS(y,labels =  rownames(des)[pos], col = colors ,top = 100, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 100 genes")
#plotMDS(y,labels =  des[pos,"CellType"], col = colors ,top = 100, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 100 genes")

#dev.off()

#pdf("mdsPlots_counts_tmm_voom_isexpr_top500.pdf", width = 12, height = 4)
#par(mfrow = c(1,3))
#plotMDS(x,labels =  des[pos,"CellType"], col = colors ,top = 500, gene.selection="common", cex = 1.2, main = "counts\nMDS on top 500 genes")
#plotMDS(x.tmm,labels =  des[pos,"CellType"], col = colors ,top = 500, gene.selection="common", cex = 1.2, main = "TMM norm counts\nMDS on top 500 genes")
#plotMDS(y,labels =  des[pos,"CellType"], col = colors ,top = 500, gene.selection="common", cex = 1.2, main = "TMM and voom norm counts\nMDS on top 500 genes")
#dev.off()

#x[1:6,1:6]
#x.tmm[1:6,1:6]
#y[1:6,1:6]

fit <- lmFit(y,design)
summary(fit)

#####
pairs <- NULL
for(i in 1:(length(colnames(design))-1)){
  pairs <- c(pairs, paste(colnames(design)[i], colnames(design)[(i+1):ncol(design)], sep = "-"))
}
pairs
print(paste(pairs, collapse = " , "))
contrast.matrix <- makeContrasts(celltype0-celltype1, levels = design)
contrast.matrix
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
attributes(fit2)

#celltypeCD4NKT - celltypeNK
tab <- topTable(fit2, coef = 1, adjust = "fdr", number = 10922)
head(tab)
head(fit2$coefficients)
tab2 <- merge(tab, match, by = 1)
head(tab2)
write.table(tab2, "diffExpr_tofi_vs_noTofi_allDonors_11KexprGenes.tab", sep = "\t", quote = FALSE, row.names = FALSE )
tab2[grep("CD4", tab2[,"GENE_NAME"]),]
tab2[grep("IFN", tab2[,"GENE_NAME"]),]
tab2[grep("IL2RA", tab2[,"GENE_NAME"]),]
tab2[grep("JAK", tab2[,"GENE_NAME"]),]
tab2[grep("TYK2", tab2[,"GENE_NAME"]),]

de <- tab2[which(tab2$adj.P.Val < 0.01),]

dim(de)

summary(-log10(de$adj.P.Val))

de$GENE_NAME[which(-log10(de$adj.P.Val) > 3)]




#gene id CD4 is ENSG00000010610.5

x.tmm[which(x.tmm$genes[,1] %in% "ENSG00000010610.5"),]
head(x.tmm)

cd4 <- (as.numeric(x.tmm$counts[which(x.tmm$genes[,1] %in% "ENSG00000010610.5"),])/(x.tmm$samples[,2]*x.tmm$samples[,3]))*1000000
#cd4 <- (as.numeric(x.tmm$counts[which(x.tmm$genes[,1] %in% "ENSG00000010610.5"),])/(x.tmm$samples[,2]*x.tmm$samples[,3]*1033))*1000000000
#gene id CD8A is ENSG00000153563.11, CD8B ENSG00000172116.17

cd8A <- (as.numeric(x.tmm$counts[which(x.tmm$genes[,1] %in% "ENSG00000153563.11"),])/(x.tmm$samples[,2]*x.tmm$samples[,3]))*1000000
cd8B <- (as.numeric(x.tmm$counts[which(x.tmm$genes[,1] %in% "ENSG00000172116.17"),])/(x.tmm$samples[,2]*x.tmm$samples[,3]))*1000000

#gene id CD3G "ENSG00000160654.5", CD3D "ENSG00000167286.5", CD3E "ENSG00000198851.5"
cd3G <- (as.numeric(x.tmm$counts[which(x.tmm$genes[,1] %in% "ENSG00000160654.5"),])/(x.tmm$samples[,2]*x.tmm$samples[,3]))*1000000
cd3D <- (as.numeric(x.tmm$counts[which(x.tmm$genes[,1] %in% "ENSG00000167286.5"),])/(x.tmm$samples[,2]*x.tmm$samples[,3]))*1000000
cd3E <- (as.numeric(x.tmm$counts[which(x.tmm$genes[,1] %in% "ENSG00000198851.5"),])/(x.tmm$samples[,2]*x.tmm$samples[,3]))*1000000

#gene id for PLZF is ENSG00000109906.9 (protein name is ZBTB16)
plzf <- (as.numeric(x.tmm$counts[which(x.tmm$genes[,1] %in% "ENSG00000109906.9"),])/(x.tmm$samples[,2]*x.tmm$samples[,3]))*1000000

#gene id for CEBPD is ENSG00000221869.4
cebpd <- (as.numeric(x.tmm$counts[which(x.tmm$genes[,1] %in% "ENSG00000221869.4"),])/(x.tmm$samples[,2]*x.tmm$samples[,3]))*1000000

#gene id for tbx21 ENSG00000073861.2
tbx21 <- (as.numeric(x.tmm$counts[which(x.tmm$genes[,1] %in% "ENSG00000073861.2"),])/(x.tmm$samples[,2]*x.tmm$samples[,3]))*1000000

#DELETED CODE FOR LOG2 EXPR BECAUSE IT WAS WRONG

###########
results <- decideTests(fit2)
rownames(results) <- fit2$genes[,1]
#check in which contrasts PLZF is diff expr

head(results)
attributes(results)
dim(results)
attributes(fit2)

resPLZF <- results["ENSG00000109906.9",]
resPLZF <- resPLZF[-which(resPLZF == 0)]
resPLZF

resCD4 <- results["ENSG00000010610.5",]
resCD4 <- resCD4[-which(resCD4 == 0)]
resCD4

resCD8A <- results["ENSG00000153563.11",]
resCD8A <- resCD8A[-which(resCD8A == 0)]
resCD8A
length(resCD8A[grep("celltypeCD8 -", names(resCD8A))])
length(resCD8A[grep("celltypeCD8NKT -", names(resCD8A))])
names(resCD8A)

resCEBPD <- results["ENSG00000221869.4",]
resCEBPD <- resCEBPD[-which(resCEBPD == 0)]
cbind(names(resCEBPD))

resTBX21 <- results["ENSG00000073861.2",]
resTBX21 <- resTBX21[-which(resTBX21 == 0)]
cbind(names(resTBX21))


####### reapeat bar plot with gene expr levels for exampel genes but without log scale
x11(width = 7, height = 10)
pdf("barplots_exprLevel_CD4_CD8_CD3_PLZF_allSamples_NOTLOG.pdf", width = 7, height = 10)
par(mar = c(5,9,4,2))
col <- rep("gray", ncol(x.tmm$counts))
col[grep("CD8", colnames(x.tmm$counts))] <- "darkgoldenrod1"
pos <- c(1,2,20,21,3,4,5,6,7,19,14,15,16,17,18,10,11,12,13,24,25,22,23,8,9)
length(pos)
colnames(x.tmm$counts)[8:9] <- c("X16.500NK", "X17.500NK")
barplot(cd8A[pos], names = colnames(x.tmm$counts)[pos], las = 1, xlab = "counts per million (TMM)", main = "CD8A", cex.axis = 1.5, cex.lab = 1.5, col = col[pos], horiz = T, cex.name = 1.2)
barplot(cd8B[pos], names = colnames(x.tmm$counts)[pos], las = 1, xlab = "counts per million (TMM)", main = "CD8B", cex.axis = 1.5, cex.lab = 1.5, col = col[pos], horiz = T, cex.name = 1.2)

col <- rep("gray", ncol(x.tmm$counts))
col[grep("CD4", colnames(x.tmm$counts))] <- "darkgoldenrod1"
barplot(cd4[pos], names = colnames(x.tmm$counts)[pos], las = 1, xlab = "counts per million (TMM)", main = "CD4", cex.axis = 1.5, cex.lab = 1.5, col = col[pos], horiz = T, cex.name = 1.2)

#CD3 should be pos in MAIT, gammas and NKT (not in NKs)
col <- rep("gray", ncol(x.tmm$counts))
col[grep("MAIT", colnames(x.tmm$counts))] <- "darkgoldenrod1"
col[grep("gamma", colnames(x.tmm$counts))] <- "darkgoldenrod1"
col[grep("CD4.NKT", colnames(x.tmm$counts))] <- "darkgoldenrod1"
col[grep("CD8.NKT", colnames(x.tmm$counts))] <- "darkgoldenrod1"
col[grep("dn.NKT", colnames(x.tmm$counts))] <- "darkgoldenrod1"
barplot(cd3G[pos], names = colnames(x.tmm$counts)[pos], las = 1, xlab = "counts per million (TMM)", main = "CD3G", cex.axis = 1.5, cex.lab = 1.5, col = col[pos], horiz = T, cex.name = 1.2)
barplot(cd3D[pos], names = colnames(x.tmm$counts)[pos], las = 1, xlab = "counts per million (TMM)", main = "CD3D", cex.axis = 1.5, cex.lab = 1.5, col = col[pos], horiz = T, cex.name = 1.2)
barplot(cd3E[pos], names = colnames(x.tmm$counts)[pos], las = 1, xlab = "counts per million (TMM)", main = "CD3E", cex.axis = 1.5, cex.lab = 1.5, col = col[pos], horiz = T, cex.name = 1.2)

#PLZF
col <- rep("gray", ncol(x.tmm$counts))
barplot(plzf[pos], names = colnames(x.tmm$counts)[pos], las = 1, xlab = "counts per million (TMM)", main = "ZBTB16 (PLZF)", cex.axis = 1.5, cex.lab = 1.5, col = col[pos], horiz = T, cex.name = 1.2)

dev.off()

####more genes
x11(width = 7, height = 10)
pdf("barplots_exprLevel_moreGenes_allSamples_NOTLOG.pdf", width = 7, height = 10)
par(mar = c(5,9,4,2))

col <- rep("gray", ncol(x.tmm$counts))
barplot(cebpd[pos], names = colnames(x.tmm$counts)[pos], las = 1, xlab = "counts per million (TMM)", main = "CEBPD", cex.axis = 1.5, cex.lab = 1.5, col = col[pos], horiz = T, cex.name = 1.2)

#tbx21
barplot(tbx21[pos], names = colnames(x.tmm$counts)[pos], las = 1, xlab = "counts per million (TMM)", main = "TBX21", cex.axis = 1.5, cex.lab = 1.5, col = col[pos], horiz = T, cex.name = 1.2)

dev.off()

