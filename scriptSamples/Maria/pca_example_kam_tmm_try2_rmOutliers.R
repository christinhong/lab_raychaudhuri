wd <- "/Users/mgutierr/Documents/work/rnaseq/access/diffExpr"
setwd(wd)

# Load convenience functions for PCA.
source("pca_kam.R")

# Dependencies.
library(ggplot2)
library(limma)
library(edgeR)
library(stringr)
# Let's start with two dataframes:
# 1. A dataframe called "dat" with rows = genes and columns = samples.
# 2. A dataframe called "meta" with covariates. It has one row for each sample.

counts <- read.table("merged_gene_counts_access_sampleNames_EFFECTIVE_LENGTH_RcolNames.txt", header = T, stringsAsFactors = F)
#x <- counts[,7:ncol(counts)]
#dim(x)
#dat <- cpm(x)

head(counts)
dim(counts)
colnames(counts)[c(7,22 )]
counts2 <- counts[,-c(7,22 )]
dim(counts2)
x <- DGEList(counts=counts2[,7:ncol(counts2)], genes=counts2[,c("Geneid","Length")])
isexpr <- rowSums(cpm(x) > 10) >= 2
length(which(isexpr))
dim(x)
x <- x[isexpr,]
dim(x)

#perform tmm norm (norm factors coded separately, and eventually they are multiplied by library size for an "effective lib size")
x.tmm <- calcNormFactors(x)

head(x.tmm$samples)
head(x.tmm)

#dat <- cpm(x.tmm)
#dat2 <- cpm(counts[isexpr,7:ncol(counts)])

#create dat matrix with log2(cpm) where cpm is calculating using TMM and summing 0.5 to counts
dat <- apply(x.tmm$counts, 1, function(gene){
      log2(((0.5 + gene)/(x.tmm$samples[,2]*x.tmm$samples[,3]))*1000000)
})

dat <- t(dat)



dat[1:4,1:4]
dat2[1:4,1:4]
#dat <- apply(x.tmm$counts,1,function(sample){
#  return(sample*x.tmm$samples[,3])
#})
#dat <- t(dat)
head(dat)
dim(dat)

cov <- read.table("/Users/mgutierr/Documents/work/rnaseq/access/covariate/master_covariate_25samples_RNA_access.txt", head = T, stringsAsFactors = FALSE)
head(cov)
rname <- sapply(cov$SampleID, function(s){
  paste("X", paste(strsplit(s, "-")[[1]], collapse = "."), sep = "")  
})
head(rname)
head(cov)
des <- cbind(rname, cov$CellType, cov$RunID, cov$NumCells, cov$percUniqMapped, cov$BAyield_ng_per_ul, cov$totalReadPairs, cov$GCpercent_R1andR2m, cov$TopOverrepresentedSeqPercent_R1andR2m, cov$FracReadsQ20less_R1andR2m, cov$InsertSizeMedian, cov$percMap2Probes)
head(des)
colnames(des) <- c("SampleID", "CellType", "RunID", "NumCells", "PercUniqMapped", "BAyield_ng_per_ul", "totalReadPairs", "GCpercent_R1andR2m", "TopOverrepresentedSeqPercent_R1andR2m", "FracReadsQ20less_R1andR2m", "InsertSizeMedian", "percMap2Probes")

pos <- sapply(colnames(dat), function(s){
  which(des[,1] %in% s)
})
pos
des[12,]
des.ord <- des[pos,]
head(des.ord)

meta <- as.data.frame(des.ord)
head(meta)
class(meta$GCpercent_R1andR2m)
for(i in 4:ncol(meta)){
  meta[,i] <- as.numeric(as.character(meta[,i]))
}

# Typically, one of my first steps is to look at the distribution of standard deviation.
row_sd <- apply(dat, 1, sd)
plot(density(row_sd))
hist(row_sd)

#or first select expressed genes and then look at sd
#isexpr <- rowSums(dat > 10) >= 2
#sum(isexpr)
#dat2 <- dat[isexpr,]
#row_sd <- apply(dat2, 1, sd)
#plot(density(row_sd))
#hist(row_sd[which(row_sd < 200)])


# It is useful look at PCA for genes with high variability.
#idx_sd <- row_sd > 4
idx_sd <- row_sd > 3
#idx_sd <- row_sd > 100
#idx_sd <- row_sd > 10
sum(idx_sd)

# I find that centering is necessary for expression data, but scaling is not always desired.
pca1 <- prcomp(dat[idx_sd, ], scale = FALSE, center = TRUE)

# Always check to see if the principal components explain a large amount of variation.
# Remember that PCA with fewer genes will always explain more variation than PCA with more genes.
x11()
dev.off()
plot_variance_explained(pca1) + theme_bw(base_size = 18)
#dev.print("barplot_varExplained_PCs_top581varGenes_log2TMMs_05offset.pdf", dev = pdf)
#dev.print("barplot_varExplained_PCs_top2215varGenes_log2TMMs_05offset.pdf", dev = pdf)
dev.print("barplot_varExplained_PCs_top2031varGenes_log2TMMs_05offset_rm100T27NKT.pdf", dev = pdf)
# We'll use this to annotate the scatter plots with variance explained.
pca1_var_labs <- sprintf(
  "PC%s, %.02g%% of variance",
  variance_explained(pca1)$Component,
  100 * variance_explained(pca1)$Variance
)

# We'll use this to annotate the plot with covariates.
pca1_r <- cbind(pca1$rotation, meta)

# For exploratory analysis, we should correlate covariates with principal components.
# This will help us decide which components to plot, and how to color points.
correlate_pcs(pca1, meta)

# I find that exploratory plotting is easiest with ggplot2.

pca1_r <- as.data.frame(pca1_r)
head(pca1_r)
class(pca1_r$PC1)
#for(i in 1:25){
#  pca1_r[,i] <- as.numeric(as.character(pca1_r[,i]))
#}

x11(width = 9, height = 7)
ggplot(data = pca1_r) +
  geom_point(size = 5, aes(PC1, PC2, color = CellType)) +
  geom_text(size = 4, hjust = 0, vjust = 0, aes(PC1, PC2, label = NumCells)) +
  theme_bw(base_size = 18) +
  labs(x = pca1_var_labs[1],
       y = pca1_var_labs[2])

sum(idx_sd)

#dev.print("PCA_top581varGenes_log2TMMs_05offset.pdf", dev = pdf)
dev.print("PCA_top2031varGenes_log2TMMs_05offset_rm100T27NKT.pdf", dev = pdf)
dev.off()
