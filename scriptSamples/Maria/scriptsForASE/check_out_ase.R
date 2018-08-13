wd <- "/Users/mgutierr/Documents/work/rnaseq/cd4timelinepilot/tofiPub/ASE"
setwd(wd)

refRatio <- read.table("summary_refRatio.txt", header = T, sep = "\t")
dim(refRatio)
head(refRatio)
tail(refRatio)
cbind(colnames(refRatio))
colnames(refRatio) <- c("SNP_ID", "noTofi", "tofi0.1", "tofi0.3")
rownames(refRatio) <- refRatio[,1]
refRatio <- refRatio[,c(2,3,4)]
head(refRatio)

bas <- read.table("summary_both_alleles_seen.txt", header = T, sep = "\t")
cbind(colnames(bas))
colnames(bas) <- c("SNP_ID", "noTofi", "tofi0.1", "tofi0.3")
rownames(bas) <- bas[,1]
bas <- bas[,c(2,3,4)]
head(bas)

cov <- read.table("summary_total_count.txt", header = T, sep = "\t")
cbind(colnames(cov))
colnames(cov) <- c("SNP_ID", "noTofi", "tofi0.1", "tofi0.3")
rownames(cov) <- cov[,1]
cov <- cov[,c(2,3,4)]
head(cov)

### Remove SNPs with mapping problems
badSnps <- scan("snps_withMappingProblems_tofiPub.txt", what = "character")
posBadSnps <- which(rownames(cov) %in% badSnps)

cov <- cov[-posBadSnps,]
refRatio <- refRatio[-posBadSnps,]
bas <- bas[-posBadSnps,]
dim(refRatio)
dim(cov)
dim(refRatio)
dim(bas)

#EAS2

eas2 <- apply(refRatio, 1, function(snp){
            numAlt <- length(which(snp < 1))
            numRef <- length(which(snp > 0))
            if(numAlt >= 2 & numRef >= 2){
              return(1)
            }
            else{
              return(0)
            }
})
head(eas2)
head(refRatio)
posEas2 <- which(eas2 == 1)
head(posEas2)
#cov <- cov[posEas2,]
#refRatio <- refRatio[posEas2,]
#bas <- bas[posEas2,]

#dim(refRatio)
#dim(cov)
#dim(bas)



#Ditribution of refRatio for each sample (all sites, and BAS sites)
i <- 2
pdf("hists_refRatio.pdf", width = 10, height = 4)
#x11(width = 10, height = 4)
par(mfrow = c(1,2))
for(i in 1:ncol(refRatio)){
  name <- colnames(refRatio)[i]
  hist(refRatio[,i], breaks = 20, main = name, xlab = "Reference Ratio", cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, col = "gray")
  abline(v = 0.5, col = "red", lwd = 2)
  hist(refRatio[which(bas[,i] == 1),i], breaks = 20, main = paste(name, "(BAS)"), xlab = "Reference Ratio", cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, col = "gray")
  abline(v = 0.5, col = "red", lwd = 2)
}
dev.off()

#Boxplot dist to 0.5
i <- 1
pdf("boxplot_distToBalance.pdf", width = 6, height = 8)
par(mfrow = c(2,1))
l <- list()
for(i in 1:ncol(refRatio)){
  name <- colnames(refRatio)[i]
  l[[name]] <- abs(0.5 - refRatio[,i])
}
boxplot(l, ylab = "Distance to 0.5", cex.lab = 1.5, cex.axis = 1.5, varwidth = TRUE, col = "gray")

l <- list()
for(i in 1:ncol(refRatio)){
  name <- colnames(refRatio)[i]
  l[[name]] <- abs(0.5 - refRatio[which(bas[,i] == 1),i])
}
boxplot(l, ylab = "Distance to 0.5", cex.lab = 1.5, cex.axis = 1.5, main = "BAS", varwidth = TRUE, col = "gray")
dev.off()


### Coverage per site

pdf("boxplot_covPerSite.pdf", width = 10, height = 8)
par(mfrow = c(2,1), mar = c(6,6,3,1))
l <- list()
for(i in 1:ncol(refRatio)){
  name <- colnames(refRatio)[i]
  l[[name]] <- cov[,i]
  name2 <- paste(name, "BAS", sep = "-")
  l[[name2]] <- cov[which(bas[,i] == 1),i]
}

boxplot(l, ylab = "", cex.lab = 1.3, cex.axis = 1.3, varwidth = TRUE, las = 2, col = "gray")
boxplot(l, ylab = "Coverage per Site", cex.lab = 1.3, cex.axis = 1.3, varwidth = TRUE, ylim = c(0, 70), las = 2, col = "gray")
dev.off()

medians <- lapply(l, function(tmpt){
  median(na.omit(tmpt))
})
class(medians)
names(medians)
head(medians)
meds <- c(medians[[1]][1],medians[[3]][1],medians[[5]][1])
head(meds)

### Proportion of BAS

pdf("barplot_propBAS.pdf", width = 10, height = 8)
par(mfrow = c(2,1), mar = c(6,6,3,1))
l <- NULL
n <- NULL
for(i in 1:ncol(refRatio)){
  n <- c(n, colnames(refRatio)[i])
  l <- c(l, length(which(bas[,i] == 1))/length(na.omit(bas[,i])))
}
barplot(l, names = n, ylim = c(0,1), ylab = "Fraction of BAS", cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.5, col = "gray")
abline(h = 0.9)
abline(h = 1)
dev.off()

l
x11(width = 5.5, height = 5)
plot(meds, l, pch = 20, cex = 1.5, xlab = "Median Coverage Per Site", ylab = "Proportion of BAS", cex.lab = 1.5, cex.axis = 1.5)
dev.print("scatterplot_medianCov_by_propBAS.pdf", dev = pdf)
