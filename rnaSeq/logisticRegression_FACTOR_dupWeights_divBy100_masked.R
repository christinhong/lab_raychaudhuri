library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(qvalue)
library(aod)

wd <- "/Users/mgutierr/Documents/work/rnaseq/cd4timelinepilot/tofiPub/ASE/dupWeights"
setwd(wd)

refRatio <- read.table("summary_refRatio_bas.txt", header = T, sep = "\t")
dim(refRatio)
head(refRatio)
tail(refRatio)
cbind(colnames(refRatio))
colnames(refRatio) <- c("SNP_ID", "noTofi", "tofi0.1", "tofi0.3")
rownames(refRatio) <- refRatio[,1]
refRatio <- refRatio[,c(2,3,4)]
head(refRatio)

cov <- read.table("summary_coverage_bas.txt", header = T, sep = "\t")
cbind(colnames(cov))
colnames(cov) <- c("SNP_ID", "noTofi", "tofi0.1", "tofi0.3")
rownames(cov) <- cov[,1]
cov <- cov[,c(2,3,4)]
head(cov)

###ref counts
ref <- read.table("summary_ref_counts_bas.txt", header = T, sep = "\t")
cbind(colnames(ref))
colnames(ref) <- c("SNP_ID", "noTofi", "tofi0.1", "tofi0.3")
rownames(ref) <- ref[,1]
ref <- ref[,c(2,3,4)]
head(ref)
tail(ref)

###alternative counts (nonref)
alt <- read.table("summary_alt_counts_bas.txt", header = T, sep = "\t")
cbind(colnames(alt))
colnames(alt) <- c("SNP_ID", "noTofi", "tofi0.1", "tofi0.3")
rownames(alt) <- alt[,1]
alt <- alt[,c(2,3,4)]
head(alt)
tail(alt)

head(ref)




### Remove SNPs with mapping problems 
badSnps <- scan("../snps_withMappingProblems_tofiPub.txt", what = "character")

posBadSnps <- which(rownames(cov) %in% badSnps)

cv <- cov[-posBadSnps,]
rr <- refRatio[-posBadSnps,]
rc <- ref[-posBadSnps,]
ac <- alt[-posBadSnps,]
dim(rr)
dim(cv)
dim(rr)
dim(rc)
dim(ac)



#take sites that have been assayed in all samples? (or 4 out 6?)
nas <- apply(rr, 1, function(snp){
  length(which(is.na(snp)))
})

head(nas)
summary(nas)
length(which(nas == 0))

rr.all <- rr[which(nas == 0),]
cv.all <- cv[which(nas == 0),]
rc.all <- rc[which(nas == 0),]
ac.all <- ac[which(nas == 0),]

mrs <- apply(cv.all, 1, function(snp){
  length(which(snp >= 1000))
})



rr.all.bs <- rr.all[which(mrs == 3),]
cv.all.bs <- cv.all[which(mrs == 3),]
rc.all.bs <- rc.all[which(mrs == 3),]
ac.all.bs <- ac.all[which(mrs == 3),]

t <- colnames(ac.all.bs)

dim(ac.all.bs)
head(ac.all.bs)
head(rc.all.bs)
####RUN LOGISTIC REGRESSION

snp <- rownames(rr.all.bs)[1]
i <- 1

i

lmt <- NULL
#dim(rc.all.bs)
#dim(ac.all.bs)
for(i in 1:nrow(rc.all.bs)){
  tab <- rbind(rc.all.bs[i,], ac.all.bs[i,])
  tab <- tab/100
  tab <- round(tab, digits = 0)
  #colt <- NULL
  #col01 <- NULL
  #colw <- NULL
  #for(i in 1:length(t)){
  #  colt <- c(colt, rep(t[i], 2))
  #  col01 <- c(col01, c(0,1))
  #  colw <- c(colw, c(tab[2,i], tab[1,i]))
  #}
  #mydata <- as.data.frame(cbind(colt,col01, colw))
  #head(mydata)
  #colnames(mydata) <- c("timepoint", "refAllele", "weights")
  #summary(mydata)
  #sapply(mydata, sd)
  #mydata$timepoint <- as.factor(mydata$timepoint)
  #mylogit <- glm(refAllele ~ timepoint, data = mydata, family = "binomial", weights = weights)
  
  x <- NULL
  y <- NULL
  for(j in 1:length(t)){
    x <- c(x, rep(t[j], (tab[1,j] + tab[2,j])))
    y <- c(y, c(rep(1, tab[1,j]), rep(0, tab[2,j])))
  }
  mydata <- as.data.frame(cbind(x,y))
  head(mydata)
  colnames(mydata) <- c("timepoint", "refAllele")
  summary(mydata)
  sapply(mydata, sd)
  mydata$timepoint <- as.factor(mydata$timepoint)
  
  xtabs(~refAllele + timepoint, data = mydata)
  
  mylogit <- glm(refAllele ~ timepoint, data = mydata, family = "binomial")
  
  
  dim(mydata)
  summary(mylogit)
  if(nrow(mydata) > 0){
    
    #summary(mylogit)
    #confint(mylogit)
    
    #chisq <- with(mylogit, (null.deviance/100) - (deviance/100))
    p <- with(mylogit, pchisq((null.deviance) - (deviance), df.null - df.residual, lower.tail = FALSE))
    #p <- with(mylogit, pchisq((null.deviance) - (deviance), df.null - df.residual, lower.tail = FALSE))
    coefs <- mylogit$coefficients
    
    #wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 2:6)
    trm2 <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 2)$result$chi2[3]
    trm3 <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 3)$result$chi2[3]
    #trm4 <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 4)$result$chi2[3]
    #trm5 <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 5)$result$chi2[3]
    #trm6 <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 6)$result$chi2[3]
    
    #trm23 <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 2:3)$result$chi2[3]
    #trm234 <- try(wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 2:4)$result$chi2[3], TRUE)
    #if(class(trm234)=="try-error"){
     # trm234 <- NA
    #}
    #trm456 <- try(wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 4:6)$result$chi2[3], TRUE)
    #if(class(trm456)=="try-error"){
    #  trm456 <- NA
    #}
    #trm56 <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 5:6)$result$chi2[3]
       
    #trm2 <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit)*100, Terms = 2)$result$chi2[3]
    #trm3 <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit)*100, Terms = 3)$result$chi2[3]
    #trm4 <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit)*100, Terms = 4)$result$chi2[3]
    #trm5 <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit)*100, Terms = 5)$result$chi2[3]
    #trm6 <- wald.test(b = coef(mylogit), Sigma = vcov(mylogit)*100, Terms = 6)$result$chi2[3]
    
    #lmt <- rbind(lmt, c(coefs, p))  
    #lmt <- rbind(lmt, c(coefs, trm2, trm3, trm4, trm5, trm6, trm23, trm234, trm456, trm56, p))
    lmt <- rbind(lmt, c(coefs, trm2, trm3, p))
  }
  else{
    lmt <- rbind(lmt, rep(NA, 6))
    #lmt <- rbind(lmt, rep(NA, 3))
    print("Something went wrong...")
  }
}
dim(lmt)
colnames(lmt) <- c("Intercept", t[2:3], paste(t[2:3], "_Pvalue", sep =""), "total_Pvalue")
#colnames(lmt)[7:16] <- c(paste(colnames(lmt)[2:6], "_Pvalue", sep = ""), "t0vst4t8", "t0vst4t8t12", "t0vst12t24t72", "t0vst24t72", "total_Pvalue")
#colnames(lmt) <- c("Intercept", "Beta", "Pvalue", "Chisq")
rownames(lmt) <- rownames(rc.all.bs)
head(lmt)
dim(rc.all.bs)
dim(lmt)

length(which(is.na(lmt[,3])))
#lmt <- lmt[-which(is.na(lmt[,3])),]
length(which(is.na(c("ole", NA, "mole"))))

write.table(lmt, "results_logReg_timepointsFACTOR_weightedDups_divBy100_mrs10all_noTofiVsOthers_maskedGenome.txt", sep = "\t", quote = FALSE)


#### check those at 5%FDR
library(qvalue)
pvalues <- lmt[,6]
head(pvalues)
#pvalues <- c(lmt[,7], lmt[,8], lmt[,9], lmt[,10], lmt[,11])
x11()
hist(pvalues, breaks = 50, main = "P-values model fit", col = "gray")
hist(lmt[,5], breaks = 50)
hist(lmt[,3], breaks = 50)
i <- 8
pdf("hists_pvalues_logRegFACTOR_divBy100_mrs10all_noTofiVsOthers_maskedG.pdf", width = 12, height = 8)
par(mfrow = c(2,3))
i <- 13
for(i in 4:ncol(lmt)){
  q <- qvalue(na.omit(lmt[,i]))
  pi1 <- 1-q$pi0
  numFDR5 <- length(which(q$qvalues <= 0.1))
  hist(lmt[,i], xlab = "P-value", main = paste(colnames(lmt)[i], "\npi1 =", round(pi1, digits = 4), ", signif 10% FDR =", numFDR5) , breaks = 50, col = "gray")
}
dev.off()

#q <- qvalue(pvalues)
q <- qvalue(lmt[,6])
pi1 <- 1-q$pi0
pi1

length(pvalues)
dim(lmt)
head(q$qvalues)
length(which(q$qvalues <= 0.5))
length(which(pvalues < 0.01))
#posFdr5  <- which(q$qvalues <= 0.2)
posFdr5 <- which(lmt[,6] < 0.01)
length(q$qvalues)
length(posFdr5)
dim(lmt)
summary(lmt[,6])
x11()
#hist(pvalues, xlab = "P-value", main = paste("P-values of t0 vs each of the others\npi1 =", round(pi1, digits = 4)) , breaks = 50, col = "gray")
#dev.print("hist_pvalues_LogRegQUANT_BASinAll_weightedDups_by100_0vsEach.pdf", dev = pdf)

exp <- runif(nrow(lmt))
head(lmt)
colnames(lmt)
x11()
chisq <- qchisq(1-lmt[,6],5)
lambda <- median(chisq)/qchisq(0.5,5)
qqplot(-log10(exp), -log10(lmt[,6]), cex.axis = 1.7, cex.lab = 1.7, xlab = "Expected", ylab = "Observed", main = paste("Multiple df test (min weighted cnts 10 in all)\nlambda =", round(lambda, digits = 3)))
abline(0,1, col = "red", lwd = 2)

dev.print("qqplot_LogRegFACTOR_BASinAll_weightedDups_divBy100_mrs10all_maskedG.pdf", dev = pdf)

###PLOT SIGNIF associations
snp2gene <- read.table("gene_snp_mapping_file_NEW_tofiPub.txt", stringsAsFactors = FALSE)
rc.all.bs[posFdr5,]
length(posFdr5)

pdf("refRatioVsTime_coverage_LogRegFACTOR_sitesBASinAll_P01_weightedDups_dividBy100_mrs10all_maskedG.pdf", width = 10, height = 4)
#x11(width = 10, height = 4)
par(mfrow = c(1,2), mar = c(4,5,4,2))
i <- posFdr5[6]
i <- 3061
lmt[grep("snp_chr18", rownames(lmt)),]
head(rc)
rc[grep("snp_chr18_77288806", rownames(rc)),]
ac[grep("snp_chr18_77288806", rownames(rc)),]

head(lmt)
#x11(width = 10, height = 4)
for(i in posFdr5){
  tab <- rbind(rc.all.bs[i,], ac.all.bs[i,])
  
  snp <- rownames(rc.all.bs)[i]
  gene <- snp2gene[which(snp2gene[,3] %in% snp),4]
  if(length(gene) == 0){
    print("ERROR: SNP not in snp2gene file")
  }else{
    sp <- strsplit(gene, ":")[[1]]
    exon <- sp[grep("-exon", sp)]
    if(length(exon) > 0){
      pc <- exon[grep("protein", exon)]
      if(length(pc) > 0){
        out <- c(sp[1], strsplit(pc[1], ";")[[1]][1])
        geneName <- out[2]
      }
      else{
        lr <- exon[grep("lincRNA", exon)]
        if(length(lr) > 0){
          out <- c(sp[1], strsplit(lr[1], ";")[[1]][1])
          geneName <- out[2]
        }
        else{
          out <- c(sp[1], strsplit(exon[1], ";")[[1]][1])
          geneName <- out[2]
        }
      }
    }
    else{
      intron <- sp[grep("-intron", sp)]
      if(length(intron) > 0){
        pc <- intron[grep("protein", intron)]
        if(length(pc) > 0){
          out <- c(sp[1], strsplit(pc[1], ";")[[1]][1])
          geneName <- out[2]
        }
        else{
          lr <- intron[grep("lincRNA", intron)]
          if(length(lr) > 0){
            out <- c(sp[1], strsplit(lr[1], ";")[[1]][1])
            geneName <- out[2]
          }
          else{
            out <- c(sp[1], strsplit(sp[4], ";")[[1]][1])
            geneName <- out[2]
          }
        }
      }
      else{
        out <- c(sp[1], strsplit(sp[4], ";")[[1]][1])
        geneName <- out[2]
      }
    }  
  }
  
  #logistic regression
  tab <- tab/100
  tab <- round(tab, digits = 0)  
  x <- NULL
  y <- NULL
  for(j in 1:length(t)){
    x <- c(x, rep(t[j], (tab[1,j] + tab[2,j])))
    y <- c(y, c(rep(1, tab[1,j]), rep(0, tab[2,j])))
  }
  
  mydata <- as.data.frame(cbind(x,y))
  head(mydata)
  colnames(mydata) <- c("timepoint", "refAllele")
  summary(mydata)
  sapply(mydata, sd)
  mydata$timepoint <- as.factor(mydata$timepoint)
  
  xtabs(~refAllele + timepoint, data = mydata)
  
  mylogit <- glm(refAllele ~ timepoint, data = mydata, family = "binomial")
  
  fit <- summary(mylogit)
  eq <- sprintf(
    "Max Beta = %.3g, P = %.3g",
    max(abs(fit$coefficients[2:nrow(fit$coefficients),1])),
    with(mylogit,pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
  )
  
  t <- c(0,1,2)
  plot(t, tab[1,]/(tab[1,]+tab[2,]), type = "o", lwd = 2, cex.lab = 1.5, cex.axis = 1.5, cex = 1.5, xlab = "Time point", ylab = "Reference ratio", pch = 19, main = paste(geneName, eq, sep = "\n"), ylim = c(0,1))
  tab <- as.matrix(tab)
  abline(h = 0.5, col = "gray", lwd = 1.5)
  rownames(tab) <- c("ref", "alt")
  barplot(tab, beside = T, col = c("dodgerblue3", "tan1" ), cex.axis = 1.5, cex.lab = 1.5, cex.name = 1.5, ylab = "Weighted counts per allele", main = snp)
    
}
dev.off()


dev.print("dqb1_missingSignifSNP_refRatioVsTime_cov_LogRegFACTOR_etc.pdf", dev = pdf)
