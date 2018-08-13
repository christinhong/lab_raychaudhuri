# Christin M. Hong
# 2017-05
# Lab of Soumya Raychaudhuri, Harvard Medical School
# Wellcome Trust Tutorial on genotype analysis (Anderson 2010)

#### Project Notes ####

# Data on Broad servers
# This script is for IDing individuals with high amounts of missing data or outlying heterozygosity rate (too low of a heterozygosity rate suggests inbreeding while too high suggests sample contamination).

#### Script ####

# import libraries
library(dplyr)

#example command R CMD BATCH '--args filename="all_dbsnp" a=3.5 b=0.1' imiss-vs-het.Rscript
args=(commandArgs(TRUE))

if(length(args)<=1){
  print("No filename arguments supplied.")
  a=3 #this is the default x standard deviation threshold for heterozygosity
  b=0.05 # this is the default missingness threshold
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

# Read in calculations of missing genotypes by individual
imiss=read.table(paste(filename,".imiss",sep=""),h=T)

# Log transform the missingness frequency (number of missing genotypes / number of total genotypes) 
imiss$logF_MISS = log10(imiss$F_MISS)

# Read in calculations of heterozygosity rates
het=read.table(paste(filename,".het",sep=""),h=T)


#### Calculate observed heterozygosity by individual ####
## Formula: ([number of non-missing genotypes] - [number of homozygous genotypes] / [number of non-missing genotypes])
## So...the frequency of heterozygote genotypes within the non-missing genotypes for each individual
het$meanHet = (het$E.HOM. - het$O.HOM.)/het$E.HOM.

# Coloring
colors  <- densCols(imiss$logF_MISS,het$meanHet)


#### Plot for QC check and save as PDF ####
pdf(paste(filename,"-imiss-vs-het.pdf",sep=""),width=10,height=7)

yrange<-c(min(het$meanHet),max(het$meanHet))*1.5

plot(imiss$logF_MISS,het$meanHet, col=colors, xlim=c(-3,0),ylim=yrange,pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=F)
     axis(2,at=seq(-0.1,0.1,length=11),tick=T)
     axis(1,at=c(-4,-3,-2,-1,0),labels=c(0.0001,0.001,0.01,0.1,1))
     lower=mean(het$meanHet)-(a*sd(het$meanHet))
     upper=mean(het$meanHet)+(a*sd(het$meanHet))
     abline(h=lower,col="RED",lty=2)
     abline(h=upper,col="RED",lty=2)
     abline(v=log10(b), col="RED", lty=2)

     fail.index<-which(het$meanHet<lower | het$meanHet>upper | imiss$logF_MISS>log10(b))
     points(imiss[fail.index,]$logF_MISS,het$meanHet[fail.index],col="orangered3",pch=1)
     
# Horizontal boundaries are +/- 3 standard deviations.  Vertical line is set at missingness rate of >3% (log10(0.03) = -1.522879).
  

#### Extract failed family IDs (FID, column 1) and individual IDs (IID, column 2) ####
# For file 'fail-imisshet-qc.txt' (one individual per line, tab delimited)

## Extracting individuals with excess missingness
fail_imiss <- filter(imiss, logF_MISS>log10(b)) # {dplyr}

count(imiss, logF_MISS>=log10(b))


# Extract family and individual IDs from failed samples
fail_imissIDs <- select(fail_imiss, FID, IID)

head(fail_imissIDs)


### Extracting individuals with outlying heterozygosity ###
head(het)

# Filtering for individuals with a heterozygosity rate less than or equal to the mean of the heterozygosity rate minus 2 standard deviations from the heterozygosity rate (so this is really a 2 SD cutoff...?)
fail_hetMin <- filter(het, meanHet<=mean(het$meanHet)-(a*sd(het$meanHet))) # {dplyr}

# Extract family and individual IDs from failed samples
fail_hetMinIDs <- select(fail_hetMin, FID, IID)

# Filtering for greater than mean + 2 SD
fail_hetMax <- filter(het, meanHet>=mean(het$meanHet)+(a*sd(het$meanHet))) # {dplyr}

# Extract family and individual IDs from failed samples
fail_hetMaxIDs <- select(fail_hetMax, FID, IID)


### Take union of failed IDs with duplicates removed ###
fail_imissHetIDs <- union(fail_imissIDs, fail_hetMinIDs, fail_hetMaxIDs)

## check results of removing the failed samples ##

# Extract samples with IIDs *not* in list of failed IIDs
imiss_qc <- imiss[!(imiss$IID %in% fail_imissHetIDs$IID), ] 
het_qc <- het[!(het$IID %in% fail_imissHetIDs$IID), ]

#### Replot qc data to check that correct samples have been removed ####
plot(imiss_qc$logF_MISS,het_qc$meanHet, col=colors, xlim=c(-3,0),ylim=yrange,pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=F, main = "Post-QC")
axis(2,at=seq(-0.1,0.1,length=11),tick=T)
axis(1,at=c(-4,-3,-2,-1,0),labels=c(0.0001,0.001,0.01,0.1,1))
lower=mean(het$meanHet)-(a*sd(het$meanHet))
upper=mean(het$meanHet)+(a*sd(het$meanHet))
abline(h=lower,col="RED",lty=2)
abline(h=upper,col="RED",lty=2)
abline(v=log10(b), col="RED", lty=2)

# looks good!

dev.off();


#### Save list of failed FIDs and IIDs to tab-delimited file ####
write.table(fail_imissHetIDs, "fail-imisshet-qc.txt", sep=" ", row.names = F, col.names = F, quote = F)


#### Moving back to bash/Broad server for using plink ####
