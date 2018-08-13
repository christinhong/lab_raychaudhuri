# Christin M. Hong
# 2015-08
# Soumya Raychaudhuri Lab, Harvard Medical School
# Differential expression analysis on mRNA from naive CD4+ T cells from human patients with and without a mutation in STAT1

#####
# Project Notes
# analyzing data from naive CD4 T cells from human patients from Hirahara-O'Shea 2015 Immunity paper.

# Controls are healthy patients.  Cases have a STAT1 mutation that seems to inhibit dimerization with STAT3, resulting in an over-accumulation of STAT1 dimers and defective IFNg, IL-22, and IL-17 production in response to IL-1b, IL-12, and IL-23.

# It's important to keep in mind that these patients were identified because they have chronic mucocutaneous candidasis (CMC), and so any differences in gene expression may be due to life with recurring fungal infections rather than a direct effect of a STAT1 variant.  (Thanks to Jason Qian for asking and bringing the question to mind.)

# Ultimately, cases have defective Th1 and Th17 responses which dramatically increase susceptibility to fungal infection.  (IFNg production from IL-18 is fine, so susceptibility to mycobacteria and virus is roughly normal.)  

# In the O'Shea paper, this is called a STAT1 "gain-of-function" due to the over-expression of STAT1-unique genes from the loss of STAT1/3 heterodimers.

# Here, "null" indicates human CD4+CD45RA(hi)CD45RO(lo) samples from PBMCs that have been stimulated with plate-coated anti-CD3 and anti-CD28 for three days.

# "IL-6" will indicate samples that also had 50 ng/mL of IL-6 in culture, which biases towards STAT3 expression.  (STAT3 KO -> severely reduced gene expression profile).

# "IL-27" samples had 50 ng/mL of IL-27 in culture, which increases STAT1 expression.  (STAT1 KO -> IL-6 and IL-27 eliciting the same atypical expression profile = loss of specificity).

# RNA was extracted from 1e6 cells by mirVana miRNA Isolation Kit.  200 ng was processed with TruSeq SR RNA Sample Prep Kit to form an mRNA library; supposedly these samples were processed like this since they were analyzed by mRNA?  Libraries were sequenced for 50 cycles (single read) with Illumina HiSeq 2000.  (NOTE: Whole genome sequencing - not biased towards immune genes.)

#! Need to add representative output to this file/analysis protocol for easier future reference

#####
# Infrastructure


# LATER ON, convert this to modular by changing paths/project-specific values to variables that can be set in this section
# For filtering 0 reads, will need to filter all samples in the set to ensure that I'm comparing the same set of genes between them...so I'll need to redo the fc analysis for these null samples with the analyses for the IL-6 and IL-27 sets.  (Do that on server.)

# split output between console and a txt file
sink(file = "Data/2015-Hirahara_OShea-stat1/2015-Hirahara_OShea-stat1 RStudio output.txt", append = TRUE, split = TRUE)

# import libraries
library(Rsubread)
library(limma)
library(edgeR)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(gtable)
library(DESeq2)
library(Biobase)
library(matrixStats)
library(mixOmics)
library(RColorBrewer)
library(tidyr) # not yet used
library(dplyr) # not yet used

wd <- "~/aa_Raychaudhuri-Lab/data/015-08_RNA_Hirahara-OShea_stat1/archive"
suffix <- "_2015-STAT1"
setwd(wd)

#####
# Far unknown reason, following code results in "unique mapping : no."  Not sure why, but I run subjunc in Subread from the CL.  (Maybe because I was using an old version of R?)

# subjunc(index = "Data/ref-hg19/hg19-index2",
#        readfile1 = "Data/2015-Hirahara_OShea-stat1/SRR1786605-S1G5-null.fastq", readfile2=NULL,
#        input_format="FASTQ", output_format="BAM",
#        nsubreads=14, 
#        TH1=1, TH2=1, nthreads=1, indels=5, phredOffset=33,
#        tieBreakQS=FALSE, tieBreakHamming=TRUE, unique=TRUE,
#        minFragLength=50, maxFragLength=100000,
#        PE_orientation="fr",
#        nTrim5=0, nTrim3=0,
#        readGroupID=NULL, readGroup=NULL,color2base=FALSE,DNAseq=FALSE,reportAllJunctions=TRUE)

#####
# generating raw count table with featureCounts

fc <- featureCounts(files = c("Data/2015-Hirahara_OShea-stat1/HC1-null.bam",
                              "Data/2015-Hirahara_OShea-stat1/HC2-null.bam", 
                              "Data/2015-Hirahara_OShea-stat1/HC3-null.bam", 
                              "Data/2015-Hirahara_OShea-stat1/HC4-null.bam", 
                              "Data/2015-Hirahara_OShea-stat1/S1G1-null.bam", 
                              "Data/2015-Hirahara_OShea-stat1/S1G2-null.bam", 
                              "Data/2015-Hirahara_OShea-stat1/S1G3-null.bam", 
                              "Data/2015-Hirahara_OShea-stat1/S1G4-null.bam", 
                              "Data/2015-Hirahara_OShea-stat1/S1G5-null.bam"),
                    annot.ext = "Data/ref-hg19/gencode-v19-chr.gtf", 
                    isGTFAnnotationFile = TRUE, 
                    useMetaFeatures = TRUE, 
                    GTF.attrType = "gene_id", 
                    allowMultiOverlap = FALSE, 
                    countMultiMappingReads = FALSE, 
                    isPairedEnd = FALSE,
                    nthreads = 4,
                    minMQS = 20)

# checking fc
# finding section names.  Expected return: [1] "counts" "annotation" "targets" "stat"  
names(fc)
# checking first 7 lines of count section ($ is operator for calling sections)
head(fc$counts)
# checking featureCount analysis stats
fc$stat


##### 
# Trying different count analysis methods
# plotting the second column (first sample) of counts
plot(fc$counts[,1])
abline(0,1) # regression line

# finding, checking, and plotting RPKM values with tools in edgeR library
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])
x_rpkm <- rpkm(x,x$genes$Length)
head(x_rpkm)

# different ways of plotting rpkm
plot(x_rpkm[,1])

hist(x_rpkm[which(x_rpkm[,1] > 50)], breaks = 100)

boxplot(x_rpkm, las = 2) # las = 2 sets labels perpendicular to axis.  I like the boxplot because I can see all 8 samples at once, but it's still too messy to get any real info

density(x_rpkm)

# "voom {limma} transforms count data to log2 counts per million (logCPM), estimates the mean-variance relationship, and uses it to compute appropriate observational-level weights for linear modeling."  
# I don't really understand what all that means, but I'll record it here for later.
fc_voom <- voom(x, plot = TRUE)
names(fc_voom)

# sample clustering with limma
# plotMDS: "Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the typical log2 fold changes between the samples."  Basically a 2-axes PCA.
plotMDS(fc_voom)

#####
# Trying out analysis methods based on Ignacio Gonzalez's online tutorial from Nov. 2014

# writing fc$counts to an external tab-delineated text file to isolate count data
write.table(x = data.frame(fc$annotation[,c("GeneID","Length")],
                           fc$counts,stringsAsFactors=FALSE),
            file = "~/aa_Raychaudhuri-Lab/data/2015-08_RNA_Hirahara-OShea_stat1/archive/2015-08-12 - fc_counts.txt",
            quote=FALSE,
            sep="\t",
            row.names = FALSE)



#####
# RE-START HERE WITH SAVED FILES


# importing above *.txt file as a table (command depends on file format; see http://rstudio-pubs-static.s3.amazonaws.com/1776_dbaebbdbde8d46e693e5cb60c768ba92.html)
fcTable = read.delim("~/aa_Raychaudhuri-Lab/data/2015-08_RNA_Hirahara-OShea_stat1/archive/2015-08-12 - fc_counts.txt")
head(fcTable) # checking output

# renaming columns
colnames(fcTable)[3:6] = paste0("HC", 1:4)
colnames(fcTable)[7:11] = paste0("S1G", 1:5)
head(fcTable)

# plotting data as a histogram
x11()
ggplot(fcTable, aes(x = HC1)) + geom_histogram(binwidth = 2000)

### Since the range of the number of counts can be large ("count values distribution is highly skewed"), people typically transform data to a log base 2 scale for easier analysis (a large part of QC is based on human interpretation of visual representations of data).  

# Log2 scale means that a difference of 1 on log2 = 2-fold change, and this scale generates a good representation of the biological range of RNA expression.  However, counts can be 0, and log2(0) is outputted as negative infinity.  Therefore, people often use pseudocounts by adding 1 to each count, as log2(1) = 0, and the addition of 1 read should have minimal impact on the values of the data.

# To convert read counts to log2: 
# making new copy of table to work on
L2fcTable <- fcTable

# converting fcTable sample values to log2 and storing in copy
## OLD CODE
## L2fcTable[,3:11] <- log(fcTable[3:11], 2)
## head(L2fcTable)
### this gives a large number of infinity counts due to having 0 reads before transformation, so standard practice is to use a psuedocount by adding +1.  This returns 0 (since log2(1) = 0), and has minimal impact on other counts.  But depending on the context, sometimes the true count is better.
L2fcTable[,3:11] <- log(fcTable[3:11] + 1, 2)
head(L2fcTable)
# much better!

# plotting log2 transformed counts.  Clearly most of the genes have 0 reads (though this could be missed by using a larger bin value.
ggplot(L2fcTable, aes(x = HC1)) + ylab(expression(log[2](count+1))) + geom_histogram(binwidth = 0.3, color = "#FFFFFF")

#####
# Structuring data for plotting

# melt {reshape} is used to reorganize data for between sample distribution.  melt transforms groups from wide to long format, as ggplot requires long-format data.
bx = melt(L2fcTable, id = c("GeneID", "Length"),
          measure.vars = c("HC1",
                           "HC2",
                           "HC3",
                           "HC4",
                           "S1G1",
                           "S1G2",
                           "S1G3",
                           "S1G4",
                           "S1G5"),
          variable.name = "Samples")
summary(bx)
head(bx)
levels(bx$Samples)

# subsetting molten list between samples and creating control and case groups
## bxSubsets <- split(bx, bx$Samples, drop = TRUE)
## head(bxSubsets[[2]]) # check
## HC <- bxSubsets[1:4]
## S1G <- bxSubsets[5:8]
# can see from this that the S1G group starts at line 231281, which is nifty, but currently not useful.


# Want to create control and case groups for graphical comparison.  After 8 hours of searching the internet, trying different lines of code, and reading a ton about data frames and levels, finally figured out how to do this!

# To add a new colunm to transformed data, turn it into a dataframe with a new column that replicates the Sample variable.  Call this new column Condition.
bx.df <- data.frame(bx, Condition = substr(bx$Samples, 1, 9))
head(bx.df)
tail(bx.df)

# for ease of undoing, make a copy of this dataframe
bx2 <- bx.df

# confirm levels in column Condition
levels(bx2$Condition)

# rename levels in Condition using a list
levels(bx2$Condition) <- list(HC = c("HC1", "HC2", "HC3", "HC4"), S1G = c("S1G1", "S1G2", "S1G3", "S1G4", "S1G5"))

# check if renaming worked
head(bx2$Condition)
head(bx2)
tail(bx2)

# Now we can graph data by condition while retaining separate replicates

##### 
# Graphing data: Between sample distribution
# can open a new window for plotting with 'dev.new()' or 'x11()'

# ggplot box plot of psuedocounts.  
x11()
ggplot(bx2, aes(x = Samples, y = value, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count + 1))) + scale_fill_manual(values = c("#619cff", "#f564e3"))
# Important to NOT put aes values in quotes, or else output won't work.

# ggplot histogram of psuedocounts by sample
ggplot(bx2, aes(x = value, color = Samples, fill = Samples)) + ylim(c(0,0.25)) + geom_density(alpha = 0.2, size = 1.25) + theme(legend.position = "top") + xlab(expression(log[2](count + 1))) + facet_wrap(~ Samples)

# ggplot histogram of psuedocounts by condition
ggplot(bx2, aes(x = value, color = Samples, fill = Samples)) + ylim(c(0,0.25)) + geom_density(alpha = 0.2, size = 1.25) + theme(legend.position = "top") + xlab(expression(log[2](count + 1))) + facet_wrap(~ Condition)

# Show multiple plots at once:
# Store plots in variables.
box <- ggplot(bx2, aes(x = Samples, y = value, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count + 1))) + scale_fill_manual(values = c("#619cff", "#f564e3"))

hist_group <- ggplot(bx2, aes(x = value, color = Samples, fill = Samples)) + ylim(c(0,0.25)) + geom_density(alpha = 0.2, size = 1.25) + theme(legend.position = "top") + xlab(expression(log[2](count + 1))) + facet_wrap(~ Condition)

# Display plots with gridExtra (for more precise control of plotting, will want to use grid.table)
x11()
grid.arrange(box, hist_group, ncol = 2)
dev.print(png,filename="Data/2015-Hirahara_OShea-stat1//Images/13-count analysis.png",width=1000);
dev.off ();

##### 
# Graphing data: Check if normalization is needed with MA plot

# MA plot = plot of log-fold change (M-values, "log of the ratio of level counts for each gene between two samples") vs. log-fold average (A-values, "log of average level counts for each gene across two samples").  So if a gene is counted 16 times in Sample 1 and 4 times in Sample 2:
# M = log2(16/4)
# A = (1/2)log2(16+4)
# repeated for every gene in the sample.

# MA plots are therefore used to visualize reproducibility between samples in an experiment.  Of course, it makes sense to limit the use of MA plots to within experimental groups, e.g. one MA plot for cases and another for controls.

# If the expression levels in two samples are similar, they will appear around M = 0 on the y-axis.
# (Note that MA plots should REALLY be called AM plots, since A is on the x axis.)

# M = log2(4/4) = log2(4) - log2(4) = log2(1) = 0
# (A = (1/2)log2(4+4) = (1/2)log2(8), will be proportional to expression level)

# To graph a basic MA plot between HC1 and HC2 from log2 normalized fc data:
head(L2fcTable)


h1 <- L2fcTable[,3] # set h1 as HC1 values
h2 <- L2fcTable[,4] # set h2 as HC2 values

h3 <- L2fcTable[,5]
h4 <- L2fcTable[,6]

s1 <- L2fcTable[,7]
s2 <- L2fcTable[,8]
s3 <- L2fcTable[,9]
s4 <- L2fcTable[,10]
s5 <- L2fcTable[,11]

# Since x and y values have already been transformed by log2:
Mh12 = h1 - h2

Ah12 = (h1 + h2)/2

ma.h12.df = data.frame(Ah12, Mh12)

# code for MA scatterplot with a loess fit (?) in red for estimating a possible trend in bias, storing in variable
ma.h12 <- ggplot(ma.h12.df, aes(x = Ah12, y = Mh12)) + geom_point(size = 1.5, alpha = 1/5) + geom_hline(color = "blue3") + stat_smooth(se = FALSE, method = "loess", color = "red3")


# To get a feel for the rest of the dataset, repeat for h3 vs h4, s1 vs s2, s2 vs s3, s3 vs s4, asd s4 vs s5.  (Paying more attention to s3 because in Hirahara's paper, the S3 samples looked particularly sparse, but there doesn't seem to be a difference in read number or mapping so far.)

Mh34 = h3 - h4
Ah34 = (h3 + h4)/2
ma.h34.df = data.frame(Ah34, Mh34)
ma.h34 <- ggplot(ma.h34.df, aes(x = Ah34, y = Mh34)) + geom_point(size = 1.5, alpha = 1/5) + geom_hline(color = "blue3") + stat_smooth(se = FALSE, method = "loess", color = "red3")

Ms12 = s1 - s2
As12 = (s1 + s2)/2
ma.s12.df = data.frame(As12, Ms12)
ma.s12 <- ggplot(ma.s12.df, aes(x = As12, y = Ms12)) + geom_point(size = 1.5, alpha = 1/5) + geom_hline(color = "blue3") + stat_smooth(se = FALSE, method = "loess", color = "red3")

Ms23 = s2 - s3
As23 = (s2 + s3)/2
ma.s23.df = data.frame(As23, Ms23)
ma.s23 <- ggplot(ma.s23.df, aes(x = As23, y = Ms23)) + geom_point(size = 1.5, alpha = 1/5) + geom_hline(color = "blue3") + stat_smooth(se = FALSE, method = "loess", color = "red3")

Ms34 = s3 - s4
As34 = (s3 + s4)/2
ma.s34.df = data.frame(As34, Ms34)
ma.s34 <- ggplot(ma.s34.df, aes(x = As34, y = Ms34)) + geom_point(size = 1.5, alpha = 1/5) + geom_hline(color = "blue3") + stat_smooth(se = FALSE, method = "loess", color = "red3")

Ms45 = s4 - s5
As45 = (s4 + s5)/2
ma.s45.df = data.frame(As45, Ms45)
ma.s45 <- ggplot(ma.s45.df, aes(x = As45, y = Ms45)) + geom_point(size = 1.5, alpha = 1/5) + geom_hline(color = "blue3") + stat_smooth(se = FALSE, method = "loess", color = "red3")

x11()
grid.arrange(ma.h12, ma.h34, ma.s12, ma.s23, ma.s34, ma.s45, ncol = 3)
ggsave("Data/2015-Hirahara_OShea-stat1/Images/2015-08-12 - Hirahara_OShea - MA plot.png", width = 8, height = 6)

## Looks like replication quality of samples is generally quite good!  But there's still some room for normalizing.  (And S3 still doesn't look particularly special...)


# Just out of curiosity...
x11()
Mh1s1 = h1 - s1
Ah1s1 = (h1 + s1)/2
ma.11.df = data.frame(Ah1s1, Mh1s1)
ggplot(ma.11.df, aes(x = Ah1s1, y = Mh1s1)) + geom_point(size = 1.5, alpha = 1/5) + geom_hline(color = "blue3") + stat_smooth(se = FALSE, method = "loess", color = "red3")
ggsave("Data/2015-Hirahara_OShea-stat1/Images/2015-08-12 - Hirahara_OShea - MA plot h1s1.svg", width = 4, height = 3)
# Pretty cool.  Can see that in general, S1 has higher gene expression than H1 (correlates with earlier histogram of counts).  OTOH, S1G samples are from patients with severe fungal infections (or possibly on treatment for infections?), so even though these are CD45RA(hi) = "naive" T cells, they may still be primed for activation by existing in an inflamed environment.  That's a factor that'll affect all aspects of comparison between these two groups, which is why biological validation is so important.  (Thanks to Jason Qian for asking if the patients had active infections over lunch and sparking that realization!)




# PCA

#http://www.r-bloggers.com/using-r-two-plots-of-principal-component-analysis/
head(L2fcTable)

L2fcTable2 <- L2fcTable[,-c(1,2)]
rownames(L2fcTable2) <- L2fcTable[,1]


# Filter out 0 counts
keep = rowSums(L2fcTable2) > 0
L2Count.filt = L2fcTable2[keep, ]

dim(L2fcTable2) # 57820 rows, 9 columns
dim(L2Count.filt) # 40057     9



L2fcTable2.ma <- data.matrix(L2fcTable2, rownames.force = NA) # Must store data.matrix
str(L2fcTable2.ma[, 1:4]) # Checking if conversion worked

group = c("HC", "HC", "HC", "HC", "S1G", "S1G", "S1G", "S1G", "S1G")


x11()
pdf("aa_RNA-02-DE_NUM_clusteringData.pdf",width = 10, height = 8, useDingbats = F);

## Color picker!!  :D
# {RColorBrewer} # http://colorbrewer2.org/ ; http://www.compbiome.com/2010/12/r-using-rcolorbrewer-to-colour-your.html
### Set the display a 2 by 2 grid
par(mfrow=c(2,2))

### Show all the colour schemes available
display.brewer.all()


# Heatmap {mixOmics}
mat.dist <- L2fcTable2.ma
mat.dist <- as.matrix(dist(t(mat.dist)))
mat.dist2 = mat.dist/(max(mat.dist))
head(mat.dist2)

hmcol = colorRampPalette(brewer.pal(9, "Blues"))(16)
cim(mat.dist2, color = rev(hmcol), symkey = F, margins = c(9,9), transpose = F)


dev.off ();


rv = rowVars(L2fcTable2.ma)
select = order(rv, decreasing = T)[1:1000] # 2nd number here = how many of the most variable genes should be used for calculating PCA
pca <- prcomp(t(L2fcTable2.ma[select, ]), scale = T)

names(pca)
head(pca$rotation)

meltedPCA <- cbind(group, melt(pca$rotation[,1:9]))

ggplot(data=meltedPCA) + geom_bar(aes(x=Var1, y=value, fill=group), stat="identity") + facet_wrap(~Var2)


scores <- data.frame(sample.groups, pca$x[,1:3])
pc1.2 <- qplot(x=PC1, y=PC2, data=scores, colour=factor(sample.groups)) +
    theme(legend.position="none")
pc1.3 <- qplot(x=PC1, y=PC3, data=scores, colour=factor(sample.groups)) +
    theme(legend.position="none")
pc2.3 <- qplot(x=PC2, y=PC3, data=scores, colour=factor(sample.groups)) +
    theme(legend.position="none")



#####
# Mosiac plot of featureCount data to demonstrate QC
head(fc$stat)
names(fc$stat)
fcStat <- (fc$stat)


# renaming columns
colnames(fcStat)[2:5] = paste0("HC", 1:4)
colnames(fcStat)[6:10] = paste0("S1G", 1:5)
head(fcStat)

## Building mosaic plot
# Invert rows and columns
fcStat.n <- fcStat$Status
fcStat.dft <- as.data.frame(t(fcStat[,-1]))
colnames(fcStat.dft) <- fcStat.n
# fcStat.dft$Samples <- factor(row.names(fcStat.dft)) # Adds a column for referencing rows by sample (see the new column 12), but I don't think I need it here...

names(fcStat.dft)
str(fcStat.dft) # check column types


## Build basic mosaic plot
x11()
mosaicplot(fcStat.dft[,c(1,2,4,5)], shade = FALSE, color=c("#2ED0D0", "#FF9700", "#45B922", "#FFD000"), las = 1, cex.axis = .9, main = "featureCount Stats", sub = "From human CD4+ CD45RA(hi) CD45RO(lo) T cells with and without STAT1 mutations.  \n Treated 3 days with anti-CD3/CD28.")
dev.print(png,filename="11-fcStat mosaic.png",width=800);
dev.off ();

## I guess this is it, but I don't really like how relative it is...  Decided to turn it into a bar graph with actual numbers for context
fc_stats_count = melt(fcStat, id = c("Status"),
                      measure.vars = c("HC1",
                                       "HC2",
                                       "HC3",
                                       "HC4",
                                       "S1G1",
                                       "S1G2",
                                       "S1G3",
                                       "S1G4",
                                       "S1G5"),
                      variable.name = "Samples")
summary(fc_stats_count)
head(fc_stats_count)
levels(fc_stats_count$Samples)


fc_stats_count.df <- data.frame(fc_stats_count, Condition = substr(fc_stats_count$Samples, 1, 9))

# confirm levels in column Condition
levels(fc_stats_count.df$Condition)

# rename levels in Condition using a list
levels(fc_stats_count.df$Condition) <- list(HC = c("HC1", "HC2", "HC3", "HC4"), S1G = c("S1G1", "S1G2", "S1G3", "S1G4", "S1G5"))

# check if renaming worked
head(fc_stats_count.df$Condition)
head(fc_stats_count.df)
tail(fc_stats_count.df)

# Making quantitative bar graphs of featureCount Stats
x11()
ggplot(fc_stats_count.df, aes(x = Samples, y = value, fill = Status)) + geom_bar(stat="identity",colour="black", position = "stack") + xlab("") + ylab(expression("Number of reads")) + guides(fill=FALSE) + ggtitle("featureCount Stats") + facet_wrap(~ Status) +  theme(axis.text.x = element_text(angle=90, vjust=1), text = element_text(size=14), axis.text=element_text(colour="#525252")) 
#vjust adjust the vertical justification of the labels, which is often useful

x11()
ggplot(fc_stats_count.df, aes(x = Samples, y = value, fill = Status)) + geom_bar(stat="identity", colour="black", position = "stack") + xlab("") + ylab(expression("Number of reads")) + scale_fill_manual(values = c("#2ED0D0","#FF9700","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#45B922","#FFFFFF","#FFFFFF","#FFD000")) + guides(fill=FALSE) + ggtitle(expression(atop("featureCount Statistics", atop("Bottom to top: Assigned, Ambiguity, No Features, and Unmapped"), ""))) +  theme(text = element_text(size=18), axis.text=element_text(colour="#525252", size=16), axis.title.y=element_text(angle=90))
dev.print(png,filename="Data/2015-Hirahara_OShea-stat1//Images/12-fcStats.png",width=800);
dev.off ();


# Much more informative, AND so much easier!  (Though to be fair, that's mainly because I can now reuse code for ggplot, whereas I had to go find and learn everything anew for mosaicplot.)

# Can clearly see difference in read counts between controls and cases.  If this were run on an immunochip, I might believe it, but from this data, I suspect that the authors ran the cases and controls separately...  Not ideal.



## Trying something fancier with {vcd} so I can label the values, but man, the defaults on mosaicplot are much easier to work with...
# head(fcStat.dft)
# fcStat.dft2 <- fcStat.dft
# fcStat.dft2$Unassigned_MultiMapping <- NULL
# fcStat.ma <- data.matrix(fcStat.dft2)
# head(fcStat.ma)

#mosaic(fcStat.ma[,1:4], shade = T, gp_labels = gpar(fontsize = 10), rot_labels = c(90, 180, 0, 0), just_labels=c("left","left","right","right"), direction = "v", margins = 3, title_margins = 2, pop = FALSE, main = "featureCount Stats", sub = "From human CD4+ CD45RA(hi) CD45RO(lo) T cells with and without STAT1 mutations.  \n Treated 3 days with anti-CD3/CD28.")





#####
## Filtering genes with < 1 counts per million (CPM) for all current samples
# Filtering by CPM adjusts for the library size/read depth (CPM = (counts / library size) * 1e6)
keep = rowSums(fc$counts) > 0 # Can't use fcTable directly because it has 2 extra columns (row ID and length), and the rowSums function only distinguishes between the first row and column.  (Get "x must be numeric error" because function can't sum the GeneID.)

filtL2Count = L2fcTable[keep, ]
dim(fc$counts) # 57820 rows, 9 columns
head(fc$counts)

dim(filtL2Count) # 40057 rows, 11 columns (row # and gene length)
head(filtL2Count)

# Analyzing filtered genes
filtL2Count.df = melt(filtL2Count, id = c("GeneID", "Length"),
                    measure.vars = c("HC1",
                                     "HC2",
                                     "HC3",
                                     "HC4",
                                     "S1G1",
                                     "S1G2",
                                     "S1G3",
                                     "S1G4",
                                     "S1G5"),
                    variable.name = "Samples")
head(filtL2Count.df)
filtL2Count.df = data.frame(filtL2Count.df, Condition = substr(filtCount.df$Samples, 1, 9))
tail(filtL2Count.df)

# rename levels in Condition using a list
levels(filtL2Count.df$Condition) <- list(HC = c("HC1", "HC2", "HC3", "HC4"), S1G = c("S1G1", "S1G2", "S1G3", "S1G4", "S1G5"))

# check if renaming worked
head(filtL2Count.df$Condition)
head(filtL2Count.df)
tail(filtL2Count.df)

ggplot(filtL2Count.df, aes(x = value, color = Samples, fill = Samples)) + geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ Condition) + theme(legend.position = "top") + xlab(expression(log[2](count + 1)))


##### 
## Normalizing for read depth/library size (after filtering and converting to log2)

##### 
# TMM (trimmed mean of M-values) {edgeR} assumes that the majority of gene aren't differentially expressed, and therefore the counts can be normalized by the expression levels of the median genes.  
# Refresher: M = log-fold change in gene expression between two data sets.  A = average of the log-fold counts of a gene from two data sets ("intensity"). 
# A trimmed mean = avg. after ignoring the upper and lower x% of data.  TMM {edgeR} is doubly trimmed because it trims by M and A.  The default is to trim M by 30% and A by 5%, but this can be adjusted (I'm not sure if that means trimming M by 15% from upper and lower, or 30%, but I can easily determine this if necessary).
# The doubly-trimmed values are then weighted in order to adjust for differences in library size. 

# default: calcNormFactors(object, method=c("TMM","RLE","upperquartile","none"), refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)

# Before h1s1:
x11()
Mh1s1 = h1 - s1
Ah1s1 = (h1 + s1)/2
ma.11.df = data.frame(Ah1s1, Mh1s1)
ggplot(ma.11.df, aes(x = Ah1s1, y = Mh1s1)) + geom_point(size = 1.5, alpha = 1/5) + geom_hline(color = "blue3") + stat_smooth(se = FALSE, method = "loess", color = "red3")

## Calculating TMM normalization factors
# Building edgeR object (rows = gene ID, columns = samples)
head(fcTable)
fcTable2 <- fcTable
# Set rownames to fcTable2$GeneID
rownames(fcTable2) <-fcTable2$GeneID 
head(fcTable2)

fcTable2$GeneID <- NULL
fcTable2$Length <- NULL
head(fcTable2)
summary(rowSums(fcTable2))

# Filtering out rows with 0 reads
keep = rowSums(fcTable2) > 0 
ffc = fcTable2[keep, ]
dim(fcTable2) # 57820 rows, 9 columns
dim(ffc) # 40057, 9


# Create DGEList object from RAW COUNTS (not log2 transformed!)
group <- c(rep("Controls", 4), rep("Cases", 5))

dge <- DGEList(counts = ffc, group = group)

fc.tmm <- calcNormFactors(dge, method = c("TMM"), refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)

# From http://davetang.org/muse/2012/01/19/the-dgelist-object-in-r/:
fc.tmm <- estimateCommonDisp(fc.tmm, verbose = T) # Result: Disp = 0.15424 , BCV = 0.3927 
fc.tmm <- estimateTagwiseDisp(fc.tmm)
fc.tgw <- exactTest(fc.tmm)


# https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html and edgeR tutorial:

#   The function topTags() takes the output from exactTest(), adjusts the raw p-values using the False Discovery Rate (FDR) correction, and returns the top differentially expressed genes.  The output is similar to that of exactTest() but with a column of adjusted p-values and sorted by increasing p-value.  The sort.by argument allows you to sort the table by p-value, concentration or fold-change if desired.  We can also use topTags() to return the original counts of the top differentially expressed genes.  By setting the n parameter to the total number of genes, we can save the entire topTags() results table

topTags(fc.tgw, n = 20, sort.by = "p.value")

# Store full topTags results table
resultsTT.tgw <- topTags(fc.tgw, n = nrow (fc.tgw$table), sort.by = "p.value")$table

# Store topTags results in csv:
write.table(resultsTT.tgw, file="Data/2015-Hirahara_OShea-stat1/2015-08-25 - edgeR TMM DE TopTags.csv", append = FALSE, quote = TRUE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = NA, qmethod = c("escape"), fileEncoding = "")



# Default: decideTestsDGE(object, adjust.method="BH", p.value=0.05)

de <- decideTestsDGE(fc.tgw, adjust.method="BH", p.value=0.05)
summary(de)
#   [,1] 
# -1     3 = -1 for down-regulated genes
# 0  40039 = 0 for non-DE genes
# 1     15 = 1 for up-regulated

# differentially expressed tags from the naive method in d1
fc.de <- rownames(fc.tmm)[as.logical(de)] 
x11()
plotSmear(fc.tgw, de.tags=fc.de) 
# plotSmear is a more sophisticated and superior way to produce an 'MA plot'. plotSmear resolves the problem of plotting tags that have a total count of zero for one of the groups by adding the 'smear' of points at low A value. The points to be smeared are identified as being equal to the minimum estimated concentration in one of the two groups. The smear is created by using random uniform numbers of width smearWidth to the left of the minimum A. plotSmear also allows easy highlighting of differentially expressed (DE) tags. 
abline(h = c(-2, 2), col = "blue")
dev.print(png,filename="Data/2015-Hirahara_OShea-stat1/Images/14-TMM.png",width=800);
dev.off ();


# Plotting the tagwise dispersion calculated from this
names(fc.tmm)
ffc$twd <- fc.tmm$tagwise.dispersion
head(ffc)
hist(ffc$twd, breaks=20, xlim=c(0,3))

#now if you wanted to save the fold changes
names(fc.tgw)
head(fc.tgw$table)
head(fc.tgw$genes)
ffc2 <- cbind(ffc, fc.tgw$table)
head(ffc2)

#I want the fdr
ffc2$PValue_fdr <- p.adjust(method="fdr",p=ffc2$PValue)
head(ffc2)

#sanity check with the decideTestsDGE() call
table(ffc2$PValue_fdr<0.05) # TRUE should be sum of -1 and 1 from summary(de).  Works here: 3 + 15 = 18

table(ffc2$PValue_fdr<0.01)


write.table(ffc2, file="Data/2015-Hirahara_OShea-stat1/2015-08-25 - edgeR TMM DE.csv", append = FALSE, quote = TRUE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = NA, qmethod = c("escape"), fileEncoding = "")
sessionInfo()


# Trying to make plot with top 20 hits labeled with gene ID, from https://support.bioconductor.org/p/56590/:
de.name <- rownames(topTags(fc.tgw, n = 20)$table)

x11()
plotSmear(fc.tgw, de.tags=de.name, cex = 1)
fc.tgw$genes <- rownames(fc.tmm$counts)
ids <- c("ENSG00000137959.11", "ENSG00000137965.6")
gene.labels <- fc.tgw$table[fc.tgw$genes %in% ids,]
abline(h = c(-2, 2), col = "blue")
points(x=gene.labels$logCPM, y=gene.labels$logFC, cex=2, col="red")
text(x = gene.labels$logCPM, y = gene.labels$logFC, labels = rownames(gene.labels), cex = 1, pos = 2)
dev.print(png,filename="Data/2015-Hirahara_OShea-stat1/Images/15-TMM.png",width=800);
dev.off ();
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
## voom {limma} seems to be an updated version of relative log expression (RLE) {edgeR}.  RLE normalizes datasets by generating a pseudo-reference sample from the geometric mean across all the genes, which it then uses to generate a weight for adjusting the library sizes.  Voom generates weights from the mean AND variance, since higher count numbers show decreased variance.

# According to limma userguide, "The limma-voom method assumes that rows with zero or very low counts have been removed.
# It is usual to apply scale normalization to RNA-seq read counts, and the TMM normalization method [29] in particular has been found to perform well in comparative studies.  To apply TMM normalization, it is convenient to create a DGEList object using the edgeR package:
# > dge <- DGEList(counts=counts)
# > dge <- calcNormFactors(dge)
# The voom transformation is then applied:
#  v <- voom(dge,design,plot=TRUE)
#The voom transformation uses the experiment design matrix, and produces an EList object.
# It is also possible to give a matrix of counts directly to voom without TMM normalization, by
# > v <- voom(counts,design,plot=TRUE)

# From https://rstudio-pubs-static.s3.amazonaws.com/85101_ad12a157e52d4fbb8361219c0c8b50c1.html, http://bioinf.wehi.edu.au/RNAseqCaseStudy/, and the limma user guide:

fc.tmmv <- calcNormFactors(dge, method = c("TMM"), refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)

x11()
voom(fc.tmmv, design = NULL, plot = T)
dev.print(png,filename="Data/2015-Hirahara_OShea-stat1/Images/16-TMMV.png",width=800);
dev.off ();

v <- voom(fc.tmmv, design = NULL, plot = T)

plotMDS(v, labels = 1:9, main = "MDS plot", cex = 2)
dev.print(png,filename="Data/2015-Hirahara_OShea-stat1/Images/17-TMMV MDS.png",width=800);
dev.off ();

fit <- lmFit(v, design = NULL) # fit linear model
fit <- eBayes(fit) # compute moderated t statistic for each gene

topTable(fit, adjust.method = "BH", sort.by = "P")
# Store topTable in csv:
tTv <- topTable(fit, adjust.method = "BH", n = nrow (fit), sort.by = "P")
write.table(tTv, file="Data/2015-Hirahara_OShea-stat1/2015-08-25 - edgeR TMM DE TopTable.csv", append = FALSE, quote = TRUE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = NA, qmethod = c("escape"), fileEncoding = "")


# Call DEGs
resultsTT.v <- decideTests(fit, adjust.method = "BH")

# save results in csv:
write.fit(fit, results = resultsTT.v, file = "Data/2015-Hirahara_OShea-stat1/2015-08-25 - voom TMM decideTests.csv", sep = ",", eol = "\n", na = "NA", dec = ".", digits = 20, adjust = "BH")

summary(resultsTT.v)
#    (Intercept) groupControls
# -1       23387             0
# 0         2309         40057
# 1        14361             0

# Umm...that seems a bit much.  Basically everything looks like a DEG.  Probably need to go back to the beginning and make target tables, etc.

x11()
vennDiagram(resultsTT.v, include = c("up", "down"))
# Venn diagram only makes sense when comparing two groups to the control, because the values are normalized to control.  :p

# From my reading, it looks like voom {limma} really shines when comparing more than 2 groups of data.  Can eventually use to compare all 6 groups (HC null, S1M null, HC IL-6, S1M IL-6, HC IL-27, S1M IL-27)


### Hmm...but the output from voom still doesn't seem right.

# I think the best thing I can do is create a v2 of this script and base it around bioinf.wehi.edu.au/RNAseqCaseStudy/.  Right now, it's getting too big and messy to work with.


# Further notes with http://genomespot.blogspot.com/2015/01/generate-rna-seq-count-matrix-with.html

# Notes on plotting with R: http://ww2.coastal.edu/kingw/statistics/R-tutorials/graphs.html

# Once I get voom and independent filtering working, make a version of DE analysis script that's cleaned up and can be submitted to cluster - set up like a CSS script, so I only need to change project-specific variables at the beginning.  Test again on cytokine null data, and then run all samples simultaneously on the cluster (maybe submit as an overnight job).  Then on to allelic bias (which is really a way to potentially link SNPs to genes!).
## Use hashtags to mark modules.


# For general tips regarding R: http://www.r-statistics.com/2014/08/simpler-r-coding-with-pipes-the-present-and-future-of-the-magrittr-package/


#####
## voom attempt 2

options(digits = 3)

# readTargets(file="Targets.txt", path=NULL, sep="\t", row.names=NULL, quote="\"",...)
targets <- readTargets(file = "Targets.csv", sep = ",", path = NULL, row.names = NULL)
group <- factor(targets$Group)
design <- model.matrix(~group)
head(targets)


#####
# Where Gordon said in 2013-01 (http://thread.gmane.org/gmane.science.biology.informatics.conductor/45469/focus=45500): 
# "I started the edgeR project way back in 2004, before RNA-Seq technology even existed.  My thought then was that non-normal methods would be essential.  And I think that it remains true that edgeR has a clear edge over other methods for the SAGE data that was available then.

# edgeR performs very well, but has not solved all our data analysis problems for RNA-Seq data, not yet anyway.

# voom has surprised me by performing much better than I expected.  It solves the biggest problem with normal-based methods by estimating the mean-variance relationship adaptively.  It holds its size correctly in almost any circumstance, scales to very large datasets, can adaptively estimate the amount of smoothing necessary, and gives immediate access to all the downstream limma pipeline including gene set testing.

# My feeling the moment is that edgeR is superior for small counts and but that voom is safer and more reliable for noisy heterogeneous data.  Only edgeR can estimate the biological coefficient of variation (as we defined this in our 2012 paper).  But we are actively working on both methods, and are open to what we find."

