# Christin M. Hong
# Started 2015-10-27
# Soumya Raychaudhuri Lab, Harvard Medical School

# Allelic-specific expression analysis on mRNA from naive CD4+ T cells from human patients with and without a mutation in STAT1


#####

# Project Notes
    # Analyzing data from naive CD4 T cells from human patients from Hirahara-O'Shea 2015 Immunity paper.

    # Controls are healthy patients.  Cases have a STAT1 mutation that seems to inhibit dimerization with STAT3, resulting in an over-accumulation of STAT1 dimers and defective IFNg, IL-22, and IL-17 production in response to IL-1b, IL-12, and IL-23.

    # Here, "null" indicates human CD4+CD45RA(hi)CD45RO(lo) samples from PBMCs that have been stimulated with plate-coated anti-CD3 and anti-CD28 for three days.

    # "IL-6" will indicate samples that were also incubated with 50 ng/mL of IL-6 in culture, which biases towards STAT3 expression.  (STAT3 KO -> severely reduced gene expression profile).

    # "IL-27" samples were treated with 50 ng/mL of IL-27 in culture, which increases STAT1 expression.  (STAT1 KO -> IL-6 and IL-27 eliciting the same atypical expression profile = loss of specificity).

    # RNA was extracted from 1e6 cells by mirVana miRNA Isolation Kit.  200 ng was processed with TruSeq SR RNA Sample Prep Kit to form an mRNA library.  Libraries were sequenced on an Illumina HiSeq 2000 for 50 cycles, single read.


    # IMPORTANT: Authors sorted by CD4+CD45RA(hi)CD45RO(lo).  This will also collect NK cells and monocytes, so these samples are almost certainly a mixed population.


#####

# OBJECTIVES OF ANALYSIS / TOC

## RNA-ASE-R01: Plots of [reference allele count]/[total count] (shown as a percentage) for 1 library.
    # Expect a normal distribution centered around 50%, as the GATK Haplotype Caller expects an allele-expression ratio of 0.5.


## RNA-ASE-R02: Plots of [reference allele count]/[total count] for all samples


## RNA-ASE-R03: Dotplots of percentage of data where both alleles are seen per sample.
    # Expect majority of called variants to have both alleles seen, but not all, since the filtering was done after merging all the libraries from a given individual and this is analyzing 1 library ay a time.

## RNA-ASE-R04: Coverage per site

#####

# INFRASTRUCTURE

# Import libraries
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(gtable)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(knitr)
library(stringr)
library(plyr)
library(qqman)

# Options
stringsAsFactors=FALSE

# Set working dir
pathWD <- "~/aa_Raychaudhuri-Lab/data/2015-08_RNA_Hirahara-OShea_stat1/ASE"
setwd(pathWD)


#####

# Downloaded final data files from ASE pipeline with the following command:
# rsync -avzs --progress ch961@erisone.partners.org:/data/srlab/cmhong/data/2015-09_RNA_Hirahara-OShea_stat1-ASE/aaa_finalASEFiles /home/christin/aa_Raychaudhuri-Lab/data/2015-08_RNA_Hirahara-OShea_stat1/ASE


#####

# CUT THIS LATER

# # Importing result from 1 library from ASE pipeline (which are saved to a tab-delimited TXT file)
# dataASE <- read.delim("aaa_finalASEFiles/ASE02f13_5SamASE_allASE01_aa_stat1-ind.1_SRR1786581.sra.ASE.txt",
#     header = TRUE, 
#     sep = "\t", 
#     quote = "\"",
#     dec = ".", 
#     fill = TRUE, 
#     comment.char = "")
# 
# head(dataASE)
# names(dataASE)

#####

# Importing the concatenated file of all results from ASE pipeline
dataAllASE <- read.delim("aaa_finalASEFiles/ASE02f17_mergedASEOutputWithIndLib.txt",
                      header = TRUE, 
                      sep = "\t", 
                      quote = "\"",
                      dec = ".", 
                      fill = TRUE, 
                      comment.char = "")

head(dataAllASE)
tail(dataAllASE)
names(dataAllASE)


#####

# DATA ANALYSES FOR SINGLE LIBRARY


##### ADJUST THIS LATER TO WORK WITH CONCATENATED FILE #####


# RNA-ASE-R01: Histogram of [reference allele count]/[total count]
dataPercentRef <- ifelse(dataASE$TOTAL_COUNT == 0, NA, 100*(dataASE$REF_COUNT/dataASE$TOTAL_COUNT))

dataPercentRef.df <- data.frame(dataPercentRef)

dataPercentRef.df$RSID <- dataASE$RSID

head(dataPercentRef.df)

x11()

pdf("ASE04f01_RefAlleleCountPercentage.pdf", width = 8, height = 8, useDingbats = F);

ggplot(dataPercentRef.df, aes(x = dataPercentRef)) + 
    geom_histogram(binwidth = 1, color = "#FFFFFF") + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8")) + 
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nHistogram of reference allele percentage from\nIndividual 1, Sample 1\n") + 
    xlab("Reference allele percentage") + ylab("Number of variants")


ggplot(dataPercentRef.df, aes(x = dataPercentRef)) + 
    geom_density(kernel = "gaussian", alpha = 0.25) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8")) +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nGaussian density plot of reference allele percentage from\nIndividual 1, Sample 1\n") + 
    xlab("Reference allele percentage") + 
    ylab("Probability of a variant showing x")

# Density plots can later be overlaid (using faceting) to quickly check batches of samples for normal distribution around 50%.

dev.off()


#####

# DATA ANALYSES FOR MULTIPLE LIBRARIES
dataAllASEv2 <- dataAllASE
tail(dataAllASEv2)


# DATA FORMATTING AND GROUPING
# Check that numeric columns are numeric (imported data under global "stringsAsFactors=FALSE" setting)
    # Expected numeric columns (should be TRUE): 4 (POS), 6-10 (Both alleles seen, ref count, nonref count, total count, and individual ID #)
sapply(dataAllASEv2, is.numeric)


# Add column for healthy vs. STAT1 mutants
dataAllASEv2$Group <- 
    ifelse(dataAllASEv2$IndID >= 1 & dataAllASEv2$IndID <= 4, 
           'HC', 
           'STAT1')



# Add column for cytokine stimulation (trickier)
    # 1) Split string to add column with just library numbers
dataAllASEv2.lib <- data.frame(do.call('rbind', strsplit(as.character(dataAllASEv2$Library), 'SRR', fixed=TRUE)))

dataAllASEv2$LibNum <- dataAllASEv2.lib$X2


    # 2) Convert LibNum to numeric
dataAllASEv2$LibNum <- sapply(dataAllASEv2$LibNum, as.numeric)

sapply(dataAllASEv2, is.numeric) # Check


    # 3) Add cytokine values
libNull <- seq.int(1, 27, 3)
libIL6 <-seq.int(2, 27, 3)
libIL27 <- seq.int(3, 27, 3)

dataAllASEv2$Treatment <- 
    ifelse(dataAllASEv2$LibNum %in% libNull, 'NULL', 
           ifelse(dataAllASEv2$LibNum %in% libIL6, 'IL6', 
                  'IL27'))

head(dataAllASEv2)
tail(dataAllASEv2)


    # 4) Check that groups are properly distributed
table(dataAllASEv2$IndID)
table(dataAllASEv2$LibNum)
table(dataAllASEv2$Treatment)
dplyr::filter(dataAllASEv2, LibNum == 6)


# Already in long format = ready for ggplot.


#####

# RNA-ASE-R02: [reference allele count]/[total count] for all samples

# Add column with the reference allele percentage seen
dataAllASEv2$PercentRefAllele <- 
    ifelse(dataAllASEv2$TOTAL_COUNT==0, NA, 
           100*(dataAllASEv2$REF_COUNT/dataAllASEv2$TOTAL_COUNT))


x11()

pdf("ASE04f02_RefAlleleCountPercentageAll.pdf",width = 8, height = 8, useDingbats = F);


# Histogram of total reference allele percentages
ggplot(dataAllASEv2, aes(x = PercentRefAllele)) + 
    geom_histogram(binwidth = 1, color = "#FFFFFF") + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8")) + 
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nHistogram of reference allele percentage,\nAll samples\n") + 
    xlab("Reference allele percentage, binwidth = 1") + 
    ylab("Number of variants")


# Histograms of reference allele percentages faceted by library
ggplot(dataAllASEv2, aes(x = PercentRefAllele)) + 
    geom_histogram(binwidth = 1, color = "#1f1f1f") + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8")) + 
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nHistograms of reference allele percentage,\nFaceted by library\n") + 
    xlab("Reference allele percentage, binwidth = 1") + 
    ylab("Number of variants") + 
    facet_wrap(~ LibNum)


# Histograms of reference allele percentages faceted by individual
ggplot(dataAllASEv2, aes(x = PercentRefAllele)) + 
    geom_histogram(binwidth = 1, color = "#1f1f1f") + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8")) + 
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nHistograms of reference allele percentage,\nFaceted by individual\n") + 
    xlab("Reference allele percentage, binwidth = 1") + 
    ylab("Number of variants") + 
    facet_wrap(~ IndID)


# Cumulative density plot
ggplot(dataAllASEv2, aes(x = PercentRefAllele)) + 
    geom_density(kernel = "gaussian", alpha = 0.05) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8"), legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nGaussian density plot of reference allele percentage,\nAll libraries\n") + 
    xlab("Reference allele percentage") + 
    ylab("Probability of a variant showing x")


# Density plots by library
ggplot(dataAllASEv2, aes(x = PercentRefAllele, color = INDIVIDUAL, fill = INDIVIDUAL)) + 
    geom_density(kernel = "gaussian", alpha = 0.05) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8"), legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nGaussian density plot of reference allele percentage,\nBy library\n") + 
    xlab("Reference allele percentage") + 
    ylab("Probability of a variant showing x") + 
    facet_wrap(~ INDIVIDUAL)


# Density plots by library, overlaid
ggplot(dataAllASEv2, aes(x = PercentRefAllele, color = INDIVIDUAL, fill = INDIVIDUAL)) + 
    geom_density(kernel = "gaussian", alpha = 0.05) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8"), legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nGaussian density plot of reference allele percentage,\nBy library\n") + 
    xlab("Reference allele percentage") + 
    ylab("Probability of a variant showing x")


# Density plots by individual
ggplot(dataAllASEv2, aes(x = PercentRefAllele, fill = IndID, color = IndID)) + 
    geom_density(kernel = "gaussian", alpha = 0.05) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8"), legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nGaussian density plot of reference allele percentage,\nBy individual\n") + 
    xlab("Reference allele percentage") + 
    ylab("Probability of a variant showing x") + 
    facet_wrap(~ IndID)


# Density plot by individual, overlaid
ggplot(dataAllASEv2, aes(x = PercentRefAllele, fill = IndID, color = IndID)) +  
    geom_density(kernel = "gaussian", alpha = 0.05, aes(group=IndID)) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8"), legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nGaussian density plot of reference allele percentage,\nBy individual\n") + 
    xlab("Reference allele percentage") + 
    ylab("Probability of a variant showing x")


# Density plot by treatment, overlaid
ggplot(dataAllASEv2, aes(x = PercentRefAllele, fill = Treatment, color = Treatment)) +  
    geom_density(kernel = "gaussian", alpha = 0.05, aes(group=Treatment)) + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8"), legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nGaussian density plot of reference allele percentage,\nBy cytokine treatment\n") + 
    xlab("Reference allele percentage") + 
    ylab("Probability of a variant showing x")


# Boxplot by library
ggplot(dataAllASEv2, aes(x = INDIVIDUAL, y = PercentRefAllele)) +
    geom_boxplot(fill="#00bbff") + 
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8"), legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nReference allele percentages by library\n") +
    xlab("Library") + 
    ylab("Reference allele (%)") +
    scale_x_discrete(labels = c(1:27))

dev.off()


#####

# RNA-ASE-R03: Percentage of data where both alleles are seen

# Get total number of variants seen for each library
dataVarBothAllelesByLib <- table(dataAllASEv2$LibNum)

# Convert to data table 
dataVarBothAllelesByLib <- data.frame(dataVarBothAllelesByLib)

# Rename columns with {plyr}
dataVarBothAllelesByLib <- rename(dataVarBothAllelesByLib, c("Var1"="LibNum", "Freq"="NumOfVar"))

head(dataVarBothAllelesByLib)

# Sum number of variants where both alleles were seen (recorded as 1, vs. 0 when only one allele was seen) by library
dataVarBothAllelesByLib.both <- aggregate(dataAllASEv2$BOTH_ALLELES_SEEN, by=list(LibNum=dataAllASEv2$LibNum), FUN=sum)

# Rename columns
dataVarBothAllelesByLib.both <- rename(dataVarBothAllelesByLib.both, c("x"="BothAllelesSeen"))

head(dataVarBothAllelesByLib.both)

# Add column to data table
dataVarBothAllelesByLib$BothAllelesSeen <- dataVarBothAllelesByLib.both$BothAllelesSeen

head(dataVarBothAllelesByLib)

# Get percentages of both alleles seen
dataVarBothAllelesByLib$PercentBothAllele <-
    ifelse(dataVarBothAllelesByLib$NumOfVar == 0, NA, 100*(dataVarBothAllelesByLib$BothAllelesSeen/dataVarBothAllelesByLib$NumOfVar))


# 3) Add cytokine values
libNull <- seq.int(1, 27, 3)
libIL6 <-seq.int(2, 27, 3)
libIL27<- seq.int(3, 27, 3)

dataVarBothAllelesByLib$Treatment <- 
    ifelse(dataVarBothAllelesByLib$LibNum %in% libNull, 'NULL', 
           ifelse(dataVarBothAllelesByLib$LibNum %in% libIL6, 'IL6', 'IL27'))

levels(dataVarBothAllelesByLib$Treatment)


# Plot percentages as scatterplot
x11()

pdf("ASE04f03_PercentAlleleSeenByLibrary.pdf",width = 8, height = 8, useDingbats = F);

head(dataVarBothAllelesByLib)

# Single color, 0-100%
ggplot(dataVarBothAllelesByLib, aes(x = LibNum, y = PercentBothAllele)) + 
    geom_point(size = 4, color = "#00bbff") +
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), panel.grid.minor = element_line(colour = "#D8D8D8"), legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nPercent of variants with both alleles seen,\nBy library (3 libraries/individual)\n") + 
    xlab("Library") + 
    ylab("Variants with both alleles seen (%)") +
    ylim(0,100)


# Colored by treatment, 99-100%
ggplot(dataVarBothAllelesByLib, aes(x = LibNum, y = PercentBothAllele)) + 
    geom_point(size = 4, aes(color = factor(Treatment))) +
    scale_color_manual(values = c("NULL"="#005e80", "IL6"="#00bbff", "IL27"="#80ddff"), name = "Cytokine treatment:", guide = guide_legend(reverse = TRUE)) +
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), 
          panel.grid.minor = element_line(colour = "#D8D8D8"), 
          legend.position = "top") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nPercent of variants with both alleles seen,\nBy library (3 libraries/individual)\n") + 
    xlab("Library") + 
    ylab("Variants with both alleles seen (%)") +
    ylim(99,100)


# Colored by treatment and labeled, 99-100%
ggplot(dataVarBothAllelesByLib, aes(x = LibNum, y = PercentBothAllele)) + 
    geom_point(size = 4, aes(color = factor(Treatment))) +
    scale_color_manual(values = c("NULL"="#005e80", "IL6"="#00bbff", "IL27"="#80ddff"), name = "Cytokine treatment:", guide = guide_legend(reverse = TRUE)) +
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), 
          panel.grid.minor = element_line(colour = "#D8D8D8"), 
          legend.position = "top") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nPercent of variants with both alleles seen,\nBy library (3 libraries/individual)\n") + 
    xlab("Library") + 
    ylab("Variants with both alleles seen (%)") +
    ylim(99,100) +
    geom_text(label = dataVarBothAllelesByLib$NumOfVar, vjust = -1, size = 3) +
    annotate("text", x = 1, y = 100, label = "Label: Total number of variants per library", size = 3.5, hjust = 0, vjust = 0)


dev.off()


#####

# RNA-ASE-R04: Mapping variants

# Source data table
dataAllASEv3 <- dataAllASEv2

# Reorder factors to plot in correct order
dataAllASEv3$CHR <- 
    factor(dataAllASEv3$CHR, levels = c(1:22,"X","M"))

    levels(dataAllASEv3$CHR)

    
dataAllASEv3$Treatment <- 
    factor(dataAllASEv3$Treatment, levels = c("NULL", "IL6", "IL27"))

    levels(dataAllASEv3$Treatment)


# Generate log2-tranformed counts for better visualization
dataAllASEv3$Log2TotalCounts <- log(dataAllASEv3$TOTAL_COUNT + 1, 2)

dataAllASEv3$Log2RefCounts <- log(dataAllASEv3$REF_COUNT + 1, 2)


# Check data table
tail(dataAllASEv3)


x11()

pdf("ASE04f04_VariantMapping.pdf",width = 12, height = 8, useDingbats = F);

# Dotplot of reads/site by chromosome from all libraries
ggplot(dataAllASEv3, aes(x = POS, y = TOTAL_COUNT)) +
    geom_point(size = 1) +
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), 
          panel.grid.minor = element_line(colour = "#D8D8D8"), 
          legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nTotal read counts per variant, all libraries\nFaceted by chromosome\n") +
    xlab("Position") + 
    ylab("Total read count") +
    facet_wrap(~ CHR, ncol = 3)


# Dotplot of log2-transformed reference reads/site by chromosome from all libraries
ggplot(dataAllASEv3, aes(x = POS, y = Log2TotalCounts)) +
    geom_point(size = 1) +
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), 
          panel.grid.minor = element_line(colour = "#D8D8D8"), 
          legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nLog2-transformed total read counts per variant, all libraries\nFaceted by chromosome\n") +
    xlab("Position") + 
    ylab(expression(log[2](count+1))) +
    facet_wrap(~ CHR, ncol = 3)


# Dotplot of log2-transformed reads/site by chromosome from all libraries
ggplot(dataAllASEv3, aes(x = POS, y = Log2RefCounts)) +
    geom_point(size = 1) +
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), 
          panel.grid.minor = element_line(colour = "#D8D8D8"), 
          legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nLog2-transformed reference read counts per variant, all libraries\nFaceted by chromosome\n") +
    xlab("Position") + 
    ylab(expression(log[2](count+1))) +
    facet_wrap(~ CHR, ncol = 3)


# Plot log2-transformed reference read counts by individual and cytokine treatment
ggplot(dataAllASEv3, aes(x = RSID, y = Log2RefCounts)) +
    geom_point(size = 1) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nLog2-transformed reference read counts, all libraries\nFaceted by individual and cytokine treatment\n") +
    xlab("Position") + 
    ylab(expression(log[2](count+1))) +
    facet_grid(Treatment ~ IndID)



# Filter for PercentRefAllele < 25th percentile or > 75th percentile (average by library) and replot
quantile(dataAllASEv3$PercentRefAllele)

dataAllASEv3.filt <- subset(dataAllASEv3, PercentRefAllele < 45.3125 | PercentRefAllele > 60.0000)
    
head(dataAllASEv3.filt)


ggplot(dataAllASEv3.filt, aes(x = POS, y = Log2RefCounts)) +
    geom_point(size = 1) +
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), 
          panel.grid.minor = element_line(colour = "#D8D8D8"), 
          legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nLog2-transformed reference read counts per variant, all libraries\nPercentRefAllele < 1st quantile or > 3rd quantile\nFaceted by chromosome\n") +
    xlab("Position") + 
    ylab(expression(log[2](count+1))) +
    facet_wrap(~ CHR, ncol = 3)


ggplot(dataAllASEv3.filt, aes(x = RSID, y = Log2RefCounts)) +
    geom_point(size = 1) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nLog2-transformed reference read counts, all libraries\nPercentRefAllele < 1st quantile or > 3rd quantile\nFaceted by individual and cytokine treatment\n") +
    xlab("Position") + 
    ylab(expression(log[2](count+1))) +
    facet_grid(Treatment ~ IndID)


# Filter for PercentRefAllele < 10th percentile or > 90th percentile (average by library) and replot
quantile(dataAllASEv3$PercentRefAllele, c(0.10, 0.90))

dataAllASEv3.filt2 <- subset(dataAllASEv3, PercentRefAllele < 36.84211 | PercentRefAllele > 68.96552)



ggplot(dataAllASEv3.filt2, aes(x = POS, y = Log2RefCounts)) +
    geom_point(size = 1) +
    theme_bw() + 
    theme(panel.grid.major = element_line(colour = "#C8C8C8"), 
          panel.grid.minor = element_line(colour = "#D8D8D8"), 
          legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nLog2-transformed reference read counts per variant, all libraries\nPercentRefAllele < 10th percentile or > 90th percentile\nFaceted by chromosome\n") +
    xlab("Position") + 
    ylab(expression(log[2](count+1))) +
    facet_wrap(~ CHR, ncol = 6)


ggplot(dataAllASEv3.filt2, aes(x = RSID, y = Log2RefCounts)) +
    geom_point(size = 1) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nLog2-transformed reference read counts, all libraries\nPercentRefAllele < 10th percentile or > 90th percentile\nFaceted by individual and cytokine treatment\n") +
    xlab("Position") + 
    ylab(expression(log[2](count+1))) +
    facet_grid(Treatment ~ IndID)


# Plot y = PercentRefAllele
head(dataAllASEv3)

ggplot(dataAllASEv3, aes(x = RSID, y = PercentRefAllele)) +
    geom_point(size = 1) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          legend.position = "none") +
    ggtitle("RNA-seq ASE from O'Shea Immunity 2015 paper:\nPosition vs. Reference allele expression\nFaceted by individual and cytokine treatment\n") +
    xlab("Position") + 
    ylab("Percent of reads mapping to reference allele") +
    facet_grid(Treatment ~ IndID)

dev.off()


#####

# RNA-ASE-R05: Looking at ASE
head(dataAllASEv3)


# Sort by PercentRefAllele by chromosome by library
dataAllASEv3.sort <- dataAllASEv3[order(dataAllASEv3$LibNum, dataAllASEv3$CHR, dataAllASEv3$PercentRefAllele),]

head(dataAllASEv3.sort)


x11()

#pdf("ASE04f05_ASE.pdf",width = 8, height = 8, useDingbats = F);


#dev.off()