# Christin M. Hong
# Lab PI: Soumya Raychaudhuri, Harvard Medical School
# Last updated: 2017-07


#### Project Notes ####

# Source: GaP Long Island registry
# Platform: QC Illumina Global Screening Array (genotype data)

# Script purpose(s):
# 1. Plot PCA results (adapted from Dr. Yang Luo).
# 2. Analyze for European ancestry by comparing with 1000 Genomes.  Individuals were determined as being European by being within 2 SD of the Mahalanobis distance for Europeans from PC1 and PC2 in 1000 Genomes. 
# 3. Analyze and plot gender, age, and number of European individuals per category.
# 4. Select candidates for Immunomics project.

# Individuals were determined as being European by being within 2 SD of the Mahalanobis distance for Europeans from PC1 and PC2 in 1000 Genomes.


#### Infrastructure ####
wd <- "~/00_SR-lab-data/immunomics/genotypeQC_GaP-GSA_2017/analysis/"
setwd(wd)


# import libraries
library(ggplot2)
library(dplyr)
library(plyr)
library(Cairo)
library(gridExtra)


CairoX11()

#### PCA of main dataset ####
pdf(
  "PDF-pca-projectData.pdf",
  width = 10,
  height = 7,
  useDingbats = F
)


evec <-
  read.table("pca-gap201707.eigenvec",
             h = F,
             stringsAsFactors = F)

df <- data.frame(evec)

eval <- scan("pca-gap201707.eigenval", as.numeric())
val <- eval / sum(eval) * 100

#pheno <- read.table("gap201707-clean-maf.fam", header = FALSE, stringsAsFactors = F)
#df$Pheno <- ifelse(pheno[match(df$V2,pheno$V2),6]==1, "Controls","Cases")


### Plot and save as PDF ###
plot(
  val,
  ylab = "Proportion of genetic variance explained (%)",
  ylim = c(0, 60),
  xlab = "Principal components",
  pch = 20,
  lwd = 2
)

ggplot(df, aes(x = V3, y = V4)) +
  geom_point() +
  xlab(paste("PC1 (", round(val[1], 2), "%)")) +
  ylab(paste("PC2 (", round(val[2], 2), "%)")) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        legend.position = c(0.9, 0.9))

ggplot(df, aes(x = V5, y = V6)) +
  geom_point() +
  xlab(paste("PC3 (", round(val[3], 2), "%)")) +
  ylab(paste("PC4 (", round(val[4], 2), "%)")) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        legend.position = c(0.9, 0.9))

dev.off()



#### PCA of data merged with 1000 Genomes ####
evec2 <-
  read.table("pca-gap201707-merged.eigenvec",
             h = F,
             stringsAsFactors = F)

dfMerged <- data.frame(evec2)

eval2 <- scan("pca-gap201707-merged.eigenval", as.numeric())
val2 <- eval2 / sum(eval2) * 100

pop <-
  read.table(
    "/Users/chong/00_SR-lab-data/aaa_resources/ref/integrated_call_samples_v3.20130502.ALL.panel",
    h = T,
    stringsAsFactors = F
  )

dfMerged$POP <- pop[match(dfMerged$V1, pop$sample),]$pop
dfMerged$REGION <- pop[match(dfMerged$V1, pop$sample),]$super_pop

dfMerged[is.na(dfMerged$POP),]$REGION <- "GAP"
dfMerged[is.na(dfMerged$POP),]$POP <- "GAP"



### Plot and save as PDF ###
pdf(
  "PDF-pca-projectMergedG1K.pdf",
  width = 10,
  height = 7,
  useDingbats = F
)

plot(
  val2,
  ylab = "Proportion of genetic variance explained (%)",
  ylim = c(0, 60),
  xlab = "Principal components",
  pch = 20,
  lwd = 2
)

ggplot(dfMerged, aes(x = V3, y = V4, colour = REGION)) +
  geom_point() +
  xlab(paste("PC1 (", round(val2[1], 2), "%)")) +
  ylab(paste("PC2 (", round(val2[2], 2), "%)")) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "dodgerblue4",
      "darkolivegreen4",
      "red",
      "goldenrod1",
      "darkorchid3",
      "pink"
    )
  )


ggplot(dfMerged, aes(x = V5, y = V6, colour = REGION)) +
  geom_point() +
  xlab(paste("PC3 (", round(val2[3], 2), "%)")) +
  ylab(paste("PC4 (", round(val2[4], 2), "%)")) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "dodgerblue4",
      "darkolivegreen4",
      "red",
      "goldenrod1",
      "darkorchid3",
      "pink"
    )
  )


dfM.1 <- dfMerged[dfMerged$REGION == "AMR" |
                    dfMerged$REGION == "GAP",]

ggplot(dfM.1, aes(x = V3, y = V4, colour = POP)) +
  geom_point() +
  xlab(paste("PC1 (", round(val2[1], 2), "%)")) +
  ylab(paste("PC2 (", round(val2[2], 2), "%)")) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "dodgerblue4",
      "darkolivegreen4",
      "darkorchid3",
      "goldenrod1",
      "red",
      "orange",
      "darkred",
      "grey"
    )
  )

dev.off()



#### Identify Europeans in GaP data by Mahalanobis distance by European PCs ####
## Method suggested by Soumya, advice on code from Kam

# Find names of regions
levels(factor(dfMerged$REGION))

# Extract 1kG European population
dfEur <- subset(dfMerged, REGION == "EUR")

# Calculate squared Mahalanobis distance based on 1kG Europeans
dfMerged.mt <- as.matrix(dfMerged[, 3:4])
dfEur.mt <- as.matrix(dfEur[, 3:4])

dfMerged$MahalanobisDistance <-
  mahalanobis(dfMerged.mt, colMeans(dfEur.mt), cov(dfEur.mt))

head(dfMerged)
hist(dfMerged$MahalanobisDistance)


# Subset 1kG European population with Mahalanobis distances
dfEur.maha <- subset(dfMerged, REGION == "EUR")

hist(dfEur.maha$MahalanobisDistance)


# Get + SD from mean for max Mahalanobis distances for Europeans
## From Kam:
# i think you want to exclude individuals that are far away from europeans, but not exclude individuals that are very close to europeans
# in other words, if Mahalanobis distance is near 0, then the person is very close to the europeans
maxEur <-
  mean(dfEur.maha$MahalanobisDistance) + 2 * sd(dfEur.maha$MahalanobisDistance)


# Filter for GaP samples within Mahal range
dfGapEur <-
  subset(dfMerged, REGION == "GAP" & MahalanobisDistance < maxEur)

hist(dfGapEur$MahalanobisDistance) # Check that filters worked properly


# Change REGION for these samples to reflect being European
dfGapEur$REGION <- gsub('GAP', 'GAP_EUR', dfGapEur$REGION)

# Number of samples
nrow(dfGapEur)


#### Plot 1000 Genomes vs. GaP ####
gg1 <- ggplot() +
  geom_point(
    data = dfMerged,
    aes(x = V3, y = V4, colour = REGION),
    alpha = 0.75,
    size = 0.75
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  xlab(paste("PC1 (", round(val2[1], 2), "%)")) +
  ylab(paste("PC2 (", round(val2[2], 2), "%)")) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "AFR" = "goldenrod1",
      "AMR" = "darkolivegreen4",
      "EAS" = "grey",
      "EUR" = "blue",
      "GAP" = "darkorchid3",
      "SAS" = "pink"
    )
  )


gg2 <- ggplot() +
  geom_point(
    data = dfMerged,
    aes(x = V3, y = V4, colour = REGION),
    alpha = 0.75,
    size = 0.75
  ) +
  geom_point(
    data = dfGapEur,
    aes(x = V3, y = V4, colour = REGION),
    alpha = 1,
    size = 1.5
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  xlab(paste("PC1 (", round(val2[1], 2), "%)")) +
  ylab(paste("PC2 (", round(val2[2], 2), "%)")) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "GAP_EUR" = "red",
      "AFR" = "goldenrod1",
      "AMR" = "darkolivegreen4",
      "EAS" = "grey",
      "EUR" = "blue",
      "GAP" = "darkorchid3",
      "SAS" = "pink"
    )
  )



#### Subset for gender and age ####
# First need to get original Excel file from the GaP and save the Sex and Year of Birth worksheet as a CSV.  I saved in the analysis folder for ease of use.

dfGapAnno <-
  read.csv("GAP July 2017_Sex_YOB.csv",
           h = TRUE,
           stringsAsFactors = TRUE)

class(dfGapAnno) # Data frame
dfGapAnno.1 <- dfGapAnno[, 3:5] # Keep only relevant columns


# Renaming columns
names(dfGapAnno.1)[names(dfGapAnno.1) == "X.Master.Query.5.Subject.ID"] <-
  "V2"


# Identify duplicates
nrow(dfGapAnno.1[duplicated(dfGapAnno.1[c("V2")]), ]) # 13 dup IDs


# saving duplicates, just in case
allDup <- function (value)
{
  duplicated(value) | duplicated(value, fromLast = TRUE)
}

dfGapAnno.dup <- dfGapAnno.1[allDup(dfGapAnno.1[c("V2")]), ]

nrow(dfGapAnno.dup)
nrow(unique(dfGapAnno.dup))


# Dedup data
dfGapAnno.dedup <- dfGapAnno.1[!duplicated(dfGapAnno.1[c("V2")]), ]

nrow(dfGapAnno.1) - nrow(dfGapAnno.dedup) # 13, cool


# Merge annotation data to GAP_EUR by ID
dfGapEurAnno <- left_join(dfGapEur, dfGapAnno.dedup, by = c("V2"))

head(dfGapEurAnno)
nrow(dfGapEurAnno)
nrow(dfGapEur)


# Make consistent gender labels
qcGap.1 <- dfGapEurAnno

levels(factor(qcGap.1$Sex))


# Remove NoSex individuals
filter(qcGap.1, Sex == "") # 3 samples with no sex

qcGap.2 <- filter(qcGap.1,!Sex == "")

nrow(qcGap.2)
levels(factor(qcGap.2$Sex))


# Consolidate and order gender levels
qcGap.2$Sex <-
  `levels<-`(factor(qcGap.2$Sex), list(
    Male = c("M", "Male"),
    Female = c("F", "Female")
  ))

levels(factor(qcGap.2$Sex))

nrow(filter(qcGap.2, Sex == "Male"))
nrow(filter(qcGap.2, Sex == "Female"))

nrow(filter(qcGap.2, Sex == "Male")) + nrow(filter(qcGap.2, Sex == "Female"))


head(qcGap.2)


## Convert birth year to numeric
qcGap.3 <- qcGap.2

qcGap.3$Birth.Year <- as.numeric(as.character(qcGap.3$Birth.Year))
# Got warning about NAs, checked data to see where they're coming from

# Checking NAs
GapNoAge <- qcGap.3[is.na(qcGap.3$Birth.Year),]
nrow(GapNoAge)

GapNoAgeOrig <-
  filter(qcGap.2, V2 %in% GapNoAge$V2) # Looks like these samples were mainly excluded due to having month/year format instead of only year, plus a few N/As.  Hmm...could correct, but since it's a small number of samples, for now I'll simply exclude.

qcGap.3 <- filter(qcGap.2,!V2 %in% GapNoAge$V2)

levels(factor(qcGap.3$Birth.Year))
nrow(qcGap.3)
head(qcGap.3)

qcGap.3$Birth.Year <- as.numeric(as.character(qcGap.3$Birth.Year))


#### Plot by age vs. gender ####
gg3 <- ggplot(data = qcGap.3, aes(x = Birth.Year, fill = Sex)) +
  geom_histogram(position = position_dodge(width = 6),
                 breaks = seq(1912, 2002, 10)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_minimal() +
  geom_text(
    stat = "bin",
    position = position_dodge(width = 6),
    breaks = seq(1912, 2002, 10),
    aes(label = ..count.., fontface = "bold"),
    vjust = -1
  ) +
  scale_x_continuous(breaks = seq(1912, 2002, 10)) +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Birth decade vs. Number of GaP Europeans, by gender", subtitle = "Within 2 SD of the Mahalanobis distance for 1kG Europeans by PC1 & PC2 \nBirth Year of 1992 = ~25 years old in 2017")


gg4 <- ggplot(data = qcGap.3, aes(x = Birth.Year, fill = Sex)) +
  geom_histogram(position = position_dodge(width = 3),
                 breaks = seq(1972, 2002, 5)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_minimal() +
  geom_text(
    stat = "bin",
    position = position_dodge(width = 3),
    breaks = seq(1972, 2002, 5),
    aes(label = ..count.., fontface = "bold"),
    vjust = -1
  ) +
  scale_x_continuous(breaks = seq(1972, 2002, 5)) +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Birth decade vs. Number of GaP Europeans, by gender", subtitle = "Within 2 SD of the Mahalanobis distance for 1kG Europeans by PC1 & PC2 \nBirth Year of 1992 = ~25 years old in 2017")



#### Save GAP_EUR plots in PDF ####
pdf(
  "999_output_GapEur-MahalDist2SD-pcaAncestryGenderAge.pdf",
  width = 11,
  height = 8.5,
  useDingbats = F
)


grid.arrange(gg1, gg2, ncol = 2)

grid.arrange(gg3, gg4, ncol = 2)

dev.off()



#### Selection of potential GaP participants ####
# From Soumya (slack): Ask for health info/phenotype data on all males and females of European descent and within 25-40 age range before proceeding with recruitment.  Will exclude those on drugs or w/ immune dz - may end up with female>male in our final numbers.


nrow(filter(qcGap.3, Sex == "Female" & Birth.Year >= 1977 & Birth.Year <= 1992))
# This is 196 because ggplot actually does Birth.Year > 1977 & Birth.Year <= 1992, whereas I'm including both end years.

nrow(filter(qcGap.3, Sex == "Male" & Birth.Year >= 1977 & Birth.Year <= 1992))


nrow(filter(qcGap.3, Sex == "Female" & Birth.Year >= 1977 & Birth.Year <= 1992)) + nrow(filter(qcGap.3, Sex == "Male" & Birth.Year >= 1977 & Birth.Year <= 1992))


# Extract FIDs and IIDs of GAP_EUR samples within 25-40 years of age
potentialGapInd <- filter(qcGap.3, Birth.Year >= 1977 & Birth.Year <= 1992)

nrow(potentialGapInd) # 311

dfGapEurIIDs <- potentialGapInd[, 2]

write.table(
  dfGapEurIIDs,
  file = "999_output_GaPEurIIDs-25to40.txt",
  sep = ",",
  row.names = F,
  col.names = F,
  quote = F
)


# Extract info on GAP_EUR samples
head(qcGap.3)
dfGapEurIIDs.All <- as.data.frame(qcGap.3[, 2])

write.table(
  dfGapEurIIDs.All,
  file = "999_output_GaPEurIIDs-all.txt",
  sep = " ",
  row.names = F,
  col.names = F,
  quote = F
)

head(dfGapEurIIDs.All)


##### Checking IID for Maria #####
head(dfMerged)

MariaD <- filter(dfMerged, V2 == "TB03073703")

filter(dfGapEur, V2 == "TB03073703")


ggplot() +
  geom_point(
    data = dfMerged,
    aes(x = V3, y = V4, colour = REGION),
    alpha = 0.75,
    size = 0.75
  ) +
  geom_point(
    data = dfGapEur,
    aes(x = V3, y = V4, colour = REGION),
    alpha = 1,
    size = 0.75
  ) +
  geom_point(
    data = MariaD,
    aes(x = V3, y = V4),
    alpha = 0.5,
    size = 3
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  xlab(paste("PC1 (", round(val2[1], 2), "%)")) +
  ylab(paste("PC2 (", round(val2[2], 2), "%)")) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "GAP_EUR" = "red",
      "AFR" = "goldenrod1",
      "AMR" = "darkolivegreen4",
      "EAS" = "grey",
      "EUR" = "blue",
      "GAP" = "darkorchid3",
      "SAS" = "pink"
    )
  )

ggsave("TB03073703.png")
