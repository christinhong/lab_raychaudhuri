# Christin M. Hong
# 2017-07
# PI: Soumya Raychaudhuri, Harvard Medical School

#### Project Notes ####

# This data was exported from FlowJo.  Samples have already been compensated and gated for live single lymphocytes.  Using scale value data from all compensated parameters.

# For installing Bioconductor packages: 
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggcyto", lib="/Users/chong/R/library")


#### Infrastructure ####

# Set working directory
wd <- "~/00_SR-lab-data/pfizer/data/"
setwd(wd)


#### Set seed for reproducibility ####
set.seed(200)


#### Set variables ####
# Number of cells/sample
#cells.sample <- 1500

# Columns with data of interest
coi <- c(8:11, 13:19, 22)

# Threshold for compensation
compMin = 0.01

# Min and max for outlier removal
outMin = 0.02
outMax = 0.999


#### Import libraries ####
# Analysis
library(Biobase)
library(flowCore)
library(ggcyto)
library(FlowSOM)
library(pheatmap)
library(dendsort)
library(Rtsne)

# Standard data formatting/tidying
library(data.table)
library(reshape2)
library(stringr)
library(tidyr)
library(plyr)
library(dplyr)
library(knitr)
library(httr)

# Standard visualization
library(Cairo) # For visualizing transparency in ggplots and PDFs
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(hexbin)
library(gtable)


#### Color palettes ####
colors.div <- brewer.pal(10, "PRGn") # For heatmap


#### START ####

#### Automate import, annotatation, and compensation adjustment of data ####
# These are QC processes to be performed by file
# Note: This data has already been gated for live single lymphocytes in FlowJo

# Exporting data from FlowJo with the following naming format: "export_## Baseline/Activated [Subset] D#_Live_[date]"

fnums <- 1:32

fnums.1 <- sprintf("%02d", fnums) # Adding leading zeros for 2-digit numbers
class(fnums.1)

flow.list <- list() # Create a list in which you intend to save your dfs


for(i in fnums.1) {
  tryCatch({
    name.pattern <- paste0("^export_", i)
  fname <- list.files(pattern = name.pattern)
  print(fname) # Getting filename
  
  fname.1 <- strsplit(fname, "_")
  fname.2 <- fname.1[[1]][2]
  fname.3 <- strsplit(fname.2, " ") # Getting data from filename
  #  print(fname.3[[1]][4])
  
  flow <- read.csv(fname, header = TRUE) # Importing file
  
  # Subset each sample to ensure that each sample has equal weight in data
  flow.1 <- sample_n(flow, cells.sample, replace = F)
  
  # Annotate
  flow.1$State <- fname.3[[1]][2]
  flow.1$Subset <- fname.3[[1]][3]
  flow.1$Donor <- fname.3[[1]][4]
  #  print(head(flow.1))
  
  ## Extract columns of interest
#  flow.2 <- select(flow.1, colKeep)
  flow.2 <- flow.1
  

  #### Adjusting for compensation - shifting data to > 0 ####
  floor <- apply(flow.2[coi], 2, quantile, probs = c(compMin), na.rm = F, type = 7) # Getting lower bound of data (see notes below)
  
  floor.df <- as.data.frame(floor)
  floor.df$names <- rownames(floor.df) # Turn rownames into a column
  floor.dfNeg <- filter(floor.df, floor < 0) # Extract negative values
  floor.dfNeg$shift <- abs(floor.dfNeg$floor) # Convert neg to positive
  floor.dfNeg$floor <- NULL # Drop original column
  floor.dfNeg.t <- t(floor.dfNeg) # transpose
  colnames(floor.dfNeg.t) <- as.character(unlist(floor.dfNeg.t[1,])) # turn first row into column name
  floor.dfNeg.t <- floor.dfNeg.t[-1, ] # drop first row
  floor.t.df <- data.frame(lapply(floor.dfNeg.t, type.convert), stringsAsFactors = T) # Convert to data.frame
  
  # Find matching columns (note that colnames function is same as names function)
  comp.match <- intersect(colnames(floor.t.df), colnames(flow.2)) 
  
  # Start Kam
  flow.3 <- flow.2
  
  for (col in comp.match) {
    flow.3[[col]] <- flow.2[[col]] + floor.t.df[[col]]
  }
  # End Kam
  
  flow.list[[i]] <- flow.3 # Add updated df to list of dfs
}
, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#### Compensation adjustment notes ####
# Compensation causes the signal from the negative population to shift below zero
# Adjusting by moving majority of data above 0 by column
# Quantile calculation assumes data is normally distributed, so using bottom 0.1% is expected to exclude outliers
# Can subtract to enhance contrast for rare populations/long-tails, but I feel that artificially inflates the importance of those markers by creating zeros
# Currently adding and log scaling instead of arcsin transforming because then will get min/max by column rather than over all data -> greater dynamic range for each marker



#### Concatenate list of data frames into one ####
flow.all <- rbind.fill(flow.list)

# Check that factors are correctly assigned and that all levels are present
factor(unique(flow.all$Donor))

names(flow.all)
head(flow.all)
tail(flow.all)


#### Log transform fluorescence data since fluorescence is on log scale ####
## Note: log in R is natural log, whereas FlowJo plots in log10, but I doubt that matters

flow1 <- flow.all

# Convert negative values to zero
flow1[flow1 < 0] <- 0

# Scale fluorescence values while retaining data annotations
flow1.1 <- flow1
flow1.1[, coi] <- log10(flow1[, coi] + 1) # Decided to use log10 for easier interpretation by immunologists



#### QC by donor ####
# Normalizing by donor to capture the full expression range of a marker within each donor.  Note that Chamith prefers to not convert to Z scores in order to maximize information in the original data, but FlowSOM prefers Z scores ("to ensure that the markers are equally weighted in importance").  Decided to NOT use Z scores because from histograms, it seems to artificially reduce distinction between positive and negative populations.
flow2 <- flow1.1

#### Converting annotations to numeric values so they can be preserved in data matrix ####
flow2$State <- revalue(flow2$State, 
                       c("Baseline"=0, 
                         "Activated"=1))

flow2$Subset <- revalue(flow2$Subset, 
                       c("Naive"=1, 
                         "Central"=2,
                         "Effector"=3,
                         "Regulatory"=4))

flow2$Donor <- gsub("D", "", flow2$Donor)

head(flow2)


#### Loop per donor ####
# This normalizes the possible range of expression to each donor.  Note that this may reduce donor-to-donor differences - some donors my be legitimately higher or lower range expressers - but it also reduces the risk of batch artifacts.  It basically results in comparing the change in expression relative to the same 0 for all donors, which I think is fair.
donors <- levels(factor(flow2$Donor))

donor.list <- list() # Create a list in which you intend to save your dfs


for (id in donors) {
  print(id)
  flow2.D <- flow2[ flow2[["Donor"]] == id , ]
  
  # Calculate minimum threshold per column.  (The 2 argument applies the function to columns, while 1 would apply to rows)
  min <- apply(flow2.D[, coi], 2, quantile, probs = c(outMin), na.rm = F, type = 7)
  
  ## Subtract minimum per column from each value in that column
  flow2.1 <- flow2.D
  flow2.1[, coi] <- sweep(flow2.D[, coi], 2, min, FUN = "-") 
  # So that's the numerator!
  
  
  # Find the value for 99.9% of the max value per column
  max <- apply(flow2.D[, coi], 2, quantile, probs = c(outMax), na.rm = F, type = 7)
  
  # Subtract min from 99th percentile to get denominator per column
  adjMax <- max - min
  
  # Final min/max normalized matrix
  flow2.2 <- flow2.1
  flow2.2[, coi] <- sweep(flow2.1[, coi], 2, adjMax, FUN = "/")
  
  #### Filtering/excluding outliers ####
  flow3 <- flow2.2
  print(nrow(flow3))
  
  # Removing samples with value < 0
  flow3.1 <- flow3[apply(flow3[, coi], 1, function(x) all(x >= 0)), ]
  print(nrow(flow3.1))
  
  # Removing samples with a value > 1
  flow3.2 <- flow3.1[apply(flow3.1[, coi], 1, function(x) all(x <= 1)), ]
  print(nrow(flow3.2))
  
  
  # Requiring 3 values to be >0.25 to ensure sample was well-stained
  flow3.3 <- flow3.2
  
  # Count number of values >0.25 per row
  flow3.3$highStain <- rowSums(flow3.3[, coi] > 0.25)
  
  # Keep only rows with >3 values that are >0.25
  flow3.4 <- filter(flow3.3, highStain > 3)
  
  print(nrow(flow3.4))
  
  
  
  #### Return to original values ####
  # Multiply values by adjMax
  flow4 <- flow3.4
  flow4[, coi] <- sweep(flow3.4[, coi], 2, adjMax, FUN = "*")
  
  # Add minimum back to column
  flow4.1 <- flow4
  flow4.1[, coi] <- sweep(flow4[, coi], 2, min, FUN = "+") 
  
  
  # Column center by subtracting the mean value for each column (improves heatmap and pca) 
  flow4.2 <- flow4.1
  
  flow4.mean <- apply(flow4.1[, coi], 2, mean)
  
  flow4.2[, coi] <- sweep(flow4.1[, coi], 2, flow4.mean, FUN = "-")
  
  
  # Divide by standard deviation (normalizes the distributions for each marker, convert to Z-score)
#  flow4.3 <- flow4.2
#  flow4.sd <- apply(flow4.1[, coi], 2, sd)
#  flow4.3[, coi] <- sweep(flow4.3[, coi], 2, flow4.sd, FUN = "/")
  
  
  print(tail(flow4.2)) # Check that annotations are still present!
  
  donor.list[[id]] <- flow4.2 # Add updated df to list of dfs
}



#### Concatenate list of data frames into one ####
flow5 <- rbind.fill(donor.list)

hist(flow5$CD45RO)


#### Write out to an external CSV for later ease of loading/appending ####
output <- paste0("flowCytometryQC_", Sys.Date( ), "_immunomics.csv")

write.table(x = data.frame(flow5, stringsAsFactors = T),
            file = output,
            quote = F,
            sep = ",",
            row.names = F)



#### FlowSOM ####

# Import template data for FlowSom
ff <- flowCore::read.FCS("20170802_Immunomics/Z_01 Baseline Naive D11_001.fcs")

fsom <- ReadInput(ff, 
                  compensate = FALSE, 
                  transform = FALSE, 
                  toTransform = NULL,
                  scale = FALSE,
                  scaled.center = FALSE,
                  scaled.scale = FALSE)

str(fsom)

# Build self-organizing map from columns of interest

# Select CoI
fsom$prettyColnames <- append(fsom$prettyColnames, c("State", "Subset", "Donor"))

fsom$prettyColnames

fsom$prettyColnames[c(9:15, 19:22, 24:26)]

fsom$prettyColnames[c(9:15, 19:22)]


# Rearrange QC'd data to match data of FlowSOM object
str(flow5)
fsom$prettyColnames

flow5.1 <- flow5[c("Time", "FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H", "SSC.W", "CD49d", "HLA.DR", "OX40", "PD.1", "CD38", "CD25", "CD45RO", "CTLA.4", "PI", "CD127", "CD62L", "ICOS", "CCR6", "CD40L", "CD69", "CD4", "State", "Subset", "Donor")]

# Replace FlowSOM object data with QC'd data
str(flow5.1)
head(fsom$data)

flow5.1 <- sapply(flow5.1, as.numeric)
head(flow5.1)
fsom$data <- flow5.1


# Can increase number of nodes by changing xdim and ydim.  E.g. Default is 10 x 10 = 100, but can set to 30 x 30 = 900 total.
# For 3000 cells per base-activated pair of samples, 200 nodes = 3000/200 = 15 cells per node.  So 200 nodes should be plenty for capturing any real signal.
fsom.1 <- BuildSOM(fsom, 
                   xdim = 20, 
                   ydim = 10,
                   colsToUse = c(9:15, 19:22)) 

str(fsom.1$map, max.level = 2)

# Build minimal spanning tree
fsom.2 <- BuildMST(fsom.1, tSNE = TRUE) 
# tSNE as TRUE is slower since it's running tSNE, and probably will be unnecessary later on.  But for now, I like seeing how the nodes cluster, and tSNE on 200 nodes is still pretty fast.


#### FlowSOM plots and PDF ####
## Save in PDF
pdf.fsom <- paste0("PDF_", Sys.Date( ), "_immunomicsAdd0Log10-flowSOM.pdf")

CairoPDF(pdf.fsom, 
         width = 11.5, 
         height = 8, 
         onefile = TRUE, 
         family = "Helvetica");


# Plotting by tree means that the distance between each node is informative, but difficult to interpret this way, especially b/c arms of tree can be on top of each other.  Possibly will be more informative after defining groups?  That said, may not be necessary for statistical analysis
PlotStars(fsom.2)

PlotMarker(fsom.2, "CD69")
PlotMarker(fsom.2, "CD40L")
PlotMarker(fsom.2, "HLA.DR")
PlotMarker(fsom.2, "CD25")
PlotMarker(fsom.2, "CD38")

PlotMarker(fsom.2, "CTLA.4")
PlotMarker(fsom.2, "PD.1")
PlotMarker(fsom.2, "ICOS")
PlotMarker(fsom.2, "OX40")
PlotMarker(fsom.2, "CCR6")

PlotMarker(fsom.2, "CD45RO")
PlotMarker(fsom.2, "State")
PlotMarker(fsom.2, "Subset")
PlotMarker(fsom.2, "Donor")


# Numbers each node for desired extraction
#PlotNumbers(UpdateNodeSize(fsom.2, reset = TRUE)) 


# Very cool - I like how the distance accurately reflects similarity - but it may be more interpretable to make each node a row and each marker a column and represent as a sorted heatmap.  Grid view sort of gets at it, but still not quite as clear as a heatmap for understanding expression patterns.
PlotStars(fsom.2, view = "grid")

PlotMarker(fsom.2, view = "grid", "CD69")
PlotMarker(fsom.2, view = "grid", "CD40L")
PlotMarker(fsom.2, view = "grid", "HLA.DR")
PlotMarker(fsom.2, view = "grid", "CD25")
PlotMarker(fsom.2, view = "grid", "CD38")

PlotMarker(fsom.2, view = "grid", "CTLA.4")
PlotMarker(fsom.2, view = "grid", "PD.1")
PlotMarker(fsom.2, view = "grid", "ICOS")
PlotMarker(fsom.2, view = "grid", "OX40")
PlotMarker(fsom.2, view = "grid", "CCR6")

PlotMarker(fsom.2, view = "grid", "CD45RO")
PlotMarker(fsom.2, view = "grid", "State")
PlotMarker(fsom.2, view = "grid", "Subset")
PlotMarker(fsom.2, view = "grid", "Donor")


# tSNE view
PlotStars(fsom.2, view = "tSNE")

PlotMarker(fsom.2, view = "tSNE", "CD69")
PlotMarker(fsom.2, view = "tSNE", "CD40L")
PlotMarker(fsom.2, view = "tSNE", "HLA.DR")
PlotMarker(fsom.2, view = "tSNE", "CD25")
PlotMarker(fsom.2, view = "tSNE", "CD38")

PlotMarker(fsom.2, view = "tSNE", "CTLA.4")
PlotMarker(fsom.2, view = "tSNE", "PD.1")
PlotMarker(fsom.2, view = "tSNE", "ICOS")
PlotMarker(fsom.2, view = "tSNE", "OX40")
PlotMarker(fsom.2, view = "tSNE", "CCR6")

PlotMarker(fsom.2, view = "tSNE", "CD45RO")
PlotMarker(fsom.2, view = "tSNE", "State")
PlotMarker(fsom.2, view = "tSNE", "Subset")
PlotMarker(fsom.2, view = "tSNE", "Donor")


# Metaclustering
fsom.2.meta10 <- metaClustering_consensus(fsom.2$map$codes, k = 10) # Generate metacluster assignment for each node by consensus hierarchical clustering.

fsom.2.meta40 <- metaClustering_consensus(fsom.2$map$codes, k = 40) # Generate metacluster assignment for each node by consensus hierarchical clustering.  Using k = 40 per recommendation in Weber 2016 (~1-2 minutes runtime).  Higher k =  greater chance of identifying rare cell types, and overly similar nodes can be merged downstream.

PlotStars(fsom.2, backgroundValues = as.factor(fsom.2.meta10))
PlotStars(fsom.2, backgroundValues = as.factor(fsom.2.meta40))

dev.off()

# Other clustering methods can be used with 
#MetaClustering(data, method, max = 20, nClus = NULL, ...)
#but flowSOM authors (van Gassen 2015) saw best results with consensus hierarchical clustering
#method: clustering method to use, given as a string. Options are metaClustering_consensus,metaClustering_hclust, metaClustering_kmeans,metaClustering_som
#max: Maximum number of clusters to try out

# Tested comparing between donors, but that requires comparing the raw FCS files, so dropped that part of FlowSOM.  Will implement that with MASC.


#### PDF 1: QC histograms, ViSNE, heatmap, PCA + k means clustering ####

## Save in PDF
pdf.1 <- paste0("PDF_", Sys.Date( ), "_immunomicsAdd0Log10.pdf")

CairoPDF(pdf.1, 
         width = 11.5, 
         height = 8, 
         onefile = TRUE, 
         family = "Helvetica");


# Compare distributions before and after QC/normalization
par(mfcol=c(3,4)) # mfcol plots by column (down) while mfrow plots by row (left to right)

hist(flow3$CCR6)
hist(flow4$CCR6)
hist(flow5$CCR6)

hist(flow3$CD69)
hist(flow4$CD69)
hist(flow5$CD69)

hist(flow3$ICOS)
hist(flow4$ICOS)
hist(flow5$ICOS)

hist(flow3$CD4)
hist(flow4$CD4)
hist(flow5$CD4)

hist(flow3$HLA.DR)
hist(flow4$HLA.DR)
hist(flow5$HLA.DR)

hist(flow3$OX40)
hist(flow4$OX40)
hist(flow5$OX40)

hist(flow3$PD.1)
hist(flow4$PD.1)
hist(flow5$PD.1)

hist(flow3$CD40L)
hist(flow4$CD40L)
hist(flow5$CD40L)

hist(flow3$CD38)
hist(flow4$CD38)
hist(flow5$CD38)

hist(flow3$CD25)
hist(flow4$CD25)
hist(flow5$CD25)

hist(flow3$CTLA.4)
hist(flow4$CTLA.4)
hist(flow5$CTLA.4)

hist(flow3$CD45RO)
hist(flow4$CD45RO)
hist(flow5$CD45RO)

par(mfrow=c(1,1))



#### Further QC - plotting as violin plots ####
#flow5$State <- ordered(flow5$State, levels = c("0", "1"))

levels(factor(flow5$Subset))

labels.state = c("Baseline", "Activated")
labels.subset = c("Naive", "Central Memory", "Effector Memory", "Regulatory")

ggplot(flow5, aes(x = Subset, y = CD4, fill = Donor)) +
  geom_violin(trim = FALSE,
              show.legend = T, 
              size = 0.5, 
              alpha = 0.25, 
              scale = "count",
              draw_quantiles = c(0.25, 0.5, 0.75),
              width = 0.75) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = labels.subset) +
  theme_light() +
  ggtitle("CD4")


ggplot(flow5, aes(x = Subset, y = CCR6, fill = Donor)) +
  geom_violin(trim = FALSE,
              show.legend = T, 
              size = 0.5, 
              alpha = 0.25, 
              scale = "count",
              draw_quantiles = c(0.25, 0.5, 0.75),
              width = 0.75) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = labels.subset) +
  theme_light() +
  ggtitle("CCR6")


ggplot(flow5, aes(x = Subset, y = CD69, fill = Donor)) +
  geom_violin(trim = FALSE,
              show.legend = T, 
              size = 0.5, 
              alpha = 0.25, 
              scale = "count",
              draw_quantiles = c(0.25, 0.5, 0.75),
              width = 0.75) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +  
  scale_x_discrete(labels = labels.subset) +
  theme_light() +
  ggtitle("CD69")


ggplot(flow5, aes(x = Subset, y = CD40L, fill = Donor)) +
  geom_violin(trim = FALSE,
              show.legend = T, 
              size = 0.5, 
              alpha = 0.25, 
              scale = "count",
              draw_quantiles = c(0.25, 0.5, 0.75),
              width = 0.75) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = labels.subset) +
  theme_light() +
  ggtitle("CD40L")


ggplot(flow5, aes(x = Subset, y = CD25, fill = Donor)) +
  geom_violin(trim = FALSE,
              show.legend = T, 
              size = 0.5, 
              alpha = 0.25, 
              scale = "count",
              draw_quantiles = c(0.25, 0.5, 0.75),
              width = 0.75) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = labels.subset) +
  theme_light() +
  ggtitle("CD25")


ggplot(flow5, aes(x = Subset, y = HLA.DR, fill = Donor)) +
  geom_violin(trim = FALSE,
              show.legend = T, 
              size = 0.5, 
              alpha = 0.25, 
              scale = "count",
              draw_quantiles = c(0.25, 0.5, 0.75),
              width = 0.75) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = labels.subset) +
  theme_light() +
  ggtitle("HLA-DR")


ggplot(flow5, aes(x = Subset, y = CD38, fill = Donor)) +
  geom_violin(trim = FALSE,
              show.legend = T, 
              size = 0.5, 
              alpha = 0.25, 
              scale = "count",
              draw_quantiles = c(0.25, 0.5, 0.75),
              width = 0.75) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = labels.subset) +
  theme_light() +
  ggtitle("CD38")


ggplot(flow5, aes(x = Subset, y = ICOS, fill = Donor)) +
  geom_violin(trim = FALSE,
              show.legend = T, 
              size = 0.5, 
              alpha = 0.25, 
              scale = "count",
              draw_quantiles = c(0.25, 0.5, 0.75),
              width = 0.75) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = labels.subset) +
  theme_light() +
  ggtitle("ICOS")


ggplot(flow5, aes(x = Subset, y = CTLA.4, fill = Donor)) +
  geom_violin(trim = FALSE,
              show.legend = T, 
              size = 0.5, 
              alpha = 0.25, 
              scale = "count",
              draw_quantiles = c(0.25, 0.5, 0.75),
              width = 0.75) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = labels.subset) +
  theme_light() +
  ggtitle("CTLA-4")


ggplot(flow5, aes(x = Subset, y = PD.1, fill = Donor)) +
  geom_violin(trim = FALSE,
              show.legend = T, 
              size = 0.5, 
              alpha = 0.25, 
              scale = "count",
              draw_quantiles = c(0.25, 0.5, 0.75),
              width = 0.75) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = labels.subset) +
  theme_light() +
  ggtitle("PD-1")


ggplot(flow5, aes(x = Subset, y = OX40, fill = Donor)) +
  geom_violin(trim = FALSE,
              show.legend = T, 
              size = 0.5, 
              alpha = 0.25, 
              scale = "count",
              draw_quantiles = c(0.25, 0.5, 0.75),
              width = 0.75) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = labels.subset) +
  theme_light() +
  ggtitle("OX40")



# Boxplot
# ggplot(flow5, aes(x = Subset, y = CD69, fill = Donor)) +
#   geom_boxplot(show.legend = T, 
#               size = 0.5, 
#               alpha = 0.25,
#               width = 0.7) +
#   guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
#   theme_light() +
#   ggtitle("CD69")



#### tSNE ####
flow5.mt <- as.matrix(flow5[,coi])

head(flow5.mt)


# Perplexity: I've used 8, 16, and 30 with no notable differences.  Going with 30 since that's Chamith's default as well.
tsne_out <- Rtsne(flow5.mt, 
                  dims = 2, 
                  perplexity = 30,
                  theta = 0.5, 
                  pca = F, 
                  verbose = T, 
                  check_duplicates = F, 
                  max_iter = 500)

names(tsne_out)
head(tsne_out$Y)

flow5$x.rtsne <- tsne_out$Y[, 1]
flow5$y.rtsne <- tsne_out$Y[, 2]

head(flow5)


# Colorless plot
gg1 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne)) +
  geom_point(color = "black", size = 0.25) +
  theme_light() +
  ggtitle("Raw plot")



# May want to turn all this graphing into some kind of loop/function to minimize code duplication

# State
gg2 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = factor(State, labels = labels.state))) +
  labs(color = "State") +
  geom_point(show.legend = T, size = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values = c("Baseline" = "deepskyblue2", 
                                "Activated" = "deeppink1")) +
  theme_light() +
  ggtitle("Baseline vs. Activated")


# Subset
gg3 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = factor(Subset, labels = labels.subset))) +
  labs(color = "State") +
  geom_point(show.legend = T, size = 0.25, alpha = 0.2) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_manual(values = c("Naive" = "#a6cee3", 
                                "Central Memory" = "#1f78b4",
                                "Effector Memory" = "#b2df8a",
                                "Regulatory" = "#33a02c")) +
  theme_light() +
  ggtitle("Subset")


# Donor
gg4 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = Donor)) +
  geom_point(show.legend = T, size = 0.25, alpha = 0.2) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_brewer(palette = "PRGn") +
  theme_light() +
  ggtitle("Donor")


# CCR6
gg5 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = CCR6)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CCR6")


# ICOS
gg6 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = ICOS)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("ICOS")


# HLA-DR
gg7 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = HLA.DR)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("HLA-DR")


# OX40
gg8 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = OX40)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("OX40")


# PD-1
gg9 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = PD.1)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("PD-1")


# CD40L
gg10 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = CD40L)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD40L")


# CD38
gg11 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = CD38)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD38")


# CD25
gg12 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = CD25)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD25")


# CD69
gg13 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = CD69)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD69")


# CTLA.4
gg14 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = CTLA.4)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CTLA-4")


# CD45RO
gg15 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = CD45RO)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD45RO")


# CD4
gg16 <- ggplot(data = flow5, aes(x = x.rtsne, y = y.rtsne, color = CD4)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD4")


grid.arrange(gg1, gg2, gg3, gg4, ncol = 2)

grid.arrange(gg5, gg6, gg7, gg8, gg9, gg10, gg11, gg12, gg13, gg14, gg15, gg16, ncol=4)



#### heatmap (what's defining these clusters?) ####
# From Kam's tutorial: http://slowkow.com/notes/heatmap-tutorial/

flow5.1 <- flow5[,coi]
nrow(flow5.1)

# Pheatmap can only handle up to 65k rows, and my R crashes before it finishes processing that, so subsampling 20k rows for analysis instead.
flow5.1 <- sample_n(flow5.1, 20000, replace = T)

## Setting coloring scheme based on quantile breaks so coloring will represent equal proportion of data
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0 , 1, length.out = n))
  breaks[!duplicated(breaks)]
}

# Quantile function dislikes data.frames -> coerce into matrix
flowHM_breaks <- quantile_breaks(as.matrix(flow5.1), n = 11)


## Sorting the dendrograms for easier interpretation
# Cluster by column
flowHM_cluster_cols <- hclust(dist(t(flow5.1)))
# Results from unsorted analysis
#plot(flowHM_cluster_cols, main = "Unsorted Dendrogram")

# Sort clustering by column
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

flowHM_cluster_cols <- sort_hclust(flowHM_cluster_cols)
# plot(flowHM_cluster_cols, main = "Sorted Dendrogram")

# Sorting for rows
flowHM_cluster_rows <- sort_hclust(hclust(dist(flow5.1)))


## Plotting sorted heatmap
pheatmap(
  mat = flow5.1,
  color = inferno(10),
  border_color = NA,
  show_colnames = T,
  show_rownames = F,
  drop_levels = T,
  fontsize = 14,
  main = "Flow cytometry analysis with \ndendrogram sorting and quantile color breaks\nSubsample of 20k cells",
  kmeans_k = NA,
  breaks = flowHM_breaks,
  cluster_cols = flowHM_cluster_cols,
  cluster_rows = flowHM_cluster_rows
)


#### k means with k = 8 (k chosen randomly) ####
autoplot(kmeans(flow5[,coi], 
                centers = 8, 
                iter.max = 1000, 
                nstart = 20, 
                algorithm="Lloyd"),
         main = "k means clustering on flow cytometry data post-QC", 
         data = flow5,
         asp = 1)


#### pca ####
## PCA by prcomp (singular value decomposition (SVD), slightly better numerical accuracy than princomp per R help = preferred function)
flow_pcaPr <- prcomp(flow5[,coi], center = F, scale. = F)
summary(flow_pcaPr)

plot(flow_pcaPr, type="lines", main="PCA of flow cytometry data post-QC \nfrom CD4 T cells activated o/n with 1 Dynabead per 1 cell \n(PCs determined by prcomp)")
title(xlab="Principle Component")

#qplot(x=PC1, y=PC2, data=flow_pcaPr) + theme(legend.position="right") + geom_point(size = 1) + scale_colour_manual(values = c("#333333", "#00ccff")) + ggtitle("PCA by prcomp of flow cytometry data") + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL)

## END
dev.off()


#### PDF 2: tSNE by subset ####

# Set seed for reproducibility
set.seed(222)

pdf.2 <- paste0("PDF_", Sys.Date( ), "_immunomics_subsets.pdf")

CairoPDF(pdf.2, 
         width = 11.5, 
         height = 8, 
         onefile = TRUE, 
         family = "Helvetica");


#### naive ####
flow5.naive <- filter(flow5, Subset == "Naive")

flow5.naive.mt <- as.matrix(flow5.naive[,1:12])

head(flow5.naive.mt)

tsne_out.naive <- Rtsne(flow5.naive.mt, 
                        dims = 2, 
                        perplexity = 16,
                        theta = 0.5, 
                        pca = F, 
                        verbose = T, 
                        check_duplicates = F, 
                        max_iter = 1000)

names(tsne_out.naive)
head(tsne_out.naive$Y)

flow5.naive$x.rtsne <- tsne_out.naive$Y[, 1]
flow5.naive$y.rtsne <- tsne_out.naive$Y[, 2]


# Colorless plot
gg17 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne)) +
  geom_point(color = "black", size = 0.25) +
  theme_light() +
  ggtitle("Raw plot")


# State
gg18 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = State)) +
  geom_point(show.legend = T, size = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values = c("Baseline" = "deepskyblue2", 
                                "Activated" = "deeppink1")) +
  theme_light() +
  ggtitle("Baseline vs. Activated")


# Subset
gg19 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = Subset)) +
  geom_point(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_manual(values = c("Naive" = "#a6cee3", 
                                "Central" = "#1f78b4",
                                "Effector" = "#b2df8a",
                                "Regulatory" = "#33a02c")) +
  theme_light() +
  ggtitle("Subset")


# Donor
gg20 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = Donor)) +
  geom_point(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_brewer(palette = "PRGn") +
  theme_light() +
  ggtitle("Donor")


# CCR6
gg21 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = CCR6)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CCR6")


# ICOS
gg22 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = ICOS)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("ICOS")


# HLA-DR
gg23 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = HLA.DR)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("HLA-DR")


# OX40
gg24 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = OX40)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("OX40")


# PD-1
gg25 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = PD.1)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("PD-1")


# CD40L
gg26 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = CD40L)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD40L")


# CD38
gg27 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = CD38)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD38")


# CD25
gg28 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = CD25)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD25")


# CD69
gg29 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = CD69)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD69")


# CTLA.4
gg30 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = CTLA.4)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CTLA-4")


# CD45RO
gg31 <-ggplot(data = flow5.naive, 
              aes(x = x.rtsne, y = y.rtsne, color = CD45RO)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD45RO")


# CD4
gg32 <- ggplot(data = flow5.naive, 
               aes(x = x.rtsne, y = y.rtsne, color = CD4)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD4")


grid.arrange(gg17, gg18, gg19, gg20, ncol = 2) 

grid.arrange(gg21, gg22, gg23, gg24, gg25, gg26, gg27, gg28, gg29, gg30, gg32, ncol=4)


#### Heatmap ####
flow5.N1 <- flow5.naive[,1:12]
nrow(flow5.N1)

# Pheatmap can only handle up to 65k rows, and my R crashes before it finishes processing that, so subsampling 20k rows for analysis instead.
flow5.N1 <- sample_n(flow5.N1, 10000, replace = T)

# Quantile function dislikes data.frames -> coerce into matrix
flowHM_breaks <- quantile_breaks(as.matrix(flow5.N1), n = 11)

## Sorting the dendrograms for easier interpretation
# Cluster by column
flowHM_cluster_cols <- hclust(dist(t(flow5.N1)))
# Results from unsorted analysis
#plot(flowHM_cluster_cols, main = "Unsorted Dendrogram")

# Sort clustering by column
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

flowHM_cluster_cols <- sort_hclust(flowHM_cluster_cols)
# plot(flowHM_cluster_cols, main = "Sorted Dendrogram")

# Sorting for rows
flowHM_cluster_rows <- sort_hclust(hclust(dist(flow5.N1)))


colors.div <- brewer.pal(10, "BrBG")


## Plotting sorted heatmap
pheatmap(
  mat = flow5.N1,
  color = colors.div,
  border_color = NA,
  show_colnames = T,
  show_rownames = F,
  drop_levels = T,
  fontsize = 14,
  main = "Flow cytometry analysis with \ndendrogram sorting and quantile color breaks\nSubsample of 20k cells",
  kmeans_k = NA,
  breaks = flowHM_breaks,
  cluster_cols = flowHM_cluster_cols,
  cluster_rows = flowHM_cluster_rows
)


#### central ####
flow5.central <- filter(flow5, Subset == "Central")

flow5.central.mt <- as.matrix(flow5.central[,1:12])

head(flow5.central.mt)

tsne_out.central <- Rtsne(flow5.central.mt, 
                          dims = 2, 
                          perplexity = 16,
                          theta = 0.5, 
                          pca = F, 
                          verbose = T, 
                          check_duplicates = F, 
                          max_iter = 1000)

names(tsne_out.central)
head(tsne_out.central$Y)

flow5.central$x.rtsne <- tsne_out.central$Y[, 1]
flow5.central$y.rtsne <- tsne_out.central$Y[, 2]


# Colorless plot
gg33 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne)) +
  geom_point(color = "black", size = 0.25) +
  theme_light() +
  ggtitle("Raw plot")


# State
gg34 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = State)) +
  geom_point(show.legend = T, size = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values = c("Baseline" = "deepskyblue2", 
                                "Activated" = "deeppink1")) +
  theme_light() +
  ggtitle("Baseline vs. Activated")


# Subset
gg35 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = Subset)) +
  geom_point(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_manual(values = c("Naive" = "#a6cee3", 
                                "Central" = "#1f78b4",
                                "Effector" = "#b2df8a",
                                "Regulatory" = "#33a02c")) +
  theme_light() +
  ggtitle("Subset")


# Donor
gg36 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = Donor)) +
  geom_point(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_brewer(palette = "PRGn") +
  theme_light() +
  ggtitle("Donor")


# CCR6
gg37 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = CCR6)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CCR6")


# ICOS
gg38 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = ICOS)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("ICOS")


# HLA-DR
gg39 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = HLA.DR)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("HLA-DR")


# OX40
gg40 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = OX40)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("OX40")


# PD-1
gg41 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = PD.1)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("PD-1")


# CD40L
gg42 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = CD40L)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD40L")


# CD38
gg43 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = CD38)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD38")


# CD25
gg44 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = CD25)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD25")


# CD69
gg45 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = CD69)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD69")


# CTLA.4
gg46 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = CTLA.4)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CTLA-4")


# CD45RO
gg47 <-ggplot(data = flow5.central, 
              aes(x = x.rtsne, y = y.rtsne, color = CD45RO)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD45RO")


# CD4
gg48 <- ggplot(data = flow5.central, 
               aes(x = x.rtsne, y = y.rtsne, color = CD4)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD4")


grid.arrange(gg33, gg34, gg35, gg36, ncol = 2)

grid.arrange(gg37, gg38, gg39, gg40, gg41, gg42, gg43, gg44, gg45, gg46, gg48, ncol=4)



#### Heatmap ####
flow5.c1 <- flow5.central[,1:12]
nrow(flow5.c1)

# Pheatmap can only handle up to 65k rows, and my R crashes before it finishes processing that, so subsampling 20k rows for analysis instead.
flow5.c1 <- sample_n(flow5.c1, 10000, replace = T)

# Quantile function dislikes data.frames -> coerce into matrix
flowHM_breaks <- quantile_breaks(as.matrix(flow5.c1), n = 11)

## Sorting the dendrograms for easier interpretation
# Cluster by column
flowHM_cluster_cols <- hclust(dist(t(flow5.c1)))

# Results from unsorted analysis
#plot(flowHM_cluster_cols, main = "Unsorted Dendrogram")

# Sort clustering by column
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

flowHM_cluster_cols <- sort_hclust(flowHM_cluster_cols)
# plot(flowHM_cluster_cols, main = "Sorted Dendrogram")

# Sorting for rows
flowHM_cluster_rows <- sort_hclust(hclust(dist(flow5.c1)))


## Plotting sorted heatmap
pheatmap(
  mat = flow5.c1,
  color = colors.div,
  border_color = NA,
  show_colnames = T,
  show_rownames = F,
  drop_levels = T,
  fontsize = 14,
  main = "Flow cytometry analysis with \ndendrogram sorting and quantile color breaks\nSubsample of 20k cells",
  kmeans_k = NA,
  breaks = flowHM_breaks,
  cluster_cols = flowHM_cluster_cols,
  cluster_rows = flowHM_cluster_rows
)



#### effector ####
flow5.effector <- filter(flow5, Subset == "Effector")

flow5.effector.mt <- as.matrix(flow5.effector[,1:12])

head(flow5.effector.mt)

tsne_out.effector <- Rtsne(flow5.effector.mt, 
                           dims = 2, 
                           perplexity = 16,
                           theta = 0.5, 
                           pca = F, 
                           verbose = T, 
                           check_duplicates = F, 
                           max_iter = 1000)

names(tsne_out.effector)
head(tsne_out.effector$Y)

flow5.effector$x.rtsne <- tsne_out.effector$Y[, 1]
flow5.effector$y.rtsne <- tsne_out.effector$Y[, 2]


# Colorless plot
gg49 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne)) +
  geom_point(color = "black", size = 0.25) +
  theme_light() +
  ggtitle("Raw plot")


# State
gg50 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = State)) +
  geom_point(show.legend = T, size = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values = c("Baseline" = "deepskyblue2", 
                                "Activated" = "deeppink1")) +
  theme_light() +
  ggtitle("Baseline vs. Activated")


# Subset
gg51 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = Subset)) +
  geom_point(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_manual(values = c("Naive" = "#a6cee3", 
                                "Central" = "#1f78b4",
                                "Effector" = "#b2df8a",
                                "Regulatory" = "#33a02c")) +
  theme_light() +
  ggtitle("Subset")


# Donor
gg52 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = Donor)) +
  geom_point(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_brewer(palette = "PRGn") +
  theme_light() +
  ggtitle("Donor")


# CCR6
gg53 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = CCR6)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CCR6")


# ICOS
gg54 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = ICOS)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("ICOS")


# HLA-DR
gg55 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = HLA.DR)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("HLA-DR")


# OX40
gg56 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = OX40)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("OX40")


# PD-1
gg57 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = PD.1)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("PD-1")


# CD40L
gg58 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = CD40L)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD40L")


# CD38
gg59 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = CD38)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD38")


# CD25
gg60 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = CD25)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD25")


# CD69
gg61 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = CD69)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD69")


# CTLA.4
gg62 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = CTLA.4)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CTLA-4")


# CD45RO
gg63 <-ggplot(data = flow5.effector, 
              aes(x = x.rtsne, y = y.rtsne, color = CD45RO)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD45RO")


# CD4
gg64 <- ggplot(data = flow5.effector, 
               aes(x = x.rtsne, y = y.rtsne, color = CD4)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD4")


grid.arrange(gg49, gg50, gg51, gg52, ncol = 2)

grid.arrange(gg53, gg54, gg55, gg56, gg57, gg58, gg59, gg60, gg61, gg62, gg64, ncol=4)



#### Heatmap ####
flow5.e1 <- flow5.effector[,1:12]
nrow(flow5.e1)

# Pheatmap can only handle up to 65k rows, and my R crashes before it finishes processing that, so subsampling 20k rows for analysis instead.
flow5.e1 <- sample_n(flow5.e1, 10000, replace = T)

# Quantile function dislikes data.frames -> coerce into matrix
flowHM_breaks <- quantile_breaks(as.matrix(flow5.e1), n = 11)

## Sorting the dendrograms for easier interpretation
# Cluster by column
flowHM_cluster_cols <- hclust(dist(t(flow5.e1)))
# Results from unsorted analysis
#plot(flowHM_cluster_cols, main = "Unsorted Dendrogram")

# Sort clustering by column
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

flowHM_cluster_cols <- sort_hclust(flowHM_cluster_cols)
# plot(flowHM_cluster_cols, main = "Sorted Dendrogram")

# Sorting for rows
flowHM_cluster_rows <- sort_hclust(hclust(dist(flow5.e1)))


# Color palette
colors.div <- brewer.pal(10, "BrBG")

## Plotting sorted heatmap
pheatmap(
  mat = flow5.e1,
  color = colors.div,
  border_color = NA,
  show_colnames = T,
  show_rownames = F,
  drop_levels = T,
  fontsize = 14,
  main = "Flow cytometry analysis with \ndendrogram sorting and quantile color breaks\nSubsample of 20k cells",
  kmeans_k = NA,
  breaks = flowHM_breaks,
  cluster_cols = flowHM_cluster_cols,
  cluster_rows = flowHM_cluster_rows
)




#### regulatory ####
flow5.reg <- filter(flow5, Subset == "Regulatory")

flow5.reg.mt <- as.matrix(flow5.reg[,1:12])

head(flow5.reg.mt)

tsne_out.reg <- Rtsne(flow5.reg.mt, 
                      dims = 2, 
                      perplexity = 16,
                      theta = 0.5, 
                      pca = F, 
                      verbose = T, 
                      check_duplicates = F, 
                      max_iter = 1000)

names(tsne_out.reg)
head(tsne_out.reg$Y)

flow5.reg$x.rtsne <- tsne_out.reg$Y[, 1]
flow5.reg$y.rtsne <- tsne_out.reg$Y[, 2]


# Colorless plot
gg65 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne)) +
  geom_point(color = "black", size = 0.25) +
  theme_light() +
  ggtitle("Raw plot")


# State
gg66 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = State)) +
  geom_point(show.legend = T, size = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values = c("Baseline" = "deepskyblue2", 
                                "Activated" = "deeppink1")) +
  theme_light() +
  ggtitle("Baseline vs. Activated")


# Subset
gg67 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = Subset)) +
  geom_point(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_manual(values = c("Naive" = "#a6cee3", 
                                "Central" = "#1f78b4",
                                "Effector" = "#b2df8a",
                                "Regulatory" = "#33a02c")) +
  theme_light() +
  ggtitle("Subset")


# Donor
gg68 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = Donor)) +
  geom_point(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_brewer(palette = "PRGn") +
  theme_light() +
  ggtitle("Donor")


# CCR6
gg69 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = CCR6)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CCR6")


# ICOS
gg70 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = ICOS)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("ICOS")


# HLA-DR
gg71 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = HLA.DR)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("HLA-DR")


# OX40
gg72 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = OX40)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("OX40")


# PD-1
gg73 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = PD.1)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("PD-1")


# CD40L
gg74 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = CD40L)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD40L")


# CD38
gg75 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = CD38)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD38")


# CD25
gg76 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = CD25)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD25")


# CD69
gg77 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = CD69)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD69")


# CTLA.4
gg78 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = CTLA.4)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CTLA-4")


# CD45RO
gg79 <-ggplot(data = flow5.reg, 
              aes(x = x.rtsne, y = y.rtsne, color = CD45RO)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD45RO")


# CD4
gg80 <- ggplot(data = flow5.reg, 
               aes(x = x.rtsne, y = y.rtsne, color = CD4)) +
  geom_point(show.legend = F, size = 0.25) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle("CD4")


grid.arrange(gg65, gg66, gg67, gg68, ncol = 2) 

grid.arrange(gg69, gg70, gg71, gg72, gg73, gg74, gg75, gg76, gg77, gg78, gg80, ncol=4)



#### Heatmap ####
flow5.R1 <- flow5.reg[,1:12]
nrow(flow5.R1)

# Pheatmap can only handle up to 65k rows, and my R crashes before it finishes processing that, so subsampling 20k rows for analysis instead.
flow5.R1 <- sample_n(flow5.R1, 10000, replace = T)

# Quantile function dislikes data.frames -> coerce into matrix
flowHM_breaks <- quantile_breaks(as.matrix(flow5.R1), n = 11)

## Sorting the dendrograms for easier interpretation
# Cluster by column
flowHM_cluster_cols <- hclust(dist(t(flow5.R1)))
# Results from unsorted analysis
#plot(flowHM_cluster_cols, main = "Unsorted Dendrogram")

# Sort clustering by column
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

flowHM_cluster_cols <- sort_hclust(flowHM_cluster_cols)
# plot(flowHM_cluster_cols, main = "Sorted Dendrogram")

# Sorting for rows
flowHM_cluster_rows <- sort_hclust(hclust(dist(flow5.R1)))


# Color palette
colors.div <- brewer.pal(10, "PRGn")

## Plotting sorted heatmap
pheatmap(
  mat = flow5.R1,
  color = colors.div,
  border_color = NA,
  show_colnames = T,
  show_rownames = F,
  drop_levels = T,
  fontsize = 14,
  main = "Flow cytometry analysis with \ndendrogram sorting and quantile color breaks\nSubsample of 20k cells",
  kmeans_k = NA,
  breaks = flowHM_breaks,
  cluster_cols = flowHM_cluster_cols,
  cluster_rows = flowHM_cluster_rows
)


dev.off()
