# Christin M. Hong
# 2017-07
# PI: Soumya Raychaudhuri, Harvard Medical School
# Soumya's homework: Analyzing flow cytometry data manually in R

#### Project Notes ####

# This data was exported from FlowJo.  Samples have already been compensated and gated for live single lymphocytes.  Using scale value data from all compensated parameters.  

#### Infrastructure ####

## Import libraries
# Analysis
library(pheatmap)
library(dendsort)
library(Rtsne)
#library(broom) # Not sure this is necessary, will check later

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


# Color palette for heatmap
colors.div <- brewer.pal(10, "PRGn")


#### Automate import, annotatation, and compensation adjustment of data ####
# These are QC processes to be performed by file

# Set working directory
wd <- "~/00_SR-lab-data/wetLab/flowCytometry/flowCytometry2017/20170905_TNF-5d-5d/"
setwd(wd)

tp <- "5d-prime_5d-activation"


# Exporting data from FlowJo with the following naming format: "export_## Baseline/Activated [Subset] D#_Live_[date]"

fnums <- 1:72

fnums.1 <- sprintf("%02d", fnums) # Adding leading zeros for 2-digit numbers
class(fnums.1)

flow.list <- list() # Create a list in which you intend to save your dfs


#### START LOOP ####

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
    
    #### Subset each sample to ensure that each sample has equal weight in data ####
    flow.1 <- sample_n(flow, 2000, replace = T)
    
    # Annotate
    flow.1$PrimeTP <- "3d"
    flow.1$ActivationTP <- "Overnight"
    flow.1$State <- fname.3[[1]][2]
    flow.1$Dose <- fname.3[[1]][3]
    flow.1$Donor <- fname.3[[1]][4]
    flow.1$Replicate <- fname.3[[1]][5]
    #  print(head(flow.1))
    
    ## Extract columns of interest
    names(flow.1)
    flow.2 <- select(flow.1, 7:8, 10:16, 18:20, 22:27)
    names(flow.2)
  
  
  # #### Adjusting for compensation - shifting 99.9 percentile of data to > 0 ####
  floor <- apply(flow.2[,1:12], 2, quantile, probs = c(0.01), na.rm = F, type = 7) # Getting lower bound of data (see notes below)

  floor.df <- as.data.frame(floor)
  floor.df$names <- rownames(floor.df) # Turn rownames into a column
  floor.dfNeg <- filter(floor.df, floor < 0) # Extract negative values
  floor.dfNeg$shift <- abs(floor.dfNeg$floor) # Convert neg to positive
  floor.dfNeg$floor <- NULL # Drop original column
  floor.dfNeg.t <- t(floor.dfNeg) # transpose
  colnames(floor.dfNeg.t) <- as.character(unlist(floor.dfNeg.t[1,])) # turn first row into column name
  floor.dfNeg.t <- floor.dfNeg.t[-1, ] # drop first row
  floor.t.df <- data.frame(lapply(floor.dfNeg.t, type.convert), stringsAsFactors = T) # Convert to data.frame


  comp.match <- intersect(colnames(floor.t.df), colnames(flow.2)) # Find matching columns (note that colnames function is same as names function)

  #### Add to shift 99.9 percentile above zero to adjust for compensation ####
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


#### Log scale data since fluorescence is on log scale ####
## Note: log in R is natural log, whereas FlowJo plots in log10, but I doubt that matters

flow1 <- flow.all

# Convert negative values to zero
flow1[flow1 < 0] <- 0

# Scale fluorescence values while retaining data annotations
flow1.1 <- flow1
flow1.1[, 1:12] <- log10(flow1[, 1:12] + 1)



#### QC by donor ####
# Normalizing by donor to capture the full expression range of a marker within each donor

# This is Chamith's preferred method - to stay as close to the original data as possible and maximize information.  Chamith recommends avoiding normalizing now because it'll reduce future ability to cluster.  He also suggests looking hard at FlowSum, which he finds to be faster than DenseVM and just as accurate (it was also recently published, which is why they haven't been using it historically).

flow2 <- flow1.1


#### Loop per donor ####
donors <- levels(factor(flow2$Donor))

donor.list <- list() # Create a list in which you intend to save your dfs


for (id in donors) {
  print(id)
  flow2.D <- flow2[ flow2[["Donor"]] == id , ]
  
  
  # Calculate minimum threshold per column.  (The 2 argument applies the function to columns, while 1 would apply to rows)
  min <- apply(flow2.D[, 1:12], 2, quantile, probs = c(0.02), na.rm = F, type = 7)
  
  ## Subtract minimum per column from each value in that column
  flow2.1 <- flow2.D
  flow2.1[, 1:12] <- sweep(flow2.D[, 1:12], 2, min, FUN = "-") 
  # So that's the numerator!
  
  
  # Find the value for 99.9% of the max value per column
  max <- apply(flow2.D[, 1:12], 2, quantile, probs = c(0.999), na.rm = F, type = 7)
  
  # Subtract min from 99th percentile to get denominator per column
  adjMax <- max - min
  
  
  # Final min/max normalized matrix
  flow2.2 <- flow2.1
  flow2.2[,1:12] <- sweep(flow2.1[,1:12], 2, adjMax, FUN = "/")
  
  
  #### Filtering/excluding outliers ####
  flow3 <- flow2.2
  print(nrow(flow3))
  
  # Removing samples with value < 0
  flow3.1 <- flow3[apply(flow3[, 1:12], 1, function(x) all(x >= 0)), ]
  print(nrow(flow3.1))
  
  # Removing samples with a value > 1
  flow3.2 <- flow3.1[apply(flow3.1[, 1:12], 1, function(x) all(x <= 1)), ]
  print(nrow(flow3.2))
  
  # Requiring 3 values to be >0.25 to ensure sample was well-stained
  flow3.3 <- flow3.2
  
  # Count number of values >0.25 per row
  flow3.3$highStain <- rowSums(flow3.3[, 1:12] > 0.25)
  
  # Keep only rows with >3 values that are >0.25
  flow3.4 <- filter(flow3.3, highStain > 3)
  
  print(nrow(flow3.4))
  
  
  
  #### Return to original values ####
  # Multiply values by adjMax
  flow4 <- flow3.4
  flow4[,1:12] <- sweep(flow3.4[,1:12], 2, adjMax, FUN = "*")
  
  # Add minimum back to column
  flow4.1 <- flow4
  flow4.1[, 1:12] <- sweep(flow4[, 1:12], 2, min, FUN = "+") 
  
  # Column center by subtracting the mean value for each column (improves heatmap and pca) 
  flow4.2 <- flow4.1
  
  flow4.mean <- apply(flow4.1[,1:12], 2, mean)
  
  flow4.2[,1:12] <- sweep(flow4.1[,1:12], 2, flow4.mean, FUN = "-")
  
  
  print(head(flow4.2))
  print(tail(flow4.2)) # Check that annotations are still present!
  
  donor.list[[id]] <- flow4.2 # Add updated df to list of dfs
}



#### Concatenate list of data frames into one ####
flow5 <- rbind.fill(donor.list)

hist(flow5$CD69)


#### Write out to an external CSV for later ease of loading/appending ####
output <- paste0("flowCytometryQC_", Sys.Date( ), tp, ".csv")

write.table(x = data.frame(flow5, stringsAsFactors = T),
            file = output,
            quote = F,
            sep = ",",
            row.names = F)


#### Data analysis ####
hist(flow5$Foxp3)
hist(flow5$CD69)
hist(flow5$CD45RA)
hist(flow5$Ki67)
hist(flow5$CD27)
hist(flow5$HLA.DR)
hist(flow5$IL.2)
hist(flow5$PD.1)
hist(flow5$IFNg)
hist(flow5$CD25)
hist(flow5$CTLA4)
hist(flow5$IL.17A)


#### Plotting values by boxplots ####
names(flow5)

flow5$State <- ordered(flow5$State, levels = c("Base", "Act"))

summarise(flow5 %>% filter(State %in% "Base", 
                               Donor %in% "LP19", 
                               Dose %in% "T0"), 
              median(Foxp3))

test <- flow5 %>% filter(State %in% "Base", 
                         Donor %in% "LP19", 
                         Dose %in% "T0")

levels(factor(test$Replicate))

# Next: Test making a new column that appends Dose to State (State-Dose) for plotting on x-axis.  
# May also need to reshape data such that the markers = rows for plotting facets.

#

d19 <- filter(flow5, Donor == "LP19")

ggplot(data = subset(d19, State %in% c("Act")), aes(x = Replicate, y = IL.2, color = Dose)) +
  geom_boxplot(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_manual(values = c("T0" = "#a6cee3", 
                                "T1" = "#1f78b4",
                                "T10" = "#b2df8a",
                                "T100" = "#33a02c")) +
  theme_light() +
  ggtitle("Dose")


ggplot(data = subset(flow5, State %in% c("Act") & Donor %in% c("LP19")), 
       aes(x = Dose, y = IFNg, color = Replicate)) +
  geom_boxplot(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_light() +
#  geom_boxplot(data = subset(flow5, State %in% c("Act") & Donor %in% c("LP20")), 
#               show.legend = T, size = 0.25, alpha = 0.25) +
  ggtitle("")


ggplot(data = subset(flow5, State %in% c("Act") & Donor %in% c("LP22")), 
       aes(x = Dose, y = IFNg, color = Replicate)) +
  geom_boxplot(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_light() +
  ggtitle("")


ggplot(data = subset(flow5, State %in% c("Base")), aes(x = Dose, y = CTLA4, color = Donor)) +
  geom_boxplot(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_light() +
  ggtitle("Donor")


gg2 <- ggplot(data = subset(flow5, State %in% c("Act")), aes(x = Dose, y = Ki67, color = Donor)) +
  geom_boxplot(show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_light() +
  ggtitle("Donor")



#### CODING k means clustering ####
# Subset data for faster runs and dev
flow5.sub <- sample_n(flow5, 20000, replace = F)

# PCA on all markers, using PC1 and PC2 for kmeans
flow_pcaPr <- prcomp(flow5.sub[, 1:12], center = F, scale. = F)


# Add PC1 and PC2 to data
kdata$PC1 <- flow_pcaPr$x[,"PC1"]
kdata$PC2 <- flow_pcaPr$x[,"PC2"]

plot(kdata$PC1, kdata$PC2)


# Set variables
k <- 7
kdata <- flow5.sub
iter <- 200 # max # of iterations
k.conv <- 1e-6 # threshold for convergence


## Initialize centroids
seeds <- sample_n(kdata, k, replace = F)


# Assign each cell (row) to the closest centroid by squared Euclidean distance [Euclidean distance = sqrt((x1 - x2)^2 + (y1 - y2)^2)]
for(i in 1:nrow(kdata)) {
  kdata[i, c("oldCentroid")] <- which.min((seeds$PC1 - kdata[i, c("PC1")])^2 + (seeds$PC2 - kdata[i, c("PC2")])^2)
} 


# Create new centroids for each cluster
k.list <- list() # Create a list for saving dfs

for(j in 1:k) {
  k.j <- kdata %>% 
    filter(oldCentroid %in% j) %>% 
    summarize(meanPC1 = mean(PC1), meanPC2 = mean(PC2))

    k.list[[j]] <- k.j # Add updated centroid to list of updated centroids
}

# Concatenate list of updated centroid dfs into one df
seeds.new <- rbind.fill(k.list)


# Re-assign data to updated centroids
for(i in 1:nrow(kdata)) {
  kdata[i, c("newCentroid")] <- 
    which.min((seeds.new$meanPC1 - kdata[i, c("PC1")])^2 + (seeds.new$meanPC2 - kdata[i, c("PC2")])^2)
} 


#### Start iteration ####
seed.history <- list() 

for(x in 1:iter) {
  # Update centroids
  k.list <- list() # Create list for storing new centroids

  for(a in 1:k) {
    k.a <- kdata %>%
      filter(newCentroid %in% a) %>%
      summarize(meanPC1 = mean(PC1), meanPC2 = mean(PC2))

    k.list[[a]] <- k.a # Collect updated centroids
  }

  seeds.new <- rbind.fill(k.list) # Create 1 df of updated centroids
  seed.history[[x]] <- seeds.new # Store centroid history

  # Re-assign data
  for(i in 1:nrow(kdata)) {
    kdata[i, c("newCentroid")] <- 
      which.min((seeds.new$meanPC1 - kdata[i, c("PC1")])^2 + (seeds.new$meanPC2 - kdata[i, c("PC2")])^2)
  } 

  print(paste0(x, " iteration(s) and ", nrow(kdata %>% filter(newCentroid %in% "1")), " rows in cluster 1."))


  # Check for centroid convergence
  if(x > 1) {
    if(sum((seed.history[[x]][c("meanPC1")] - seed.history[[x-1]][c("meanPC1")])^2) < k.conv & 
      sum((seed.history[[x]][c("meanPC2")] - seed.history[[x-1]][c("meanPC2")])^2) < k.conv) {
      print(paste0("Success!  Sum of differences between current and previous centroids is < ", k.conv))
      break
    }
  }
}


# Check and save in PDF
pdf.kmeans <- paste0("PDF_", Sys.Date( ), "_kmeans.pdf")

CairoPDF(pdf.kmeans, 
         width = 11.5, 
         height = 8, 
         onefile = TRUE, 
         family = "Helvetica");


gg1 <- ggplot() +  
  geom_point(data = kdata, aes(x = PC1, 
                               y = PC2, 
                               color = factor(newCentroid)), 
             show.legend = T, size = 0.25, alpha = 0.25) +
  geom_point(data = seed.history[[1]], aes(x = meanPC1, y = meanPC2),
             show.legend = F, size = 2, alpha = 1) +
    guides(colour = guide_legend(title = "Clusters", 
                               override.aes = list(size = 2, alpha = 1))) +
  theme_light() +
  ggtitle("Centroids: Iteration 1")


gg2 <- ggplot() +  
  geom_point(data = kdata, aes(x = PC1, 
                               y = PC2, 
                               color = factor(newCentroid)), 
             show.legend = T, size = 0.25, alpha = 0.25) +
  geom_point(data = seed.history[[10]], aes(x = meanPC1, y = meanPC2),
             show.legend = F, size = 2, alpha = 1) +
  guides(colour = guide_legend(title = "Clusters", 
                               override.aes = list(size = 2, alpha = 1))) +
  theme_light() +
  ggtitle("Centroids: Iteration 10")


gg3 <- ggplot() +  
  geom_point(data = kdata, aes(x = PC1, 
                               y = PC2, 
                               color = factor(newCentroid)), 
             show.legend = T, size = 0.25, alpha = 0.25) +
  geom_point(data = seed.history[[x]], aes(x = meanPC1, y = meanPC2),
             show.legend = F, size = 2, alpha = 1) +
  guides(colour = guide_legend(title = "Clusters", 
                               override.aes = list(size = 2, alpha = 1))) +
  theme_light() +
  ggtitle("Centroids: Iteration Final")


gg4 <- ggplot(data = kdata, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = factor(newCentroid)), show.legend = T, size = 0.25, alpha = 0.25) +
  guides(colour = guide_legend(title = "Clusters", 
                               override.aes = list(size = 2, alpha = 1))) +
  theme_light() +
  ggtitle("")

grid.arrange(gg1, gg2, gg3, gg4, ncol = 2)

autoplot(kmeans(kdata[, 1:12], 
                centers = k, 
                iter.max = iter, 
                nstart = 20, 
                algorithm="Lloyd"),
         main = "k means clustering on flow cytometry data post-QC\nPCA by prcomp and kmeans by R", 
         data = kdata,
         size = 0.5,
         alpha = 0.5,
         asp = 1)


dev.off();


# More elegant, R-friendly method: https://stackoverflow.com/questions/27082378/how-to-compute-distances-between-centroids-and-data-matrix-for-kmeans-algorithm

