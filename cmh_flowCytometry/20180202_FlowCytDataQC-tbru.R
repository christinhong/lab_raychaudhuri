#### Background ####

# Christin M. Hong
# Last updated: 2018-02-05
# PI: Soumya Raychaudhuri, Harvard Medical School

# This data was exported as CSVs from FlowJo.
# Baseline is gated on Lymphocytes (magnetic gate) -> Singlets (magnetic) -> Time (2 - 180 seconds) vs. Live (-1000 - 1000) -> CD4 (magnetic gate).

# Activated is gated on Lymphocytes (magnetic) -> Singlets (magnetic) -> Time (2 - 180 seconds) vs. Live (-1000 - 1000) -> CD4 (magnetic gate).
# (Currently excluding blasting cells due to irregular frequency of blasts.)

# Currently clustering on uncompensated scale values and coloring on compensated channel values, as coloring is just for visualization, not statistics.  See more info on channel vs. scale in the Notes section of this script.

# This script is still in progress - expect inefficiencies.


#### Infrastructure ####
#### Import libraries ####
# Analysis
library(Biobase)
library(flowCore)
library(FlowSOM)
library(pheatmap)
library(dendsort)
library(Rtsne)
library(cytofkit)
library(sva)
library(pamr)
library(limma)


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


#### Set seed for reproducibility ####
set.seed(200)


#### Set variables ####
# Set working directory
wd <-
  "~/00_SR-lab-data/data_wet/flowCytometry/flowCytometry2018/20180130_batchCtrlChecks/pilot-tbru-batches/"
setwd(wd)


#### * File name * ####
out.name <-
  paste0(Sys.Date(),
         "_fcQC_scaleUncomp_keep001-999_compAdjusted_transform-log10_colorChanComp_batchByDate-withinStates")

# Columns with data of interest
coi <- c(8:11, 13:19, 22)


# Min and max thresholds for outlier removal
# Note: PCA and SSW/SST suggest that top and bottom outlier removal matters.  Flow cytometry outliers seem independent of each other, as if there's technical variation in the cytometer measurements rather than issues with staining.  To minimize data loss while maximizing QC, best to set minimal thresholds.
outMin = 0.001 # Default: 0.001.  So for 10 markers, will drop ~1% of data.
outMax = 0.999 # Default: 0.999


# Number of cells/sample for analysis
cells.final <- 500



#### START ####

#### Automate import, annotatation, and compensation adjustment of data ####
# These are QC processes to be performed by file
# Exporting data from FlowJo with the following naming format: "export_## Baseline/Activated [Subset] D#_Live_[date]"

fnums <- 1:32

fnums.1 <-
  sprintf("%02d", fnums) # Adding leading zeros for 2-digit numbers

class(fnums.1)

flow.list <-
  list() # Create a list in which you intend to save your dfs


# Annotate and QC while saving printed status messages
capture.output(for (i in fnums.1) {
  name.pattern <- paste0("^export_scaleUncomp_", i)
  fname <- list.files(pattern = name.pattern) # Will return list of names
  length(fname)
  
  # Since the filenames for each batch differ by date...
  for (a in fname) {
    tryCatch( {
      print(a) # Getting filename being processed
      fname.1 <- strsplit(a, "_")
      fname.2 <- fname.1[[1]][3]
      fname.3 <- strsplit(fname.2, " ") # Getting data from filename
      fdate <- fname.1[[1]][5]
    
      # Import scale uncompensated files
      flow <- read.csv(a, header = TRUE)
    
      # Annotate
      flow$State <- fname.3[[1]][2]
      flow$Subset <- fname.3[[1]][3]
      flow$Donor <- fname.3[[1]][4]
      flow$Batch <- fdate
      flow$Cell <- rownames(flow) # Turn rownames into a column
    #  print(tail(flow))
    
  
    
      #### Import and annotate matching compensated file ####
      fname2 <- gsub("scaleUncomp", "chanComp", a, fixed = TRUE)

      # Import matching compensated file
      flow.comp <- read.csv(fname2, header = TRUE)
      
      # Annotate
      flow.comp$Cell <- rownames(flow.comp) # Assuming that FlowJo outputs cells in same order no matter export format
      

      #### Merge on cell number ####
      flow1 <-
        left_join(flow,
                  flow.comp,
                  by = c("Cell" = "Cell"),
                  copy = FALSE,
                  c("_scUncomp", "_chComp")) # Annotating with _ instead of . for ease of string-splitting down the line.
    
    
      #### Remove outliers ####
      flow2 <- flow1
    
      # Calculate minimum threshold for outliers.  (The 2 argument applies the function to columns, while 1 would apply to rows)
      min <-
        apply(
          flow2[, coi],
          2,
          quantile,
          probs = c(outMin),
          na.rm = F,
          type = 7
        )
    
      # Calculate max threshold for outliers
      max <-
        apply(
          flow2[, coi],
          2,
          quantile,
          probs = c(outMax),
          na.rm = F,
          type = 7
        )
    
      # Print status message
      print(paste0("Number of cells pre-filtering: ", nrow(flow2)))
    
  
      # Keep only values that fall within min and max
      flow3 <-
        flow2[apply(flow2[, coi], 1, function(x)
          all(x > min & x < max)),]
    
      # Print status message
      print(paste0("Number of cells post-outlier removal: ", nrow(flow3)))
    
      
    
      #### Adjust for negatives by adding back min(x) ####
      # Get lowest remaining value in columns of interest
      floor <- apply(flow3[, coi], 2, min)
    
      # Convert floor data into a workable format.  For consistency, will add floor to all markers, not just markers with negative min values.
      floor.df <- as.data.frame(floor)
      floor.df$floor <- abs(floor.df$floor) # Convert neg to positive
      floor.t.mat <- t(floor.df) # transpose
      floor.t.df <-
        as.data.frame(floor.t.mat, stringsAsFactors = T) # Convert from matrix to data.frame
    
      # Select matching columns (note that colnames function is same as names function)
      comp.match <- intersect(colnames(floor.t.df), colnames(flow3))
    
  
      # Add absolute value of xmin
      flow4 <- flow3
    
      # Start Kam
      for (col in comp.match) {
        flow4[[col]] <- flow3[[col]] + floor.t.df[[col]]
      }
      # End Kam
    
    
      #### Subset each sample so that each sample has equal weight in analysis ####
      flow5 <- flow4
      flow5 <- sample_n(flow4, cells.final, replace = F)
    
      flow.list[[a]] <- flow5 # Add updated df to list of dfs
      
      
      }   # closing bracket for tryCatch
    , error = function(e) {
      cat("ERROR :", conditionMessage(e), "\n")
      }) # closing parenthesis for tryCatch (if there's an error in processing, e.g. a missing sample or not enough cells, will skip and continue to next line in list instead of aborting the loop completely)
    
    } # closing bracket for fname list loop
  } # closing bracket for full loop
  , file = paste0("fcQC-by-sample_", out.name, ".txt")) # closing parenthesis for capture.output, which saves output to TXT file



#### Concatenate list of data frames into one ####
flow.all <- rbind.fill(flow.list)


# Check that factors are correctly assigned and that all levels are present
factor(unique(flow.all$Donor))

names(flow.all)
head(flow.all)
tail(flow.all)

flow7 <- flow.all


#### Log10 transform ####
# Convert any negative values to zero
flow7[flow7 < 0] <- 0

# Log10 transform
flow8 <- flow7
flow8[, coi] <- log10(flow7[, coi] + 1)
# If not converting to Z scores, can use log10 for easier interpretation, as value = # of zeros/decimals (e.g. log10(100) = 2)


# #### Inverse hyperbolic sin transform ####
# flow8 <- flow7
# flow8[, coi] <- asinh(flow7[, coi])



#### Converting annotations to numeric values so they can be preserved in data matrices ####
flow8$State <- revalue(flow8$State,
                       c("Baseline" = 0,
                         "Activated" = 1))

levels(factor(flow8$State))
labels.state = c("Baseline", "Activated")


flow8$Subset <- revalue(flow8$Subset,
                       c("Naive"=1,
                         "Central"=2,
                         "Effector"=3))

levels(factor(flow8$Subset))
labels.subset = c("Naive", "Central Memory", "Effector Memory")

flow8$Donor <- gsub("D", "", flow8$Donor)

head(flow8)



# #### OUTDATED - Standardize entire dataset ####
# Using entire dataset to standardize to capture the full possible range of expression (not sure this is ideal - may want to standardize by batch)
# 
# # Column center by subtracting the mean value for each column (improves heatmap and pca)
# flow8.mean <- apply(flow8[, coi], 2, mean)
# flow8.1 <- flow8
# 
# flow8.2 <- flow8.1
# flow8.2[, coi] <- sweep(flow8.1[, coi], 2, flow8.mean, FUN = "-")
# 
# # Divide by standard deviation (normalizes the distributions for each marker, convert to Z-score)
# flow8.3 <- flow8.2
# flow8.sd <- apply(flow8.2[, coi], 2, sd)
# flow8.3[, coi] <- sweep(flow8.3[, coi], 2, flow8.sd, FUN = "/")
# 
# flow9 <- flow8.3
# 
# hist(flow9$PD.1_scUncomp)


#### Adjust for batches using ComBat {sva} ####
# See https://www.rdocumentation.org/packages/sva/versions/3.20.0/topics/ComBat
# ComBat assumes data has been normalized (e.g. for RNA-seq, converted to FPKM), then does its own standardization of the normalized data (so no need to convert to Z-scores in advance).
flow10 <- flow8
head(flow10)


#### ComBat on Baseline samples ####
## For this project, State is 100% correlated with batch -> will adjust for batch within State
flow10.base <- filter(flow10, State == 0)
nrow(flow10.base)

# Split data into pheno vs. expression
names(flow10.base)

# Phenotype info dataframe: Samples are rows.  For each sample, has sample number, outcome, batch, and cancer status. 
flow10.b.pheno <- flow10.base[-coi]
head(flow10.b.pheno)

# Expression matrix: Rows are genes/markers, columns are samples.
flow10.b.e <- t(as.matrix(flow10.base[coi]))
head(flow10.b.e[, 1:3])
class(flow10.b.e)


# Create full model matrix with the variable of interest - but no variable of interest for this (though later may treat state like case/control)
# mod.flow10 <- model.matrix(~as.factor(State), data = flow10.pheno)

# Create null matrix (not adjusting for anything else, as batch is seen across all subsets and don't ever want to regress out donor since this'll eventually be linked back to genotype).
mod0.b.flow10 <- model.matrix(~1, data = flow10.b.pheno)

# Get batch column (date run)
flow10.b.batch <- flow10.b.pheno$Batch

# For ComBat, model may include covariate of interest - see https://support.bioconductor.org/p/63082/.  But that isn't in the tutorial, so will use null for now.

# Apply ComBat to standardize data across genes
flow10.b.combat_edata <- ComBat(dat = flow10.b.e, 
                       batch = flow10.b.batch,
                       mod = mod0.b.flow10,
                       par.prior = TRUE,
                       prior.plots = TRUE)

#### UNDER CONSIDERATION ####
# Running with par.prior = TRUE to see prior.plots/have some idea of how ComBat it working, but some markers having highly skewed distributions.  May be better to run with par.prior = FALSE.
# Note: QC plots look a little strange - will probably want to read more to interpret.  See https://support.bioconductor.org/p/71198/


#### ComBat Options ####
##  "par.prior = FALSE": non-parametric Bayes adjustments - will take longer to run

## "prior.plots = TRUE": produces prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric estimate.

## "mean.only = TRUE": only adjusts the mean of the batch effects across batches (default adjusts the mean and variance).  This option is recommended for cases where milder batch effects are expected (so no need to adjust the variance), ***OR in cases where the variances are expected to be different across batches due to the biology.***

## "ref.batch": allows users to select one batch as a reference to which other batches will be adjusted.  Specifically, the means and variances of the non-reference batches will be adjusted to make the mean/variance of the reference batch. This is a useful feature for cases where one batch is larger or better quality. In addition, this will be useful in biomarker situations where the researcher wants to fix the traning set/model and then adjust test sets to the reference/training batch. This avoids test-set bias in such studies.


# Perform significance analysis between model matrix and null matrix with the ComBat adjusted values (very slow on my flow cytometry data and, I think, not relevant right now -> skipping)
# pValuesComBat = f.pvalue(flow10.combat_edata, mod.flow10, mod0.flow10)
# qValuesComBat = p.adjust(pValuesComBat, method="BH")


# Re-format data
head(flow10.b.combat_edata)

flow10.b.e.combat <- as.data.frame(t(flow10.b.combat_edata))
head(flow10.b.e.combat)
class(flow10.b.e.combat)

# Replace previous values with date/batch adjusted values
flow11 <- flow8
flow11.base <- filter(flow11, State == 0)

flow11.base[, names(flow10.b.e.combat)] <- flow10.b.e.combat

names(flow8[coi])
names(flow11.base[coi]) # Cool, looks like columns of interest are still where they should be

head(flow10.base$CCR6_scUncomp)
head(flow10.b.e.combat$CCR6_scUncomp)
head(flow11.base$CCR6_scUncomp)
# Data looks properly replaced in flow10, nice.  Now for activated...


#### ComBat on Activated samples ####
flow10.act <- filter(flow10, State == 1)
names(flow10.act)

# Phenotype info dataframe: Samples are rows.  For each sample, has sample number, outcome, batch, and cancer status. 
flow10.a.pheno <- flow10.act[-coi]
head(flow10.a.pheno)

# Expression matrix: Rows are genes/markers, columns are samples.
flow10.a.e <- as.matrix(t(flow10.act[coi]))
head(flow10.a.e)
class(flow10.a.e)

# Create full model matrix with the variable of interest - but no variable of interest for this (though later may treat state like case/control)
# mod.flow10 <- model.matrix(~as.factor(State), data = flow10.pheno)

# Create null matrix (not adjusting for anything else, as batch is seen across all subsets and don't ever want to regress out donor since this'll eventually be linked back to genotype).
mod0.a.flow10 <- model.matrix(~1, data = flow10.a.pheno)

# Get batch column (date run)
flow10.a.batch <- flow10.a.pheno$Batch

# For ComBat, model may include covariate of interest - see: https://support.bioconductor.org/p/63082/.  But that isn't in the tutorial, so will use null for now.

# Apply ComBat to standardize data across genes
flow10.a.combat_edata <- ComBat(dat = flow10.a.e, 
                                batch = flow10.a.batch,
                                mod = mod0.a.flow10,
                                par.prior = TRUE,
                                prior.plots = TRUE)


# Re-format data
head(flow10.a.combat_edata)

flow10.a.e.combat <- as.data.frame(t(flow10.a.combat_edata))
head(flow10.a.e.combat)
class(flow10.a.e.combat)

# Replace previous values with date/batch adjusted values
flow11 <- flow8
flow11.act <- filter(flow11, State == 1)

flow11.act[, names(flow10.a.e.combat)] <- flow10.a.e.combat

names(flow8[coi])
names(flow11.act[coi]) # Cool, looks like columns of interest are still where they should be

head(flow10.act$CCR6_scUncomp)
head(flow10.a.e.combat$CCR6_scUncomp)
head(flow11.act$CCR6_scUncomp)
# Data looks properly replaced in flow10, nice.  Now for activated...


#### Rejoin baseline and activated data frames ####
flow11 <- rbind(flow11.base, flow11.act)

names(flow8[coi])
names(flow11[coi]) # Still looks good, nice.



#### Save data frame as CSV for ease of import later ####
write.table(
  x = data.frame(flow11, stringsAsFactors = T),
  file = paste0(out.name, ".csv"),
  quote = F,
  sep = ",",
  row.names = F
)


#### Start SSW/SST analysis ####
df.t <- flow11

#### Calculate SST ####
# Subtract mean and square each value
mean.df.t <- apply(df.t[, coi], 2, mean)

df.s.t <- df.t
df.s.t[, coi] <- sweep(df.t[, coi], 2, mean.df.t, FUN = "-")

df.ss.t <- df.s.t
df.ss.t[, coi] <- df.s.t[, coi] ^ 2

# Sum by marker/column
sst <- apply(df.ss.t[, coi], 2, sum)


# Format and merge
df.sst <- as.data.frame(sst)
df.sst.2 <- cbind(marker = rownames(df.sst), 
                  data.frame(df.sst, row.names = NULL))
df.sst.2$state <- "all"
df.sst.2$subset <- "all"
df.sst.2$donor <- "all"



#### Calculating SSW ####
# Get names for each unique combination
ssw.all <- unique(df.t[c("State", "Subset", "Donor")])

# Get variables of interest for subsetting
lvls.state <- levels(factor(df.t$State))
lvls.subset <- levels(factor(df.t$Subset))
lvls.donor <- levels(factor(df.t$Donor))

class(lvls.state) # Expecting "character"



#### Calculate sum of SSW/SST by sample (by state, subset, and donor) ####
# Auto subset and calculate SSW
list.ssw <- list()

for (i in 1:nrow(ssw.all)) {
  a <- ssw.all$State[[i]]
  b <- ssw.all$Subset[[i]]
  c <- ssw.all$Donor[[i]]
  print(paste0("State: ", a))
  print(paste0("Subset: ", b))
  print(paste0("Donor: ", c))
  
  df <- df.t %>%
    filter(State == a & Subset == b & Donor == c)
  
  # Subtract mean and square each value
  mean.df <- apply(df[, coi], 2, mean)
  
  df.s <- df
  df.s[, coi] <- sweep(df[, coi], 2, mean.df, FUN = "-")
  
  df.ss <- df.s
  df.ss[, coi] <- df.s[, coi] ^ 2
  
  # Sum by marker/column
  ssw <- apply(df.ss[, coi], 2, sum)
  
  # Format and merge
  df.ssw <- as.data.frame(ssw)
  df.ssw.2 <- cbind(marker = rownames(df.ssw), 
                    data.frame(df.ssw, row.names = NULL))
  df.ssw.2$state <- a
  df.ssw.2$subset <- b
  df.ssw.2$donor <- c
  
  # Store new SSW
  list.ssw[[i]] <- df.ssw.2
}


# Merge SSW dfs
df.ssw.m <- rbind.fill(list.ssw)
head(df.ssw.m)

# Check levels/unique combinations
levels(factor(df.ssw.m$state))

# Calculate sum of SSW of all conditions (state, subset, and donor) over total SST
df.ssw.sum <- df.ssw.m %>%
  group_by(marker) %>%
  summarise(sum(ssw))

# Join SSW and SST by matching on values in column for marker
df.ssw.sst <-
  left_join(df.ssw.sum, df.sst.2, by = 'marker')


# Calculate SSW/SST
df.ssw.sst$ssw.sst <-
  df.ssw.sst$`sum(ssw)` / df.ssw.sst$sst



#### Calculate sum of SSW over SST of baseline and activated within SST ####
# Auto subset and calculate SSW
list.ssw.state <- list()

for (j in lvls.state) {
  print(j)
  
  df <- df.t %>%
    filter(State == j)
  
  # Subtract mean and square each value
  mean.df <- apply(df[, coi], 2, mean)
  
  df.s <- df
  df.s[, coi] <- sweep(df[, coi], 2, mean.df, FUN = "-")
  
  df.ss <- df.s
  df.ss[, coi] <- df.s[, coi] ^ 2
  
  # Sum by marker/column
  ssw.state <- apply(df.ss[, coi], 2, sum)
  
  # Format and merge
  df.ssw.state <- as.data.frame(ssw.state)
  df.ssw.state.2 <- cbind(marker = rownames(df.ssw.state), 
                          data.frame(df.ssw.state, row.names = NULL))
  df.ssw.state.2$state <- j
  
  # Store new SSW
  list.ssw.state[[j]] <- df.ssw.state.2
}


# Merge SSW dfs
df.ssw.state.m <- rbind.fill(list.ssw.state)
nrow(df.ssw.state.m) # Expecting 24 rows for 12 markers at baseline and activated

# Calculate sum of SSW of all conditions (state, subset, and donor) over total SST
df.ssw.state.sum <- df.ssw.state.m %>%
  group_by(marker) %>%
  summarise(sum(ssw.state))

# Join SSW and SST by matching on values in column for marker
df.ssw.state.sst <-
  left_join(df.ssw.state.sum, df.sst.2, by = 'marker')

# Calculate SSW/SST
df.ssw.state.sst$ssw.sst <-
  df.ssw.state.sst$`sum(ssw.state)` / df.ssw.state.sst$sst



#### Calculate sum of SSW/SST by donor ####
# Auto subset and calculate SSW
list.ssw.donor <- list()

for (k in lvls.donor) {
  print(k)
  
  df <- df.t %>%
    filter(Donor == k)
  
  # Subtract mean and square each value
  mean.df <- apply(df[, coi], 2, mean)
  
  df.s <- df
  df.s[, coi] <- sweep(df[, coi], 2, mean.df, FUN = "-")
  
  df.ss <- df.s
  df.ss[, coi] <- df.s[, coi] ^ 2
  
  # Sum by marker/column
  ssw.donor <- apply(df.ss[, coi], 2, sum)
  
  # Format and merge
  df.ssw.donor <- as.data.frame(ssw.donor)
  df.ssw.donor.2 <- cbind(marker = rownames(df.ssw.donor), 
                          data.frame(df.ssw.donor, row.names = NULL))
  df.ssw.donor.2$donor <- k
  
  # Store new SSW
  list.ssw.donor[[k]] <- df.ssw.donor.2
}


# Merge SSW dfs
df.ssw.donor.m <- rbind.fill(list.ssw.donor)
nrow(df.ssw.donor.m) # Expecting 24 rows for 12 markers at baseline and activated

# Calculate sum of SSW of all conditions (state, subset, and donor) over total SST
df.ssw.donor.sum <- df.ssw.donor.m %>%
  group_by(marker) %>%
  summarise(sum(ssw.donor))

# Join SSW and SST by matching on values in column for marker
df.ssw.donor.sst <-
  left_join(df.ssw.donor.sum, df.sst.2, by = 'marker')

# Calculate SSW/SST
df.ssw.donor.sst$ssw.sst <-
  df.ssw.donor.sst$`sum(ssw.donor)` / df.ssw.donor.sst$sst

#### End SSW/SST analysis ####


#### Copy as data matrix ####
flow12 <- flow11

names(flow12)

flow12.mt <- as.matrix(flow12[coi])

head(flow12.mt)



#### Start PDF ####
pdf.name <- paste0("PDF_", out.name, ".pdf")

CairoPDF(
  pdf.name,
  width = 11.5,
  height = 8,
  onefile = TRUE,
  family = "Helvetica"
)



#### Histograms ####
par(mfrow = c(3, 4)) # mfcol plots by column (down) while mfrow plots by row (left to right)

hist(flow12$CD4_scUncomp)
hist(flow12$CD45RO_scUncomp)
hist(flow12$CD69_scUncomp)
hist(flow12$CD40L_scUncomp)
hist(flow12$CD25_scUncomp)
hist(flow12$HLA.DR_scUncomp)
hist(flow12$CD38_scUncomp)
hist(flow12$ICOS_scUncomp)
hist(flow12$CTLA.4_scUncomp)
hist(flow12$PD.1_scUncomp)
hist(flow12$OX40_scUncomp)
hist(flow12$CCR6_scUncomp)

hist(flow12$CD4_chComp)
hist(flow12$CD45RO_chComp)
hist(flow12$CD69_chComp)
hist(flow12$CD40L_chComp)
hist(flow12$CD25_chComp)
hist(flow12$HLA.DR_chComp)
hist(flow12$CD38_chComp)
hist(flow12$ICOS_chComp)
hist(flow12$CTLA.4_chComp)
hist(flow12$PD.1_chComp)
hist(flow12$OX40_chComp)
hist(flow12$CCR6_chComp)

par(mfrow = c(1, 1))



#### SSW/SST of donors ####
ggplot(df.ssw.donor.sst, aes(x = marker, y = ssw.sst)) +
  geom_col() +
  geom_text(
    aes(
      x = marker,
      y = ssw.sst,
      label = round(ssw.sst, digits = 3)
    ),
    size = 4,
    position = position_dodge(width = 1),
    hjust = -0.5,
    inherit.aes = TRUE
  ) +
  geom_hline(yintercept = 1) +
  #  scale_y_continuous(limits = c(0, 1.15)) +
  coord_flip() +
  ggtitle(paste0(
    out.name,
    "\nSSW/SST of donors / Total \nSum of SSW/SST = ",
    round(sum(df.ssw.donor.sst$ssw.sst), 3)
  )) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18, face = "bold")
  )


#### SSW/SST of state/total ####
ggplot(df.ssw.state.sst, aes(x = marker, y = ssw.sst)) +
  geom_col() +
  geom_text(
    aes(
      x = marker,
      y = ssw.sst,
      label = round(ssw.sst, digits = 3)
    ),
    size = 4,
    position = position_dodge(width = 1),
    hjust = -0.5,
    inherit.aes = TRUE
  ) +
  geom_hline(yintercept = 1) +
  coord_flip() +
  ggtitle(paste0(
    out.name,
    "\nSSW/SST of Baseline & Activated / Total \nSum of SSW/SST = ",
    round(sum(df.ssw.state.sst$ssw.sst), 3)
  )) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18, face = "bold")
  )


#### SSW/SST of all samples ####
ggplot(df.ssw.sst, aes(x = marker, y = ssw.sst)) +
  geom_col() +
  geom_text(
    aes(
      x = marker,
      y = ssw.sst,
      label = round(ssw.sst, digits = 3)
    ),
    size = 4,
    position = position_dodge(width = 1),
    hjust = -0.5,
    inherit.aes = TRUE
  ) +
  geom_hline(yintercept = 1) +
  #  scale_y_continuous(limits = c(0, 1.15)) +
  coord_flip() +
  ggtitle(paste0(
    out.name,
    "\nSSW/SST of all samples / Total \nSum of SSW/SST = ",
    round(sum(df.ssw.sst$ssw.sst), 3)
  )) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18, face = "bold")
  )



#### PCA ####
autoplot(
  kmeans(
    flow12[, coi],
    centers = 7,
    iter.max = 1000,
    nstart = 20,
    algorithm = "Lloyd"
  ),
  main = paste0("k means clustering on flow cytometry data post-QC \n",
                pdf.name),
  data = flow12,
  asp = 1
)

flow_pcaPr <- prcomp(flow12[, coi], center = F, scale. = F)
summary(flow_pcaPr)

plot(flow_pcaPr, type = "lines", main = "PCA of flow cytometry data post-QC \nCD4 T cells activated o/n with 1 Dynabead per 1 cell \nPCs determined by prcomp")
title(xlab = "Principle Component")



#### Phenograph and SNE ####

# Phenograph
cluster_PhenoGraph <-
  cytof_cluster(xdata = flow12.mt, method = "Rphenograph") # Default for Rphenograph k = 30

flow12.mt.tsne <-
  cytof_dimReduction(data = flow12.mt, method = "tsne")

# cluster_FlowSOM <- cytof_cluster(xdata = flow12.mt, method = "FlowSOM", FlowSOM_k = 30)


# Concatenate data
flow12.all <- cbind(flow12, flow12.mt.tsne,
                   PhenoGraph = cluster_PhenoGraph)

# PhenoGraph
cytof_clusterPlot(
  data = flow12.all,
  xlab = "tsne_1",
  ylab = "tsne_2",
  cluster = "PhenoGraph",
  sampleLabel = FALSE,
  title = pdf.name
)


# Heatmap of PhenoGraph clusters
names(flow12.all) # Getting compensated and PhenoGraph column numbers

coi.comp.pheno <- c(35:38, 40:46, 49, 53)

PhenoGraph_cluster_median <- aggregate(. ~ PhenoGraph, 
                                       flow12.all[, c(coi.comp.pheno)], 
                                       median)

cytof_heatmap(PhenoGraph_cluster_median[, 2:13], baseName = "PhenoGraph Cluster Median, compensated")


PhenoGraph_cluster_median2 <- aggregate(. ~ PhenoGraph, 
                                       flow12.all[, c(coi, 53)], 
                                       median)

cytof_heatmap(PhenoGraph_cluster_median2[, 2:13], baseName = "PhenoGraph Cluster Median, uncompensated")



#### tSNE annotations ####
ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = factor(State, labels = labels.state))) +
  labs(color = "State") +
  geom_point(show.legend = T, size = 0.25, alpha = 0.4, position = position_jitter()) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_manual(values = c("Baseline" = "deepskyblue2", 
                                "Activated" = "deeppink1")) +
  theme_light() +
  ggtitle(paste0(out.name, "\nBaseline vs. Activated"))


ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = factor(Donor))) +
  labs(color = "Donor") +
  geom_point(show.legend = T, size = 0.25, alpha = 0.4, position = position_jitter()) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_light() +
  ggtitle(paste0(out.name, "\nBy Donor"))


ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = factor(Subset, labels = labels.subset))) +
  labs(color = "Subset") +
  geom_point(show.legend = T, size = 0.25, alpha = 0.4, position = position_jitter()) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_light() +
  scale_color_manual(values = c("Naive" = "#a6cee3", 
                                "Central Memory" = "#1f78b4",
                                "Effector Memory" = "#b2df8a")) +
  ggtitle(paste0(out.name, "\nBy Subset"))
# For Reg: "Regulatory" = "#33a02c")) +



ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = factor(Batch))) +
  labs(color = "Batch") +
  geom_point(show.legend = T, size = 0.25, alpha = 0.4, position = position_jitter()) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_light() +
  ggtitle(paste0(out.name, "\nBy Batch (date)"))



#### Colored on compensated channel values ####
gg1.2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                                      color = CD4_chComp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCD4, colored on compensated"))


gg2.2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                                      color = CD45RO_chComp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCD45RO, colored on compensated"))


gg3.2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                                      color = CD69_chComp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +  
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCD69, colored on compensated"))


gg4.2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                                      color = CD40L_chComp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCD40L, colored on compensated"))


gg5.2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                                      color = CD25_chComp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCD25, colored on compensated"))


gg6.2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                                      color = HLA.DR_chComp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nHLA-DR, colored on compensated"))


gg7.2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                                      color = CD38_chComp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCD38, colored on compensated"))


gg8.2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                                      color = ICOS_chComp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nICOS, colored on compensated"))


gg9.2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                                      color = CTLA.4_chComp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCTLA-4, colored on compensated"))


gg10.2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                                       color = PD.1_chComp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nPD-1, colored on compensated"))


gg11.2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                                       color = OX40_chComp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nOX40, colored on compensated"))


gg12.2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                                       color = CCR6_chComp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCCR6, colored on compensated"))


grid.arrange(gg1.2, gg2.2, gg3.2, gg4.2, gg5.2, gg6.2, gg7.2, gg8.2, gg9.2, gg10.2, gg11.2, gg12.2, ncol=4)



#### Colored on uncompensated ####
gg1 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = CD4_scUncomp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCD4, colored on uncompensated"))


gg2 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = CD45RO_scUncomp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCD45RO, colored on uncompensated"))


gg3 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = CD69_scUncomp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCD69, colored on uncompensated"))


gg4 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = CD40L_scUncomp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCD40L, colored on uncompensated"))


gg5 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = CD25_scUncomp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCD25, colored on uncompensated"))


gg6 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = HLA.DR_scUncomp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nHLA-DR, colored on uncompensated"))


gg7 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = CD38_scUncomp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCD38, colored on uncompensated"))


gg8 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = ICOS_scUncomp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nICOS, colored on uncompensated"))


gg9 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = CTLA.4_scUncomp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCTLA-4, colored on uncompensated"))


gg10 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = PD.1_scUncomp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nPD-1, colored on uncompensated"))


gg11 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = OX40_scUncomp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nOX40, colored on uncompensated"))


gg12 <- ggplot(data = flow12.all, aes(x = tsne_1, y = tsne_2,
                             color = CCR6_scUncomp)) +
  geom_point(show.legend = F, size = 0.25, alpha = 0.2, position = position_jitter()) +
  scale_color_viridis(option = "inferno", guide = F) +
  theme_light() +
  ggtitle(paste0(out.name, "\nCCR6, colored on uncompensated"))


grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9, gg10, gg11, gg12, ncol=4)



#### heatmap (what's defining these clusters?) ####
# From Kam's tutorial: http://slowkow.com/notes/heatmap-tutorial/

# Pheatmap can only handle up to 65k rows, and my R crashes before it finishes processing that, so subsampling 20k rows for analysis instead.
flow12.1 <- sample_n(flow12[, coi], 20000, replace = T)

## Setting coloring scheme based on quantile breaks so coloring will represent equal proportion of data
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0 , 1, length.out = n))
  breaks[!duplicated(breaks)]
}

# Quantile function dislikes data.frames -> coerce into matrix
flowHM_breaks <- quantile_breaks(as.matrix(flow12.1), n = 11)


## Sorting the dendrograms for easier interpretation
# Cluster by column
flowHM_cluster_cols <- hclust(dist(t(flow12.1)))
# Results from unsorted analysis
#plot(flowHM_cluster_cols, main = "Unsorted Dendrogram")

# Sort clustering by column
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

flowHM_cluster_cols <- sort_hclust(flowHM_cluster_cols)
# plot(flowHM_cluster_cols, main = "Sorted Dendrogram")

# Sorting for rows
flowHM_cluster_rows <- sort_hclust(hclust(dist(flow12.1)))


## Plotting sorted heatmap
pheatmap(
  mat = flow12.1,
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


dev.off()
#### End PDF ####



#### Notes ####

# For installing Bioconductor packages:
#source("https://bioconductor.org/biocLite.R")
#biocLite("flowCore", lib="/Users/chong/R/library")


#### On exporting data from FlowJo as scale vs. channel ####
# From http://flowjo.typepad.com/the_daily_dongle/2006/07/channels_vs_sca.html:
## In analog data, the scale values are already in log. In Digital data, scale is in lin, and channels are in log. Another difference in digital data is that FlowJo's channel information is on a log-transformed 14 bit scale, whereas the original information is in 18 bits linear dynamic range. You could definintely see the advantage there - the linear scale option contains 64 times more information. However, before any meaningful fluorescence analysis, the scale data will have to somehow be log-scaled. Why? according to this paper's first reference,
## "Logarithmic data scales, which show log-normal distributions as symmetrical peaks, are widely used and accepted as those facilitating analysis of fluorescence measurements in biological systems."
## So the short answer is, scale contains more information about the data but it may not be transformed in logarithmic fashion. The longer answer would probably involve a philosophical debate about how most of the 18 bit resolution goes to waste because most of biological data has wide CV's.
## Everything in FlowJo (statistics, gating, clustering) is done in channel.

# Another answer from http://flowjo.typepad.com/the_daily_dongle/2006/06/qa_scales_chann.html:
## I am using FlowJo to export .fcs files to .txt files. I noticed that there are two ways of exporting data: by channel and by scale. Which one should I use?
## channels correspond to values which have been scaled to the way you're used to look at your data. This includes lin/log scaling and the scale you see in your graph window axis labels.  scale corresponds to the raw output of the cytometer. these values will also allow you to analyze the data, but they are not going to yield distributions you're used to look at.  So When you export your data files, use "channels" options.



#### Error message when sample doesn't exist ####
# When import is looking for a file that doesn't exist (e.g. a missing Treg sample), it'll print the following error into the output:
# "character(0)
# ERROR : subscript out of bounds"
# May later figure out how to print a more informative error message.


#### Adjusting for negatives ####
# When the BD Fortessa collects data, there's noise around the measurement of fluorescence intensity, which is why the populations form log-normal-esque peaks.  The Fortessa does corrections for this noise, and the noise + corrections causes some of the lower values to be labeled as negative.  (These are clearly artificial, as the fluorophores aren't negatively fluorescent.)
# Compesation further spreads signal from the negative population below zero.
# When transforming, it's necessary to adjust the negatives to maintain their relationship with the positive values.  The problem isn't the existence of negative values themselves, but what negative values mean relative to the positive values.
# With financial data, a positive value (e.g. income) is meaningfully different from a negative value (e.g. debt).  Then the fact that the base asinh transform has an inflection point at 0 is harmless - maybe even helpful - because 0 is a meaningful boundary in the actual data.
# On the other hand, with raw flow cytometry data, negative values are _not_ meaningfully different from the positive ones.  Technically, they should all be positive (there's no such thing as negative fluorescence), but they aren't for reasons that would take too long to write here.  But in this case, using asinh on the raw data creates an artificial distinction by separating the values around 0, which eventually leads to non-meaningful visualization and false clusters.
# Adjusting by adding x(min) for each column (after outlier QC) is a simply way to correct for artificial negatives.  This also decreases the SSW/SST of the known cell types, which indicates better separation (higher SSB), which is preferred for accurate clustering.
