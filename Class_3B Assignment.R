#### Assignment ####
#-------------

# Work through the full preprocessing workflow with your own dataset

# 1. Perform quality control before and after normalization and 
# check whether any arrays are flagged as outliers. 
# note down how many you found before and after normalization

# 2. Normalize the data and then apply filtering to remove low-intensity probes 
# and note how many transcripts remain. 

# 3. Use the phenotype information to define your target groups and re-label them (e.g normal vs cancer)

#### Solution ####

# Install Bioconductor 
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))

# Install dplyr
install.packages("dplyr")

# Install Hmisc (required for arrayQualityMetrics)
install.packages("Hmisc")

# Load Required Libraries
library(GEOquery)  
library(affy)                 
library(arrayQualityMetrics)  
library(dplyr)                

#### Download Series Matrix Files
gse_data <- getGEO("GSE79973", GSEMatrix = TRUE)

# Extract phenotype data
phenotype_data <-  pData(gse_data$GSE79973_series_matrix.txt)

#### Download Raw Data ####

# Raw data requires full preprocessing (e.g., RMA normalization, QC)

# CEL files are large. 
#It's recommended to download raw data directly from NCBI GEO

# Fetch GEO supplementry files
getGEOSuppFiles("GSE79973", baseDir = "Raw_Data", makeDirectory = TRUE)


# Untar CEL files if compressed as .tar
untar("Raw_Data/GSE79973_RAW.tar", exdir = "Raw_Data/CEL_Files")

# Read CEL files into R
raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_Files")

raw_data   


# Quality Control (QC) before normalization
arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)


#### Number of outliers before normalization are 5 i.e. sample 1,7,8,10,20


# RMA (Robust Multi-array Average) Normalization
normalized_data <- rma(raw_data)

# QC after data normalization
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)


#### Number of outliers after normalization are 1 i.e. sample 8


# Extract normalized expression values into a data frame
processed_data <- as.data.frame(exprs(normalized_data))

dim(processed_data)


#### Filtering to remove low-intensity probes ####

# Calculate median intensity per probe across samples
row_median <- rowMedians(as.matrix(processed_data))

row_median

# plot median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

threshold <- 3.5 
abline(v = threshold, col = "black", lwd = 2) 

# Select probes above threshold
indx <- row_median > threshold 
filtered_data <- processed_data[indx, ] 

# Rename filtered expression data with sample metadata
colnames(filtered_data) <- rownames(phenotype_data)

# Overwrite processed data with filtered dataset
processed_data <- filtered_data 

# Number of transcripts remaining are 41601.

#### Phenotype Data Preparation ####

class(phenotype_data$source_name_ch1)

# Define experimental groups (normal vs cancer)
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("gastric mucosa", "gastric adenocarcinoma"),
                 label = c("normal", "cancer"))

class(groups)
levels(groups)

