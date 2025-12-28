# ================== AI and Omics Internship Final Project ==================== #
 


# === Bioconductor provides R packages for analyzing omics data === #

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")



# === Install Bioconductor packages === #

BiocManager::install(c("GEOquery","affy","arrayQualityMetrics", "oligo",
                       "pd.hugene.1.0.st.v1", "hugene10sttranscriptcluster.db",
                       "R.utils"))



# === Load Required packages === #

library(GEOquery) # Download GEO datasets (series matrix, raw CEL files)
library(affy)
library(arrayQualityMetrics)  # QC reports for microarray data
library(oligo)  # Pre-processing of Affymetrix ST microarray data (RMA normalization)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
library(R.utils)
library(dplyr) # Data manipulation



# === Download Raw Data (CEL files) === #
# CEL files contain raw probe-level intensity values

getGEOSuppFiles("GSE28735", baseDir = "Raw_Data", makeDirectory = TRUE)



# === Untar CEL files === #

untar("Raw_Data/GSE28735/GSE28735_RAW.tar", 
      exdir = "Raw_Data/GSE28735/CEL_Files")



# === Unzip files === #

gz_files <- list.files("Raw_Data/GSE28735/CEL_Files", 
                       pattern="\\.gz$", full.names=TRUE)
lapply(gz_files, gunzip, remove=FALSE)



# === Read CEL files === #

raw_data <- read.celfiles(list.files("Raw_Data/GSE28735/CEL_Files", 
                                     pattern="\\.CEL$", full.names=TRUE))



# === Quality Control (QC) Before Pre-processing === #

arrayQualityMetrics(expressionset = Raw_Data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)



# === RMA (Robust Multi-array Average) Normalization === #

Normalized_Data = rma(Raw_Data)



# === QC after data normalization === #

arrayQualityMetrics(expressionset = Normalized_Data,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)



# === Extract normalized expression values into a data frame === #

Processed_Data = as.data.frame(exprs(Normalized_Data))
dim(Processed_Data) # -> Dimensions --> number of probes × number of samples


 
# === Filter Low-Variance Transcripts (“soft” intensity based filtering) === #
# Calculate median intensity per probe across samples 

Row_Median = rowMedians(as.matrix(Processed_Data))




# == Visualize distribution of Probe Median intensities === #

hist(Row_Median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")



#  === Set a threshold to remove low variance probes === #

Threshold = 3.5 



# === Select probes above threshold === #

Index = Row_Median > Threshold 
Filtered_Data = Processed_Data[Index, ] 



# === Overwrite processed data with filtered data === #

Processed_Data = Filtered_Data



# === Extract probe IDs from processed microarray data === #

Probe_IDs = rownames(Processed_Data)




# === Map probe IDs to gene symbols using the platform annotation database === #

Gene_Symbols = mapIds(
  hugene10sttranscriptcluster.db, # -> Database used for mapping
  keys = Probe_IDs,        # -> Input probe IDs
  keytype = "PROBEID",     # -> Probe ID key type
  column = "SYMBOL",       # -> Desired annotation column (gene symbols)
  multiVals = "first"      # -> Return first match if multiple exist
)



# === Convert mapping to a data frame and rename columns === #

Gene_Map_df = Gene_Symbols %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::rename(SYMBOL = 2)



# === Handle multiple probes mapping to a single gene === #
# Summarize number of probes per gene symbol

Duplicate_Summary = Gene_Map_df %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))

# === Identify genes associated with multiple probes === #

Duplicate_Genes = Duplicate_Summary %>%
  filter(probes_per_gene > 1)
sum(Duplicate_Genes$probes_per_gene)



# === Merge annotation mapping with expression data === #
# -> Verify if probe IDs in mapping correspond to expression data

all(Gene_Map_df$PROBEID == row.names(Processed_Data))



# === Merge annotation (SYMBOL) with expression matrix === #

Processed_Data_df = Processed_Data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = Gene_Symbols[PROBEID]) %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)



# === Remove probes without valid gene symbol annotation === #

Processed_Data_df = Processed_Data_df %>%
  dplyr::filter(!is.na(SYMBOL))



# === Select only numeric expression columns === #

Expressed_only = Processed_Data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)


# === Collapse multiple probes per gene using average expression === #
# === limma::avereps() -> computes the average for probes representing the same gene === #

Avg_data = limma::avereps(Expressed_only, ID = Processed_Data_df$SYMBOL)

# === Feature Selection with Boruta === #

install.packages("Boruta")
library(Boruta)
Boruta_Results  = Boruta

gse = getGEO("GSE28735", GSEMatrix = TRUE)
Expressed_Data = exprs(GSE_Data[[1]])
Expr = t(Expressed_Data)
Express = as.data.frame(Expr)
Expressed_Data$Class = as.factor(Phenotype_Data$characteristics_ch1)
Phenotype_Data = pData(GSE_Data[[1]])
levels(Expressed_Data$Class)

Train = as.data.frame(Expressed_Data)
class(Train)
is.factor(Train$Class)

set.seed(123)
Boruta_Output = Boruta(Class ~ ., data = train, doTrace = 0)

# === Select only confirmed features === #

Confirmed_Gene <- getSelectedAttributes(Boruta_Output, withTentative = FALSE)
length(Confirmed_Gene)



final_Data  = colnames(train)[colnames(train) %in% Confirmed_Gene]
length(final_Data)



set.seed(123)
boruta_output <- Boruta(Class ~ ., data = train, doTrace = 0)

Control <- rfeControl(functions = rfFuncs,
                      method = "cv",
                      verbose = FALSE)

set.seed(123)
rfe_Output = rfe(train[, -ncol(train)], train$Class,
                  sizes = c(5, 10, 20, 30, 40, 50),
                  rfeControl = Control)

rfe_genes = predictors(rfe_Output)
length(rfe_genes)

# === Spliting the Dataset into training === #

set.seed(123)
Indx = sample(seq_len(nrow(final_Data)),
              size = 0.7*
      nrow(final_Data))
Train_ML = final_data[idx, ]
Test_ML = final_Data[-idx, ]


# Training a Random Forest classification for the training dataset === #

install.packages("randomforest")
library(randomForest)

RF_Model = randomForest(
  class~ .,
  data = Train_ML,
  ntree = 500
  )

# === Predicting Class Labels for the test dataset === #

Predict_RF = predict(RF_Model, newdata = Test_ML)

# === Confusion Matrix and Calculating Accuracy === #

Confusion_Matrix = table(Predicted = Predict_RF, Actual = Test_ML$Class)
Accuracy = mean(Predict_RF == Test_ML$Class)


varImpPlot(RF_Model, n.var = 20)