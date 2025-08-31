#### Assignment 2 #####
# In this assignment you will work with the results of differential gene expression (DGE) analysis. 
#The analysis produces two key measures for each gene:

# log2FoldChange (log2FC): 
# Indicates the magnitude and direction of change in gene expression. 
# Positive values suggest higher expression(upregulated gene) in the experimental condition compared to control. 
# Negative values suggest lower expression (downregulated gene). 
# The absolute value reflects the strength of the change.

# Adjusted p-value (padj): 
# Represents the statistical significance of the observed difference, corrected for multiple testing. 
# A smaller value indicates stronger evidence that the observed difference is not due to chance.

# Write a function classify_gene() 

# that takes:
#   - logFC (log2FoldChange)
#   - padj  (adjusted p-value)

# and returns:
#   - "Upregulated" if log2FC > 1 and padj < 0.05
#   - "Downregulated" if log2FC < -1 and padj < 0.05
#   - "Not_Significant" otherwise


# Then:
#   - Apply it in a for-loop to process both datasets (DEGs_data_1.csv, DEGs_data_2.csv)
#   - Replace missing padj values with 1
#   - Add a new column 'status'
#   - Save processed files into Results folder
#   - Print summary counts of significant, upregulated, and downregulated genes
#   - Use table() for summaries

# Data Availability
# The input files are available in the GitHub repository:
#      DEGs_Data_1.csv
#      DEGs_Data_2.csv

# Each file contains three columns: 
# Gene_Id	
# padj	
# logFC

#### Solution #####

# Function to classify gene
classify_gene <- function(logFC, padj) {
  if (logFC > 1 & padj < 0.05) {return("Upregulated")} 
  else if (logFC < -1 & padj < 0.05) {return("Downregulated")} 
  else {return("Not_Significant")}
}

# Define input and output folders
input_dir <- "data" 
output_dir <- "results"

# Create output folder if it does not already exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# List of files to process
files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

# Prepare empty list to store results
result_list <- list()

# For-loop to process both datasets. 
for (file_name in files_to_process) {
  cat("\nProcessing:", file_name, "\n")
  input_file_path <- file.path(input_dir, file_name)
  
  # Import dataset
  degs_data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing values...\n")
  
  # Handle missing padj values
  if ("padj" %in% names(degs_data)) {
    missing_count <- sum(is.na(degs_data$padj))
    cat("Missing values in 'padj':", missing_count, "\n")
    degs_data$padj[is.na(degs_data$padj)] <- 1
  }
  
  # Classify genes
  degs_data$status <- mapply(classify_gene, degs_data$logFC, degs_data$padj)
  cat("Gene has been classified successfully.\n")
  
  # Save results in R
  result_list[[file_name]] <- degs_data 
  
  # Save results in Results folder
  output_file_path <- file.path(output_dir, paste0("Gene_Expression_", file_name))
  write.csv(degs_data, output_file_path, row.names = FALSE)
  cat("Results saved to:", output_file_path, "\n")
  
  # Print summary with table()
  summary_counts <- table(degs_data$status)
  print(paste("Summary for", file_name))
  print(summary_counts)
}
