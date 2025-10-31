#### Assignment ####

# 1. Set Working Directory
# Create a new folder on your computer "AI_Omics_Internship_2025".

# 2. Create Project Folder

# In RStudio, create a new project named "Module_I" in your "AI_Omics_Internship_2025" folder.

# Inside the project directory, create the following subfolders using R code:
# raw_data, clean_data, scripts, results or Tasks, plots etc

# ---------------------------------------------------------------------------
# 3. Download "patient_info.csv" dataset from GitHub repository

# load the dataset into your R environment

# Inspect the structure of the dataset using appropriate R functions

# Identify variables with incorrect or inconsistent data types.

# Convert variables to appropriate data types where needed

# Create a new variable for smoking status as a binary factor:

# 1 for "Yes", 0 for "No"

# Save the cleaned dataset in your clean_data folder with the name patient_info_clean.csv
# Save your R script in your script folder with name "class_Ib"
# Upload "class_Ib" R script into your GitHub repository

####Solution####

#Check working directory
getwd()

#Set working directory
setwd("/Users/adil/Desktop/AI_Omics_Internship_2025/Module_I")

#Creating subfolders inside project directory
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("tasks")
dir.create("plots")

#Loading "patient_info.csv" into R
data <- read.csv(file.choose())

# Check structure of  data
str(data)

# Convert 'smoker' to factor
data$smoker_fac <- as.factor(data$smoker)
str(data)

# Convert smoker factor to numeric using ifelse statement (Smoker = 1, Non-Smoker = 0)
data$smoking_status <- ifelse(data$smoker_fac == "Yes", 1, 0)
str(data)

# Save the entire R workspace
save.image(file = "class_Ib.RData")