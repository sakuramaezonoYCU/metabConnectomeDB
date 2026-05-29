#!/Library/Frameworks/R.framework/Versions/4.6/Resources/bin/Rscript
# Script Purpose: Download raw TCGA files for our 5 cancer types directly to input/TCGA/
# We avoid merging in R to prevent memory crashes and timeout issues.

rm(list=ls()) # Clean up environment
set.seed(123) # set random state
options(timeout = 3600) # Increase timeout to 1 hour for large TCGA files

### R packages
required_packages <- c("tidyverse", "UCSCXenaTools")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")

library(tidyverse)
library(UCSCXenaTools)

dir_input <- "input/TCGA"
dir.create(dir_input, recursive = TRUE, showWarnings = FALSE)

# We only need the 5 cancer types currently analyzed in this project:
# Breast (BRCA), Colorectal (COAD, READ), Lung (LUAD, LUSC), Melanoma (SKCM), Ovarian (OV)
tcga_cancer_codes <- c(
  "BRCA", "COAD", "READ", "LUAD", "LUSC", "SKCM", "OV"
)

#### Step 1: Download TCGA data (gene exp matrix, metadata, survival data) ####

downloadTCGARaw <- function(cancer) {
  message("--------------------------------------------------")
  message("Processing TCGA-", cancer, "...")
  df_todo_cancer <- UCSCXenaTools::XenaGenerate(subset = XenaHostNames %in% c("tcgaHub", "gdcHub")) |>
    UCSCXenaTools::XenaFilter(filterDatasets = cancer) 
  
  # 1. Download FPKM Gene Expression
  exp_filter <- paste0("TCGA-", cancer, ".star_fpkm.tsv")
  exp_query <- df_todo_cancer |>
    UCSCXenaTools::XenaFilter(filterDatasets = exp_filter)
    
  # Helper to safely check if query has rows
  has_datasets <- function(q) {
    if (is.data.frame(q)) return(nrow(q) > 0)
    if (isS4(q)) return(length(q@datasets) > 0)
    return(FALSE)
  }
  
  if (has_datasets(exp_query)) {
    message("Downloading FPKM expression data for ", cancer)
    UCSCXenaTools::XenaQuery(exp_query) |> 
      UCSCXenaTools::XenaDownload(destdir = dir_input, force = FALSE)
  } else {
    message("WARNING: FPKM dataset not found for ", cancer)
  }
  
  # 2. Download Survival Data
  surv_filter1 <- paste0("TCGA-", cancer, ".survival.tsv")
  surv_filter2 <- paste0("survival/", cancer, "_survival.txt")
  
  surv_query1 <- df_todo_cancer |> UCSCXenaTools::XenaFilter(filterDatasets = surv_filter1)
  surv_query2 <- df_todo_cancer |> UCSCXenaTools::XenaFilter(filterDatasets = surv_filter2)
  
  if (has_datasets(surv_query1)) {
    message("Downloading Survival data (GDC format) for ", cancer)
    UCSCXenaTools::XenaQuery(surv_query1) |> 
      UCSCXenaTools::XenaDownload(destdir = dir_input, force = FALSE)
  } else if (has_datasets(surv_query2)) {
    message("Downloading Survival data (TCGA Hub format) for ", cancer)
    UCSCXenaTools::XenaQuery(surv_query2) |> 
      UCSCXenaTools::XenaDownload(destdir = dir_input, force = FALSE)
  } else {
    message("WARNING: Survival dataset not found for ", cancer)
  }
}

# Apply to all cancers
lapply(tcga_cancer_codes, downloadTCGARaw)

message("--------------------------------------------------")
message("All downloads completed! Files are saved in ", dir_input)
