#!/usr/bin/env Rscript

# This script downloads the complete METABRIC dataset (all genes) 
# using the cBioDataPack function which downloads curated tarballs 
# from Bioconductor's ExperimentHub (AWS) rather than the REST API or GitHub.

if (!require("cBioPortalData", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
  }
  BiocManager::install("cBioPortalData")
}

suppressPackageStartupMessages(library(cBioPortalData))

cat("Initializing cBioPortal API...\n")
# We don't need the cbio connection for cBioDataPack, but it's fine.

# Ensure the data directory exists
dir.create("input/metabric", recursive = TRUE, showWarnings = FALSE)

cat("Fetching complete MultiAssayExperiment for brca_metabric via cBioDataPack (ExperimentHub)...\n")
cat("Note: Downloading the full tarball might take a few minutes. Please wait.\n")

# Use cBioDataPack which is designed for downloading full studies (all genes).
# This avoids the REST API restrictions and pagination issues.
mae <- cBioDataPack("brca_metabric", ask = FALSE)

cat("Extracting clinical data...\n")
clin <- as.data.frame(colData(mae))

if ("patientId" %in% colnames(clin)) {
    colnames(clin)[colnames(clin) == "patientId"] <- "PATIENT_ID"
} else if ("PATIENT_ID" %in% toupper(colnames(clin))) {
    colnames(clin)[toupper(colnames(clin)) == "PATIENT_ID"] <- "PATIENT_ID"
}

write.csv(clin, "input/metabric/clinical.csv", row.names = FALSE)
cat("Successfully saved clinical data to input/metabric/clinical.csv\n")

cat("Extracting mRNA expression data...\n")
# The profile is typically named "brca_metabric_mrna"
if ("mrna_illumina_microarray" %in% names(mae)) {
    expr_se <- mae[["mrna_illumina_microarray"]]
    expr_matrix <- assay(expr_se)
    expr_df <- as.data.frame(expr_matrix)
    
    expr_df$Hugo_Symbol <- rownames(expr_df)
    expr_df <- expr_df[, c("Hugo_Symbol", setdiff(colnames(expr_df), "Hugo_Symbol"))]
    
    write.csv(expr_df, "input/metabric/expression.csv", row.names = FALSE)
    cat("Successfully saved expression data to input/metabric/expression.csv\n")
} else {
    cat("Available molecular profiles:\n")
    print(names(mae))
    stop("Could not find mRNA profile 'mrna_illumina_microarray' in the downloaded MAE.")
}

cat("\n✅ Data download complete. You can now run the Python script:\n")
cat("python scripts/generate_ml_prognostic_classifier_notebook.py\n")
