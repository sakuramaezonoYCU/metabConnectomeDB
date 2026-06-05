#!/usr/bin/env Rscript

# This script downloads a complete cBioPortal dataset (all genes) 
# using the cBioDataPack function which downloads curated tarballs 
# from Bioconductor's ExperimentHub (AWS) rather than the REST API or GitHub.

# Example way to use the script:
# Rscript scripts/download_cbioportal_data.R brca_aurora_2023 brca aurora
# Rscript scripts/download_cbioportal_data.R brca_mbcproject_2022 brca mbcproject

args = commandArgs(trailingOnly=TRUE)
study_id <- if (length(args) > 0) args[1] else "brca_metabric"
cancer_prefix <- if (length(args) > 1) args[2] else "brca"
database <- if (length(args) > 2) args[3] else "metabric"

if (!require("cBioPortalData", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
  }
  BiocManager::install("cBioPortalData")
}

suppressPackageStartupMessages(library(cBioPortalData))

cat("Initializing cBioPortal API...\n")

# Ensure the data directory exists
out_dir <- file.path("input", database)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Fetching complete MultiAssayExperiment for %s via cBioDataPack (ExperimentHub)...\n", study_id))
cat("Note: Downloading the full tarball might take a few minutes. Please wait.\n")

# Use cBioDataPack which is designed for downloading full studies (all genes).
mae <- cBioDataPack(study_id, ask = FALSE)

cat("Extracting clinical data...\n")
clin <- as.data.frame(colData(mae))

if ("patientId" %in% colnames(clin)) {
    colnames(clin)[colnames(clin) == "patientId"] <- "PATIENT_ID"
} else if ("PATIENT_ID" %in% toupper(colnames(clin))) {
    colnames(clin)[toupper(colnames(clin)) == "PATIENT_ID"] <- "PATIENT_ID"
}

clin_out <- file.path(out_dir, sprintf("%s_clinical.csv", cancer_prefix))
write.csv(clin, clin_out, row.names = FALSE)
cat(sprintf("Successfully saved clinical data to %s\n", clin_out))

cat("Extracting mRNA expression data...\n")
# Auto-detect the mRNA expression profile name
mae_names <- names(mae)
expr_profile <- NULL

# Common naming conventions for mRNA expression in cBioPortal
possible_patterns <- c("mrna_illumina_microarray", "rna_seq_v2_mrna", "mrna_U133_microarray", "mrna", "rna")

for (pat in possible_patterns) {
    matches <- grep(pat, mae_names, ignore.case=TRUE, value=TRUE)
    if (length(matches) > 0) {
        expr_profile <- matches[1]
        break
    }
}

if (!is.null(expr_profile)) {
    cat(sprintf("Found mRNA profile: %s\n", expr_profile))
    expr_se <- mae[[expr_profile]]
    expr_matrix <- assay(expr_se)
    expr_df <- as.data.frame(expr_matrix)
    
    # Check if rowData has Hugo_Symbol (often the case when rownames are Entrez IDs)
    rd <- as.data.frame(rowData(expr_se))
    if ("Hugo_Symbol" %in% colnames(rd)) {
        expr_df$Hugo_Symbol <- rd$Hugo_Symbol
    } else if ("hugo_gene_symbol" %in% colnames(rd)) {
        expr_df$Hugo_Symbol <- rd$hugo_gene_symbol
    } else {
        expr_df$Hugo_Symbol <- rownames(expr_df)
    }
    
    expr_df <- expr_df[, c("Hugo_Symbol", setdiff(colnames(expr_df), "Hugo_Symbol"))]
    
    expr_out <- file.path(out_dir, sprintf("%s_expression.csv", cancer_prefix))
    write.csv(expr_df, expr_out, row.names = FALSE)
    cat(sprintf("Successfully saved expression data to %s\n", expr_out))
} else {
    cat("Available molecular profiles:\n")
    print(names(mae))
    stop("Could not find any suitable mRNA expression profile in the downloaded MAE.")
}

cat("\n✅ Data download complete. You can now run the Python script:\n")
cat(sprintf("python scripts/generate_ml_prognostic_classifier_notebook.py --database %s --cancer %s --study-id %s\n", database, cancer_prefix, study_id))
