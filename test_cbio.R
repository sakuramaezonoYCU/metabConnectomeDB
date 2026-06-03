suppressPackageStartupMessages(library(cBioPortalData))
cbio <- cBioPortal()

clin <- clinicalData(cbio, "brca_metabric")
samp_ids <- clin$sampleId
if(is.null(samp_ids)) samp_ids <- clin$patientId

cat(sprintf("Found %d samples\n", length(samp_ids)))
cat("Trying to fetch molecular data without entrezGeneIds...\n")

res <- tryCatch({
  molecularData(cbio, molecularProfileIds = "brca_metabric_mrna", sampleIds = samp_ids[1:2])
}, error = function(e) {
  cat("Error without entrezGeneIds:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(res)) {
  cat(sprintf("Rows returned: %d\n", nrow(res)))
}

cat("Trying cBioPortalData with a large number of genes...\n")
gt <- geneTable(cbio)
all_genes <- gt$hugoGeneSymbol[1:100] # just 100 genes for testing
cat(sprintf("Fetching %d genes using cBioPortalData...\n", length(all_genes)))
res2 <- tryCatch({
  cBioPortalData(cbio, studyId = "brca_metabric", genes = all_genes, ask=FALSE)
}, error = function(e) {
  cat("Error with cBioPortalData:", conditionMessage(e), "\n")
})
if (inherits(res2, "MultiAssayExperiment")) {
  cat("Successfully fetched MAE with genes subset!\n")
}
