library(cBioPortalData)
library(ExperimentHub)
eh <- ExperimentHub(localHub=TRUE)
mae <- cBioDataPack("brca_aurora_2023", ask=FALSE, localHub=TRUE)
expr_se <- mae[["mrna_seq_v2_rsem"]]
write.csv(as.data.frame(rowData(expr_se)), "aurora_rowdata.csv")
