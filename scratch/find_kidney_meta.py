import cellxgene_census

diseases = ["nonpapillary renal cell carcinoma", "clear cell renal carcinoma", "chromophobe renal cell carcinoma", "Wilms tumor"]
primary_tissues = ["kidney", "cortex of kidney", "renal medulla", "renal papilla", "renal pelvis"]

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    df = cellxgene_census.get_obs(census, "homo_sapiens", column_names=["disease", "tissue_general", "is_primary_data"])
    df = df[df["is_primary_data"] == True]
    df_cancer = df[df["disease"].isin(diseases)]
    df_meta = df_cancer[~df_cancer["tissue_general"].isin(primary_tissues)]
    top_3_meta = df_meta["tissue_general"].value_counts().head(3).index.tolist()
    print("Top 3 metastatic tissues:", top_3_meta)
