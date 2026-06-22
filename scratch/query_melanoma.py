import cellxgene_census

print("Loading CellxGene Census 2025-11-08...")
with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    df = cellxgene_census.get_obs(census, "homo_sapiens", column_names=["disease", "tissue_general", "is_primary_data"])
    
    # Check all diseases that contain 'melanoma'
    melanoma_df = df[df["disease"].str.contains("melanoma", case=False, na=False)]
    
    print("\n--- All Melanoma Cells (Primary + Duplicate) ---")
    print(melanoma_df["disease"].value_counts())
    
    print("\n--- Primary Only (is_primary_data == True) ---")
    primary_df = melanoma_df[melanoma_df["is_primary_data"] == True]
    print(primary_df["disease"].value_counts())
    
    print("\n--- Primary Tissues for 'melanoma' AND 'metastatic melanoma' ---")
    selected_df = primary_df[primary_df["disease"].isin(["melanoma", "metastatic melanoma"])]
    print(selected_df["tissue_general"].value_counts())
