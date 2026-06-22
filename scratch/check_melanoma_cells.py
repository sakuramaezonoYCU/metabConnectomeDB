import cellxgene_census

print("Querying CellxGene Census for Melanoma...")

diseases = ["melanoma", "metastatic melanoma", "skin melanoma", "cutaneous melanoma", "malignant melanoma"]

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    df = cellxgene_census.get_obs(census, "homo_sapiens", column_names=["disease", "tissue_general", "is_primary_data"])
    df = df[df["is_primary_data"] == True]
    
    df_cancer = df[df["disease"].isin(diseases)]
    
    print(f"Total primary cells matching {diseases}: {len(df_cancer)}")
    
    tissues = df_cancer["tissue_general"].value_counts()
    print("\nCells per tissue:")
    print(tissues.head(10))
