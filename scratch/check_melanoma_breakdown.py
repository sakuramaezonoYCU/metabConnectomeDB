import cellxgene_census

print("Querying CellxGene Census for Melanoma breakdown...")

diseases = ["melanoma", "metastatic melanoma", "skin melanoma", "cutaneous melanoma", "malignant melanoma", "uveal melanoma"]

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    df = cellxgene_census.get_obs(census, "homo_sapiens", column_names=["disease", "is_primary_data"])
    df = df[df["is_primary_data"] == True]
    
    print("\nBreakdown of cell counts by disease keyword:")
    total = 0
    for disease in diseases:
        count = len(df[df["disease"] == disease])
        print(f"  {disease.ljust(25)} : {count} cells")
        total += count
    
    print(f"\nTotal: {total} cells")
