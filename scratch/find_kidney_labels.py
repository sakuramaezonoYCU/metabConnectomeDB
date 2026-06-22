import cellxgene_census
import pandas as pd

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    obs = census["census_data"]["homo_sapiens"].obs.read(column_names=["disease", "tissue", "tissue_general"]).concat().to_pandas()
    
    kidney_data = obs[obs["tissue_general"].str.lower().str.contains("kidney", na=False) | obs["tissue"].str.lower().str.contains("kidney", na=False)]
    print("Top tissues related to kidney:")
    print(kidney_data["tissue"].value_counts().head(10))
    print("\nTop diseases related to kidney:")
    print(kidney_data["disease"].value_counts().head(20))
