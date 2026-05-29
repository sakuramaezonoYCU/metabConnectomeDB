import cellxgene_census
import pandas as pd

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    print("1. Querying exactly as run_all_cancers.py does (no is_primary_data filter, using 'tissue'):")
    df = cellxgene_census.get_obs(census, "homo_sapiens", column_names=["disease", "tissue", "is_primary_data", "tissue_general"])
    
    df_cancer = df[df["disease"] == "colorectal cancer"]
    df_meta_run = df_cancer[~df_cancer["tissue"].isin(["colon", "large intestine"])]
    print("\nrun_all_cancers.py Top 'tissue' counts:")
    print(df_meta_run["tissue"].value_counts().head(10))
    
    print("\n2. Querying as the Notebook does (is_primary_data == True, using 'tissue_general'):")
    df_cancer_nb = df[(df["disease"] == "colorectal cancer") & (df["is_primary_data"] == True)]
    df_meta_nb = df_cancer_nb[~df_cancer_nb["tissue_general"].isin(["colon", "large intestine"])]
    print("\nNotebook Top 'tissue_general' counts:")
    print(df_meta_nb["tissue_general"].value_counts().head(10))
