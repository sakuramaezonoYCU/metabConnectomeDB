import cellxgene_census

with cellxgene_census.open_soma() as census:
    obs_df = cellxgene_census.get_obs(census, "homo_sapiens", column_names=["disease"])
df = obs_df.to_pandas()
cancer_terms = [d for d in df["disease"].unique() if "cancer" in d.lower() or "carcinoma" in d.lower() or "melanoma" in d.lower() or "leukemia" in d.lower()]
df_cancer = df[df["disease"].isin(cancer_terms)]
print("Top 5 Cancers by Cell Count in Census:")
print(df_cancer["disease"].value_counts().head(5))
