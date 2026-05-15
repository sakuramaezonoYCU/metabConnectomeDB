import pandas as pd, numpy as np

df1 = pd.DataFrame({'Metabolite_Name': ['A', 'B'], 'HMDB_ID': [1, 2], 'Receptor': ['R1', 'R2'], 'PMID': ['111', '222']})
df2 = pd.DataFrame({'Metabolite_Name': ['A', 'B'], 'Concentration': [0.5, 0.6], 'PMID': [np.nan, '333']})

keys = [k for k in ['Metabolite_Name', 'HMDB_ID', 'Gene_Name', 'Receptor_Gene_Symbol'] if k in df1.columns and k in df2.columns]
res = pd.merge(df1, df2, on=keys, how='outer')

for c in res.columns:
    if c.endswith('_x'):
        orig = c[:-2]
        y = orig + '_y'
        if y in res.columns:
            res[orig] = res[c].fillna(res[y])
            res.drop(columns=[c, y], inplace=True)

print(res)
