import pandas as pd
import re

md_file = 'output/AI_summary_and_insights.md'
with open(md_file, 'r') as f:
    content = f.read()

# 1. Replace 21-gene table (now with Sensor_Type, Diseases and GEPIA Link)
df_sig = pd.read_csv('output/ai_summary_tables/21_gene_directed_signature_annotation_Br500k_Co100k_Lu500k_Me100k_Ov100k.csv')

table_md = """> [!NOTE]
> **Data Provenance**
> - **Source File:** `output/ai_summary_tables/21_gene_directed_signature_annotation_Br500k_Co100k_Lu500k_Me100k_Ov100k.csv`
> - **Scripts:** `scripts/generate_ai_summary_tables.py`, `scripts/fetch_uniprot_roles.py`, `scripts/fetch_opentargets.py`, `scripts/update_md.py`

"""
table_md += "| Gene | Source_Database | Sensor_Type | Key Metabolite(s) | Top OpenTargets Diseases | MRCLinkDB_Disease | GEPIA Link |\n|:---|:---|:---|:---|:---|:---|:---|\n"
for _, row in df_sig.iterrows():
    gene = row['Gene']
    source_db = str(row['Source_Database']) if pd.notna(row.get('Source_Database')) else ""
    sensor = str(row['Sensor_Type']) if pd.notna(row['Sensor_Type']) else ""
    met = str(row['Key Metabolite(s)'])
    
    # We purposefully exclude Biological Role, Rhea_Reaction, and PMID(s) from the Markdown because it's too long/not useful
    diseases = str(row['Top OpenTargets Diseases']) if pd.notna(row.get('Top OpenTargets Diseases')) else 'Not fetched'
    mrclink = str(row['MRCLinkDB_Disease']) if pd.notna(row.get('MRCLinkDB_Disease')) else ""
    link = str(row['GEPIA Link']) if pd.notna(row.get('GEPIA Link')) else ""
    
    # Optional: format link as markdown link
    if link:
        link = f"[GEPIA]({link})"
        
    table_md += f"| **{gene}** | {source_db} | {sensor} | {met} | {diseases} | {mrclink} | {link} |\n"

# Replace the existing markdown table
# The original table has columns: Gene, Direction, Key Metabolite(s), Biological Role
# Or it has: Gene, Sensor_Type, Key Metabolite(s), Rhea_Reaction, PMID(s), Biological Role, GEPIA Link
# We want to replace it up to the first blank line.
content = re.sub(
    r'\| Gene \|.*?(?=\n\n)',
    table_md.strip(),
    content,
    count=1,
    flags=re.DOTALL
)

with open(md_file, 'w') as f:
    f.write(content)

print("Updated AI_summary_and_insights.md with Sensor_Type table.")
