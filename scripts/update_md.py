import pandas as pd
import re
import sys
import os

if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX

md_file = 'output/AI_summary_and_insights.md'
with open(md_file, 'r') as f:
    content = f.read()

# Helper function to convert csv to markdown table
def csv_to_md_table(csv_path):
    df = pd.read_csv(csv_path)
    headers = "| " + " | ".join(df.columns) + " |"
    separators = "| " + " | ".join([":---" if i==0 else ":---:" for i in range(len(df.columns))]) + " |"
    rows = []
    for _, row in df.iterrows():
        rows.append("| " + " | ".join([str(x) for x in row.values]) + " |")
    return "\n".join([headers, separators] + rows)

# 1. Conserved signature table
df_sig = pd.read_csv(f'output/ai_summary_tables/conserved_gene_directed_signature_annotation{ANALYSIS_SUFFIX}.csv')
table_md = f"""> [!NOTE]
> **Data Provenance**
> - **Source File:** `output/ai_summary_tables/conserved_gene_directed_signature_annotation{ANALYSIS_SUFFIX}.csv`
> - **Scripts:** `scripts/generate_ai_summary_tables.py`, `scripts/fetch_uniprot_roles.py`, `scripts/fetch_opentargets.py`, `scripts/update_md.py`

"""
table_md += "| Gene | Source_Database | Sensor_Type | Key Metabolite(s) | Top OpenTargets Diseases | MRCLinkDB_Disease | GEPIA Link |\n|:---|:---|:---|:---|:---|:---|:---|\n"
for _, row in df_sig.iterrows():
    gene = row['Gene']
    source_db = str(row['Source_Database']) if pd.notna(row.get('Source_Database')) else ""
    sensor = str(row['Sensor_Type']) if pd.notna(row['Sensor_Type']) else ""
    met = str(row['Key Metabolite(s)'])
    diseases = str(row['Top OpenTargets Diseases']) if pd.notna(row.get('Top OpenTargets Diseases')) else 'Not fetched'
    mrclink = str(row['MRCLinkDB_Disease']) if pd.notna(row.get('MRCLinkDB_Disease')) else ""
    link = str(row['GEPIA Link']) if pd.notna(row.get('GEPIA Link')) else ""
    if link: link = f"[GEPIA]({link})"
    table_md += f"| **{gene}** | {source_db} | {sensor} | {met} | {diseases} | {mrclink} | {link} |\n"

# Replace conserved gene table
content = re.sub(r'\| Gene \|.*?(?=\n\n)', table_md.strip(), content, count=1, flags=re.DOTALL)

# 2. Version 6/7 Dataset Overview
overview_500k = csv_to_md_table('output/ai_summary_tables/dataset_overview_500k.csv')
overview_md = overview_500k + """

**Column Definitions:**
- **Total Primary TME Cells**: All cells in the primary tumor microenvironment, excluding the actual malignant (cancer) cells. This includes immune cells, fibroblasts, and endothelial cells located at the site of the primary tumor.
- **Primary Malignant Cells**: The actual cancer cells located at the primary tumor site.
- **Total Metastatic TME Cells**: All cells in the metastatic tumor microenvironment, excluding the disseminated malignant cells. This includes the immune and stromal cells located at the distant metastatic site.
- **Metastatic Malignant Cells**: The disseminated cancer cells located at a distant metastatic site."""

# Replace the Expanded Dataset Overview table
content = re.sub(
    r'\| Dataset \| Total Primary.*?(?=\n\n### 1\. STAT3)',
    overview_md.strip(),
    content,
    count=1,
    flags=re.DOTALL
)

# 3. Version 6/7 Subclone Summary
subclone_defs = """
**Column Definitions:**
- **Primary Cells Scored**: The absolute number of malignant cells from the primary tumor that successfully passed quality control and received a conserved Metastatic Signature Score.
- **Score Distribution**: The shape of the score distribution across the primary tumor cells, computed mathematically using skewness (Left-skewed: < -0.5, Right-skewed: > +0.5, Symmetric otherwise).
- **Pre-Metastatic Subclone (%)**: The percentage of primary tumor cells whose signature score is greater than the mean plus one standard deviation (> +1 SD). This represents the highly metastatic "tail" or subclone already present in the primary tumor prior to dissemination."""

subclone_500k = csv_to_md_table('output/ai_summary_tables/subclone_summary_500k.csv')
subclone_500k_md = "#### Pre-Metastatic Subclone Resolution (Latest Version)\n\n" + subclone_500k + "\n\n" + subclone_defs

# Let's insert it right before the "### 1. STAT3" section if it's not already there
if "Pre-Metastatic Subclone Resolution (Latest Version)" not in content:
    content = re.sub(
        r'(\*\*Column Definitions:\*\*.*?\n\n)(### 1\. STAT3)',
        r'\g<1>' + subclone_500k_md + r'\n\n\g<2>',
        content,
        count=1,
        flags=re.DOTALL
    )
else:
    # update existing
    content = re.sub(
        r'#### Pre-Metastatic Subclone Resolution \(Latest Version\)\n\n\| Cancer \|.*?(?=\n\n### 1\. STAT3)',
        subclone_500k_md.strip(),
        content,
        count=1,
        flags=re.DOTALL
    )

# 4. Version 5 Q5 Table (100k)
subclone_100k = csv_to_md_table('output/ai_summary_tables/subclone_summary_100k.csv')
subclone_100k_md = subclone_100k + "\n\n" + subclone_defs

content = re.sub(
    r'(#### Q5: Can the .*?pan-cancer signature predict metastatic potential from primary tumor biopsies\?.*?)\| Cancer \| Primary Cells Scored.*?(?=\n\n- \*\*Lung cancer\*\*)',
    r'\g<1>' + subclone_100k_md.strip() + r'\n\n',
    content,
    count=1,
    flags=re.DOTALL
)

with open(md_file, 'w') as f:
    f.write(content)
print("Updated AI_summary_and_insights.md with all dynamic tables.")
