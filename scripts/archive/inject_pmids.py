import os
import json
import pandas as pd
import sys
sys.path.append('scripts')
from pan_cancer_config import CANCERS_TO_RUN, CANCER_PO2_CSV_MAPPING

csv_path = "input/pO2_guide_24588669.csv"
po2_df = pd.read_csv(csv_path)

md_file = "pipeline_execution_checklist.md"
with open(md_file, "r") as f:
    content = f.read()

citations = []
for cancer in CANCERS_TO_RUN:
    csv_name = CANCER_PO2_CSV_MAPPING.get(cancer, "")
    if csv_name:
        row = po2_df[po2_df['Tumour type'] == csv_name]
        if not row.empty:
            pmid_str = str(row['Reference'].values[0]).replace('"', '')
            pmids = pmid_str.split(',')
            links = []
            for p in pmids:
                p = p.strip()
                if p and p != 'nan':
                    links.append(f"[PMID {p}](https://pubmed.ncbi.nlm.nih.gov/{p}/)")
            if links:
                citations.append(f"- **{cancer.capitalize()}**: {', '.join(links)}")

if citations:
    citation_text = "\n\n## References & PMIDs\n" + "\n".join(citations) + "\n"
    if "## References & PMIDs" not in content:
        with open(md_file, "a") as f:
            f.write(citation_text)
        print("Injected PMIDs into checklist.")
    else:
        print("PMIDs already present.")
