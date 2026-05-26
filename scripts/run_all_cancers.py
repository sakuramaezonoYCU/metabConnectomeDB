import subprocess
import cellxgene_census
import sys
import os

# Set default CAP for CellxGene downloads (number of cells). 100k = 100000 cells.
# Users can override by setting the environment variable CELLXGENE_CAP before running this script.
DEFAULT_CAP = "100000"
os.environ.setdefault('CELLXGENE_CAP', DEFAULT_CAP)

cancer_to_primary = {
    'breast cancer': ['breast', 'mammary gland'],
    'colorectal cancer': ['colon', 'large intestine'],
    'melanoma': ['skin'],
    'lung adenocarcinoma': ['lung'],
    'ovarian cancer': ['ovary']
}

print("Querying CellxGene Census for the Top 3 Metastatic Tissues...")

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    df = cellxgene_census.get_obs(census, "homo_sapiens", column_names=["disease", "tissue"])

# For testing, we will only run the first cancer
test_mode = False
cancers_to_run = list(cancer_to_primary.keys())
if test_mode:
    cancers_to_run = [cancers_to_run[0]]
    print(f"TEST MODE: Only running pipeline for {cancers_to_run[0]}")

processes = []
script_path = os.path.join("scripts", "run_cancer_pipeline.py")

for cancer in cancers_to_run:
    primary_tissues = cancer_to_primary[cancer]
    
    # Filter df for this cancer
    df_cancer = df[df["disease"] == cancer]
    
    # Filter OUT primary tissues to find metastasis
    df_meta = df_cancer[~df_cancer["tissue"].isin(primary_tissues)]
    
    # Get Primary + top 3 metastatic tissues by cell count
    top_3_meta = df_meta["tissue"].value_counts().head(4).index.tolist()
    
    print(f"\n--- {cancer.upper()} ---")
    print(f"Primary Tissues: {primary_tissues}")
    print(f"Top 3 Metastatic Tissues: {top_3_meta}")
    
    # Combine them for the pipeline
    all_tissues = primary_tissues + top_3_meta
    tissue_filter_str = ",".join(all_tissues)
    primary_tissue_str = ",".join(primary_tissues)
    
    print(f"Launching headless pipeline for: {cancer}...")
    # Launch as a subprocess. We pass: script, disease, tissue_filter, primary_tissues
    cmd = [sys.executable, script_path, cancer, tissue_filter_str, primary_tissue_str]
    # Ensure CAP env var is passed
    env = os.environ.copy()
    env.setdefault('CELLXGENE_CAP', DEFAULT_CAP)
    p = subprocess.Popen(cmd, env=env)
    processes.append(p)

# Wait for all processes to finish
exit_codes = [p.wait() for p in processes]

if all(code == 0 for code in exit_codes):
    print("\n✅ All cancer pipelines completed successfully.")
else:
    print("\n❌ One or more cancer pipelines failed.")
    sys.exit(1)
