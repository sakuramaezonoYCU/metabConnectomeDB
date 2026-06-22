import subprocess
import cellxgene_census
import sys
import os

# Set default CAP for CellxGene downloads (number of cells). 100k = 100000 cells.
# Users can override by setting the environment variable CELLXGENE_CAP before running this script.
DEFAULT_CAP = "500000"
os.environ.setdefault('CELLXGENE_CAP', DEFAULT_CAP)

import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import CANCERS_TO_RUN, CANCER_PRIMARY_TISSUE, CANCER_DISEASE_QUERIES

print("Querying CellxGene Census for the Top 3 Metastatic Tissues...")

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    # Query is_primary_data to avoid duplicate/redundant cells, and use tissue_general
    df = cellxgene_census.get_obs(census, "homo_sapiens", column_names=["disease", "tissue_general", "is_primary_data"])
    df = df[df["is_primary_data"] == True]

# For testing, we will only run the first cancer
test_mode = False
if test_mode:
    CANCERS_TO_RUN = [CANCERS_TO_RUN[0]]
    print(f"TEST MODE: Only running pipeline for {CANCERS_TO_RUN[0]}")

import time

processes = []
exit_codes = []
active_processes = []
max_workers = 3 # Safe limit for 64GB RAM (each cancer can peak at 15-20GB)
script_path = os.path.join("scripts", "run_cancer_pipeline.py")

for cancer_key in CANCERS_TO_RUN:
    primary_tissues = CANCER_PRIMARY_TISSUE[cancer_key]
    disease_query = CANCER_DISEASE_QUERIES[cancer_key]
    
    # Filter df for this cancer (supports multiple comma-separated diseases or logical splits based on '||')
    diseases = [d.strip() for d in disease_query.replace('||', ',').split(',')]
    df_cancer = df[df["disease"].isin(diseases)]
    
    # Filter OUT primary tissues to find metastasis
    df_meta = df_cancer[~df_cancer["tissue_general"].isin(primary_tissues)]
    
    print(f"\n--- {cancer_key.upper()} ---")
    output_dir = os.path.join("output", f"{cancer_key}_results")
    if os.path.exists(output_dir) and os.listdir(output_dir):
        print(f"Output directory '{output_dir}' already exists. Re-running to recover missing HTML reports.")

    # Get top 3 metastatic tissues by cell count
    top_3_meta = df_meta["tissue_general"].value_counts().head(3).index.tolist()
    print(f"Primary Tissues: {primary_tissues}")
    print(f"Top 3 Metastatic Tissues: {top_3_meta}")
    
    # Combine them for the pipeline
    all_tissues = primary_tissues + top_3_meta
    tissue_filter_str = ",".join(all_tissues)
    primary_tissue_str = ",".join(primary_tissues)
    
    # Wait if we hit max workers
    while len(active_processes) >= max_workers:
        for i in range(len(active_processes)-1, -1, -1):
            p, c_key, log_f = active_processes[i]
            ret = p.poll()
            if ret is not None:
                exit_codes.append(ret)
                print(f"\n[{c_key}] Finished with exit code {ret}.")
                log_f.close()
                active_processes.pop(i)
        if len(active_processes) >= max_workers:
            time.sleep(5)

    print(f"Launching parallel pipeline for: {cancer_key} (Output redirected to output/{cancer_key}_pipeline.log)...")
    cmd = [sys.executable, script_path, disease_query, tissue_filter_str, primary_tissue_str]
    env = os.environ.copy()
    env.setdefault('CELLXGENE_CAP', DEFAULT_CAP)
    
    # Redirect output to a log file to avoid terminal garble
    log_file = open(f"output/{cancer_key}_pipeline.log", "w")
    p = subprocess.Popen(cmd, env=env, stdout=log_file, stderr=subprocess.STDOUT)
    active_processes.append((p, cancer_key, log_file))
    
# Wait for remaining processes
for p, c_key, log_f in active_processes:
    ret = p.wait()
    exit_codes.append(ret)
    print(f"\n[{c_key}] Finished with exit code {ret}.")
    log_f.close()

if all(code == 0 for code in exit_codes):
    print("\n✅ All cancer pipelines completed successfully.")
else:
    print("\n❌ One or more cancer pipelines failed.")
    sys.exit(1)
