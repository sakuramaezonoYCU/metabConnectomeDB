import os
import sys
import json
import subprocess
from execute_pancancer_notebooks import execute_and_export

script_dir = os.path.dirname(os.path.abspath(__file__))
ml_gen_script = os.path.join(script_dir, "generate_ml_prognostic_classifier_notebook.py")
ml_nb = os.path.join(script_dir, "ml_prognostic_classifier.ipynb")

config_path = os.path.join(os.path.dirname(script_dir), "input", "pipeline.config.json")
tcga_mapping = {}
cancers_to_run = []
if os.path.exists(config_path):
    with open(config_path, "r") as f:
        cfg = json.load(f)
        tcga_mapping = cfg.get("PHASE_6_REPORTING", {}).get("TCGA_MAPPING", {})
        cancers_to_run = cfg.get("PHASE_2_SINGLE_CELL_INTEGRATION", {}).get("CANCERS_TO_RUN", [])

# 1. Per-Cancer Reports
for cancer_name in cancers_to_run:
    tcga_prefixes = tcga_mapping.get(cancer_name, [])
    if isinstance(tcga_prefixes, str):
        tcga_prefixes = [tcga_prefixes]
    
    valid_prefixes = []
    for prefix in tcga_prefixes:
        prefix = prefix.lower()
        expr_file = os.path.join(os.path.dirname(script_dir), "input", "TCGA", f"TCGA-{prefix.upper()}.star_fpkm.tsv.gz")
        if os.path.exists(expr_file):
            valid_prefixes.append(prefix)
            
    if not valid_prefixes:
        continue
        
    display_prefixes = ', '.join([p.upper() for p in valid_prefixes])
    print(f"Generating ML Classifier for {cancer_name.capitalize()} ({display_prefixes})...")
    subprocess.run([sys.executable, ml_gen_script, "--database", "tcga", "--cancer"] + valid_prefixes, check=True)
    
    ml_html = os.path.join(os.path.dirname(script_dir), "output", f"{cancer_name}_ml_prognostic_classifier_report.html")
    execute_and_export(ml_nb, ml_html, f"{cancer_name.capitalize()} ML Prognostic Classifier")

# 2. Pan-Cancer Report
print("Generating Pan-Cancer ML Classifier...")
subprocess.run([sys.executable, ml_gen_script, "--database", "tcga", "--cancer", "all"], check=True)
pancancer_html = os.path.join(os.path.dirname(script_dir), "output", "pancancer_ml_prognostic_classifier_report.html")
execute_and_export(ml_nb, pancancer_html, "Pan-Cancer ML Prognostic Classifier")
