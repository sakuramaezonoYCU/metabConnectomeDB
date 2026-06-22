import os
import sys
import nbformat as nbf
import subprocess
import argparse


# Example case use:


def main():
    parser = argparse.ArgumentParser(
        description="Generate ML Prognostic Classifier Notebook",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with the 12-gene STAT3 Core Axis (default):
  python generate_ml_prognostic_classifier_notebook.py

  # Run with the 21-Gene Directed Metastatic Signature:
  python scripts/generate_ml_prognostic_classifier_notebook.py \
  --database aurora \
  --cancer brca \
  --study-id brca_aurora_2023 \
  --signature-name "Directed Metastatic Signature" \
  --genes GLS SGMS1 SPTLC1 GBE1 SLC16A7 AUH FZD6 NR1D2 CD46 MTMR1 ESRRG ITGA4 SLC11A2 ERAP1 C1GALT1 ADAM10 TRPM8 SLC22A1 AMDHD1 EPOR PDE3B

"""
    )
    parser.add_argument("--genes", nargs="+", 
                        default=['ADAM10', 'C1GALT1', 'ESRRG', 'FZD6', 'GBE1', 'GLS', 'ITGA4', 'PDE3B', 'SGMS1', 'SLC11A2', 'SLC16A7', 'SLC22A1'], 
                        help="List of genes for the signature")
    parser.add_argument("--signature-name", default="STAT3 Core Axis", help="Name of the signature")
    parser.add_argument("--database", default="metabric", help="Data source to use (e.g. metabric, broad, tcga)")
    parser.add_argument("--cancer", default="brca", help="Cancer prefix (e.g., brca, luad)")
    parser.add_argument("--study-id", default="brca_metabric", help="cBioPortal study ID (for fallback downloading)")
    parser.add_argument("--profile-expr", default="mrna_illumina_microarray", help="cBioPortal mRNA profile name")
    args = parser.parse_args()
    
    genes_list = args.genes
    signature_name = args.signature_name
    database = args.database
    cancer = args.cancer
    study_id = args.study_id
    profile_expr = args.profile_expr
    num_genes = len(genes_list)
    genes_str = ", ".join(genes_list)

    # 1. Create a new notebook
    nb = nbf.v4.new_notebook()
    
    # Notebook metadata
    nb['metadata'] = {
        'kernelspec': {
            'display_name': 'Python 3',
            'language': 'python',
            'name': 'python3'
        },
        'language_info': {
            'name': 'python'
        }
    }
    
    cells = []
    
    # --- Cell 1: Title and Purpose ---
    cell_1_md = f"""\
# ML Prognostic Classifier using the {num_genes}-gene {signature_name} ({database.upper()} - {cancer.upper()})

**Purpose:** We have identified a {num_genes}-gene {signature_name} signature: `{genes_str}`. This notebook builds and validates a Machine Learning Prognostic Classifier using this signature on the {cancer.upper()} cohort from {database.upper()}.

**Methodology:**
1. Load {database.upper()} clinical and mRNA expression data.
2. Train baseline Cox Proportional Hazards, Random Forest, and MLP Neural Network models.
3. Evaluate models using C-index, ROC-AUC, Kaplan-Meier curves, and feature importance.
"""
    cells.append(nbf.v4.new_markdown_cell(cell_1_md))

    # --- Cell 2: Imports and Setup ---
    cells.append(nbf.v4.new_code_cell("""\
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import joblib

# ML & Survival Modeling
from lifelines import CoxPHFitter, KaplanMeierFitter
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, roc_curve

# If scikit-survival is installed, try importing RandomSurvivalForest
try:
    from sksurv.ensemble import RandomSurvivalForest
    SKSURV_AVAILABLE = True
except ImportError:
    SKSURV_AVAILABLE = False
    print("scikit-survival is not installed. Will fallback to standard Random Forest Classifier for binary 5-year OS status.")

import warnings
warnings.filterwarnings('ignore')

# Set aesthetic style
sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 300
"""))

    # --- Cell 3: Data Loading & Preprocessing ---
    cells.append(nbf.v4.new_markdown_cell(f"""\
## 1. Data Fetching and Preprocessing

**Methodology:** We load the {database.upper()} clinical and expression data.

**Interpretation Guide:** Patients with missing survival times or status are dropped. We convert `OS_STATUS` to a boolean event indicator (1 = Death, 0 = Censored).
"""))

    cell_3_code = f"""\
genes = {repr(genes_list)}
database = "{database}"
cancer = "{cancer}"
study_id = "{study_id}"
profile_expr = "{profile_expr}"

if database.lower() == "tcga":
    data_dir = '../input/TCGA'
    clin_path = os.path.join(data_dir, f'TCGA-{{cancer.upper()}}.survival.tsv.gz')
    expr_path = os.path.join(data_dir, f'TCGA-{{cancer.upper()}}.star_fpkm.tsv.gz')
    
    if not os.path.exists(clin_path) or not os.path.exists(expr_path):
        raise FileNotFoundError(f"TCGA data files not found in {{data_dir}}. Expected {{clin_path}} and {{expr_path}}")
        
    clin_df = pd.read_csv(clin_path, sep='\\t')
    # clinical preprocessing: sample -> PATIENT_ID, OS.time -> OS_MONTHS, OS -> event
    if '_PATIENT' in clin_df.columns:
        clin_df.rename(columns={{'_PATIENT': 'PATIENT_ID'}}, inplace=True)
    elif 'sample' in clin_df.columns:
        clin_df.rename(columns={{'sample': 'PATIENT_ID'}}, inplace=True)
        
    if 'OS.time' in clin_df.columns:
        clin_df.rename(columns={{'OS.time': 'OS_MONTHS'}}, inplace=True)
        # TCGA OS.time is often in days
        clin_df['OS_MONTHS'] = clin_df['OS_MONTHS'] / 30.44
        
    if 'OS' in clin_df.columns:
        clin_df.rename(columns={{'OS': 'event'}}, inplace=True)
        
    # We don't have OS_STATUS string, just event 0/1
    clin_df['OS_STATUS'] = clin_df['event'].apply(lambda x: 'DECEASED' if x == 1 else 'LIVING')
    clin_df = clin_df[['PATIENT_ID', 'OS_MONTHS', 'OS_STATUS', 'event']].dropna()

    expr_df = pd.read_csv(expr_path, sep='\\t')
    # TCGA expression uses Ensembl_ID. We need to map Ensembl IDs to Hugo_Symbol.
    if 'Ensembl_ID' in expr_df.columns:
        expr_df['Ensembl_ID_clean'] = expr_df['Ensembl_ID'].apply(lambda x: str(x).split('.')[0])
        
        # Load HGNC mapping
        import json
        hgnc_path = '../input/hgnc_approved_genes.json'
        if os.path.exists(hgnc_path):
            with open(hgnc_path, 'r') as f:
                hgnc_data = json.load(f)
            
            ensembl_to_hugo = {{}}
            for doc in hgnc_data.get('response', {{}}).get('docs', []):
                if 'ensembl_gene_id' in doc and 'symbol' in doc:
                    ensembl_to_hugo[doc['ensembl_gene_id']] = doc['symbol']
                    
            expr_df['Hugo_Symbol'] = expr_df['Ensembl_ID_clean'].map(ensembl_to_hugo)
            # Fill unmapped with Ensembl ID
            expr_df['Hugo_Symbol'] = expr_df['Hugo_Symbol'].fillna(expr_df['Ensembl_ID_clean'])
        else:
            print("Warning: hgnc_approved_genes.json not found. Using Ensembl IDs as symbols.")
            expr_df['Hugo_Symbol'] = expr_df['Ensembl_ID_clean']
            
        expr_df = expr_df.drop(columns=['Ensembl_ID', 'Ensembl_ID_clean']).set_index('Hugo_Symbol').T
    else:
        expr_df = expr_df.set_index(expr_df.columns[0]).T

    expr_df.index.name = 'PATIENT_ID'
    expr_df = expr_df.reset_index()
    # TCGA samples often have -01A (sample type) appended, but clinical might be up to patient ID.
    # Align IDs: take first 12 chars for patient ID (e.g. TCGA-BH-A0BD)
    expr_df['PATIENT_ID'] = expr_df['PATIENT_ID'].apply(lambda x: str(x)[:12])
    clin_df['PATIENT_ID'] = clin_df['PATIENT_ID'].apply(lambda x: str(x)[:12])

else:
    data_dir_from_root = f'input/{{database}}'
    data_dir_from_scripts = f'../input/{{database}}'
    
    if os.path.exists('input'):
        data_dir = data_dir_from_root
    else:
        data_dir = data_dir_from_scripts
        
    os.makedirs(data_dir, exist_ok=True)
    
    # Check output directory as well for later
    if os.path.exists('output'):
        out_base = 'output'
    else:
        out_base = '../output'

    clin_url = f'https://raw.githubusercontent.com/cBioPortal/datahub/master/public/{{study_id}}/data_clinical_patient.txt'
    clin_path = os.path.join(data_dir, f'{{cancer}}_clinical.csv')
    expr_path = os.path.join(data_dir, f'{{cancer}}_expression.csv')

    # Download data if needed
    try:
        if not os.path.exists(clin_path) or not os.path.exists(expr_path) or os.path.getsize(expr_path) < 1000:
            import urllib.request
            import tarfile
            
            print(f"\\nDownloading {{study_id}} tarball from cBioPortal AWS Datahub...")
            tar_url = f"https://cbioportal-datahub.s3.amazonaws.com/{{study_id}}.tar.gz"
            tar_path = os.path.join(data_dir, f"{{study_id}}.tar.gz")
            urllib.request.urlretrieve(tar_url, tar_path)
            
            print("Extracting clinical and expression data...")
            with tarfile.open(tar_path, "r:gz") as tar:
                # Extract Clinical Patient
                clin_member = [m for m in tar.getmembers() if 'data_clinical_patient.txt' in m.name][0]
                clin_df = pd.read_csv(tar.extractfile(clin_member), sep='\\t', skiprows=4)
                
                # Extract Clinical Sample (for SAMPLE_ID mapping)
                sample_members = [m for m in tar.getmembers() if 'data_clinical_sample.txt' in m.name]
                if sample_members:
                    sample_df = pd.read_csv(tar.extractfile(sample_members[0]), sep='\\t', skiprows=4)
                    clin_df = pd.merge(clin_df, sample_df, on='PATIENT_ID', how='inner')
                clin_df.to_csv(clin_path, index=False)
                
                # Extract Expression
                expr_members = [m for m in tar.getmembers() if f'data_{{profile_expr}}.txt' in m.name]
                if not expr_members:
                    # Fallback to any mrna file if the specific profile is not found
                    expr_members = [m for m in tar.getmembers() if 'mrna' in m.name.lower() and m.name.endswith('.txt')]
                if not expr_members:
                    raise ValueError(f"Could not find {{profile_expr}} or any mRNA expression file in the tarball.")
                
                expr_df = pd.read_csv(tar.extractfile(expr_members[0]), sep='\\t')
                expr_df.to_csv(expr_path, index=False)
                
            os.remove(tar_path)
            print("Download and extraction complete.")
        else:
            clin_df = pd.read_csv(clin_path)
            expr_df = pd.read_csv(expr_path)
    except Exception as e:
        print(f"\\n🚨 ERROR: Failed to download {{database}} data automatically.")
        print(f"Exception details: {{e}}")
        raise e

    # Preprocessing Clinical
    # Try to find Survival Time column
    time_col = None
    for col in ['OS_MONTHS', 'OS_DAYS', 'PFS_MONTHS', 'DFS_MONTHS']:
        if col in clin_df.columns:
            time_col = col
            break
            
    # Try to find Survival Status column
    status_col = None
    for col in ['OS_STATUS', 'VITAL_STATUS', 'DFS_STATUS', 'PFS_STATUS']:
        if col in clin_df.columns:
            status_col = col
            break
            
    if time_col is None or status_col is None:
        raise ValueError(f"Could not find required survival columns in the clinical data for {{database}}. (Needs a Time and a Status column). Available columns: {{clin_df.columns.tolist()}}")

    clin_cols = ['PATIENT_ID', time_col, status_col]
    if 'SAMPLE_ID' in clin_df.columns:
        clin_cols.append('SAMPLE_ID')
        
    clin_df = clin_df[clin_cols].dropna(subset=[time_col, status_col])
    
    # Standardize to OS_MONTHS
    if 'DAYS' in time_col:
        clin_df['OS_MONTHS'] = clin_df[time_col] / 30.44
    else:
        clin_df['OS_MONTHS'] = clin_df[time_col]
        
    # Standardize to event (1 = death/progression, 0 = alive/censored)
    def parse_event(x):
        val = str(x).upper()
        if any(keyword in val for keyword in ['DECEASED', 'PROGRESSION', 'RECURRENCE']):
            return 1
        if val in ['1', '1.0', 'TRUE', 'YES']:
            return 1
        return 0
        
    clin_df['event'] = clin_df[status_col].apply(parse_event)

    # Preprocessing Expression
    expr_df = expr_df.drop(columns=['Entrez_Gene_Id'], errors='ignore').set_index('Hugo_Symbol').T
    
    merge_key = 'SAMPLE_ID' if 'SAMPLE_ID' in clin_df.columns else 'PATIENT_ID'
    expr_df.index.name = merge_key
    expr_df = expr_df.reset_index()

# Merge
df = pd.merge(clin_df, expr_df, on=merge_key, how='inner')
# Drop duplicate patients if any
df = df.drop_duplicates(subset=['PATIENT_ID'])
print(f"Final dataset shape: {{df.shape}}")
if df.shape[0] == 0:
    raise ValueError("After merging clinical and expression data, the dataset is empty. Check if PATIENT_ID/SAMPLE_ID match between files.")
df.head()
"""
    cells.append(nbf.v4.new_code_cell(cell_3_code))

    # --- Cell 4: Train / Test Split ---
    cells.append(nbf.v4.new_code_cell("""\
# Features & Target
# Ensure all signature genes are present (some might be missing in METABRIC array, check first)
present_genes = [g for g in genes if g in df.columns]
print(f"Found {len(present_genes)} out of {len(genes)} signature genes in expression data.")

X = df[present_genes]
# Target for Cox (time, event)
y_surv = df[['OS_MONTHS', 'event']]
# Target for Classifier: 5-year survival status (60 months)
# 1 = lived > 60 months, 0 = died <= 60 months
df['5yr_survival'] = np.where((df['OS_MONTHS'] >= 60), 1, 
                              np.where((df['event'] == 1) & (df['OS_MONTHS'] < 60), 0, np.nan))
                              
# We will use train_test_split on the main dataset
X_train, X_test, y_train, y_test = train_test_split(df[present_genes], df, test_size=0.2, random_state=42)

out_dir = f'{out_base}/ml_prognostic_results/{database}/{cancer}/' + str(num_genes) + '-gene'
os.makedirs(out_dir, exist_ok=True)
"""))

    # --- Cell 5: Cox Proportional Hazards Model ---
    cells.append(nbf.v4.new_markdown_cell(f"""\
## 2. Cox Proportional Hazards Model (Baseline)

**Methodology:** A multivariate Cox model to assess the linear risk effects of the gene signature.

**Interpretation Guide:** Hazard Ratios (HR) > 1 imply higher expression increases risk (worse prognosis). HR < 1 implies protective effect.
"""))
    cell_5_code = f"""\
cph_data = pd.concat([X_train, y_train[['OS_MONTHS', 'event']]], axis=1)
cph = CoxPHFitter(penalizer=0.1) # Add small penalization for stability
cph.fit(cph_data, duration_col='OS_MONTHS', event_col='event')
cph.print_summary()

plt.figure(figsize=(8, 6))
cph.plot()
plt.title(f'CoxPH Hazard Ratios for {num_genes}-gene Signature')
plt.tight_layout()
plt.savefig(f"{{out_dir}}/cox_hazard_ratios.png")
plt.show()

# Predictions (Hazard)
cph_risk_train = cph.predict_partial_hazard(X_train)
cph_risk_test = cph.predict_partial_hazard(X_test)

# Concordance index
print(f"CoxPH Test Concordance Index: {{cph.score(pd.concat([X_test, y_test[['OS_MONTHS', 'event']]], axis=1), scoring_method='concordance_index'):.3f}}")
"""
    cells.append(nbf.v4.new_code_cell(cell_5_code))

    # --- Cell 6: Random Forest & Neural Net (MLP) ---
    cells.append(nbf.v4.new_markdown_cell("""\
## 3. Machine Learning Models (Random Forest & Neural Network)

**Methodology:** Since standard classifiers don't natively handle censoring, we train them to predict 5-year overall survival (binary classification). Patients censored before 5 years are excluded from this specific training subset.

**Interpretation:** These non-linear models capture complex epistatic gene interactions that linear Cox models might miss.
"""))
    cells.append(nbf.v4.new_code_cell("""\
# Filter valid cases for binary classification (5-year survival)
train_bin = y_train.dropna(subset=['5yr_survival'])
test_bin = y_test.dropna(subset=['5yr_survival'])
X_train_bin = X_train.loc[train_bin.index]
X_test_bin = X_test.loc[test_bin.index]
y_train_bin = train_bin['5yr_survival']
y_test_bin = test_bin['5yr_survival']

# Random Forest
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_train_bin, y_train_bin)

# Neural Net (MLP)
mlp = MLPClassifier(hidden_layer_sizes=(32, 16), max_iter=500, random_state=42)
mlp.fit(X_train_bin, y_train_bin)

joblib.dump(rf, f"{out_dir}/rf_model.pkl")
joblib.dump(mlp, f"{out_dir}/mlp_model.pkl")

# Generate Risk scores (we use 1 - prob(survival) as risk score)
rf_risk_test = 1 - rf.predict_proba(X_test_bin)[:, 1]
mlp_risk_test = 1 - mlp.predict_proba(X_test_bin)[:, 1]

# Plot ROC Curves
plt.figure(figsize=(8, 6))

fpr_rf, tpr_rf, _ = roc_curve(1 - y_test_bin, rf_risk_test) # predicting 'event' (death) within 5yr
auc_rf = roc_auc_score(1 - y_test_bin, rf_risk_test)

fpr_mlp, tpr_mlp, _ = roc_curve(1 - y_test_bin, mlp_risk_test)
auc_mlp = roc_auc_score(1 - y_test_bin, mlp_risk_test)

plt.plot(fpr_rf, tpr_rf, label=f'Random Forest (AUC = {auc_rf:.3f})')
plt.plot(fpr_mlp, tpr_mlp, label=f'MLP Neural Net (AUC = {auc_mlp:.3f})')
plt.plot([0,1], [0,1], 'k--', label='Random Chance')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve: 5-Year Survival Prediction')
plt.legend()
plt.tight_layout()
plt.savefig(f"{out_dir}/roc_curves.png")
plt.show()
"""))

    # --- Cell 7: Feature Importance ---
    cells.append(nbf.v4.new_markdown_cell("""\
## 4. Feature Importance

**Methodology:** Extract Random Forest Gini importance to see which genes are the strongest predictive drivers.
"""))
    cell_7_code = f"""\
fi = pd.Series(rf.feature_importances_, index=present_genes).sort_values(ascending=False)

plt.figure(figsize=(10, 5))
sns.barplot(x=fi.values, y=fi.index, palette='viridis')
plt.title(f'Random Forest Feature Importance ({num_genes}-gene signature)')
plt.xlabel('Importance')
plt.tight_layout()
plt.savefig(f"{{out_dir}}/rf_feature_importance.png")
plt.show()
"""
    cells.append(nbf.v4.new_code_cell(cell_7_code))

    # --- Cell 8: Kaplan-Meier Risk Stratification ---
    cells.append(nbf.v4.new_markdown_cell("""\
## 5. Kaplan-Meier Risk Stratification

**Methodology:** We stratify the test cohort into 'High Risk' and 'Low Risk' groups based on the median risk score predicted by the baseline Cox model. We then plot Kaplan-Meier survival curves.

**Interpretation:** A significant separation (visual gap, log-rank test) demonstrates the clinical utility of the prognostic classifier.
"""))
    cell_8_code = f"""\
# Stratify based on median risk in the test set (using Cox model)
median_risk = np.median(cph_risk_test)
high_risk_mask = cph_risk_test > median_risk

test_high = y_test[high_risk_mask]
test_low = y_test[~high_risk_mask]

kmf = KaplanMeierFitter()
plt.figure(figsize=(8, 6))

kmf.fit(test_low['OS_MONTHS'], event_observed=test_low['event'], label=f'Low Risk (n={{len(test_low)}})')
ax = kmf.plot()

kmf.fit(test_high['OS_MONTHS'], event_observed=test_high['event'], label=f'High Risk (n={{len(test_high)}})')
kmf.plot(ax=ax)

# Try importing logrank test
try:
    from lifelines.statistics import multivariate_logrank_test
    res = multivariate_logrank_test(y_test['OS_MONTHS'], high_risk_mask, y_test['event'])
    plt.text(0.05, 0.05, f"Log-rank p-value: {{res.p_value:.2e}}", transform=ax.transAxes, fontsize=12)
except Exception as e:
    pass

plt.title(f'Kaplan-Meier Survival Curve: High vs Low Risk ({num_genes}-gene Signature)')
plt.xlabel('Time (Months)')
plt.ylabel('Survival Probability')
plt.tight_layout()
plt.savefig(f"{{out_dir}}/km_survival_curve.png")
plt.show()
"""
    cells.append(nbf.v4.new_code_cell(cell_8_code))

    # --- Cell 9: Auto-HTML Export ---
    cell_9_code = f"""\\
import subprocess
import sys
import os

try:
    from pan_cancer_config import ANALYSIS_SUFFIX
except ImportError:
    ANALYSIS_SUFFIX = ''

notebook_filename = 'ml_prognostic_classifier.ipynb'
output_base = 'ml_prognostic_classifier_report' + ANALYSIS_SUFFIX

jupyter_bin = os.path.join(os.path.dirname(sys.executable), 'jupyter')
if not os.path.exists(jupyter_bin): jupyter_bin = 'jupyter'

cmd_html = [jupyter_bin, "nbconvert", "--to", "html", notebook_filename, "--output-dir", out_dir, "--output", output_base]
res_html = subprocess.run(cmd_html, capture_output=True, text=True)

if res_html.returncode == 0:
    print(f"🎉 SUCCESS: Notebook successfully exported to '{{os.path.join(out_dir, output_base)}}.html'")
else:
    print("❌ HTML export failed.")
    print(res_html.stderr)
"""
    cells.append(nbf.v4.new_code_cell(cell_9_code))

    nb['cells'] = cells
    
    # Save the notebook
    notebook_filename = 'scripts/ml_prognostic_classifier.ipynb'
    with open(notebook_filename, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)
    
    print(f"Notebook successfully created at: {notebook_filename}")
if __name__ == "__main__":
    main()
