import os
import sys
import nbformat as nbf
import subprocess
import argparse
import pandas as pd
import json
import urllib.request
import tarfile

def standardize_clinical(clin_df, database):
    """
    Standardizes clinical dataframe to ensure PATIENT_ID, OS_MONTHS, and event exist.
    """
    if database.lower() == 'tcga':
        if '_PATIENT' in clin_df.columns:
            clin_df.rename(columns={'_PATIENT': 'PATIENT_ID'}, inplace=True)
        elif 'sample' in clin_df.columns:
            clin_df.rename(columns={'sample': 'PATIENT_ID'}, inplace=True)
        if 'OS.time' in clin_df.columns:
            clin_df.rename(columns={'OS.time': 'OS_MONTHS'}, inplace=True)
            clin_df['OS_MONTHS'] = clin_df['OS_MONTHS'] / 30.44
        if 'OS' in clin_df.columns:
            clin_df.rename(columns={'OS': 'event'}, inplace=True)
            
    # Try to find Survival Time column
    time_col = None
    for col in ['OS_MONTHS', 'OS_DAYS', 'PFS_MONTHS', 'DFS_MONTHS']:
        if col in clin_df.columns:
            time_col = col
            break
            
    # Try to find Survival Status column
    status_col = None
    for col in ['OS_STATUS', 'VITAL_STATUS', 'DFS_STATUS', 'PFS_STATUS', 'event']:
        if col in clin_df.columns:
            status_col = col
            break
            
    if time_col is None or status_col is None:
        raise ValueError(f"Could not find required survival columns in the clinical data. Available columns: {clin_df.columns.tolist()}")

    clin_cols = ['PATIENT_ID', time_col, status_col]
    if 'SAMPLE_ID' in clin_df.columns:
        clin_cols.append('SAMPLE_ID')
        
    # Drop rows missing time or status
    clin_df = clin_df.dropna(subset=[time_col, status_col]).copy()
    
    # Standardize to OS_MONTHS
    if 'DAYS' in time_col:
        clin_df['OS_MONTHS'] = clin_df[time_col] / 30.44
    else:
        clin_df['OS_MONTHS'] = clin_df[time_col]
        
    # Standardize to event (1 = death/progression, 0 = alive/censored)
    if status_col != 'event':
        def parse_event(x):
            val = str(x).upper()
            if any(keyword in val for keyword in ['DECEASED', 'PROGRESSION', 'RECURRENCE']):
                return 1
            if val in ['1', '1.0', 'TRUE', 'YES']:
                return 1
            return 0
            
        clin_df['event'] = clin_df[status_col].apply(parse_event)
        
    if database.lower() == 'tcga':
        clin_df['PATIENT_ID'] = clin_df['PATIENT_ID'].apply(lambda x: str(x)[:12])
        
    cols_to_keep = ['PATIENT_ID', 'OS_MONTHS', 'event']
    if 'SAMPLE_ID' in clin_df.columns:
        cols_to_keep.append('SAMPLE_ID')
        
    return clin_df[cols_to_keep]

def prepare_data(database, cancers, study_ids, profile_expr):
    """
    Downloads and prepares CSV files for the requested cancers.
    """
    if os.path.exists('input'):
        data_dir = f'input/{database}'
        tcga_dir = 'input/TCGA'
    else:
        data_dir = f'../input/{database}'
        tcga_dir = '../input/TCGA'
        
    os.makedirs(data_dir, exist_ok=True)
    
    for idx, cancer in enumerate(cancers):
        clin_path = os.path.join(data_dir, f'{cancer}_clinical.csv')
        expr_path = os.path.join(data_dir, f'{cancer}_expression.csv')
        
        if os.path.exists(clin_path) and os.path.exists(expr_path) and os.path.getsize(expr_path) > 1000:
            print(f"[{cancer.upper()}] Preprocessed data already exists.")
            continue
            
        print(f"[{cancer.upper()}] Preprocessing {database} data...")
        
        if database.lower() == 'tcga':
            tcga_clin = os.path.join(tcga_dir, f'TCGA-{cancer.upper()}.survival.tsv.gz')
            tcga_expr = os.path.join(tcga_dir, f'TCGA-{cancer.upper()}.star_fpkm.tsv.gz')
            
            if not os.path.exists(tcga_clin) or not os.path.exists(tcga_expr):
                raise FileNotFoundError(f"TCGA files missing. Expected: {tcga_clin} and {tcga_expr}")
            
            # Clinical
            clin_df = pd.read_csv(tcga_clin, sep='\t')
            clin_df = standardize_clinical(clin_df, 'tcga')
            clin_df.to_csv(clin_path, index=False)
            
            # Expression
            expr_df = pd.read_csv(tcga_expr, sep='\t')
            if 'Ensembl_ID' in expr_df.columns:
                expr_df['Ensembl_ID_clean'] = expr_df['Ensembl_ID'].apply(lambda x: str(x).split('.')[0])
                
                hgnc_path = os.path.join(tcga_dir, '..', 'hgnc_approved_genes.json')
                if os.path.exists(hgnc_path):
                    with open(hgnc_path, 'r') as f:
                        hgnc_data = json.load(f)
                    ensembl_to_hugo = {}
                    for doc in hgnc_data.get('response', {}).get('docs', []):
                        if 'ensembl_gene_id' in doc and 'symbol' in doc:
                            ensembl_to_hugo[doc['ensembl_gene_id']] = doc['symbol']
                    expr_df['Hugo_Symbol'] = expr_df['Ensembl_ID_clean'].map(ensembl_to_hugo)
                    expr_df['Hugo_Symbol'] = expr_df['Hugo_Symbol'].fillna(expr_df['Ensembl_ID_clean'])
                else:
                    expr_df['Hugo_Symbol'] = expr_df['Ensembl_ID_clean']
                    
                expr_df = expr_df.drop(columns=['Ensembl_ID', 'Ensembl_ID_clean'])
            
            # Set Hugo_Symbol as the first column for standard merging later
            if 'Hugo_Symbol' in expr_df.columns:
                cols = ['Hugo_Symbol'] + [c for c in expr_df.columns if c != 'Hugo_Symbol']
                expr_df = expr_df[cols]
                
            expr_df.to_csv(expr_path, index=False)
            
        else:
            # cBioPortal
            if idx >= len(study_ids):
                raise ValueError(f"Missing study_id for cancer: {cancer}. Provide matching study-ids.")
            study_id = study_ids[idx]
            
            tar_url = f"https://cbioportal-datahub.s3.amazonaws.com/{study_id}.tar.gz"
            tar_path_local = os.path.join(data_dir, f"{study_id}.tar.gz")
            
            print(f"Downloading {study_id} tarball from AWS Datahub...")
            urllib.request.urlretrieve(tar_url, tar_path_local)
            
            with tarfile.open(tar_path_local, "r:gz") as tar:
                # Clinical Patient
                clin_member = [m for m in tar.getmembers() if 'data_clinical_patient.txt' in m.name][0]
                clin_df = pd.read_csv(tar.extractfile(clin_member), sep='\t', skiprows=4)
                
                sample_members = [m for m in tar.getmembers() if 'data_clinical_sample.txt' in m.name]
                if sample_members:
                    sample_df = pd.read_csv(tar.extractfile(sample_members[0]), sep='\t', skiprows=4)
                    clin_df = pd.merge(clin_df, sample_df, on='PATIENT_ID', how='inner')
                
                clin_df = standardize_clinical(clin_df, 'cbioportal')
                clin_df.to_csv(clin_path, index=False)
                
                # Expression
                expr_members = [m for m in tar.getmembers() if f'data_{profile_expr}.txt' in m.name]
                if not expr_members:
                    expr_members = [m for m in tar.getmembers() if 'mrna' in m.name.lower() and m.name.endswith('.txt')]
                if not expr_members:
                    raise ValueError(f"No mRNA expression file found in {study_id} tarball.")
                
                expr_df = pd.read_csv(tar.extractfile(expr_members[0]), sep='\t')
                expr_df.to_csv(expr_path, index=False)
                
            os.remove(tar_path_local)

def generate_notebook(database, cancers, study_ids, profile_expr, genes_list, signature_name, genes_str, num_genes):
    nb = nbf.v4.new_notebook()
    
    nb['metadata'] = {
        'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
        'language_info': {'name': 'python'}
    }
    
    cells = []
    
    cancers_str = ", ".join([c.upper() for c in cancers])
    
    cells.append(nbf.v4.new_markdown_cell(f"""\
# ML Prognostic Classifier using the {num_genes}-gene {signature_name} ({database.upper()} - {cancers_str})

### Goal
Build and validate a Machine Learning Prognostic Classifier using the `{genes_str}` signature across {cancers_str} cohorts from {database.upper()}.

### Purpose
To determine whether the multi-gene metabolic signature possesses independent non-linear prognostic utility for predicting 5-year overall survival, providing translational insight beyond traditional linear survival models.

### Interpretation
- **Cox Proportional Hazards C-index**: Establishes baseline linear survival predictive power.
- **Random Forest & MLP ROC-AUC**: Evaluates non-linear and interactive predictive power for 5-year binary overall survival.
- **Kaplan-Meier High vs Low Risk**: Stratifies patients based on machine learning risk scores to demonstrate clinical significance.

### Inputs/Parameters
- **Database:** `{database.upper()}`
- **Cancers:** `{cancers_str}`
- **Gene Signature:** `{signature_name}` (`{genes_str}`)

### Outputs
- **Plots & Models:** Saved to `output/ml_prognostic_results/{database}/{"_".join(cancers)}/{num_genes}-gene/`
"""))

    cells.append(nbf.v4.new_code_cell("""\
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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

sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 300
"""))

    cells.append(nbf.v4.new_markdown_cell(f"""\
## 1. Data Loading

**Methodology:** We load the cleaned {database.upper()} clinical and expression data that was pre-processed by our generator script. We also transpose and merge the expression data into the clinical cohort, yielding a unified pan-cancer dataset (if multiple cancers were specified).
"""))

    cell_3_code = f"""\
genes = {repr(genes_list)}
database = "{database}"
cancers = {repr(cancers)}

data_dir = '../input/{database}'
if not os.path.exists(data_dir):
    data_dir = 'input/{database}'

clin_dfs = []
expr_dfs = []

for cancer in cancers:
    clin_path = os.path.join(data_dir, f'{{cancer}}_clinical.csv')
    expr_path = os.path.join(data_dir, f'{{cancer}}_expression.csv')
    
    c_df = pd.read_csv(clin_path, low_memory=False)
    e_df = pd.read_csv(expr_path, low_memory=False)
    
    c_df['CANCER_TYPE'] = cancer.upper()
    
    # Preprocessing Expression for ML format
    e_df = e_df.drop(columns=['Entrez_Gene_Id'], errors='ignore')
    if 'Hugo_Symbol' in e_df.columns:
        e_df = e_df.set_index('Hugo_Symbol').T
    else:
        e_df = e_df.set_index(e_df.columns[0]).T
        
    e_df.index.name = 'SAMPLE_ID' if 'SAMPLE_ID' in c_df.columns else 'PATIENT_ID'
    e_df = e_df.reset_index()
    
    if database.lower() == 'tcga':
        # Align TCGA IDs
        e_df[e_df.columns[0]] = e_df[e_df.columns[0]].apply(lambda x: str(x)[:12])
        
    clin_dfs.append(c_df)
    expr_dfs.append(e_df)

# Pan-Cancer Merge
clin_df = pd.concat(clin_dfs, ignore_index=True)
expr_df = pd.concat(expr_dfs, ignore_index=True)

# Merge based on available ID type
merge_key = 'SAMPLE_ID' if 'SAMPLE_ID' in clin_df.columns else 'PATIENT_ID'
df = pd.merge(clin_df, expr_df, on=merge_key, how='inner')

# Drop duplicate patients if any
df = df.drop_duplicates(subset=['PATIENT_ID'])
print(f"Final dataset shape: {{df.shape}}")
if df.shape[0] == 0:
    raise ValueError("After merging clinical and expression data, the dataset is empty. Check if PATIENT_ID/SAMPLE_ID match between files.")
df.head()
"""
    cells.append(nbf.v4.new_code_cell(cell_3_code))

    cells.append(nbf.v4.new_code_cell("""\
# Features & Target
present_genes = [g for g in genes if g in df.columns]
missing_genes = set(genes) - set(present_genes)
if len(present_genes) == 0:
    raise ValueError("None of the signature genes were found in the expression dataset.")
if len(missing_genes) > 0:
    print(f"Warning: The following genes were missing from the dataset and will be ignored: {missing_genes}")

num_genes = len(present_genes)
print(f"Proceeding with {num_genes} genes: {present_genes}")

# Output Directory
out_dir = f'../output/ml_prognostic_results/{database}/{"_".join(cancers)}/' + str(num_genes) + '-gene'
os.makedirs(out_dir, exist_ok=True)

df['5yr_survival'] = np.where((df['OS_MONTHS'] >= 60), 1, 
                              np.where((df['event'] == 1) & (df['OS_MONTHS'] < 60), 0, np.nan))

# train_test_split
X_train, X_test, y_train, y_test = train_test_split(df[present_genes], df, test_size=0.2, random_state=42)

print(f"Train size: {X_train.shape[0]}")
print(f"Test size: {X_test.shape[0]}")
"""))

    cells.append(nbf.v4.new_markdown_cell("""\
## 2. Model 1: Baseline Cox Proportional Hazards Model
"""))

    cells.append(nbf.v4.new_code_cell("""\
cph_train = pd.concat([X_train, y_train[['OS_MONTHS', 'event']]], axis=1)
cph = CoxPHFitter(penalizer=0.1)
cph.fit(cph_train, duration_col='OS_MONTHS', event_col='event')

cph.print_summary()
c_index_cph = cph.concordance_index_
print(f"\\nCox PH C-index on training set: {c_index_cph:.3f}")
"""))

    cells.append(nbf.v4.new_code_cell("""\
cph_test = pd.concat([X_test, y_test[['OS_MONTHS', 'event']]], axis=1)
c_index_test = cph.score(cph_test, scoring_method='concordance_index')
print(f"Cox PH C-index on TEST set: {c_index_test:.3f}")

plt.figure(figsize=(10, 6))
cph.plot()
plt.title(f'Cox PH Hazard Ratios ({num_genes}-gene Signature)')
plt.tight_layout()
plt.savefig(f"{out_dir}/cox_hazard_ratios.png")
plt.show()
"""))

    cells.append(nbf.v4.new_markdown_cell("""\
## 3. Machine Learning Models (Random Forest & Neural Network)

**Methodology:**
This section evaluates the multi-gene signature's non-linear predictive power using two standard approaches:
1. **Random Forest Classifier**
2. **Multi-Layer Perceptron (Neural Network)**

Patients are labeled based on their **5-year Overall Survival (OS)** status (1 = survived > 5 years, 0 = died < 5 years). The models evaluate performance using the **ROC-AUC** metric, explicitly plot the overlapping ROC curves, and save the trained models for future inference.
"""))

    cells.append(nbf.v4.new_code_cell("""\
from sklearn.neural_network import MLPClassifier
import joblib

# Filter valid cases for binary classification (5-year survival)
train_bin = y_train.dropna(subset=['5yr_survival'])
test_bin = y_test.dropna(subset=['5yr_survival'])
X_train_bin = X_train.loc[train_bin.index]
X_test_bin = X_test.loc[test_bin.index]
y_train_bin = train_bin['5yr_survival']
y_test_bin = test_bin['5yr_survival']

# Random Forest
rf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42, class_weight='balanced')
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

plt.plot(fpr_rf, tpr_rf, color='darkorange', lw=2, label=f'Random Forest (AUC = {auc_rf:.3f})')
plt.plot(fpr_mlp, tpr_mlp, color='green', lw=2, label=f'MLP Neural Net (AUC = {auc_mlp:.3f})')
plt.plot([0,1], [0,1], 'k--', lw=2, label='Random Chance')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title(f'ROC Curve: 5-Year Survival Prediction ({num_genes}-gene Signature)')
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(f"{out_dir}/roc_curves.png")
plt.show()

# Plot Feature Importances from Random Forest
feature_imp_rf = pd.Series(rf.feature_importances_, index=X_train.columns).sort_values(ascending=False)
plt.figure(figsize=(10, 6))
sns.barplot(x=feature_imp_rf.values, y=feature_imp_rf.index, palette='viridis')
plt.title(f'RF Gini Importances ({num_genes}-gene Signature)')
plt.tight_layout()
plt.savefig(f"{out_dir}/rf_feature_importances.png")
plt.show()
"""))

    cells.append(nbf.v4.new_markdown_cell("""\
## 4. Kaplan-Meier Survival Analysis (High vs Low Risk)
"""))

    cells.append(nbf.v4.new_code_cell("""\
cph_risk_test = cph.predict_partial_hazard(X_test)
median_risk = np.median(cph_risk_test)
high_risk_mask = cph_risk_test > median_risk

test_high = y_test[high_risk_mask]
test_low = y_test[~high_risk_mask]

kmf = KaplanMeierFitter()
plt.figure(figsize=(8, 6))

kmf.fit(test_low['OS_MONTHS'], event_observed=test_low['event'], label=f'Low Risk (n={len(test_low)})')
ax = kmf.plot()

kmf.fit(test_high['OS_MONTHS'], event_observed=test_high['event'], label=f'High Risk (n={len(test_high)})')
kmf.plot(ax=ax)

try:
    from lifelines.statistics import multivariate_logrank_test
    res = multivariate_logrank_test(y_test['OS_MONTHS'], high_risk_mask, y_test['event'])
    plt.text(0.05, 0.05, f"Log-rank p-value: {res.p_value:.2e}", transform=ax.transAxes, fontsize=12)
except Exception as e:
    pass

plt.title(f'Kaplan-Meier Survival Curve: High vs Low Risk ({num_genes}-gene Signature)')
plt.xlabel('Time (Months)')
plt.ylabel('Survival Probability')
plt.tight_layout()
plt.savefig(f"{out_dir}/km_survival_curve.png")
plt.show()
"""))

    cells.append(nbf.v4.new_code_cell(f"""\
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
"""))

    nb['cells'] = cells
    notebook_filename = 'scripts/ml_prognostic_classifier.ipynb'
    with open(notebook_filename, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)
    
    print(f"Notebook successfully created at: {notebook_filename}")

def main():
    parser = argparse.ArgumentParser(
        description="Generate ML Prognostic Classifier Notebook with Pan-Cancer Support",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # TCGA single cancer:
  python scripts/generate_ml_prognostic_classifier_notebook.py --database tcga --cancer luad

  # TCGA pan-cancer (BRCA and LUAD combined):
  python scripts/generate_ml_prognostic_classifier_notebook.py --database tcga --cancer brca luad

  # cBioPortal pan-cancer (BRCA AURORA and LUAD BROAD combined):
  python scripts/generate_ml_prognostic_classifier_notebook.py --database cbioportal --cancer brca luad --study-id brca_aurora_2023 luad_broad
"""
    )
    parser.add_argument("--genes", nargs="+", 
                        default=['ADAM10', 'C1GALT1', 'ESRRG', 'FZD6', 'GBE1', 'GLS', 'ITGA4', 'PDE3B', 'SGMS1', 'SLC11A2', 'SLC16A7', 'SLC22A1'], 
                        help="List of genes for the signature")
    parser.add_argument("--signature-name", default="STAT3 Core Axis", help="Name of the signature")
    parser.add_argument("--database", default="tcga", help="Data source to use (e.g. cbioportal, tcga, metabric)")
    parser.add_argument("--cancer", nargs="+", default=["brca"], help="List of cancer prefixes (e.g., brca, luad)")
    parser.add_argument("--study-id", nargs="+", default=["brca_metabric"], help="List of cBioPortal study IDs matching the cancers (for fallback downloading)")
    parser.add_argument("--profile-expr", default="mrna_illumina_microarray", help="cBioPortal mRNA profile name")
    args = parser.parse_args()
    
    genes_list = args.genes
    signature_name = args.signature_name
    database = args.database
    cancers = args.cancer
    
    if len(cancers) == 1 and cancers[0].lower() == 'all':
        import glob
        if database.lower() == 'tcga':
            tcga_dir = '../input/TCGA' if not os.path.exists('input') else 'input/TCGA'
            if os.path.exists(tcga_dir):
                files = glob.glob(os.path.join(tcga_dir, 'TCGA-*.survival.tsv.gz'))
                cancers = [os.path.basename(f).split('-')[1].split('.')[0].lower() for f in files]
                print(f"Auto-detected {len(cancers)} cancers for TCGA: {', '.join(cancers)}")
            else:
                raise ValueError(f"TCGA directory {tcga_dir} not found.")
        else:
            data_dir = f'../input/{database}' if not os.path.exists('input') else f'input/{database}'
            if os.path.exists(data_dir):
                files = glob.glob(os.path.join(data_dir, '*_clinical.csv'))
                cancers = [os.path.basename(f).replace('_clinical.csv', '').lower() for f in files]
                print(f"Auto-detected {len(cancers)} cancers for {database}: {', '.join(cancers)}")
            else:
                raise ValueError(f"Data directory {data_dir} not found for {database}.")
                
    study_ids = args.study_id
    profile_expr = args.profile_expr
    num_genes = len(genes_list)
    genes_str = ", ".join(genes_list)

    # 1. Prepare and standardize all datasets before writing notebook
    prepare_data(database, cancers, study_ids, profile_expr)

    # 2. Generate the notebook
    generate_notebook(database, cancers, study_ids, profile_expr, genes_list, signature_name, genes_str, num_genes)

if __name__ == "__main__":
    main()
