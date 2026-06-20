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
    
    try:
        from pan_cancer_config import ANALYSIS_SUFFIX
    except ImportError:
        ANALYSIS_SUFFIX = ''
        
    cells = []
    cancers_str = ", ".join([c.upper() for c in cancers])
    
    cells.append(nbf.v4.new_markdown_cell(f"""\
# ML Prognostic Classifier using Pan-Cancer Signatures ({database.upper()} - {cancers_str})

### Goal
Build and validate Machine Learning Prognostic Classifiers evaluating **each 4-cancer combination separately**.

### Interpretation
- **Cox Proportional Hazards C-index**: Establishes baseline linear survival predictive power.
- **Random Forest & MLP ROC-AUC**: Evaluates non-linear and interactive predictive power for 5-year binary overall survival.
- **Kaplan-Meier High vs Low Risk**: Stratifies patients based on machine learning risk scores.
"""))

    cells.append(nbf.v4.new_code_cell("""\
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
import glob

# ML & Survival Modeling
from lifelines import CoxPHFitter, KaplanMeierFitter
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, roc_curve

import warnings
warnings.filterwarnings('ignore')

sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 300
"""))

    cells.append(nbf.v4.new_markdown_cell(f"""\
## 1. Data Loading

We load the cleaned {database.upper()} clinical and expression data and merge them into a unified dataset.
"""))

    cells.append(nbf.v4.new_code_cell(f"""\
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

merge_key = 'SAMPLE_ID' if 'SAMPLE_ID' in clin_df.columns else 'PATIENT_ID'
df = pd.merge(clin_df, expr_df, on=merge_key, how='inner')
df = df.drop_duplicates(subset=['PATIENT_ID'])
print(f"Final dataset shape: {{df.shape}}")

df['5yr_survival'] = np.where((df['OS_MONTHS'] >= 60), 1, 
                              np.where((df['event'] == 1) & (df['OS_MONTHS'] < 60), 0, np.nan))

# Glob all signatures to evaluate
try:
    from pan_cancer_config import ANALYSIS_SUFFIX
except ImportError:
    ANALYSIS_SUFFIX = ''

meta_dir = '../output/pan_cancer_meta_results'
if not os.path.exists(meta_dir):
    meta_dir = 'output/pan_cancer_meta_results'
    
signature_files = glob.glob(os.path.join(meta_dir, f"pan_cancer_signature_*{ANALYSIS_SUFFIX}.csv"))
if not signature_files:
    strict_sig = os.path.join(meta_dir, f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv")
    if os.path.exists(strict_sig):
        signature_files = [strict_sig]

genes_dict = {{}}
if not signature_files:
    print("No signatures found. Falling back to dynamic genes list.")
    if 'scripts' not in sys.path and '.' not in sys.path:
        sys.path.append('scripts')
        sys.path.append('.')
    from dynamic_genes import get_dynamic_genes
    genes_dict["Dynamic_All"] = get_dynamic_genes(base_dir='..')
else:
    for sig_file in signature_files:
        sig_name = os.path.basename(sig_file).replace('.csv', '').replace('pan_cancer_signature_', '')
        df_sig = pd.read_csv(sig_file)
        genes_dict[sig_name] = df_sig.iloc[:, 0].tolist()

print(f"\\nFound {{len(genes_dict)}} combinations to evaluate.")
"""))

    # Now loop over signatures
    cells.append(nbf.v4.new_code_cell("""\
# Loop through signatures and evaluate models
out_base_dir = f'../output/ml_prognostic_results/{database}/{"_".join(cancers)}'
os.makedirs(out_base_dir, exist_ok=True)

all_metrics = []

for sig_name, genes in genes_dict.items():
    print("\\n" + "="*80)
    print(f"Evaluating Signature: {sig_name}")
    print("="*80)
    
    present_genes = [g for g in genes if g in df.columns]
    if len(present_genes) < 2:
        print(f"Skipping {sig_name}: Not enough valid genes ({len(present_genes)}) found in expression dataset.")
        continue
        
    num_genes = len(present_genes)
    print(f"Proceeding with {num_genes} genes.")
    
    out_dir = os.path.join(out_base_dir, sig_name)
    os.makedirs(out_dir, exist_ok=True)
    
    X_train, X_test, y_train, y_test = train_test_split(df[present_genes], df, test_size=0.2, random_state=42)
    
    # 1. Drop zero-variance
    variances = X_train.var()
    zero_var_cols = variances[variances < 1e-5].index
    if len(zero_var_cols) > 0:
        X_train = X_train.drop(columns=zero_var_cols)
        X_test = X_test.drop(columns=zero_var_cols)

    # 2. Drop collinear
    corr_matrix = X_train.corr().abs()
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    to_drop = [column for column in upper.columns if any(upper[column] > 0.95)]
    if len(to_drop) > 0:
        X_train = X_train.drop(columns=to_drop)
        X_test = X_test.drop(columns=to_drop)
        
    cph_train = pd.concat([X_train, y_train[['OS_MONTHS', 'event']]], axis=1)
    
    # Find penalizer
    from lifelines.utils import k_fold_cross_validation
    penalizers = [0.01, 0.05, 0.1, 0.5, 1.0]
    best_score = -np.inf
    best_penalizer = 0.1
    for p in penalizers:
        cph_cv = CoxPHFitter(penalizer=p, l1_ratio=0.1)
        try:
            scores = k_fold_cross_validation(cph_cv, cph_train, duration_col='OS_MONTHS', event_col='event', k=5, scoring_method="concordance_index")
            if np.mean(scores) > best_score:
                best_score = np.mean(scores)
                best_penalizer = p
        except:
            pass
            
    if best_score == -np.inf: best_penalizer = 1.0
    
    print(f"Train size: {len(X_train)}")
    print(f"Test size: {len(X_test)}")
    print(f"Selected Optimal Penalizer: {best_penalizer}")
    if best_score != -np.inf:
        print(f"(CV C-index: {best_score:.3f})")
    
    try:
        cph = CoxPHFitter(penalizer=best_penalizer, l1_ratio=0.1)
        cph.fit(cph_train, duration_col='OS_MONTHS', event_col='event')
    except Exception as e:
        print(f"Skipping {sig_name}: CoxPH failed to converge ({e}).")
        continue
        
    cph_test = pd.concat([X_test, y_test[['OS_MONTHS', 'event']]], axis=1)
    c_index_test = cph.score(cph_test, scoring_method='concordance_index')
    print(f"Cox PH C-index on TEST set: {c_index_test:.3f}")
    
    all_metrics.append({
        "Signature": sig_name,
        "Train_Size": len(X_train),
        "Test_Size": len(X_test),
        "Optimal_Penalizer": best_penalizer,
        "CV_C_Index": best_score if best_score != -np.inf else np.nan,
        "Test_C_Index": c_index_test
    })
    
    fig_height = max(8.0, X_train.shape[1] * 0.5)
    fig, ax = plt.subplots(figsize=(15, fig_height))
    cph.plot(ax=ax)
    plt.title(f'Cox PH Hazard Ratios - {sig_name}')
    plt.tight_layout()
    plt.savefig(f"{out_dir}/cox_hazard_ratios.png")
    plt.show()
    
    # ML Models
    train_bin = y_train.dropna(subset=['5yr_survival'])
    test_bin = y_test.dropna(subset=['5yr_survival'])
    X_train_bin = X_train.loc[train_bin.index]
    X_test_bin = X_test.loc[test_bin.index]
    y_train_bin = train_bin['5yr_survival']
    y_test_bin = test_bin['5yr_survival']
    
    rf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42, class_weight='balanced')
    rf.fit(X_train_bin, y_train_bin)
    
    mlp = MLPClassifier(hidden_layer_sizes=(32, 16), max_iter=500, random_state=42)
    mlp.fit(X_train_bin, y_train_bin)
    
    joblib.dump(rf, f"{out_dir}/rf_model.pkl")
    joblib.dump(mlp, f"{out_dir}/mlp_model.pkl")
    
    rf_risk_test = 1 - rf.predict_proba(X_test_bin)[:, 1]
    mlp_risk_test = 1 - mlp.predict_proba(X_test_bin)[:, 1]
    
    plt.figure(figsize=(8, 6))
    fpr_rf, tpr_rf, _ = roc_curve(1 - y_test_bin, rf_risk_test)
    auc_rf = roc_auc_score(1 - y_test_bin, rf_risk_test)
    fpr_mlp, tpr_mlp, _ = roc_curve(1 - y_test_bin, mlp_risk_test)
    auc_mlp = roc_auc_score(1 - y_test_bin, mlp_risk_test)
    plt.plot(fpr_rf, tpr_rf, color='darkorange', lw=2, label=f'RF (AUC = {auc_rf:.3f})')
    plt.plot(fpr_mlp, tpr_mlp, color='green', lw=2, label=f'MLP (AUC = {auc_mlp:.3f})')
    plt.plot([0,1], [0,1], 'k--', lw=2, label='Random Chance')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curve: 5-Year Survival - {sig_name}')
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/roc_curves.png")
    plt.show()
    
    feature_imp_rf = pd.Series(rf.feature_importances_, index=X_train.columns).sort_values(ascending=False)
    fig_height = max(8.0, len(feature_imp_rf) * 0.5)
    fig, ax = plt.subplots(figsize=(10, fig_height))
    sns.barplot(x=feature_imp_rf.values, y=feature_imp_rf.index, palette='viridis', ax=ax)
    plt.title(f'RF Gini Importances - {sig_name}')
    plt.tight_layout()
    plt.savefig(f"{out_dir}/rf_feature_importances.png")
    plt.show()
    
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
    except:
        pass
        
    plt.title(f'Kaplan-Meier Survival Curve - {sig_name}')
    plt.xlabel('Time (Months)')
    plt.ylabel('Survival Probability')
    plt.tight_layout()
    plt.savefig(f"{out_dir}/km_survival_curve.png")
    plt.show()

if all_metrics:
    metrics_df = pd.DataFrame(all_metrics)
    metrics_df.to_csv(f"{out_base_dir}/ml_metrics.csv", index=False)
    print(f"\\nSaved ML metrics summary to {out_base_dir}/ml_metrics.csv")

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
out_dir = f'../output/ml_prognostic_results/{database}/{"_".join(cancers)}'
output_base = 'ml_prognostic_classifier_report_' + "_".join(cancers) + ANALYSIS_SUFFIX

jupyter_bin = os.path.join(os.path.dirname(sys.executable), 'jupyter')
if not os.path.exists(jupyter_bin): jupyter_bin = 'jupyter'

# Generate HTML in the output root directly so it's easy to find
cmd_html = [jupyter_bin, "nbconvert", "--to", "html", notebook_filename, "--output-dir", "../output", "--output", output_base]
res_html = subprocess.run(cmd_html, capture_output=True, text=True)

if res_html.returncode == 0:
    print(f"🎉 SUCCESS: Notebook successfully exported to '{{os.path.join('../output', output_base)}}.html'")
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
    parser.add_argument("--database", default="tcga", help="Data source to use (e.g. cbioportal, tcga, metabric)")
    parser.add_argument("--cancer", nargs="+", default=["brca"], help="List of cancer prefixes (e.g., brca, luad)")
    parser.add_argument("--study-id", nargs="+", default=["brca_metabric"], help="List of cBioPortal study IDs matching the cancers (for fallback downloading)")
    parser.add_argument("--profile-expr", default="mrna_illumina_microarray", help="cBioPortal mRNA profile name")
    args = parser.parse_args()
    
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

    # 1. Prepare and standardize all datasets before writing notebook
    prepare_data(database, cancers, study_ids, profile_expr)

    # 2. Generate the notebook
    generate_notebook(database, cancers, study_ids, profile_expr, [], "Dynamic_Signature", "Dynamic", 0)

if __name__ == "__main__":
    main()
