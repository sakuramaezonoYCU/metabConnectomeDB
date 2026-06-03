import os
import sys
import nbformat as nbf
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Generate ML Prognostic Classifier Notebook",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with the 12-gene STAT3 Core Axis (default):
  python generate_ml_prognostic_classifier_notebook.py

  # Run with the 21-Gene Directed Metastatic Signature:
  python generate_ml_prognostic_classifier_notebook.py \\
    --signature-name "Directed Metastatic Signature" \\
    --genes GLS SGMS1 SPTLC1 GBE1 SLC16A7 AUH FZD6 NR1D2 CD46 MTMR1 ESRRG ITGA4 SLC11A2 ERAP1 C1GALT1 ADAM10 TRPM8 SLC22A1 AMDHD1 EPOR PDE3B
"""
    )
    parser.add_argument("--genes", nargs="+", 
                        default=['ADAM10', 'C1GALT1', 'ESRRG', 'FZD6', 'GBE1', 'GLS', 'ITGA4', 'PDE3B', 'SGMS1', 'SLC11A2', 'SLC16A7', 'SLC22A1'], 
                        help="List of genes for the signature")
    parser.add_argument("--signature-name", default="STAT3 Core Axis", help="Name of the signature")
    args = parser.parse_args()
    
    genes_list = args.genes
    signature_name = args.signature_name
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
# ML Prognostic Classifier using the {num_genes}-gene {signature_name}

**Purpose:** We have identified a {num_genes}-gene {signature_name} signature: `{genes_str}`. This notebook builds and validates a Machine Learning Prognostic Classifier using this signature on the METABRIC breast cancer cohort.

**Methodology:**
1. Download METABRIC clinical and mRNA expression data from cBioPortal.
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
    cells.append(nbf.v4.new_markdown_cell("""\
## 1. Data Fetching and Preprocessing

**Methodology:** We fetch the METABRIC clinical and expression data directly from the cBioPortal datahub. Data is saved locally in `../input/metabric/` so we don't redownload it repeatedly.

**Interpretation Guide:** Patients with missing survival times or status are dropped. We convert `OS_STATUS` to a boolean event indicator (1 = Death, 0 = Censored).
"""))

    cell_3_code = f"""\
genes = {repr(genes_list)}

data_dir = '../input/metabric'
os.makedirs(data_dir, exist_ok=True)

clin_url = 'https://raw.githubusercontent.com/cBioPortal/datahub/master/public/brca_metabric/data_clinical_patient.txt'
expr_url = 'https://raw.githubusercontent.com/cBioPortal/datahub/master/public/brca_metabric/data_mrna_illumina_microarray.txt'

clin_path = os.path.join(data_dir, 'clinical.csv')
expr_path = os.path.join(data_dir, 'expression.csv')

# Download clinical data if needed
try:
    if not os.path.exists(clin_path) or os.path.getsize(clin_path) < 1000:
        clin_df = pd.read_csv(clin_url, sep='\\t', skiprows=4)
        clin_df.to_csv(clin_path, index=False)
    else:
        clin_df = pd.read_csv(clin_path)

    # Download expression data if needed
    if not os.path.exists(expr_path) or os.path.getsize(expr_path) < 1000:
        expr_df = pd.read_csv(expr_url, sep='\\t')
        expr_df.to_csv(expr_path, index=False)
    else:
        expr_df = pd.read_csv(expr_path)
except Exception as e:
    print(f"\\n🚨 ERROR: Failed to download METABRIC data automatically.")
    print(f"This is likely due to a network proxy blocking raw.githubusercontent.com.")
    print(f"Please manually download the following two files:")
    print(f"1. {{clin_url}}")
    print(f"   --> Save as: {{os.path.abspath(clin_path)}}")
    print(f"2. {{expr_url}}")
    print(f"   --> Save as: {{os.path.abspath(expr_path)}}")
    print(f"\\nOnce saved, please re-run this notebook. Exception details: {{e}}")
    raise e

# Preprocessing Clinical
# We need OS_MONTHS and OS_STATUS
clin_df = clin_df[['PATIENT_ID', 'OS_MONTHS', 'OS_STATUS']].dropna()
clin_df['event'] = clin_df['OS_STATUS'].apply(lambda x: 1 if 'DECEASED' in str(x).upper() else 0)

# Preprocessing Expression
# Set index to Hugo_Symbol, drop Entrez, transpose so rows are patients
expr_df = expr_df.drop(columns=['Entrez_Gene_Id'], errors='ignore').set_index('Hugo_Symbol').T
expr_df.index.name = 'PATIENT_ID'
expr_df = expr_df.reset_index()

# Merge
df = pd.merge(clin_df, expr_df, on='PATIENT_ID', how='inner')
print(f"Final dataset shape: {{df.shape}}")
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

out_dir = '../output/ml_prognostic_results/""" + str(num_genes) + """-gene'
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
