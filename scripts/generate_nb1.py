import nbformat as nbf
import os

def create_notebook():
    nb = nbf.v4.new_notebook()
    
    nb.cells.append(nbf.v4.new_markdown_cell("""# Deep-Dive: Conserved Metastatic Metabolic Signature
## Exploring STAT3, Oxygen Gradients, Directional CCC, and Clinical Prognosis
This notebook addresses **Priority 1** and **Priority 2** next steps from the AI Summary (Version 5).
Specifically, we cover:
1. **STAT3 Regulatory Network Reconstruction** (Step 1)
2. **Intratumoural Oxygen Gradient Simulation** (Step 2)
3. **Directionality-Aware Metabolic Communication Scoring** (Step 4)
4. **Conserved Signature Score TCGA Clinical Validation** (Step 5)
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import os
import sys
import pandas as pd
import subprocess

# Add scripts dir to path
sys.path.append(os.path.abspath("."))
"""))

    nb.cells.append(nbf.v4.new_markdown_cell("""### Step 1: STAT3 Regulatory Network Reconstruction
**Purpose**: Map all strictly conserved pan-cancer genes within the STAT3 transcriptional network using ChEA 2022 ChIP-Seq data.
**Interpretation**: Validates STAT3 as an upstream regulatory master switch for the metastatic metabolic signature.
"""))
    
    nb.cells.append(nbf.v4.new_code_cell("""from compute_stat3_network import compute_stat3_network
compute_stat3_network()
# Load output
stat3_out = os.path.join("..", "output", "deepdive_conserved_metabGeneSig", "stat3_network", "stat3_u87_targets_strictly_conserved.csv")
if os.path.exists(stat3_out):
    df = pd.read_csv(stat3_out)
    display(df)
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import networkx as nx
import matplotlib.pyplot as plt

if os.path.exists(stat3_out):
    G = nx.from_pandas_edgelist(df, 'Source', 'Target')
    plt.figure(figsize=(8,6))
    pos = nx.spring_layout(G, seed=42)
    nx.draw_networkx_nodes(G, pos, nodelist=['STAT3'], node_color='red', node_size=1000)
    nx.draw_networkx_nodes(G, pos, nodelist=[n for n in G.nodes if n != 'STAT3'], node_color='lightblue', node_size=500)
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    nx.draw_networkx_labels(G, pos)
    plt.title("STAT3 Transcriptional Network of Metastatic Metabolic Genes")
    plt.axis('off')
    plt.show()
"""))

    nb.cells.append(nbf.v4.new_markdown_cell("""### Step 2: Intratumoural Oxygen Gradient Simulation
**Purpose**: Reconstruct oxygen gradients using hypoxia signature proxy scores (VEGFA, SLC2A1, BNIP3) from established HIF-1α literature (Semenza et al.) and project the conserved signature score onto this gradient.
**Interpretation**: Determines if the pre-metastatic subclone maps to the hypoxic core of the primary tumor.

> [!NOTE]
> **Data Resolution Fallback Protocol**
> If specific cell type filtering (e.g., restricting only to malignant epithelial cells) yields too few cells in a primary tumor to build a robust statistical gradient, the simulation automatically falls back to analyzing all available primary cells across the tumor microenvironment. This ensures empirical modeling of the full tumor oxygen gradient instead of failing abruptly, while still accurately projecting the metabolic signature onto the hypoxia continuum.
"""))

    nb.cells.append(nbf.v4.new_code_cell("""from simulate_oxygen_gradient import simulate_oxygen_gradient
simulate_oxygen_gradient()
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import seaborn as sns
import glob

oxy_dir = os.path.join("..", "output", "deepdive_conserved_metabGeneSig", "oxygen_gradient")
oxy_files = glob.glob(os.path.join(oxy_dir, "*_primary_oxygen_gradient_scores.csv"))

if oxy_files:
    # Set up the matplotlib figure
    n_files = len(oxy_files)
    fig, axes = plt.subplots(1, n_files, figsize=(6 * n_files, 5), sharey=True)
    if n_files == 1:
        axes = [axes]
        
    for ax, oxy_file in zip(axes, oxy_files):
        cancer_prefix = os.path.basename(oxy_file).split('_')[0].capitalize()
        oxy_df = pd.read_csv(oxy_file)
        
        sns.regplot(data=oxy_df, x='hypoxia_score', y='metastatic_score', 
                    scatter_kws={'alpha':0.3, 'color': 'gray'}, 
                    line_kws={'color': 'red'}, ax=ax)
        ax.set_title(f'Hypoxia vs Metastatic Score ({cancer_prefix} Primary)')
        ax.set_xlabel('Hypoxia Score (VEGFA, SLC2A1, BNIP3)')
        if ax == axes[0]:
            ax.set_ylabel('Metastatic Signature Score')
        else:
            ax.set_ylabel('')
            
    plt.tight_layout()
    plt.show()
else:
    print("No oxygen gradient results found.")
"""))

    nb.cells.append(nbf.v4.new_markdown_cell("""### Step 4: Directionality-Aware Metabolic Communication Scoring
**Purpose**: Separate enzymes into 'producing' and 'consuming' classes based on MetalinksDB reaction stoichiometry and score metabolic cross-talk directionally.
"""))

    nb.cells.append(nbf.v4.new_code_cell("""from compute_directional_ccc import compute_directional_ccc
compute_directional_ccc()
"""))

    nb.cells.append(nbf.v4.new_code_cell("""ccc_file = os.path.join("..", "output", "deepdive_conserved_metabGeneSig", "directional_ccc", "metalinks_direction_classes.csv")
if os.path.exists(ccc_file):
    ccc_df = pd.read_csv(ccc_file)
    plt.figure(figsize=(6,4))
    sns.countplot(data=ccc_df, x='Direction_Class', palette='Set2')
    plt.title('Directional Metabolic Communication Classes')
    plt.ylabel('Number of Genes')
    plt.show()
"""))

    nb.cells.append(nbf.v4.new_markdown_cell("""### Step 5: TCGA Survival Validation
**Purpose**: Validate the Conserved Metastatic Metabolic Score on TCGA primary tumor data against distant metastasis-free survival.
"""))

    nb.cells.append(nbf.v4.new_code_cell("""from validate_tcga_signature import validate_tcga_signature
validate_tcga_signature()
"""))

    nb.cells.append(nbf.v4.new_markdown_cell("""### Step 5.1: Empirical Validation (Permutation Test)
**Purpose**: To verify that our conserved signature's predictive power is statistically significantly better than randomly chosen genes, we generate 100 random signatures of the same size, compute their Hazard Ratios, and plot our true signature against this null distribution.
"""))

    nb.cells.append(nbf.v4.new_code_cell("""from compute_permutation_null import compute_permutation_null
# Compute the null distribution (this reads the TCGA file and iterates Cox models)
compute_permutation_null(n_permutations=100)
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import matplotlib.pyplot as plt
import seaborn as sns

null_file = os.path.join("..", "output", "deepdive_conserved_metabGeneSig", "tcga_validation", "null_distribution_metrics.csv")
true_file = os.path.join("..", "output", "deepdive_conserved_metabGeneSig", "tcga_validation", "true_signature_metrics.csv")

if os.path.exists(null_file) and os.path.exists(true_file):
    null_df = pd.read_csv(null_file)
    true_df = pd.read_csv(true_file)
    
    # Plot true HR vs null distribution for each cancer
    cancers = true_df['TCGA_Cohort'].unique()
    fig, axes = plt.subplots(len(cancers), 1, figsize=(8, 4 * len(cancers)))
    if len(cancers) == 1:
        axes = [axes]
        
    for ax, cancer in zip(axes, cancers):
        cancer_null = null_df[null_df['TCGA_Cohort'] == cancer]['Hazard_Ratio']
        cancer_true = true_df[true_df['TCGA_Cohort'] == cancer]['Hazard_Ratio'].iloc[0]
        
        sns.histplot(cancer_null, bins=20, kde=True, color='lightgray', ax=ax, label=f'Random Gene Signatures (Null)')
        ax.axvline(cancer_true, color='red', linestyle='--', linewidth=2, label=f'True Signature (HR={cancer_true:.2f})')
        
        ax.set_title(f'TCGA-{cancer}: True Signature vs. Null Distribution')
        ax.set_xlabel('Hazard Ratio (High vs Low Risk)')
        ax.set_ylabel('Frequency')
        ax.legend()
        
    plt.tight_layout()
    plt.show()
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import subprocess
import sys
import os

notebook_filename = 'deepdive_conserved_metabGeneSig.ipynb'
output_base = 'deepdive_conserved_metabGeneSig'
output_dir = os.path.join('..', 'output', 'deepdive_conserved_metabGeneSig')
os.makedirs(output_dir, exist_ok=True)

cmd_html = [sys.executable, "-m", "jupyter", "nbconvert", "--to", "html", "--execute", notebook_filename, "--output-dir", output_dir, "--output", output_base]
res_html = subprocess.run(cmd_html, capture_output=True, text=True)

if res_html.returncode == 0:
    print(f"🎉 SUCCESS: Notebook successfully exported to '{os.path.join(output_dir, output_base)}.html'")
else:
    print("❌ HTML export failed.")
    print(res_html.stderr)
"""))

    out_path = os.path.join(os.path.dirname(__file__), "deepdive_conserved_metabGeneSig.ipynb")
    if os.path.exists(out_path):
        print(f"⚠️ {out_path} already exists. Skipping overwrite to protect user edits.")
    else:
        with open(out_path, 'w') as f:
            nbf.write(nb, f)
        print(f"Created {out_path}")

if __name__ == "__main__":
    create_notebook()
