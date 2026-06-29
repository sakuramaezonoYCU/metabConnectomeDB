import nbformat as nbf
import os
import base64
import subprocess
import sys
import argparse

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')

def create_notebook(signature_name, genes):
    import pandas as pd
    META_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'massspec_metabolomics', signature_name)
    match_file = os.path.join(META_RESULTS_DIR, 'metabolite_match_table.csv')
    num_matches = 0
    if os.path.exists(match_file):
        try:
            df = pd.read_csv(match_file)
            num_matches = len(df)
        except Exception:
            pass
            raise
            
    nb = nbf.v4.new_notebook()
    
    # 1. Title and Goal
    nb.cells.append(nbf.v4.new_markdown_cell(f"""# Pan-Cancer Mass Spec Metabolomics Integration

### Goal
Integrate the {signature_name} metabolic signature with orthogonal pan-cancer mass spectrometry metabolomics data (PMID: 29396322).

### Purpose
To show the association between the {signature_name} and actual metabolite abundances across multiple cancer types. While the primary signature is derived from single-cell transcriptomics, verifying these metabolic axes at the metabolomic layer provides orthogonal validation.

### Interpretation
- **Metabolite Signature Score Shifts**: Significant PCA score differences between Tumor and Normal samples across the 7 major cancer types (KIRC, BLCA, BRCA, OV, PAAD, PRAD, LGG) indicates that the transcriptomic signature reliably translates to actual metabolic pathway dysregulation.
- **Cross-Cohort Consistency**: If the metabolomic shifts in BRCA (Breast) and OV (Ovarian) align with our scRNA-seq findings, it strongly supports the robustness of the metabolic network across independent cohorts and multi-omic layers.
- **Biomarker Utility**: Detecting consistent metabolite changes linked to the signature opens avenues for blood-based or tissue-based mass spectrometry diagnostics.

### Inputs
- **Mass Spec Data:** `input/massSpecDataMetabolicData_7panCancer_PMID29396322.csv`
- **Metabolite DB:** `output/human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv`
- **Gene Signature:** {signature_name} ({len(genes)} genes)

### Outputs
- **Final Notebook Report:** `output/massspec_metabolomics/massspec_metabolomics_integration_{signature_name}.html`
"""))

    print_config_code = """
print('--- INJECTED PIPELINE CONFIGURATION ---')
from pan_cancer_config import CANCERS_TO_RUN
print(f'CANCERS_TO_RUN: {CANCERS_TO_RUN}')
"""
    nb.cells.append(nbf.v4.new_code_cell(print_config_code))

    # 2. Config Code
    nb.cells.append(nbf.v4.new_code_cell(f"""import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import Image, display

BASE_DIR = os.path.dirname(os.path.abspath('.'))
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
META_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'massspec_metabolomics', '{signature_name}')

TARGET_GENES = {genes}

print(f"==========================================")
print(f"PARAMETERS FOR THIS RUN:")
print(f"Signature Name: {signature_name}")
print(f"Target Genes Investigated (n={{len(TARGET_GENES)}}):")
print(f"{{TARGET_GENES}}")
print(f"==========================================")
"""))

    # Helper function to embed images in Markdown
    def append_image_markdown(title, description, image_filename, attachment_name):
        image_path = os.path.join(META_RESULTS_DIR, image_filename)
        if os.path.exists(image_path):
            with open(image_path, 'rb') as f:
                png_data = base64.b64encode(f.read()).decode('utf-8')
            
            md_content = f"### {title}\n{description}\n\n![{title}](attachment:{attachment_name})"
            cell = nbf.v4.new_markdown_cell(md_content)
            cell.attachments = {attachment_name: {"image/png": png_data}}
            nb.cells.append(cell)
        else:
            nb.cells.append(nbf.v4.new_markdown_cell(f"### {title}\n*Plot not found: {image_filename}*"))

    if num_matches == 0:
        nb.cells.append(nbf.v4.new_markdown_cell("### ❌ 0 matched metabolites\n\nNo mass spec metabolites could be matched to the genes in this signature. Please try other signatures."))
    else:
        # 3. PCA Score Boxplot
        append_image_markdown(
            "Metabolite PCA Signature Score",
            "PCA score derived from the matched metabolites, showing the metabolic axis activity across Tumor vs Normal samples.",
            "metabolite_pca_boxplot.png",
            "pca_boxplot.png"
        )
    
        # 4. Pan-Cancer Heatmap
        append_image_markdown(
            "Pan-Cancer Median Tumor Metabolite Abundance (Z-score)",
            "Median abundance of matched metabolites across cancer types.",
            "pan_cancer_metabolite_heatmap.png",
            "heatmap.png"
        )
    
        # 5. Network Plot
        append_image_markdown(
            "Gene-Metabolite Bipartite Network",
            "Green nodes are metabolites detected in the mass spec data; gray are missing.",
            "gene_metabolite_network.png",
            "network.png"
        )
    
        # 6. Cross-Cohort Comparison
        append_image_markdown(
            "Cross-Cohort Comparison: scRNA vs Mass Spec",
            "Comparison of primary vs metastatic scRNA expression (y-axis) against normal vs tumor mass spec abundance (x-axis).",
            "cross_cohort_comparison_scatter.png",
            "cross_scatter.png"
        )




    scripts_dir = os.path.join(BASE_DIR, 'scripts')
    notebook_filename = os.path.join(scripts_dir, f'massspec_metabolomics_integration.ipynb')
    with open(notebook_filename, 'w') as f:
        nbf.write(nb, f)
        
    print(f"Notebook generated: {notebook_filename}")
    
    # Optional: Execute and convert to HTML
    try:
        env = os.environ.copy()
        env['JUPYTER_DATA_DIR'] = os.path.abspath('./.jupyter_data_cache')
        
        # Run nbconvert from the scripts directory
        subprocess.run([sys.executable, '-m', 'jupyter', 'nbconvert', '--execute', '--to', 'html', os.path.basename(notebook_filename), '--output-dir', '../output/massspec_metabolomics', '--output', f'massspec_metabolomics_integration'], check=True, env=env, cwd=scripts_dir)
        print(f"Notebook exported to HTML: ../output/massspec_metabolomics/massspec_metabolomics_integration.html")
    except subprocess.CalledProcessError as e:
        print(f"Failed to execute and export notebook: {e}")
        sys.exit(1)

if __name__ == '__main__':
    import pandas as pd
    parser = argparse.ArgumentParser(description="Generate Mass Spec Notebook")
    parser.add_argument('--signature_csv', required=True, help="Path to the signature CSV file containing a 'Gene' or 'Strictly_Conserved_Gene' column")
    args = parser.parse_args()
    
    if not os.path.exists(args.signature_csv):
        raise FileNotFoundError(f"CRITICAL ERROR: Signature CSV {args.signature_csv} does not exist.")

    df_sig = pd.read_csv(args.signature_csv)
    if 'Gene' in df_sig.columns:
        genes = df_sig['Gene'].dropna().unique().tolist()
    elif 'Target' in df_sig.columns:
        genes = df_sig['Target'].dropna().unique().tolist()
    elif 'Strictly_Conserved_Gene' in df_sig.columns:
        genes = df_sig['Strictly_Conserved_Gene'].dropna().unique().tolist()
    else:
        raise ValueError(f"CRITICAL ERROR: Could not find 'Gene', 'Target', or 'Strictly_Conserved_Gene' column in {args.signature_csv}")

    signature_name = os.path.basename(args.signature_csv).replace('.csv', '')
    
    print(f"Generating Notebook for {signature_name}...")
    create_notebook(signature_name, genes)
    
    # We DO NOT delete the intermediate meta_results_dir here, because the generated notebook 
    # needs the images! Instead of cleaning up, we just let the notebook reference them.
    # The previous code removed the directory, causing broken image links in the notebook.
    
    print(f"\n✅ All done! The self-contained HTML notebook is in: output/massspec_metabolomics/")
