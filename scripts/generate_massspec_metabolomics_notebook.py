import nbformat as nbf
import os
import base64
import subprocess
import sys
import argparse

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')

def create_notebook(signature_name, genes):
    META_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'massspec_metabolomics', signature_name)
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
- **Mass Spec Data:** `input/databases/massSpecDataMetabolicData_7panCancer_PMID29396322.csv`
- **Metabolite DB:** `input/databases/human_database_merge_unique_metab_with_HMDB_Info.csv`
- **Gene Signature:** {signature_name} ({len(genes)} genes)

### Outputs
- **Final Notebook Report:** `output/massspec_metabolomics/massspec_metabolomics_integration_{signature_name}.html`
"""))

    # 2. Config Code
    nb.cells.append(nbf.v4.new_code_cell(f"""import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import Image, display

BASE_DIR = os.path.dirname(os.path.abspath('.'))
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
META_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'massspec_metabolomics', '{signature_name}')
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

    # 7. Export code
    nb.cells.append(nbf.v4.new_code_cell(f"""import subprocess
import sys

notebook_filename = 'massspec_metabolomics_integration_{signature_name}.ipynb'
output_base = 'massspec_metabolomics_integration_{signature_name}'
output_dir = os.path.join(OUTPUT_DIR, 'massspec_metabolomics')
os.makedirs(output_dir, exist_ok=True)

jupyter_bin = os.path.join(os.path.dirname(sys.executable), 'jupyter')
if not os.path.exists(jupyter_bin): jupyter_bin = 'jupyter'

cmd_html = [jupyter_bin, "nbconvert", "--to", "html", "--execute", notebook_filename, "--output-dir", output_dir, "--output", output_base]
res_html = subprocess.run(cmd_html, capture_output=True, text=True)

if res_html.returncode == 0:
    print(f"🎉 SUCCESS: Notebook exported to '{{os.path.join(output_dir, output_base)}}.html'")
else:
    print("❌ HTML export failed.")
"""))

    scripts_dir = os.path.join(BASE_DIR, 'scripts')
    nb.cells.append(nbf.v4.new_code_cell(f"""# ==========================================
# HTML EXPORT
# ==========================================
import subprocess
import sys
import os

notebook_filename = 'massspec_metabolomics_integration_{signature_name}.ipynb'
output_base = 'massspec_metabolomics_integration_{signature_name}'
output_dir = os.path.join('..', 'output', 'massspec_metabolomics')
os.makedirs(output_dir, exist_ok=True)

jupyter_bin = os.path.join(os.path.dirname(sys.executable), 'jupyter')
if not os.path.exists(jupyter_bin): jupyter_bin = 'jupyter'

cmd_html = [jupyter_bin, "nbconvert", "--to", "html", notebook_filename, "--output-dir", output_dir, "--output", output_base]
res_html = subprocess.run(cmd_html, capture_output=True, text=True)

if res_html.returncode == 0:
    print(f"🎉 SUCCESS: Notebook successfully exported to '{{os.path.join(output_dir, output_base)}}.html'")
else:
    print("❌ HTML export failed.")
    print(res_html.stderr)
"""))

    notebook_filename = os.path.join(scripts_dir, f'massspec_metabolomics_integration_{signature_name}.ipynb')
    with open(notebook_filename, 'w') as f:
        nbf.write(nb, f)
        
    print(f"Notebook generated: {notebook_filename}")
    
    # Optional: Execute and convert to HTML
    try:
        env = os.environ.copy()
        env['JUPYTER_DATA_DIR'] = os.path.abspath('./.jupyter_data_cache')
        
        # Run nbconvert from the scripts directory
        subprocess.run([sys.executable, '-m', 'jupyter', 'nbconvert', '--execute', '--to', 'html', os.path.basename(notebook_filename), '--output-dir', '../output/massspec_metabolomics', '--output', f'massspec_metabolomics_integration_{signature_name}'], check=True, env=env, cwd=scripts_dir)
        print(f"Notebook exported to HTML: ../output/massspec_metabolomics/massspec_metabolomics_integration_{signature_name}.html")
    except subprocess.CalledProcessError as e:
        print(f"Failed to execute and export notebook: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate Mass Spec Notebook")
    parser.add_argument('--signature-name', default="Directed Metastatic Signature", help="Name of the signature")
    parser.add_argument('--genes', nargs='+', default="GLS SGMS1 SPTLC1 GBE1 SLC16A7 AUH FZD6 NR1D2 CD46 MTMR1 ESRRG ITGA4 SLC11A2 ERAP1 C1GALT1 ADAM10 TRPM8 SLC22A1 AMDHD1 EPOR PDE3B".split(), help="List of genes in the signature")
    args = parser.parse_args()
    args = parser.parse_args()
    
    print(f"Running full mass spec pipeline for {args.signature_name}...")
    
    # 1. Run core analysis
    try:
        print("Step 1: Running massspec_metabolomics_analysis.py...")
        analysis_script = os.path.join(BASE_DIR, 'scripts', 'massspec_metabolomics_analysis.py')
        subprocess.run([sys.executable, analysis_script, '--signature-name', args.signature_name, '--genes'] + args.genes, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running core analysis: {e}")
        sys.exit(1)
        
    # 2. Run cross-cohort comparison
    try:
        print("Step 2: Running massspec_cross_cohort_comparison.py...")
        cross_script = os.path.join(BASE_DIR, 'scripts', 'massspec_cross_cohort_comparison.py')
        subprocess.run([sys.executable, cross_script, '--signature-name', args.signature_name], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running cross-cohort comparison: {e}")
        sys.exit(1)
        
    print("Step 3: Generating Notebook...")
    create_notebook(args.signature_name, args.genes)
    
    print("Step 4: Cleaning up intermediate raw files to keep output clean...")
    import shutil
    meta_results_dir = os.path.join(OUTPUT_DIR, 'massspec_metabolomics', args.signature_name)
    if os.path.exists(meta_results_dir):
        shutil.rmtree(meta_results_dir)
        print(f"Removed intermediate data folder: {meta_results_dir}")
        
    print(f"\n✅ All done! The self-contained HTML notebook is in: output/massspec_metabolomics/")
