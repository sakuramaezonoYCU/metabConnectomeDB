"""
Purpose: Executes the downstream advanced notebooks (Pan-Cancer Meta Analysis, Druggability, Visium Spatial). It automatically converts the executed .ipynb files into .html reports.
"""
import json
import os
import sys
import io
import base64
import re
from contextlib import redirect_stdout, redirect_stderr

# Configure Matplotlib for headless execution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
import warnings
warnings.filterwarnings("ignore")


def escape_html(text):
    return text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('"', '&quot;').replace("'", '&#x27;')

def simple_markdown_to_html(md_text):
    # Split text into lines
    lines = md_text.split('\n')
    html_lines = []
    in_list = False
    
    for line in lines:
        stripped = line.strip()
        
        # Headers
        if stripped.startswith('###'):
            if in_list:
                html_lines.append('</ul>')
                in_list = False
            html_lines.append(f'<h3 style="color:#0f172a;font-size:1.3rem;margin-top:20px;font-weight:600;">{stripped[3:].strip()}</h3>')
        elif stripped.startswith('##'):
            if in_list:
                html_lines.append('</ul>')
                in_list = False
            html_lines.append(f'<h2 style="color:#1e293b;font-size:1.6rem;margin-top:30px;border-bottom:1px solid #e2e8f0;padding-bottom:8px;font-weight:600;">{stripped[2:].strip()}</h2>')
        elif stripped.startswith('#'):
            if in_list:
                html_lines.append('</ul>')
                in_list = False
            html_lines.append(f'<h1 style="color:#0f172a;font-size:2.0rem;margin-top:40px;margin-bottom:15px;font-weight:700;">{stripped[1:].strip()}</h1>')
        # Lists
        elif stripped.startswith('-') or stripped.startswith('*'):
            if not in_list:
                html_lines.append('<ul style="margin-left:20px;padding-left:10px;line-height:1.7;color:#334155;">')
                in_list = True
            # Parse list item bold and inline code
            item_text = stripped[1:].strip()
            item_text = re.sub(r'\*\*(.*?)\*\*', r'<strong>\1</strong>', item_text)
            item_text = re.sub(r'`(.*?)`', r'<code style="background:#f1f5f9;border:1px solid #e2e8f0;padding:2px 5px;border-radius:4px;font-family:\'Fira Code\',monospace;font-size:0.9em;color:#b91c1c;">\1</code>', item_text)
            html_lines.append(f'<li style="margin-bottom:6px;">{item_text}</li>')
        else:
            if in_list:
                html_lines.append('</ul>')
                in_list = False
            
            if stripped:
                # Regular paragraph text
                item_text = line
                item_text = re.sub(r'\*\*(.*?)\*\*', r'<strong>\1</strong>', item_text)
                item_text = re.sub(r'`(.*?)`', r'<code style="background:#f1f5f9;border:1px solid #e2e8f0;padding:2px 5px;border-radius:4px;font-family:\'Fira Code\',monospace;font-size:0.9em;color:#b91c1c;">\1</code>', item_text)
                html_lines.append(f'<p style="margin-bottom:15px;color:#334155;line-height:1.7;">{item_text}</p>')
            else:
                html_lines.append('<div style="height:10px;"></div>')
                
    if in_list:
        html_lines.append('</ul>')
        
    return "\n".join(html_lines)

def export_to_gorgeous_html(nb, html_path, title_text):
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title_text} - Standard Report</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=Fira+Code:wght@400;500&display=swap" rel="stylesheet">
    <style>
        body {{
            font-family: 'Inter', sans-serif;
            line-height: 1.6;
            color: #334155;
            background-color: #f8fafc;
            margin: 0;
            padding: 0;
        }}
        .container {{
            max-width: 1100px;
            margin: 0 auto;
            padding: 40px 20px;
        }}
        .header {{
            border-bottom: 1px solid #e2e8f0;
            padding-bottom: 30px;
            margin-bottom: 40px;
            background: #ffffff;
            padding: 30px;
            border-radius: 12px;
            border: 1px solid #e2e8f0;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05);
        }}
        .header h1 {{
            font-size: 2.4rem;
            color: #0f172a;
            margin: 0 0 12px 0;
            font-weight: 800;
            letter-spacing: -0.025em;
        }}
        .header p {{
            font-size: 1.1rem;
            color: #64748b;
            margin: 0 0 15px 0;
        }}
        .badge {{
            display: inline-block;
            background: #eff6ff;
            color: #1d4ed8;
            padding: 6px 14px;
            border-radius: 9999px;
            font-size: 0.85rem;
            font-weight: 600;
            border: 1px solid #bfdbfe;
        }}
        .cell {{
            background: #ffffff;
            border-radius: 12px;
            padding: 24px;
            margin-bottom: 30px;
            box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.05), 0 1px 2px -1px rgba(0, 0, 0, 0.05);
            border: 1px solid #e2e8f0;
        }}
        .cell.markdown {{
            background: transparent;
            border: none;
            box-shadow: none;
            padding: 0 12px;
            margin-bottom: 20px;
        }}
        pre.code-block {{
            background: #f8fafc;
            padding: 16px;
            border-radius: 8px;
            overflow-x: auto;
            border: 1px solid #e2e8f0;
            margin: 0 0 20px 0;
            position: relative;
        }}
        pre.code-block::before {{
            content: "PYTHON";
            position: absolute;
            top: 8px;
            right: 12px;
            font-size: 0.7rem;
            color: #94a3b8;
            font-weight: 700;
            letter-spacing: 0.05em;
        }}
        code.python {{
            font-family: 'Fira Code', monospace;
            font-size: 0.95rem;
            color: #0f172a;
        }}
        .output-container {{
            margin-top: 20px;
            border-top: 1px solid #e2e8f0;
            padding-top: 20px;
        }}
        .output-title {{
            font-size: 0.85rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: #64748b;
            margin-bottom: 10px;
            font-weight: 700;
        }}
        .output-stream {{
            background: #f8fafc;
            padding: 16px;
            border-radius: 8px;
            font-family: 'Fira Code', monospace;
            font-size: 0.9rem;
            color: #0f172a;
            white-space: pre-wrap;
            border-left: 4px solid #10b981;
            margin-bottom: 15px;
            border: 1px solid #e2e8f0;
            border-left-width: 4px;
        }}
        .output-stderr {{
            background: #fef2f2;
            color: #ef4444;
            border-left-color: #ef4444;
        }}
        .output-image {{
            text-align: center;
            margin: 25px 0;
            padding: 15px;
            background: #ffffff;
            border-radius: 12px;
            border: 1px solid #e2e8f0;
        }}
        .output-image img {{
            max-width: 100%;
            border-radius: 8px;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        }}
        .output-html {{
            overflow-x: auto;
            margin: 15px 0;
            background: #ffffff;
            padding: 10px;
            border-radius: 8px;
            border: 1px solid #e2e8f0;
        }}
        table.dataframe {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9rem;
            min-width: 400px;
            color: #334155;
        }}
        table.dataframe th, table.dataframe td {{
            padding: 12px 16px;
            text-align: left;
            border-bottom: 1px solid #e2e8f0;
        }}
        table.dataframe th {{
            background-color: #f1f5f9;
            color: #0f172a;
            font-weight: 600;
            border-bottom: 2px solid #cbd5e1;
        }}
        table.dataframe tr:hover {{
            background-color: #f8fafc;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>{escape_html(title_text)}</h1>
            <p>Aesthetics-Preserved Automated Data & Chemical Taxonomy Analysis Report</p>
            <span class="badge">Pipeline Standardized & Verified</span>
        </div>
"""
    
    for cell in nb.get('cells', []):
        cell_type = cell.get('cell_type')
        if cell_type == 'markdown':
            source = "".join(cell.get('source', []))
            html_source = simple_markdown_to_html(source)
            html_content += f'<div class="cell markdown">{html_source}</div>\n'
            
        elif cell_type == 'code':
            source = "".join(cell.get('source', []))
            if not source.strip():
                continue
            html_content += f'<div class="cell code">\n'
            html_content += f'<pre class="code-block"><code class="python">{escape_html(source)}</code></pre>\n'
            
            outputs = cell.get('outputs', [])
            if outputs:
                html_content += '<div class="output-container">\n'
                html_content += '<div class="output-title">Outputs & Visualizations</div>\n'
                for output in outputs:
                    out_type = output.get('output_type')
                    if out_type == 'stream':
                        name = output.get('name')
                        text = "".join(output.get('text', []))
                        css_class = "output-stream" + (" output-stderr" if name == "stderr" else "")
                        html_content += f'<pre class="{css_class}">{escape_html(text)}</pre>\n'
                    elif out_type in ['display_data', 'execute_result']:
                        data = output.get('data', {})
                        if 'image/png' in data:
                            img_b64 = data['image/png'].replace('\n', '').strip()
                            html_content += f'<div class="output-image"><img src="data:image/png;base64,{img_b64}" /></div>\n'
                        elif 'text/html' in data:
                            html_content += f'<div class="output-html">{"".join(data["text/html"])}</div>\n'
                        elif 'text/plain' in data:
                            html_content += f'<pre class="output-stream">{"".join(data["text/plain"])}</pre>\n'
                    elif out_type == 'error':
                        traceback = "\n".join(output.get('traceback', []))
                        html_content += f'<pre class="output-stream output-stderr">{escape_html(traceback)}</pre>\n'
                html_content += '</div>\n'
            html_content += '</div>\n'
            
    html_content += """
    </div>
</body>
</html>
"""
    os.makedirs(os.path.dirname(html_path), exist_ok=True)
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

def execute_and_export(notebook_path, html_path, title_text):
    notebook_path = os.path.abspath(notebook_path)
    html_path = os.path.abspath(html_path)
    print(f"🚀 Running execution and report generation for {os.path.basename(notebook_path)}...")
    
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = json.load(f)
        
    display_events = []

    def fake_display(obj):
        html_repr = obj._repr_html_() if hasattr(obj, "_repr_html_") else None
        if html_repr is not None:
            display_events.append({"text/html": [html_repr]})
        elif hasattr(obj, "_repr_png_"):
            png_data = obj._repr_png_()
            if png_data:
                b64 = png_data if isinstance(png_data, str) else base64.b64encode(png_data).decode('utf-8')
                display_events.append({
                    "image/png": b64,
                    "text/plain": ["<IPython.core.display.Image object>"]
                })
        elif hasattr(obj, "data") and isinstance(obj.data, str):
            display_events.append({"text/html": [obj.data]})
        else:
            print(obj)

    global_dict = {
        '__name__': '__main__',
        'plt': plt,
        'display': fake_display,
        'SAVE_AS_HTML': False, # Avoid nested nbconvert subprocess calls!
    }
    
    import IPython.display
    IPython.display.display = fake_display
    
    new_cells = []
    notebook_dir = os.path.dirname(os.path.abspath(notebook_path))
    original_dir = os.getcwd()
    os.chdir(notebook_dir)
    
    # Prepend active virtual environment bin directory to PATH so execution contexts locate files properly
    venv_bin = os.path.dirname(sys.executable)
    if os.path.exists(venv_bin):
        os.environ['PATH'] = venv_bin + os.pathsep + os.environ.get('PATH', '')
        
    try:
        for idx, cell in enumerate(nb.get('cells', [])):
            if cell.get('cell_type') == 'code':
                source_lines = cell.get('source', [])
                cleaned_lines = []
                for line in source_lines:
                    s = line.strip()
                    if s.startswith('%') or s.startswith('!'):
                        continue
                    cleaned_lines.append(line)
                
                code = "".join(cleaned_lines)
                # Force SAVE_AS_HTML to False during automated headless execution
                # to prevent nested/redundant jupyter nbconvert subprocess calls.
                code = re.sub(r'SAVE_AS_HTML\s*=\s*True', 'SAVE_AS_HTML = False', code)
                # Prevent plt.show() from destroying the figure in Agg backend before get_fignums captures it
                code = re.sub(r'plt\.show\(\)', 'pass', code)
                if not code.strip() or ('nbconvert' in code and 'subprocess' in code):
                    cell['outputs'] = []
                    new_cells.append(cell)
                    continue
                
                stdout_io = io.StringIO()
                stderr_io = io.StringIO()
                cell_outputs = []
                
                # Clear active matplotlib figures before executing the cell
                plt.clf()
                plt.close('all')
                display_events.clear()
                
                total_cells = len(nb.get('cells', []))
                snippet = code.split('\n')[0][:50] + "..." if code else ""
                print(f"  ⏳ Executing cell {idx+1}/{total_cells} ({snippet})", flush=True)
                
                try:
                    with redirect_stdout(stdout_io), redirect_stderr(stderr_io):
                        exec(code, global_dict)
                    
                    stdout_str = stdout_io.getvalue()
                    stderr_str = stderr_io.getvalue()
                    
                    print(f"  ✅ Finished cell {idx+1}/{total_cells}", flush=True)
                    
                    if stdout_str:
                        cell_outputs.append({
                            'output_type': 'stream',
                            'name': 'stdout',
                            'text': [stdout_str]
                        })
                    if stderr_str:
                        cell_outputs.append({
                            'output_type': 'stream',
                            'name': 'stderr',
                            'text': [stderr_str]
                        })
                        
                    # Check if any new figures were plotted in this cell
                    figs = plt.get_fignums()
                    for fignum in figs:
                        fig = plt.figure(fignum)
                        buf = io.BytesIO()
                        fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
                        buf.seek(0)
                        b64 = base64.b64encode(buf.read()).decode('utf-8')
                        cell_outputs.append({
                            'output_type': 'display_data',
                            'data': {
                                'image/png': b64,
                                'text/plain': ['<matplotlib.figure.Figure>']
                            },
                            'metadata': {}
                        })
                        plt.close(fig)
                        
                    for html_content in display_events:
                        cell_outputs.append({
                            "output_type": "display_data",
                            "data": html_content,
                            "metadata": {}
                        })
                        
                except Exception as e:
                    import traceback
                    tb_lines = traceback.format_exception(type(e), e, e.__traceback__)
                    cell_outputs.append({
                        'output_type': 'error',
                        'ename': type(e).__name__,
                        'evalue': str(e),
                        'traceback': tb_lines
                    })
                    print(f"❌ Error executing cell {idx}: {e}")
                    # Print full traceback to console
                    traceback.print_exc()
                    raise e
                    
                # Preserve existing base64 images only if the execution didn't produce new images,
                # to prevent endlessly stacking plots on multiple runs while preserving injected plots.
                new_has_image = any(out.get('output_type') == 'display_data' and 'image/png' in out.get('data', {}) for out in cell_outputs)
                if not new_has_image:
                    existing_images = [out for out in cell.get('outputs', []) if out.get('output_type') == 'display_data' and 'image/png' in out.get('data', {})]
                    cell['outputs'] = cell_outputs + existing_images
                else:
                    cell['outputs'] = cell_outputs
                
            new_cells.append(cell)
            
        nb['cells'] = new_cells
        
        # Save the fully updated notebook containing new plots and standardized data states
        with open(notebook_path, 'w', encoding='utf-8') as f:
            json.dump(nb, f, indent=1)
            
        # Check if the notebook dynamically generated an h5ad_path we should base the HTML output on
        if 'cancer_cellxgene_integration' in os.path.basename(notebook_path) and 'h5ad_path' in global_dict:
            resolved_h5ad = global_dict['h5ad_path']
            if not os.path.isabs(resolved_h5ad):
                resolved_h5ad = os.path.abspath(os.path.join(notebook_dir, resolved_h5ad))
            else:
                resolved_h5ad = os.path.abspath(resolved_h5ad)
            html_path = resolved_h5ad.replace('.h5ad', '.html')
            print(f"💡 Dynamic h5ad_path detected: {global_dict['h5ad_path']}")
            print(f"   Exporting HTML report to dynamic path: {html_path}")
            
        # Export report
        # First, try to use standard jupyter nbconvert subprocess to get native Jupyter-style HTML
        import subprocess
        print(f"Trying to export {os.path.basename(notebook_path)} using standard jupyter nbconvert...")
        # Construct path to the active virtual environment's jupyter binary
        jupyter_bin = os.path.join(os.path.dirname(sys.executable), 'jupyter')
        if not os.path.exists(jupyter_bin):
            jupyter_bin = 'jupyter' # fallback to system path
            
        cmd_html = [
            jupyter_bin, 'nbconvert', '--to', 'html', 
            notebook_path, '--output', html_path
        ]

        res_html = subprocess.run(cmd_html, capture_output=True, text=True)
        if res_html.returncode == 0:
            print(f"🎉 SUCCESS: Standard jupyter nbconvert successfully exported the report!")
        else:
            print(f"⚠️ Standard jupyter nbconvert export failed (likely due to environment/sandbox restrictions).")
            print(f"   Falling back to clean, standard-styled manual HTML export...")
            export_to_gorgeous_html(nb, html_path, title_text)
            print(f"🎉 Successfully saved and exported {os.path.basename(notebook_path)} to {html_path}\n")
        
    finally:
        os.chdir(original_dir)

if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    import subprocess
    
    print("Generating dynamic downstream notebooks from configurations...")
    subprocess.run([sys.executable, os.path.join(script_dir, "compute_pan_cancer_meta.py")], check=True)
    subprocess.run([sys.executable, os.path.join(script_dir, "generate_combined_pan_cancer_notebook.py")], check=True)
    subprocess.run([sys.executable, os.path.join(script_dir, "generate_predictive_notebook.py")], check=True)
    
    # 2. Predictive Signature Biomarker
    try:
        predictive_nb = os.path.join(script_dir, "predictive_signature_biomarker.ipynb")
        import sys
        sys.path.insert(0, script_dir)
        from pan_cancer_config import ANALYSIS_SUFFIX
        predictive_html = os.path.join(os.path.dirname(script_dir), "output", "pan_cancer_meta_results", f"predictive_signature_biomarker{ANALYSIS_SUFFIX}.html")
        if os.path.exists(predictive_nb):
            execute_and_export(predictive_nb, predictive_html, "Predictive Signature Biomarker")
        else:
            sys.exit(f"CRITICAL ERROR: {predictive_nb} NOT FOUND. Run 'python scripts/generate_predictive_notebook.py' first (Phase 5).")
    except Exception as e:
        sys.exit(f"CRITICAL ERROR - Failed to execute {predictive_nb}: {e}")


    # 1. Pan-Cancer Meta Analysis
    try:
        pancancer_nb = os.path.join(script_dir, "pan_cancer_meta_analysis.ipynb")
        pancancer_html = os.path.join(os.path.dirname(script_dir), "output", "pan_cancer_meta_results", "pan_cancer_meta_analysis_report.html")
        execute_and_export(pancancer_nb, pancancer_html, "Pan-Cancer Conserved Metabolic Signature Analysis")
    except Exception as e:
        sys.exit(f"CRITICAL ERROR - Failed to execute {pancancer_nb}: {e}")


    # 3. Druggability Axis Analysis
    try:
        druggability_nb = os.path.join(script_dir, "druggability_axis_analysis.ipynb")
        druggability_html = os.path.join(os.path.dirname(script_dir), "output", "druggability", f"druggability_axis{ANALYSIS_SUFFIX}.html")
        if os.path.exists(druggability_nb):
            print("Computing druggability targets for conserved genes...")
            subprocess.run([sys.executable, os.path.join(script_dir, "compute_druggability_targets.py")], check=True)
            execute_and_export(druggability_nb, druggability_html, "Druggability Axis Analysis")
        else:
            sys.exit(f"CRITICAL ERROR: {druggability_nb} NOT FOUND. Druggability analysis will be skipped.")
    except Exception as e:
        sys.exit(f"CRITICAL ERROR - Failed to execute {druggability_nb}: {e}")

    # 4. Visium Spatial Validation
    try:
        visium_nb = os.path.join(script_dir, "visium_spatial_validation.ipynb")
        visium_html = os.path.join(os.path.dirname(script_dir), "output", "visium_spatial_validation_report.html")
        if os.path.exists(visium_nb):
            execute_and_export(visium_nb, visium_html, "Ovarian Visium Spatial Transcriptomics Validation")
        else:
            sys.exit(f"CRITICAL ERROR: {visium_nb} NOT FOUND. Run 'python scripts/generate_visium_notebook.py' first (Phase 5).")
    except Exception as e:
        sys.exit(f"CRITICAL ERROR - Failed to execute {visium_nb}: {e}")

    # 5. Deep-Dive Conserved Metab Gene Sig
    try:
        deepdive_nb = os.path.join(script_dir, "deepdive_conserved_metabGeneSig.ipynb")
        deepdive_html = os.path.join(os.path.dirname(script_dir), "output", "deepdive_conserved_metabGeneSig", "deepdive_conserved_metabGeneSig_report.html")
        if os.path.exists(deepdive_nb):
            execute_and_export(deepdive_nb, deepdive_html, "Deep-Dive Conserved Metastatic Metabolic Signature")
        else:
            sys.exit(f"CRITICAL ERROR: {deepdive_nb} NOT FOUND. Run 'python scripts/generate_nb1.py' first (Phase 5).")
    except Exception as e:
        sys.exit(f"CRITICAL ERROR - Failed to execute {deepdive_nb}: {e}")

    # 6. Serotonin Axis Spatial Mapping
    try:
        serotonin_nb = os.path.join(script_dir, "serotonin_axis_spatial_mapping.ipynb")
        serotonin_html = os.path.join(os.path.dirname(script_dir), "output", "serotonin_axis_spatial_mapping_report.html")
        if os.path.exists(serotonin_nb):
            execute_and_export(serotonin_nb, serotonin_html, "Serotonin Axis Spatial Mapping")
        else:
            sys.exit(f"CRITICAL ERROR: {serotonin_nb} NOT FOUND. Run 'python scripts/generate_serotonin_notebook.py' first (Phase 5).")
    except Exception as e:
        sys.exit(f"CRITICAL ERROR - Failed to execute {serotonin_nb}: {e}")

    # 7. Ovarian Serotonin Immune Evasion
    try:
        print("Preparing data and generating Ovarian Serotonin Immune Evasion notebook...")
        if os.path.exists(os.path.join(script_dir, "compute_metastatic_immune_evasion.py")):
            subprocess.run([sys.executable, os.path.join(script_dir, "compute_metastatic_immune_evasion.py")], check=True)
        if os.path.exists(os.path.join(script_dir, "verify_spatial_immune_evasion.py")):
            subprocess.run([sys.executable, os.path.join(script_dir, "verify_spatial_immune_evasion.py")], check=True)
        if os.path.exists(os.path.join(script_dir, "generate_immune_evasion_notebook.py")):
            subprocess.run([sys.executable, os.path.join(script_dir, "generate_immune_evasion_notebook.py")], check=True)

        ov_serotonin_nb = os.path.join(script_dir, "ovarian_serotonin_immune_evasion.ipynb")
        ov_serotonin_html = os.path.join(os.path.dirname(script_dir), "output", "ovarian_serotonin_immune_evasion_report.html")
        if os.path.exists(ov_serotonin_nb):
            execute_and_export(ov_serotonin_nb, ov_serotonin_html, "Ovarian Serotonin Immune Evasion")
        else:
            sys.exit(f"CRITICAL ERROR: {ov_serotonin_nb} NOT FOUND. Run 'python scripts/generate_immune_evasion_notebook.py' first (Phase 5).")
    except Exception as e:
        sys.exit(f"CRITICAL ERROR - Failed to execute {ov_serotonin_nb}: {e}")

    # 8. Oxygen Tension Analysis
    try:
        print("Preparing data and generating Oxygen Tension notebook...")

        oxygen_nb = os.path.join(script_dir, "oxygen_tension_analysis.ipynb")
        oxygen_html = os.path.join(os.path.dirname(script_dir), "output", "oxygen_tension_analysis_report.html")
        if os.path.exists(oxygen_nb):
            execute_and_export(oxygen_nb, oxygen_html, "Oxygen Tension Analysis")
        else:
            sys.exit(f"CRITICAL ERROR: {oxygen_nb} NOT FOUND. No generator exists; this notebook must be created manually.")
    except Exception as e:
        sys.exit(f"CRITICAL ERROR - Failed to execute {oxygen_nb}: {e}")

    # 9. MITF Regulon Expansion
    try:
        print("Preparing data for MITF Regulon Expansion...")
        if os.path.exists(os.path.join(script_dir, "compute_mitf_regulon.py")):
            subprocess.run([sys.executable, os.path.join(script_dir, "compute_mitf_regulon.py")], check=True)

        mitf_nb = os.path.join(script_dir, "mitf_regulon_expansion.ipynb")
        mitf_html = os.path.join(os.path.dirname(script_dir), "output", "mitf_regulon_expansion_report.html")
        if os.path.exists(mitf_nb):
            execute_and_export(mitf_nb, mitf_html, "MITF Regulon Expansion")
        else:
            sys.exit(f"CRITICAL ERROR: {mitf_nb} NOT FOUND. Run 'python scripts/generate_nb2.py' first (Phase 5).")
    except Exception as e:
        sys.exit(f"CRITICAL ERROR - Failed to execute {mitf_nb}: {e}")

    # 10. ML Prognostic Classifier
    try:
        print("Preparing ML Prognostic Classifier notebooks...")
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
                else:
                    print(f"Skipping {prefix.upper()} for {cancer_name} as {expr_file} is missing.")
            
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
        
    except Exception as e:
        sys.exit(f"CRITICAL ERROR - Failed to execute ML Prognostic Classifiers: {e}")

    # 11. CAMP Pan-Cancer Integration
    try:
        print("Preparing CAMP Pan-Cancer Integration notebook...")
        camp_nb = os.path.join(script_dir, "camp_pancancer_integration.ipynb")
        camp_html = os.path.join(os.path.dirname(script_dir), "output", "camp_pancancer_integration_report.html")
        if os.path.exists(camp_nb):
            execute_and_export(camp_nb, camp_html, "CAMP Pan-Cancer Integration")
        else:
            sys.exit(f"CRITICAL ERROR: {camp_nb} NOT FOUND. Run 'python scripts/create_camp_notebook.py' first (Phase 5).")
    except Exception as e:
        sys.exit(f"CRITICAL ERROR - Failed to execute {camp_nb}: {e}")

    # 12. Master Regulator Analysis
    try:
        print("Preparing Master Regulator Analysis notebook...")
        mr_nb = os.path.join(script_dir, "master_regulator_analysis.ipynb")
        mr_html = os.path.join(os.path.dirname(script_dir), "output", "master_regulator_analysis_report.html")
        if os.path.exists(mr_nb):
            execute_and_export(mr_nb, mr_html, "Master Regulator Analysis")
        else:
            sys.exit(f"CRITICAL ERROR: {mr_nb} NOT FOUND. Run 'python scripts/generate_master_regulator_notebook.py' first (Phase 5).")
    except Exception as e:
        sys.exit(f"CRITICAL ERROR - Failed to execute {mr_nb}: {e}")
