import json
import os
import sys
import io
import base64
import re
import glob
from contextlib import redirect_stdout, redirect_stderr

# Configure Matplotlib for headless execution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
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

def execute_and_export(notebook_path, html_path, title_text, inject_globals=None, skip_if_exists=False):
    notebook_path = os.path.abspath(notebook_path)
    html_path = os.path.abspath(html_path)
    
    if skip_if_exists and os.path.exists(html_path):
        print(f"⏭️ Skipping {os.path.basename(notebook_path)} as {os.path.basename(html_path)} already exists.")
        
        # If this is the cellxgene notebook, we need to recover the h5ad_path
        if 'cancer_cellxgene_integration' in notebook_path:
            # First try parsing the HTML
            try:
                with open(html_path, 'r', encoding='utf-8') as hf:
                    html_content = hf.read()
                    # Try to find the cache path in the parameters block
                    match = re.search(r'Cache path:\s*([^<\n\s]+(?:\.h5ad))', html_content)
                    if not match:
                        # Fallback to the saving message
                        match = re.search(r'Saving downloaded AnnData to:\s*([^<\n\s]+(?:\.h5ad))', html_content)
                    if match:
                        return match.group(1).strip()
            except Exception as e:
                pass
                
            # Fallback to finding the .h5ad file in the output directory
            # Sort by modification time to get the most recent one
            output_dir = os.path.dirname(html_path)
            h5ads = glob.glob(os.path.join(output_dir, '*.h5ad'))
            if h5ads:
                h5ads.sort(key=os.path.getmtime, reverse=True)
                return h5ads[0]
                
        # For other notebooks, just return whatever h5ad_path was passed in
        return inject_globals.get('h5ad_path', None) if inject_globals else None

    print(f"🚀 Running execution and report generation for {os.path.basename(notebook_path)}...")
    
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = json.load(f)
        
    try:
        from IPython.display import display
    except ImportError:
        display = print

    def fake_display(obj):
        if hasattr(obj, "_repr_html_"):
            display_events.append({"text/html": [obj._repr_html_()]})
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
    if inject_globals:
        global_dict.update(inject_globals)
        
    # Monkey-patch matplotlib.text.Text to prevent upsetplot 0-dimensional array crash
    patch_code = """
import matplotlib.text
original_text_init = matplotlib.text.Text.__init__
original_set_position = matplotlib.text.Text.set_position

def patched_text_init(self, x=0, y=0, text='', **kwargs):
    try: x = float(x)
    except Exception:
        try: x = x.item()
        except Exception: pass
    try: y = float(y)
    except Exception:
        try: y = y.item()
        except Exception: pass
    original_text_init(self, x=x, y=y, text=text, **kwargs)

def patched_set_position(self, pos):
    x, y = pos
    try: x = float(x)
    except Exception:
        try: x = x.item()
        except Exception: pass
    try: y = float(y)
    except Exception:
        try: y = y.item()
        except Exception: pass
    original_set_position(self, (x, y))

matplotlib.text.Text.__init__ = patched_text_init
matplotlib.text.Text.set_position = patched_set_position


# Headless bypass for Plotly fig.show() to prevent browser blocks
try:
    import plotly.graph_objects as go
    go.Figure.show = lambda *args, **kwargs: None
except Exception:
    pass
"""
    exec(patch_code, global_dict)
    
    # ── Sampling fix ─────────────────────────────────────────────────────────
    # The notebook's original stratified sampler divides CAP equally across all
    # tissues regardless of their available cell counts, causing primary tissues
    # (which typically have fewer cells) to be overrepresented.
    # We replace it at execution time with a proportional sampler.
    BIASED_SAMPLER   = "obs_df.groupby('tissue_general', group_keys=False).apply(lambda x: x.sample(n=min(len(x), CAP // len(obs_df['tissue_general'].unique())), random_state=42))"
    BALANCED_SAMPLER = (
        "obs_df.groupby('tissue_general', group_keys=False).apply("
        "lambda x: x.sample("
        "n=min(len(x), max(1, round(CAP * len(x) / len(obs_df)))), "
        "random_state=42))"
    )
    # ─────────────────────────────────────────────────────────────────────────
    
    new_cells = []
    notebook_dir = os.path.dirname(os.path.abspath(notebook_path))
    original_dir = os.getcwd()
    
    # Prepend a Markdown cell showing the injected globals so the HTML report clearly shows the configuration
    # Prepend a Markdown cell showing the injected globals so the HTML report clearly shows the configuration
    if inject_globals:
        md_text = "### ⚙️ Pipeline Injected Configuration\nThis notebook was automatically executed by the pipeline with the following parameters:\n\n"
        for k, v in inject_globals.items():
            if k not in ['plt', 'display', '__name__', 'SAVE_AS_HTML']:
                md_text += f"- **{k}**: `{v}`\n"
        new_cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [md_text]
        })

    os.chdir(notebook_dir)
    
    # Prepend active virtual environment bin directory to PATH so execution contexts locate files properly
    venv_bin = os.path.dirname(sys.executable)
    if os.path.exists(venv_bin):
        os.environ['PATH'] = venv_bin + os.pathsep + os.environ.get('PATH', '')
        
    try:
        import IPython.display
        
        # Monkey-patch IPython.display's members so imported display/HTML work perfectly
        class MockHTML:
            def __init__(self, data):
                self.data = data
            def _repr_html_(self):
                return self.data
            def __repr__(self):
                return self.data
        IPython.display.HTML = MockHTML
        
        original_ipython_display = IPython.display.display
        IPython.display.display = global_dict['display']
        
        try:
            from tqdm import tqdm
            cell_iterator = tqdm(enumerate(nb.get('cells', [])), total=len(nb.get('cells', [])), desc=f"Executing {os.path.basename(notebook_path)}", leave=False)
        except ImportError:
            cell_iterator = enumerate(nb.get('cells', []))
            
        for idx, cell in cell_iterator:
            if cell.get('cell_type') == 'code':
                source_lines = cell.get('source', [])
                
                code = "".join(source_lines)
                original_code = "".join(source_lines)
                
                # Strip out IPython magics (e.g. %matplotlib inline, !pip install)
                # because exec() does not support them and will throw SyntaxError.
                code = re.sub(r'^\s*[%!].*$', '', code, flags=re.MULTILINE)
                
                # Force SAVE_AS_HTML to False during automated headless execution
                # to prevent nested/redundant jupyter nbconvert subprocess calls.
                code = re.sub(r'SAVE_AS_HTML\s*=\s*True', 'SAVE_AS_HTML = False', code)
                original_code = re.sub(r'SAVE_AS_HTML\s*=\s*True', 'SAVE_AS_HTML = False', original_code)
                
                # Prevent plt.show() from closing the figure before we can capture it
                code = re.sub(r'plt\.show\(.*?\)', 'pass', code)
                
                # Dynamically inject the DISEASE_FILTER and TISSUE_FILTER if they exist in global_dict
                if 'DISEASE_FILTER' in global_dict:
                    df_json = json.dumps(global_dict['DISEASE_FILTER'])
                    code = re.sub(r'^DISEASE_FILTER\s*=\s*\[[^\]]*\]', f'DISEASE_FILTER = {df_json}', code, flags=re.MULTILINE)
                    original_code = re.sub(r'^DISEASE_FILTER\s*=\s*\[[^\]]*\]', f'DISEASE_FILTER = {df_json}', original_code, flags=re.MULTILINE)
                
                if 'TISSUE_FILTER' in global_dict:
                    if global_dict['TISSUE_FILTER'] is None:
                        code = re.sub(r'^TISSUE_FILTER\s*=\s*\[[^\]]*\]|^TISSUE_FILTER\s*=\s*None', 'TISSUE_FILTER = None', code, flags=re.MULTILINE)
                        original_code = re.sub(r'^TISSUE_FILTER\s*=\s*\[[^\]]*\]|^TISSUE_FILTER\s*=\s*None', 'TISSUE_FILTER = None', original_code, flags=re.MULTILINE)
                    else:
                        tf_json = json.dumps(global_dict['TISSUE_FILTER'])
                        code = re.sub(r'^TISSUE_FILTER\s*=\s*\[[^\]]*\]|^TISSUE_FILTER\s*=\s*None', f'TISSUE_FILTER = {tf_json}', code, flags=re.MULTILINE)
                        original_code = re.sub(r'^TISSUE_FILTER\s*=\s*\[[^\]]*\]|^TISSUE_FILTER\s*=\s*None', f'TISSUE_FILTER = {tf_json}', original_code, flags=re.MULTILINE)
                
                if 'h5ad_path' in global_dict:
                    h5_json = json.dumps(global_dict['h5ad_path'])
                    code = re.sub(r'^h5ad_path\s*=.*$', f'h5ad_path = {h5_json}', code, flags=re.MULTILINE)
                    original_code = re.sub(r'^h5ad_path\s*=.*$', f'h5ad_path = {h5_json}', original_code, flags=re.MULTILINE)

                    # If this cell uses the old separate primary/meta h5ad loading pattern,
                    # replace it entirely with unified h5ad loading from h5ad_path.
                    # This avoids needing to correctly reconstruct PRIMARY_PREFIX/META_PREFIX filenames.
                    if 'sc.read_h5ad(primary_h5ad_file)' in code or 'sc.read_h5ad(meta_h5ad_file)' in code:
                        pt = global_dict.get('PRIMARY_TISSUES', None)
                        if pt is None:
                            raise RuntimeError("CRITICAL ERROR: PRIMARY_TISSUES must be injected into the globals dict.")
                        pt_json = json.dumps(pt)
                        code = f"""# Pipeline mode: load unified h5ad and split by PRIMARY_TISSUES
import scanpy as sc, numpy as np
adata = sc.read_h5ad(h5ad_path)
primary_list = globals().get('PRIMARY_TISSUES', {pt_json})
primary_name = primary_list[0].capitalize() if primary_list else 'Primary'
adata.obs['site'] = np.where(adata.obs['tissue_general'].isin(primary_list), primary_name, 'Metastasis')
adata_meta = adata[adata.obs['site'] == 'Metastasis'].copy()
adata_prim = adata[adata.obs['site'] == primary_name].copy()
print(f"Loaded combined dataset: {{adata.shape}}")
print(adata.obs['site'].value_counts())
"""
                    
                if 'PRIMARY_TISSUES' in global_dict:
                    pt_json = json.dumps(global_dict['PRIMARY_TISSUES'])
                    code = re.sub(r'^PRIMARY_TISSUES\s*=.*$', f'PRIMARY_TISSUES = {pt_json}', code, flags=re.MULTILINE)
                    original_code = re.sub(r'^PRIMARY_TISSUES\s*=.*$', f'PRIMARY_TISSUES = {pt_json}', original_code, flags=re.MULTILINE)
                    
                if 'CANCER_TYPE_NAME' in global_dict:
                    ctn_json = json.dumps(global_dict['CANCER_TYPE_NAME'])
                    code = re.sub(r'^CANCER_TYPE_NAME\s*=.*$', f'CANCER_TYPE_NAME = {ctn_json}', code, flags=re.MULTILINE)
                    original_code = re.sub(r'^CANCER_TYPE_NAME\s*=.*$', f'CANCER_TYPE_NAME = {ctn_json}', original_code, flags=re.MULTILINE)

                if 'PRIMARY_PREFIX' in global_dict:
                    pp_json = json.dumps(global_dict['PRIMARY_PREFIX'])
                    code = re.sub(r'^PRIMARY_PREFIX\s*=.*$', f'PRIMARY_PREFIX = {pp_json}', code, flags=re.MULTILINE)
                    original_code = re.sub(r'^PRIMARY_PREFIX\s*=.*$', f'PRIMARY_PREFIX = {pp_json}', original_code, flags=re.MULTILINE)
                    
                if 'META_PREFIX' in global_dict:
                    mp_json = json.dumps(global_dict['META_PREFIX'])
                    code = re.sub(r'^META_PREFIX\s*=.*$', f'META_PREFIX = {mp_json}', code, flags=re.MULTILINE)
                    original_code = re.sub(r'^META_PREFIX\s*=.*$', f'META_PREFIX = {mp_json}', original_code, flags=re.MULTILINE)

                if 'OUTPUT_DIR' in global_dict:
                    od_json = json.dumps(global_dict['OUTPUT_DIR'])
                    code = re.sub(r'^OUTPUT_DIR\s*=.*$', f'OUTPUT_DIR = {od_json}', code, flags=re.MULTILINE)
                    code = re.sub(r'^output_dir\s*=.*$', f'output_dir = {od_json}', code, flags=re.MULTILINE)
                    original_code = re.sub(r'^OUTPUT_DIR\s*=.*$', f'OUTPUT_DIR = {od_json}', original_code, flags=re.MULTILINE)
                    original_code = re.sub(r'^output_dir\s*=.*$', f'output_dir = {od_json}', original_code, flags=re.MULTILINE)
                    
                    # Intercept unique target pairs CSV and point it to the global output directory
                    global_output_dir = os.path.dirname(global_dict['OUTPUT_DIR'])
                    global_csv_path = os.path.join(global_output_dir, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')
                    hmdb_fname = 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv'
                    # Catch all variants: lowercase output_dir, uppercase OUTPUT_DIR, and BASE_DIR-relative
                    for expr in [
                        f"os.path.join(output_dir, '{hmdb_fname}')",
                        f'os.path.join(output_dir, "{hmdb_fname}")',
                        f"os.path.join(OUTPUT_DIR, '{hmdb_fname}')",
                        f'os.path.join(OUTPUT_DIR, "{hmdb_fname}")',
                        f"os.path.join(BASE_DIR, 'output', '{hmdb_fname}')",
                        f'os.path.join(BASE_DIR, "output", "{hmdb_fname}")',
                    ]:
                        code = code.replace(expr, f"'{global_csv_path}'")
                    # Brute-force: replace any remaining assignment to metab_db_path
                    code = re.sub(
                        r"metab_db_path\s*=\s*os\.path\.join\([^)]+human_database[^)]+\)",
                        f"metab_db_path = '{global_csv_path}'",
                        code
                    )

                if 'CAP' in global_dict:
                    if global_dict['CAP'] is None:
                        code = re.sub(r'^CAP\s*=\s*.*$', 'CAP = None', code, flags=re.MULTILINE)
                        original_code = re.sub(r'^CAP\s*=\s*.*$', 'CAP = None', original_code, flags=re.MULTILINE)
                    else:
                        code = re.sub(r'^CAP\s*=\s*.*$', f"CAP = {global_dict['CAP']}", code, flags=re.MULTILINE)
                        original_code = re.sub(r'^CAP\s*=\s*.*$', f"CAP = {global_dict['CAP']}", original_code, flags=re.MULTILINE)

                if 'cap_str' in global_dict:
                    cs_json = json.dumps(global_dict['cap_str'])
                    code = re.sub(r'^cap_str\s*=.*$', f'cap_str = {cs_json}', code, flags=re.MULTILINE)
                    original_code = re.sub(r'^cap_str\s*=.*$', f'cap_str = {cs_json}', original_code, flags=re.MULTILINE)

                # Apply proportional sampling fix: replaces equal-per-tissue splitter with
                # a proportional splitter that respects each tissue's true cell availability.
                if BIASED_SAMPLER in code:
                    code = code.replace(BIASED_SAMPLER, BALANCED_SAMPLER)
                    original_code = original_code.replace(BIASED_SAMPLER, BALANCED_SAMPLER)

                # Update the cell's source code so the HTML display shows the injected parameters!
                cell['source'] = [original_code]
                
                if not code.strip() or ('nbconvert' in code and 'subprocess' in code):
                    cell['outputs'] = []
                    new_cells.append(cell)
                    continue
                
                stdout_io = io.StringIO()
                stderr_io = io.StringIO()
                cell_outputs = []
                display_events = []
                global_dict["display_events"] = display_events
                
                # Clear active matplotlib figures before executing the cell
                plt.clf()
                plt.close('all')
                
                try:
                    with redirect_stdout(stdout_io), redirect_stderr(stderr_io):
                        exec(code, global_dict)
                    
                    stdout_str = stdout_io.getvalue()
                    stderr_str = stderr_io.getvalue()
                    
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
                        
                    for html_content in display_events:
                        cell_outputs.append({
                            "output_type": "display_data",
                            "data": html_content,
                            "metadata": {}
                        })
                    # Check if any new figures were plotted in this cell
                    figs = plt.get_fignums()
                    for fignum in figs:
                        fig = plt.figure(fignum)
                        buf = io.BytesIO()
                        try:
                            fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
                        except Exception as e1:
                            buf.seek(0)
                            buf.truncate()
                            try:
                                fig.savefig(buf, format='png', dpi=150)
                            except Exception as e2:
                                print(f"Warning: Failed to save figure {fignum} to PNG. Exception: {e2}")
                                plt.close(fig)
                                continue
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
                    
                cell['outputs'] = cell_outputs
                
            new_cells.append(cell)
            
        nb['cells'] = new_cells
        
        # Save the fully updated notebook to a TEMPORARY output path so we don't pollute the source notebooks
        # We save it in the same output directory as the html report
        temp_notebook_path = os.path.join(os.path.dirname(html_path), os.path.basename(notebook_path))
        with open(temp_notebook_path, 'w', encoding='utf-8') as f:
            json.dump(nb, f, indent=1)
            
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
            temp_notebook_path, '--output', html_path
        ]

        res_html = subprocess.run(cmd_html, capture_output=True, text=True)
        if res_html.returncode == 0:
            print(f"🎉 SUCCESS: Standard jupyter nbconvert successfully exported the report!")
        else:
            print(f"⚠️ Standard jupyter nbconvert export failed (likely due to environment/sandbox restrictions).")
            print(f"   Falling back to clean, standard-styled manual HTML export...")
            export_to_gorgeous_html(nb, html_path, title_text)
            print(f"🎉 Successfully saved and exported {os.path.basename(notebook_path)} to {html_path}\n")
            
        # Clean up the temporary executed notebook to avoid cluttering the output directory
        if os.path.exists(temp_notebook_path):
            try:
                os.remove(temp_notebook_path)
            except:
                pass

        # Try to restore original display function
        try:
            import IPython.display
            IPython.display.display = original_ipython_display
        except:
            pass
        
        return global_dict.get('h5ad_path', None)
        
    finally:
        os.chdir(original_dir)

if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    if len(sys.argv) < 5:
        print("Usage: python run_cancer_pipeline.py \"<CANCER_KEY>\" \"<DISEASE_FILTER>\" \"<TISSUE_FILTER>\" \"<PRIMARY_TISSUES>\"")
        print("Example: python run_cancer_pipeline.py colorectal \"colorectal cancer\" \"colon,liver,lung\" \"colon,large intestine\"")
        sys.exit(1)
        
    cancer_key = sys.argv[1]
    disease_filter_str = sys.argv[2]
    tissue_filter_str = sys.argv[3]
    primary_tissues_str = sys.argv[4]
    
    disease_filter = [d.strip() for d in disease_filter_str.split(',') if d.strip()]
    if tissue_filter_str.lower() == "all" or tissue_filter_str.lower() == "none" or tissue_filter_str == "":
        tissue_filter = None
    else:
        tissue_filter = [t.strip() for t in tissue_filter_str.split(',') if t.strip()]
        
    primary_tissues = [t.strip() for t in primary_tissues_str.split(',') if t.strip()]
        
    print(f"Initializing pipeline run for Disease: {disease_filter} | Tissue: {tissue_filter} | Primary: {primary_tissues}")
    
    cap_val = os.environ.get('CELLXGENE_CAP', None)
    if cap_val is None or cap_val.lower() == 'none' or cap_val == 'all':
        cap_str = "all"
    else:
        try:
            cap_int = int(cap_val)
            if cap_int >= 1000:
                cap_str = f"{cap_int//1000}k"
            else:
                cap_str = f"{cap_int}"
        except:
            cap_str = str(cap_val)
    
    cancer_name_safe = f"{cancer_key}_results"
    cancer_output_dir = os.path.abspath(os.path.join(os.path.dirname(script_dir), "output", cancer_name_safe))
    os.makedirs(cancer_output_dir, exist_ok=True)
    
    cellxgene_nb = os.path.join(script_dir, "cancer_cellxgene_integration.ipynb")
    # HTML name includes tissue slug so runs with the same cap but different tissue sets don't collide
    # (placeholder slug used here; actual slug is refined after h5ad is known below)
    cellxgene_html = os.path.join(cancer_output_dir, f"cancer_cellxgene_integration_{cap_str}.html")
    
    try:
        cap_int = int(cap_val)
    except:
        cap_int = None
        
    inject_globals = {
        'DISEASE_FILTER': disease_filter,
        'TISSUE_FILTER': tissue_filter,
        'OUTPUT_DIR': cancer_output_dir,
        'CAP': cap_int,
        'cancer_key': cancer_key
    }
    
    generated_h5ad = execute_and_export(cellxgene_nb, cellxgene_html, f"CellxGene Integration: {disease_filter_str}", inject_globals, skip_if_exists=True)
    
    if generated_h5ad:
        print(f"Running downstream scripts for generated h5ad: {generated_h5ad}")
        
        print(f"Ensuring LIANA target network is generated via patch script...")
        import subprocess
        subprocess.run([sys.executable, os.path.join(script_dir, "patch_liana_csvs.py")])
        
        import re
        disease_slug = re.sub(r'[^a-z0-9]+', '-', disease_filter_str.lower()).strip('-')
        primary_slug = "_".join(re.sub(r'[^a-z0-9]+', '-', t.lower()).strip('-') for t in primary_tissues)
        met_tissues = [t for t in tissue_filter if t not in primary_tissues] if tissue_filter else []
        met_slug = "_".join(re.sub(r'[^a-z0-9]+', '-', t.lower()).strip('-') for t in met_tissues) if met_tissues else "all"
        # primary_prefix = f"{primary_slug}_{cap_str}_whole_transcriptome_2025-11-08"
        # meta_prefix = f"{met_slug}_{cap_str}_whole_transcriptome_2025-11-08"
        primary_prefix = f"{primary_slug}_{cap_str}_whole_transcriptome_2025-11-08"
        meta_prefix = f"{met_slug}_{cap_str}_whole_transcriptome_2025-11-08"
        all_tissues = tissue_filter if tissue_filter else primary_tissues
        all_slug = "_".join(re.sub(r'[^a-z0-9]+', '-', t.lower()).strip('-') for t in all_tissues)
        combined_h5ad = generated_h5ad
        print(f"DEBUG tissue_filter={tissue_filter} primary_tissues={primary_tissues}")
        print(f"DEBUG combined_h5ad={combined_h5ad}")
        # Build the tissue slug from all_tissues to make output filenames unique per tissue combination
        pvsm_nb = os.path.join(script_dir, "primary_vs_metastasis_comparison.ipynb")
        pvsm_html = os.path.join(cancer_output_dir, f"primary_vs_metastasis_{all_slug}_{cap_str}.html")
        project_root = os.path.dirname(script_dir)
        metab_db_path = os.path.join(project_root, 'output', 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')
        # Update the cellxgene HTML to include the tissue slug too (now that we have all_slug)
        cellxgene_html = os.path.join(cancer_output_dir, f"cancer_cellxgene_integration_{all_slug}_{cap_str}.html")
        pvsm_globals = {
            'h5ad_path': combined_h5ad,
            'PRIMARY_TISSUES': primary_tissues,
            'CANCER_TYPE_NAME': cancer_name_safe,
            'OUTPUT_DIR': cancer_output_dir,
            'PRIMARY_PREFIX': primary_prefix,
            'META_PREFIX': meta_prefix,
            'metab_db_path': metab_db_path,
            'cap_str': cap_str,  # Passed so notebooks can include it in output filenames
            'cancer_key': cancer_key
        }
        execute_and_export(pvsm_nb, pvsm_html, f"Primary vs Metastasis: {disease_filter_str}", pvsm_globals, skip_if_exists=True)
        
        orphan_nb = os.path.join(script_dir, "orphan_metabolic_immune_evasion.ipynb")
        orphan_html = os.path.join(cancer_output_dir, f"orphan_immune_{all_slug}_{cap_str}.html")
        # Orphan CSV is uniquely named per tissue-combination and cell count to prevent overwrites
        orphan_csv = os.path.join(cancer_output_dir, f"immune_evasion_orphan_metabolic_candidates_{all_slug}_{cap_str}.csv")
        orphan_globals = {
            'h5ad_path': pvsm_globals['h5ad_path'],
            'PRIMARY_TISSUES': primary_tissues,
            'output_csv': orphan_csv,
            'CANCER_TYPE_NAME': cancer_name_safe,
            'OUTPUT_DIR': cancer_output_dir,
            'PRIMARY_PREFIX': primary_prefix,
            'META_PREFIX': meta_prefix,
            'cancer_key': cancer_key
        }
        execute_and_export(orphan_nb, orphan_html, f"Orphan Metabolic Targets: {disease_filter_str}", orphan_globals, skip_if_exists=True)
