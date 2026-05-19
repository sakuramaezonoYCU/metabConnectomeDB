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
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category=UserWarning, message=".*FigureCanvasAgg is non-interactive.*")


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
            html_lines.append(f'<h3 style="color:#93c5fd;font-size:1.4rem;margin-top:20px;font-weight:600;">{stripped[3:].strip()}</h3>')
        elif stripped.startswith('##'):
            if in_list:
                html_lines.append('</ul>')
                in_list = False
            html_lines.append(f'<h2 style="color:#60a5fa;font-size:1.8rem;margin-top:30px;border-bottom:1px solid #334155;padding-bottom:8px;font-weight:600;">{stripped[2:].strip()}</h2>')
        elif stripped.startswith('#'):
            if in_list:
                html_lines.append('</ul>')
                in_list = False
            html_lines.append(f'<h1 style="color:#3b82f6;font-size:2.2rem;margin-top:40px;margin-bottom:15px;font-weight:700;">{stripped[1:].strip()}</h1>')
        # Lists
        elif stripped.startswith('-') or stripped.startswith('*'):
            if not in_list:
                html_lines.append('<ul style="margin-left:20px;padding-left:10px;line-height:1.7;color:#d1d5db;">')
                in_list = True
            # Parse list item bold and inline code
            item_text = stripped[1:].strip()
            item_text = re.sub(r'\*\*(.*?)\*\*', r'<strong>\1</strong>', item_text)
            item_text = re.sub(r'`(.*?)`', r'<code style="background:#1e293b;padding:2px 5px;border-radius:4px;font-family:\'Fira Code\',monospace;font-size:0.9em;color:#f43f5e;">\1</code>', item_text)
            html_lines.append(f'<li style="margin-bottom:6px;">{item_text}</li>')
        else:
            if in_list:
                html_lines.append('</ul>')
                in_list = False
            
            if stripped:
                # Regular paragraph text
                item_text = line
                item_text = re.sub(r'\*\*(.*?)\*\*', r'<strong>\1</strong>', item_text)
                item_text = re.sub(r'`(.*?)`', r'<code style="background:#1e293b;padding:2px 5px;border-radius:4px;font-family:\'Fira Code\',monospace;font-size:0.9em;color:#f43f5e;">\1</code>', item_text)
                html_lines.append(f'<p style="margin-bottom:15px;color:#d1d5db;line-height:1.7;">{item_text}</p>')
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
    <title>{title_text} - Premium Report</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=Fira+Code:wght@400;500&display=swap" rel="stylesheet">
    <style>
        body {{
            font-family: 'Inter', sans-serif;
            line-height: 1.6;
            color: #e2e8f0;
            background-color: #0f172a;
            margin: 0;
            padding: 0;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 40px 20px;
        }}
        .header {{
            border-bottom: 1px solid #334155;
            padding-bottom: 30px;
            margin-bottom: 40px;
            background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
            padding: 30px;
            border-radius: 16px;
            border: 1px solid #334155;
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.3);
        }}
        .header h1 {{
            font-size: 2.8rem;
            color: #3b82f6;
            margin: 0 0 12px 0;
            font-weight: 800;
            letter-spacing: -0.025em;
        }}
        .header p {{
            font-size: 1.2rem;
            color: #94a3b8;
            margin: 0 0 15px 0;
        }}
        .badge {{
            display: inline-block;
            background: rgba(59, 130, 246, 0.2);
            color: #60a5fa;
            padding: 6px 14px;
            border-radius: 9999px;
            font-size: 0.85rem;
            font-weight: 600;
            border: 1px solid rgba(59, 130, 246, 0.3);
        }}
        .cell {{
            background: #1e293b;
            border-radius: 16px;
            padding: 28px;
            margin-bottom: 35px;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);
            border: 1px solid #334155;
            transition: transform 0.2s ease, box-shadow 0.2s ease;
        }}
        .cell:hover {{
            box-shadow: 0 20px 25px -5px rgba(0, 0, 0, 0.2), 0 10px 10px -5px rgba(0, 0, 0, 0.04);
        }}
        .cell.markdown {{
            background: transparent;
            border: none;
            box-shadow: none;
            padding: 0 12px;
            margin-bottom: 25px;
        }}
        .cell.markdown:hover {{
            transform: none;
            box-shadow: none;
        }}
        pre.code-block {{
            background: #090d16;
            padding: 20px;
            border-radius: 10px;
            overflow-x: auto;
            border: 1px solid #1e293b;
            margin: 0 0 20px 0;
            position: relative;
        }}
        pre.code-block::before {{
            content: "PYTHON";
            position: absolute;
            top: 8px;
            right: 12px;
            font-size: 0.7rem;
            color: #475569;
            font-weight: 700;
            letter-spacing: 0.05em;
        }}
        code.python {{
            font-family: 'Fira Code', monospace;
            font-size: 0.95rem;
            color: #38bdf8;
        }}
        .output-container {{
            margin-top: 20px;
            border-top: 1px solid #334155;
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
            background: #020617;
            padding: 18px;
            border-radius: 8px;
            font-family: 'Fira Code', monospace;
            font-size: 0.9rem;
            color: #10b981;
            white-space: pre-wrap;
            border-left: 4px solid #10b981;
            margin-bottom: 15px;
            box-shadow: inset 0 2px 4px 0 rgba(0, 0, 0, 0.06);
        }}
        .output-stderr {{
            background: #1c0e15;
            color: #ef4444;
            border-left-color: #ef4444;
        }}
        .output-image {{
            text-align: center;
            margin: 25px 0;
            padding: 15px;
            background: #0f172a;
            border-radius: 12px;
            border: 1px solid #1e293b;
        }}
        .output-image img {{
            max-width: 100%;
            border-radius: 8px;
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.4);
        }}
        .output-html {{
            overflow-x: auto;
            margin: 15px 0;
            background: #0f172a;
            padding: 10px;
            border-radius: 8px;
            border: 1px solid #334155;
        }}
        table.dataframe {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9rem;
            min-width: 400px;
            color: #cbd5e1;
        }}
        table.dataframe th, table.dataframe td {{
            padding: 12px 16px;
            text-align: left;
            border-bottom: 1px solid #1e293b;
        }}
        table.dataframe th {{
            background-color: #1e293b;
            color: #60a5fa;
            font-weight: 600;
            border-bottom: 2px solid #334155;
        }}
        table.dataframe tr:hover {{
            background-color: #1e293b;
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
    print(f"🚀 Running execution and report generation for {os.path.basename(notebook_path)}...")
    
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = json.load(f)
        
    global_dict = {
        '__name__': '__main__',
        'plt': plt,
        'SAVE_AS_HTML': False, # Avoid nested nbconvert subprocess calls!
    }
    
    new_cells = []
    notebook_dir = os.path.dirname(os.path.abspath(notebook_path))
    original_dir = os.getcwd()
    os.chdir(notebook_dir)
    
    # Prepend venv/bin to PATH so execution contexts locate files properly
    venv_bin = os.path.join(os.path.dirname(notebook_dir), 'venv', 'bin')
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
                if not code.strip():
                    cell['outputs'] = []
                    new_cells.append(cell)
                    continue
                
                stdout_io = io.StringIO()
                stderr_io = io.StringIO()
                cell_outputs = []
                
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
        
        # Save the fully updated notebook containing new plots and standardized data states
        with open(notebook_path, 'w', encoding='utf-8') as f:
            json.dump(nb, f, indent=1)
            
        # Export report to gorgeous, style-preserved premium HTML
        export_to_gorgeous_html(nb, html_path, title_text)
        print(f"🎉 Successfully saved and exported {os.path.basename(notebook_path)} to {html_path}\n")
        
    finally:
        os.chdir(original_dir)

if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    unique_metab_nb = os.path.join(script_dir, "unique_metab_data_exploration.ipynb")
    unique_metab_html = os.path.join(os.path.dirname(script_dir), "output", "unique_metab_data_exploration_full_report.html")
    execute_and_export(unique_metab_nb, unique_metab_html, "Metabolite Catalog Exploration")
    
    pair_analysis_nb = os.path.join(script_dir, "metab_targetPair_analysis.ipynb")
    pair_analysis_html = os.path.join(os.path.dirname(script_dir), "output", "metab_targetPair_analysis_full_report.html")
    execute_and_export(pair_analysis_nb, pair_analysis_html, "Metabolite-Target Interaction Pair Analysis")
