import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from nbconvert import HTMLExporter
import os

notebook_filename = 'oxygen_tension_analysis.ipynb'
output_filename = 'oxygen_tension_analysis_100k.html'

print(f"Executing and exporting {notebook_filename}...")
with open(notebook_filename, 'r', encoding='utf-8') as f:
    nb = nbformat.read(f, as_version=4)

ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
try:
    ep.preprocess(nb, {'metadata': {'path': './'}})
except Exception as e:
    print(f"Error executing notebook: {e}")

html_exporter = HTMLExporter()
html_exporter.theme = 'light'
body, resources = html_exporter.from_notebook_node(nb)

with open(output_filename, 'w', encoding='utf-8') as f:
    f.write(body)

print(f"Successfully exported to {output_filename}")
