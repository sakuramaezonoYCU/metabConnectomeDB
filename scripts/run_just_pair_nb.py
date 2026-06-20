import os
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from execute_and_export_notebooks import execute_and_export

pair_analysis_nb = os.path.join(script_dir, "metab_targetPair_analysis.ipynb")
pair_analysis_html = os.path.join(os.path.dirname(script_dir), "output", "metab_targetPair_analysis_full_report.html")

print("Executing target pair notebook...")
execute_and_export(pair_analysis_nb, pair_analysis_html, "Metabolite-Target Interaction Pair Analysis")
print("Done.")
