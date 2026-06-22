import sys
import os

sys.path.append(os.path.join(os.getcwd(), 'scripts'))
from execute_and_export_notebooks import execute_and_export

if __name__ == '__main__':
    # Fix matplotlib cache dir issue
    os.environ['MPLCONFIGDIR'] = os.path.join(os.getcwd(), 'matplotlib_cache')
    os.makedirs(os.environ['MPLCONFIGDIR'], exist_ok=True)
    
    # Change working directory to scripts/ as notebooks expect to be run from there
    os.chdir('scripts')
    
    print("Executing target pair analysis notebook...")
    execute_and_export('metab_targetPair_analysis.ipynb', '../output/metab_targetPair_analysis_full_report.html', "Metabolite Target Pair Analysis")
    print("Done.")
