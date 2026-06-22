import sys
import os

sys.path.append(os.getcwd())
from execute_and_export_notebooks import execute_and_export

if __name__ == '__main__':
    # Fix matplotlib cache dir issue
    os.environ['MPLCONFIGDIR'] = os.path.join(os.getcwd(), 'matplotlib_cache')
    os.makedirs(os.environ['MPLCONFIGDIR'], exist_ok=True)
    
    print("Executing unique metab exploration notebook...")
    execute_and_export('unique_metab_data_exploration.ipynb', '../output/unique_metab_data_exploration_full_report.html', "Unique Metab Data Exploration")
    print("Done.")
