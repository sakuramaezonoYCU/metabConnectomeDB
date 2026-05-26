import subprocess
print("Running pipeline for colorectal cancer...")
cmd = [
    '/Users/sakuramaezono/venvs/metabConnectomeDB/bin/python', 'scripts/run_cancer_pipeline.py',
    'colorectal cancer', 'colon,large intestine,liver', 'colon,large intestine'
]
import os
os.environ['CELLXGENE_CAP'] = '500' # Very fast
res = subprocess.run(cmd, capture_output=True, text=True)
print(res.stdout)
if res.returncode != 0:
    print(res.stderr)
