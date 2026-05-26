import subprocess
import os

print("Running pipeline TEST for lung adenocarcinoma...")
cmd = [
    '/Users/sakuramaezono/venvs/metabConnectomeDB/bin/python', 'scripts/run_cancer_pipeline.py',
    'lung adenocarcinoma', 'lung,brain,bone,liver,adrenal gland', 'lung'
]

# Set a cap of 5000 cells so this test finishes in ~1-2 minutes!
os.environ['CELLXGENE_CAP'] = '5000'

res = subprocess.run(cmd, capture_output=True, text=True)
print("STDOUT:")
print(res.stdout)

if res.returncode != 0:
    print("STDERR:")
    print(res.stderr)
else:
    print("SUCCESS: Pipeline completed for lung adenocarcinoma (Test Mode)!")
