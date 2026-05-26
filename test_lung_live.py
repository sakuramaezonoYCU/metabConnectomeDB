import subprocess
import os
import sys

print("Running pipeline TEST live...")
cmd = [
    '/Users/sakuramaezono/venvs/metabConnectomeDB/bin/python', 'scripts/run_cancer_pipeline.py',
    'lung adenocarcinoma', 'lung,brain,bone,liver,adrenal gland', 'lung'
]

os.environ['CELLXGENE_CAP'] = '5000'

process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
for line in process.stdout:
    print(line, end="")

process.wait()
if process.returncode != 0:
    print("FAILED")
    sys.exit(process.returncode)
else:
    print("SUCCESS")
