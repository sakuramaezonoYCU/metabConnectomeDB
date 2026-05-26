import subprocess
import os
import sys
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")

cancers = [
    {
        "disease": "lung adenocarcinoma", 
        "tissues": "lung,brain,skeletal system,liver,adrenal gland",
        "primary_tissues": "lung"
    },
    {
        "disease": "breast cancer", 
        "tissues": "breast,liver,skeletal system,brain,axilla",
        "primary_tissues": "breast,mammary gland"
    },
    {
        "disease": "colorectal cancer", 
        "tissues": "colon,large intestine,intestine,rectum,liver,lung,peritoneum",
        "primary_tissues": "colon,large intestine,intestine,rectum"
    },
    {
        "disease": "ovarian cancer", 
        "tissues": "ovary,fallopian tube,peritoneum,liver,lung,omentum",
        "primary_tissues": "ovary,fallopian tube"
    },
    {
        "disease": "brain cancer", 
        "tissues": "brain,spinal cord",
        "primary_tissues": "brain"
    },
    {
        "disease": "melanoma", 
        "tissues": "skin,lung,liver,brain,skeletal system,lymph node",
        "primary_tissues": "skin"
    }
]

python_bin = "/Users/sakuramaezono/venvs/metabConnectomeDB/bin/python"
script_path = "scripts/run_cancer_pipeline.py"

def run_all():
    logging.info("Starting pipeline for all cancers with CAP=None (Test mode OFF)")
    for c in cancers:
        disease = c["disease"]
        tissues = c["tissues"]
        primary = c["primary_tissues"]
        logging.info(f"--- Running pipeline for {disease} ---")
        cmd = [python_bin, script_path, disease, tissues, primary]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, env=os.environ)
            if result.returncode == 0:
                logging.info(f"SUCCESS: {disease}")
            else:
                logging.error(f"FAILED: {disease}")
                logging.error(result.stderr)
        except Exception as e:
            logging.error(f"EXCEPTION for {disease}: {e}")
            
if __name__ == "__main__":
    run_all()
