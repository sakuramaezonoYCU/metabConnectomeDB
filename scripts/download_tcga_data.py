import os
import requests

def download_tcga_cohort(cohort="TCGA-BRCA"):
    """
    Downloads clinical and RNA-seq data for a specified TCGA cohort from UCSC Xena.
    Defaulting to BRCA as a representative large cohort for validating the conserved metastatic score.
    output_dir = os.path.join(os.path.dirname(__file__), "..", "input", "TCGA")
    os.makedirs(output_dir, exist_ok=True)
    
    # TCGA Pan-Cancer Clinical Data (includes PFI, OS, DFI etc.)
    clinical_url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp"
    clinical_path = os.path.join(output_dir, "TCGA_Survival_SupplementalTable.tsv")
    
    # RNA-seq data for the specific cohort (HTSeq - FPKM or similar)
    # UCSC Xena Hub standard URL format for TCGA cohorts
    rnaseq_url = f"https://gdc.xenahubs.net/download/{cohort}.htseq_fpkm.tsv.gz"
    rnaseq_path = os.path.join(output_dir, f"{cohort}.htseq_fpkm.tsv.gz")
    
    # Download Clinical
    if not os.path.exists(clinical_path):
        print(f"Downloading TCGA Clinical Data from {clinical_url}...")
        download_file(clinical_url, clinical_path)
    else:
        print(f"TCGA Clinical Data already exists at {clinical_path}.")
        
    # Download RNA-seq
    if not os.path.exists(rnaseq_path):
        print(f"Downloading RNA-seq data for {cohort} from {rnaseq_url}...")
        print("Note: This file is large and may take a few minutes.")
        download_file(rnaseq_url, rnaseq_path)
    else:
        print(f"RNA-seq data for {cohort} already exists at {rnaseq_path}.")

def download_file(url, path):
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        with open(path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192*4):
                f.write(chunk)
                
        print(f"Successfully downloaded to {path}")
    except Exception as e:
        print(f"Error downloading {url}: {e}")

if __name__ == "__main__":
    # We will use BRCA to validate the metastatic subclone hypothesis since it has high sample count
    download_tcga_cohort("TCGA-BRCA")
