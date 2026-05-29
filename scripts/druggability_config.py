# druggability_config.py
# Configuration and parameters for druggability analysis of gene axes

import os

# ---------------------------------------------------------
# PATH CONFIGURATION
# ---------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "druggability")

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------------------------------------
# INPUT PARAMETERS
# ---------------------------------------------------------

# The gene axis to be investigated
TARGET_GENES = [
    "GLS",
    "SGMS1",
    "SLC16A7",
    "SPTLC1"
]

# Analysis suffix used for file naming and tracking
# Change this parameter to create distinct output HTMLs
import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX

# ---------------------------------------------------------
# API CONFIGURATION
# ---------------------------------------------------------

DGIDB_API_URL = "https://dgidb.org/api/v2/interactions.json"
OPENTARGETS_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"
DISEASES_API_URL = "https://diseases.jensenlab.org/api/v2/entity"

# ---------------------------------------------------------
# DEPMAP CONFIGURATION
# ---------------------------------------------------------

# Path to the local DepMap CRISPRGeneEffect.csv file.
DEPMAP_DATA_PATH = os.path.join(BASE_DIR, "input", "CRISPRGeneEffect.csv")

# ---------------------------------------------------------
# OUTPUT CONFIGURATION
# ---------------------------------------------------------

# Generate filename base for outputs and HTML notebook
OUTPUT_BASENAME = f"druggability_axis{ANALYSIS_SUFFIX}"

print(f"[CONFIG] Loaded configurations. Output will be prefixed with: {OUTPUT_BASENAME}")
