import os

# Base paths
BASE_DIR = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB"
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')

# Cancers to analyze
CANCERS = ['breast', 'colorectal', 'lung', 'melanoma', 'ovarian']

# Enrichr settings
ENRICHR_LIBRARIES = [
    'ChEA_2022',
    'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
    'TRRUST_Transcription_Factors_2019'
]

# Output files
NR1D2_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'nr1d2_results')
os.makedirs(NR1D2_RESULTS_DIR, exist_ok=True)
