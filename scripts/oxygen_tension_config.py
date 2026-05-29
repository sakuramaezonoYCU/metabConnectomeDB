import os

# Known generic oxygen tensions of dominant metastatic sites (approximate percentage O2)
OXYGEN_TENSION_MAP = {
    "Breast": 5.4,       # Liver dominant
    "Colorectal": 5.4,   # Liver dominant
    "Lung": 1.3,         # Pleural / Brain (Hypoxic pleura dominant)
    "Melanoma": 4.4,     # Brain dominant
    "Ovarian": 5.5       # Peritoneal / Omentum (~5-6%)
}

# Core OXPHOS and Glycolysis metabolic gene sets (KEGG/Reactome derived)
GLYCOLYSIS_GENES = [
    "HK1", "HK2", "HK3", "GPI", "PFKM", "PFKL", "PFKP", 
    "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "PGK1", "PGK2", 
    "PGAM1", "PGAM2", "ENO1", "ENO2", "ENO3", "PKM", "PKLR", "LDHA", "LDHB"
]

OXPHOS_GENES = [
    "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA11", "NDUFA13",
    "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB7", "NDUFB8", "NDUFB10", "NDUFB11", 
    "NDUFC1", "NDUFC2", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS5", "NDUFS6", "NDUFS7",
    "SDHA", "SDHB", "SDHC", "SDHD",
    "UQCRB", "UQCRC1", "UQCRC2", "UQCRH", "UQCRQ", "UQCR10", "UQCR11",
    "COX4I1", "COX4I2", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX7A2", "COX7C", "COX8A",
    "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E", "ATP5PB", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5ME", "ATP5MG", "ATP5PF", "ATP5PO"
]

CANCERS = ["breast", "colorectal", "lung", "melanoma", "ovarian"]

# Paths
BASE_DIR = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB"
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "oxygen_tension")
os.makedirs(OUTPUT_DIR, exist_ok=True)
