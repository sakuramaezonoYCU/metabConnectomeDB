"""
serotonin_config.py
===================
Centralized configurations and gene lists for the Serotonin-HTR7-TAM-IP4-HR repair axis.
References: 
- Li et al. 2026 (Cancer Cell)
- Reactome (R-HSA-1483257, R-HSA-5663220, R-HSA-1181150)
- HGNC (SLC6A4, TPH1/2, MAOA/B)
- IUPHAR/BPS (5-HT receptors)
"""

# -----------------------------------------------------------------------------
# Module 1: The HTR7+ TAM Signature (Li et al. 2026, Reactome)
# -----------------------------------------------------------------------------
HTR7_TAM_SIGNATURE = [
    "HTR7",       # Defining receptor (5-HT7)
    "CD163",      # Core M2/TAM marker
    "MRC1",       # CD206, Core M2/TAM marker
    "ARG1",       # Immunosuppressive enzyme
    "IL10",       # Immunosuppressive cytokine
    "VEGFA",      # Angiogenesis driver
    "GNAS"        # Gs-alpha subunit (downstream of HTR7)
]

# -----------------------------------------------------------------------------
# Module 2: The Full Serotonin Axis
# -----------------------------------------------------------------------------
SEROTONIN_SYNTHESIS = [
    "TPH1", "TPH2", "DDC"
]
SEROTONIN_TRANSPORT = [
    "SLC6A4"      # SERT
]
SEROTONIN_DEGRADATION = [
    "MAOA", "MAOB", "ALDH1A1", "ALDH2" # Degradation to 5-HIAA
]
# Complete receptor family (IUPHAR)
SEROTONIN_RECEPTORS = [
    # 5-HT1 (Gi/o)
    "HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F",
    # 5-HT2 (Gq/11)
    "HTR2A", "HTR2B", "HTR2C",
    # 5-HT3 (Ligand-gated ion channel)
    "HTR3A", "HTR3B", "HTR3C", "HTR3D", "HTR3E",
    # 5-HT4 (Gs)
    "HTR4",
    # 5-HT5 (Gi/o)
    "HTR5A",
    # 5-HT6 (Gs)
    "HTR6",
    # 5-HT7 (Gs) - The key TAM receptor
    "HTR7"
]

# -----------------------------------------------------------------------------
# Module 3: Paracrine Signaling (Ligand/Receptor)
# -----------------------------------------------------------------------------
PARACRINE_PAIRS = [
    ("TPH1", "HTR7"),
    ("SLC6A4", "HTR7")
]

# -----------------------------------------------------------------------------
# Module 4: Inositol Phosphate (IP4) Pathway (Reactome R-HSA-1483257)
# -----------------------------------------------------------------------------
INOSITOL_PATHWAY = [
    "ITPKB", "IPMK", "INPP5D", "INPPL1", "PIK3CA", "PIK3R1", "PTEN", "PLCB1", "PLCB2"
]

# -----------------------------------------------------------------------------
# Module 5: Homologous Recombination (HR) Repair (Reactome R-HSA-5663220)
# -----------------------------------------------------------------------------
HR_REPAIR_GENES = [
    "RAD51", "BRCA1", "BRCA2", "PALB2", "BARD1", "RAD51B", "RAD51C", "RAD51D", 
    "XRCC2", "XRCC3", "MRE11", "RAD50", "NBN"
]

# -----------------------------------------------------------------------------
# Module 6: Ovarian Metastasis Immune Evasion Target Signatures
# -----------------------------------------------------------------------------
# strict Database Constraint: Every gene listed here is an explicit metabolic target 
# found in the MetabConnectomeDB unique pairs dataset.

EXHAUSTION_TARGETS = [
    'PDCD1', 'CTLA4', 'HAVCR2', 'TIGIT'
]
SUPPRESSIVE_LIGAND_TARGETS = [
    'CD274', 'PDCD1LG2'
]
TREG_TARGETS = [
    'IL2RA'
]

# -----------------------------------------------------------------------------
# Module 7: EV Machinery (Biogenesis, Cargo, Secretion) (Reactome R-HSA-1181150)
# -----------------------------------------------------------------------------
EV_BIOGENESIS = [
    "RAB27A", "RAB27B", "SYYT7", "TSG101", "ALIX", "PDCD6IP", "CD9", "CD63", "CD81"
]
EV_CARGO_SORTING = [
    "VPS4A", "VPS4B", "CHMP4A", "CHMP4B", "CHMP4C"
]
SEROTONIN_DEP_SECRETION = [
    "VAMP3", "SNAP23", "STX4"
]

def get_all_genes():
    """Returns a flat list of all genes tracked in this module."""
    all_genes = (
        HTR7_TAM_SIGNATURE + 
        SEROTONIN_SYNTHESIS + 
        SEROTONIN_TRANSPORT + 
        SEROTONIN_DEGRADATION + 
        SEROTONIN_RECEPTORS + 
        INOSITOL_PATHWAY + 
        HR_REPAIR_GENES + 
        EXHAUSTION_TARGETS +
        SUPPRESSIVE_LIGAND_TARGETS +
        TREG_TARGETS +
        EV_BIOGENESIS + 
        EV_CARGO_SORTING + 
        SEROTONIN_DEP_SECRETION
    )
    return list(set(all_genes))
