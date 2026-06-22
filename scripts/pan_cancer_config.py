"""
pan_cancer_config.py
====================
Central configuration for the pan-cancer meta-analysis pipeline.

To adjust which cell-count run is used for a cancer:
  1. Change the value in CANCER_CAP below (e.g. '100k' → '100k').
  2. Re-run compute_pan_cancer_meta.py  →  generates new upset plot + strictly conserved gene list.
  3. Re-run generate_predictive_notebook.py  →  re-scores primary tumors.

Nothing else needs to be touched.
"""

import os
import json
import re


def _project_root():
    """Returns the project root (one level above scripts/)."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Load JSON Configuration from input folder
CONFIG_PATH = os.path.join(_project_root(), 'input', 'pipeline.config.json')
try:
    with open(CONFIG_PATH, 'r') as _f:
        _cfg = json.load(_f)
except Exception as e:
    print(f"Warning: Could not parse {CONFIG_PATH}: {e}")
    _cfg = {}

# ── Per-Phase Config Extraction ─────────────────────────────────────────────────
_p1 = _cfg.get("PHASE_1_DATABASE_REPORTING", {})
_p2 = _cfg.get("PHASE_2_SINGLE_CELL_INTEGRATION", {})
_p3 = _cfg.get("PHASE_3_METRICS", {})
_p45 = _cfg.get("PHASE_4_5_META_VALIDATION", {})
_p6 = _cfg.get("PHASE_6_REPORTING", {})

# ── KEGG Pathways Configuration ────────────────────────────────────────────────
KEGG_PATHWAYS = _cfg.get("KEGG_PATHWAYS", {})

# ── Per-cancer cell-count cap ─────────────────────────────────────────────────
# Set the cap that was used for each cancer's primary_vs_metastasis DE CSV.
CANCER_CAP = _p2.get("CANCER_CAP", {})
DE_TESTING_METHOD = _p2.get("DE_TESTING_METHOD", "t-test")

# ── Cancers to execute in Phase 2 ─────────────────────────────────────────────
CANCERS_TO_RUN = _p2.get("CANCERS_TO_RUN", list(CANCER_CAP.keys()))

# ── Primary tissue label (used for filtering in predictive scoring) ───────────
CANCER_PRIMARY_TISSUE = _p2.get("CANCER_PRIMARY_TISSUE", {})

# ── Disease queries for CellxGene census downloading ──────────────────────────
CANCER_DISEASE_QUERIES = _p2.get("CANCER_DISEASE_QUERIES", {})

# ── Masking settings for tumor cells ──────────────────────────────────────────
STRICT_MASK_CANCERS = _p3.get("STRICT_MASK_CANCERS", [])

# ── Color mapping for plots ───────────────────────────────────────────────────
CANCER_COLORS = _p6.get("CANCER_COLORS", {})

# ── CSV mapping for oxygen tension ────────────────────────────────────────────
CANCER_PO2_CSV_MAPPING = _p6.get("CANCER_PO2_CSV_MAPPING", {})

# ── TCGA abbreviation mapping (used for cross-cohort validation) ──────────────
TCGA_MAPPING = _p6.get("TCGA_MAPPING", {})

# ── CellxGene census version (used in h5ad filenames) ────────────────────────
CENSUS_VERSION = _p2.get("CENSUS_VERSION", "2025-11-08")

# ── Normalization Function ──────────────────────────────────────────────────────
def normalize_cancer_name(cancer_name: str) -> str:
    """
    Normalizes a cancer name for consistent dictionary lookups.
    Uses the JSON configuration to dynamically map variations (like 'ovary' to 'ovarian').
    """
    if not isinstance(cancer_name, str):
        return cancer_name
        
    # Base normalization
    normalized = re.sub(r'\s+', '', cancer_name.strip()).lower()
    
    # 1. Check if exact match in canonical keys
    if normalized in CANCER_CAP:
        return normalized
        
    # 2. Reverse lookup against CANCER_PRIMARY_TISSUE synonyms
    for canonical, synonyms in CANCER_PRIMARY_TISSUE.items():
        for syn in synonyms:
            syn_norm = re.sub(r'\s+', '', syn.strip()).lower()
            if normalized == syn_norm:
                return canonical
                
    # 3. Fuzzy match (e.g., matching the first 4 characters like 'ovar' from 'ovary' and 'ovarian')
    for canonical in CANCER_CAP.keys():
        if len(normalized) >= 4 and len(canonical) >= 4:
            if normalized.startswith(canonical[:4]) and canonical.startswith(normalized[:4]):
                return canonical
                
    return normalized

# ── API and Advanced Analysis Configuration ─────────────────────────────────────
DGIDB_API_URL = _p45.get("DGIDB_API_URL")
OPENTARGETS_API_URL = _p45.get("OPENTARGETS_API_URL")
DISEASES_API_URL = _p45.get("DISEASES_API_URL")
DEPMAP_DATA_PATH = os.path.join(_project_root(), _p45.get("DEPMAP_DATA_PATH", "input/CRISPRGeneEffect.csv"))

# ── Analysis label for output HTML report names ───────────────────────────────
# Auto-generated from CANCER_CAP — reflects which cancers use which cap.
def _build_analysis_suffix():
    caps = set(CANCER_CAP.values())
    cancer_count = len(CANCER_CAP)
    if len(caps) == 1:
        cap = next(iter(caps))
        return f'_{cancer_count}MetCan_{cap}'
    # Mixed caps: list each cancer's cap
    parts = [f"{c[:2].capitalize()}{v}" for c, v in CANCER_CAP.items()]
    return '_' + '_'.join(parts)

ANALYSIS_SUFFIX = _build_analysis_suffix()  # e.g. '_5MetCan_mixed' or '_5MetCan_100k'

# ── Derived paths (do not edit — computed from values above) ──────────────────

def get_de_csv_path(cancer: str) -> str:
    """Full path to the primary_vs_metastasis DE metabolic targets CSV for a cancer."""
    root = _project_root()
    cancer = normalize_cancer_name(cancer)
    cap  = CANCER_CAP[cancer]
    return os.path.join(
        root, 'output', f'{cancer}_results',
        f'primary_vs_metastasis_{cancer}_results_DE_metabolic_targets.csv'
    )

def get_h5ad_path(cancer: str) -> str:
    """Full path to the h5ad file for a cancer."""
    import glob
    root = _project_root()
    cancer = normalize_cancer_name(cancer)
    cap  = CANCER_CAP[cancer]
    pattern = os.path.join(
        root, 'output', f'{cancer}_results',
        f'*_{cap}_whole_transcriptome_{CENSUS_VERSION}.h5ad'
    )
    matches = glob.glob(pattern)
    if matches:
        # Sort by modification time (descending) to ensure we always grab the latest generated file
        # in case a newer census version generated a slightly different tissue slug.
        matches.sort(key=os.path.getmtime, reverse=True)
        return matches[0]
    return pattern

# Predictive subclone detection thresholds
SKEW_THRESHOLD = _p45.get("SKEW_THRESHOLD", 0.5)
SUBCLONE_SD_MULTIPLIER = _p45.get("SUBCLONE_SD_MULTIPLIER", 1.0)
