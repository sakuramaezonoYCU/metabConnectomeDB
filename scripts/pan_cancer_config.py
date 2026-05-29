"""
pan_cancer_config.py
====================
Central configuration for the pan-cancer meta-analysis pipeline.

To adjust which cell-count run is used for a cancer:
  1. Change the value in CANCER_CAP below (e.g. '100k' → '500k').
  2. Re-run compute_pan_cancer_meta.py  →  generates new upset plot + 23-gene list.
  3. Re-run generate_predictive_notebook.py  →  re-scores primary tumors.

Nothing else needs to be touched.
"""

import os

# ── Per-cancer cell-count cap ─────────────────────────────────────────────────
# Set the cap that was used for each cancer's primary_vs_metastasis DE CSV.
# Breast and Lung have been re-run at 500k; the others remain at 100k.
CANCER_CAP = {
    'breast':     '500k',   # ← change to '100k' to use the original run
    'colorectal': '100k',
    'lung':       '500k',   # ← change to '100k' to use the original run
    'melanoma':   '100k',
    'ovarian':    '100k',
}

# ── Tissue slugs (must match the h5ad filename prefix) ───────────────────────
# Format: {primary_tissue}_{meta_tissue1}_{meta_tissue2}_{meta_tissue3}
CANCER_TISSUE_SLUG = {
    'breast':     'breast_mammary-gland_liver_axilla_chest-wall',
    'colorectal': 'colon_large-intestine_liver_intestine_lung',
    'lung':       'lung_lymph-node_brain_pleural-fluid',
    'melanoma':   'skin-of-body_brain_abdomen_paracolic-gutter',
    'ovarian':    'ovary_abdomen_omentum_uterus',
}

# ── Primary tissue label (used for filtering in predictive scoring) ───────────
CANCER_PRIMARY_TISSUE = {
    'breast':     'breast',
    'colorectal': 'colon',
    'lung':       'lung',
    'melanoma':   'skin of body',
    'ovarian':    'ovary',
}

# ── CellxGene census version (used in h5ad filenames) ────────────────────────
CENSUS_VERSION = '2025-11-08'

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
def _project_root():
    """Returns the project root (one level above scripts/)."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def get_de_csv_path(cancer: str) -> str:
    """Full path to the primary_vs_metastasis DE metabolic targets CSV for a cancer."""
    root = _project_root()
    cap  = CANCER_CAP[cancer]
    return os.path.join(
        root, 'output', f'{cancer}_results',
        f'primary_vs_metastasis_{cancer}_results_DE_metabolic_targets_{cap}.csv'
    )

def get_h5ad_path(cancer: str) -> str:
    """Full path to the h5ad file for a cancer."""
    root = _project_root()
    slug = CANCER_TISSUE_SLUG[cancer]
    cap  = CANCER_CAP[cancer]
    return os.path.join(
        root, 'output', f'{cancer}_results',
        f'{slug}_{cap}_whole_transcriptome_{CENSUS_VERSION}.h5ad'
    )
