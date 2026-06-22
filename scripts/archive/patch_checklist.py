import os
import shutil

src = "/Users/sakuramaezono/.gemini/antigravity-ide/brain/f6037e62-7fba-4e29-a100-df190cd45ddf/pipeline_execution_checklist.md"
dest = "pipeline_execution_checklist.md"
guide_to_delete = "Validation_Phase_Architecture_Guide.md"

if os.path.exists(guide_to_delete):
    os.remove(guide_to_delete)

with open(src, "r", encoding="utf-8") as f:
    content = f.read()

# Add Phase 6 and shift Phase 6 to Phase 7
old_phase_6 = """## Phase 6: Dynamic AI Insights Report Generation
*This final phase physically scrapes the HTML reports and CSVs to dynamically inject the ACTUAL results into the Markdown document.*

- `[ ]` **9. Build Dynamic AI Summary Document**
  - **Command:** `python scripts/tmp_build_md.py`
  - **Purpose:** Scrapes the exact numbers directly from the generated `*_full_report.html` files, ensuring `AI_summary_and_insights.md` NEVER contains mocked or hardcoded interpretations."""

new_phase_6 = """## Phase 6: Dynamic Gene Signature Validation
*This phase is orchestrator-led and validates all dynamically identified 4-cancer combinations generated in Phase 4.*

- `[ ]` **9. Validate Derived Signatures**
  - **Command:** `python scripts/run_validation_phase.py`
  - **Purpose:** The master execution script that automatically identifies all output combination signatures and routes them through the following tests:
    - `massspec_metabolomics_analysis.py`: Verifies signature genes using mass-spectrometry clinical cohorts.
    - `validate_tcga_signature.py`: Runs Cox Proportional Hazard regressions on TCGA survival datasets.
    - `verify_spatial.py`: Applies spatial enrichment scoring and calculates Moran's I on high-resolution Visium slides.
    - `generate_predictive_notebook.py`: Scores primary tumor cells directly to identify left/right skewed pre-metastatic subclones.

## Phase 7: Dynamic AI Insights Report Generation
*This final phase physically scrapes the HTML reports and CSVs to dynamically inject the ACTUAL results into the Markdown document.*

- `[ ]` **10. Build Dynamic AI Summary Document**
  - **Command:** `python scripts/tmp_build_md.py --phase all`
  - **Purpose:** Scrapes the exact numbers directly from the generated `*_full_report.html` files and CSVs (including Phase 6 validation results), ensuring `AI_summary_and_insights.md` NEVER contains mocked or hardcoded interpretations."""

content = content.replace(old_phase_6, new_phase_6)

# also fix Step 6 checkbox stuff
content = content.replace("""- `[ ]` **6. Pre-Metastatic Subclone Resolution**

- [x] Ensure `predictive_signature_biomarker.ipynb` parses subsets correctly and outputs proper results tables.
- [/] Ensure full pipeline execution succeeds without errors. *(Waiting for background pipeline to complete after fixing JSON string escape error in druggability)*""", """- `[ ]` **6. Pre-Metastatic Subclone Resolution**
  - **Command:** `python scripts/generate_predictive_notebook.py`
  - **Purpose:** Parses subset combination signatures and scores primary tumor cells directly.""")

with open(dest, "w", encoding="utf-8") as f:
    f.write(content)
