import os
import glob
import subprocess
import sys

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from pan_cancer_config import ANALYSIS_SUFFIX

META_RESULTS_DIR = os.path.join(BASE_DIR, 'output', 'pan_cancer_meta_results')

def run_validation_phase():
    print(f"==================================================")
    print(f"      PHASE 5: DYNAMIC SIGNATURE VALIDATION       ")
    print(f"==================================================")
    
    signatures = []
    
    # Strictly conserved signature
    strict_sig = os.path.join(META_RESULTS_DIR, f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv")
    if os.path.exists(strict_sig):
        signatures.append(strict_sig)
        
    # 4-cancer combinations
    combo_sigs = glob.glob(os.path.join(META_RESULTS_DIR, f"pan_cancer_signature_*{ANALYSIS_SUFFIX}.csv"))
    signatures.extend(combo_sigs)
    
    if not signatures:
        print(f"No signatures found in {META_RESULTS_DIR}.")
        sys.exit(1)
        
    print(f"Found {len(signatures)} signatures to validate.")
    
    # 1. Run TCGA, MassSpec, and Spatial validation for each signature
    for sig in signatures:
        sig_name = os.path.basename(sig)
        print(f"\n---> Validating {sig_name} <---")
        
        # TCGA
        try:
            subprocess.run(["python", "scripts/validate_tcga_signature.py", "--signature_csv", sig], check=True, cwd=BASE_DIR)
        except subprocess.CalledProcessError as e:
            print(f"CRITICAL ERROR: validate_tcga_signature.py failed on {sig_name}. Halting.")
            sys.exit(1)
            
        # MassSpec
        try:
            subprocess.run(["python", "scripts/massspec_metabolomics_analysis.py", "--signature_csv", sig], check=True, cwd=BASE_DIR)
            
            # Run cross-cohort comparison based on the analysis output
            sig_basename = sig_name.replace('.csv', '')
            subprocess.run(["python", "scripts/massspec_cross_cohort_comparison.py", "--signature-name", sig_basename], check=True, cwd=BASE_DIR)
            
            # Generate the notebook for this signature
            subprocess.run(["python", "scripts/generate_massspec_metabolomics_notebook.py", "--signature_csv", sig], check=True, cwd=BASE_DIR)
        except subprocess.CalledProcessError as e:
            print(f"CRITICAL ERROR: massspec analysis/comparison/notebook failed on {sig_name}. Halting.")
            sys.exit(1)
            
        # Spatial
        try:
            subprocess.run(["python", "scripts/verify_spatial.py", "--signature_csv", sig], check=True, cwd=BASE_DIR)
        except subprocess.CalledProcessError as e:
            print(f"CRITICAL ERROR: verify_spatial.py failed on {sig_name}. Halting.")
            sys.exit(1)
            
    # 3. Generate ALL downstream notebooks
    # Each block includes required compute steps followed by the generator.
    # This ensures execute_pancancer_notebooks.py will NEVER silently skip a missing .ipynb.
    print(f"\n{'='*60}")
    print(f"  GENERATING ALL DOWNSTREAM NOTEBOOKS")
    print(f"{'='*60}")

    generators = [
        # (description, compute_steps[], generator_script)
        (
            "Pan-Cancer Meta-Analysis Notebook",
            [],
            "scripts/generate_combined_pan_cancer_notebook.py"
        ),
        (
            "Predictive Signature Biomarker Notebook",
            [],
            "scripts/generate_predictive_notebook.py"
        ),
        (
            "ML Prognostic Classifier Notebook",
            [],
            None  # Special: handled below with --cancer all
        ),
        (
            "CAMP Pan-Cancer Integration Notebook",
            [],
            "scripts/create_camp_notebook.py"
        ),
        (
            "Master Regulator Analysis Notebook",
            [],
            "scripts/generate_master_regulator_notebook.py"
        ),
        (
            "Serotonin Axis Spatial Mapping Notebook",
            ["scripts/compute_serotonin_htr7_axis.py", "scripts/compute_serotonin_spatial.py"],
            "scripts/generate_serotonin_notebook.py"
        ),
        (
            "Deep-Dive Conserved Metab Gene Sig Notebook",
            [],
            "scripts/generate_nb1.py"
        ),
        (
            "MITF Regulon Expansion Notebook",
            [],
            "scripts/generate_nb2.py"
        ),
        (
            "Ovarian Serotonin Immune Evasion Notebook",
            ["scripts/compute_metastatic_immune_evasion.py", "scripts/verify_spatial_immune_evasion.py"],
            "scripts/generate_immune_evasion_notebook.py"
        ),
        (
            "Visium Spatial Validation Notebook",
            [],
            "scripts/generate_visium_notebook.py"
        ),
    ]

    for desc, compute_steps, generator in generators:
        print(f"\n---> {desc} <---")
        try:
            # Run compute dependencies first
            for compute_script in compute_steps:
                compute_path = os.path.join(BASE_DIR, compute_script)
                if os.path.exists(compute_path):
                    print(f"  Computing: {compute_script}")
                    subprocess.run(["python", compute_script], check=True, cwd=BASE_DIR)
                else:
                    print(f"  ⚠️  Compute script not found: {compute_script} (skipping)")

            # Run the generator
            if generator is not None:
                generator_path = os.path.join(BASE_DIR, generator)
                if os.path.exists(generator_path):
                    subprocess.run(["python", generator], check=True, cwd=BASE_DIR)
                else:
                    print(f"  ❌ GENERATOR NOT FOUND: {generator}")
                    sys.exit(1)
        except subprocess.CalledProcessError as e:
            print(f"  ❌ CRITICAL ERROR: {desc} failed: {e}")
            sys.exit(1)

    # ML Prognostic Classifier (special: needs --cancer all flag)
    print(f"\n---> ML Prognostic Classifier Notebook <---")
    try:
        subprocess.run(["python", "scripts/generate_ml_prognostic_classifier_notebook.py", "--cancer", "all"], check=True, cwd=BASE_DIR)
    except subprocess.CalledProcessError as e:
        print(f"  ❌ CRITICAL ERROR: ML Prognostic Classifier generation failed: {e}")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"  Phase 5 COMPLETE. All notebooks generated.")
    print(f"  Next: python scripts/execute_pancancer_notebooks.py")
    print(f"{'='*60}")

if __name__ == "__main__":
    run_validation_phase()

