import os
import pandas as pd
import json

def get_chea_mitf_targets(chea_path):
    targets = set()
    if os.path.exists(chea_path):
        with open(chea_path, 'r') as f:
            for line in f:
                if line.startswith("MITF"):
                    parts = line.strip().split("\t")
                    # First two parts are signature name, background. Genes start at index 2
                    genes = [g.split(",")[0] for g in parts[2:] if g]
                    targets.update(genes)
    return targets

def compute_mitf_regulon():
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    chea_path = os.path.join(BASE_DIR, "input", "databases", "chea", "ChEA_2022.txt")
    mitf_targets = get_chea_mitf_targets(chea_path)
    print(f"Total unique MITF targets across all ChEA datasets: {len(mitf_targets)}")
    
    # Load the 1,669 metabolic target universe.
    db_merge_path = os.path.join(BASE_DIR, "output", "human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv")
    if os.path.exists(db_merge_path):
        df = pd.read_csv(db_merge_path)
        if 'Target' in df.columns:
            universe = set(df['Target'].dropna().unique())
            print(f"Metabolic target universe size: {len(universe)}")
            overlap = universe.intersection(mitf_targets)
            print(f"Overlap size: {len(overlap)}")
            
            from pan_cancer_config import ANALYSIS_SUFFIX
            # Save results
            out_dir = os.path.join(BASE_DIR, "output", "mitf_regulon")
            os.makedirs(out_dir, exist_ok=True)
            overlap_df = df[df['Target'].isin(overlap)].copy()
            overlap_df.to_csv(os.path.join(out_dir, f"mitf_metabolic_regulon_pairs{ANALYSIS_SUFFIX}.csv"), index=False)
            print("Saved MITF metabolic regulon pairs.")
            
            pd.DataFrame({"Target": list(overlap)}).to_csv(os.path.join(out_dir, f"mitf_metabolic_regulon_genes{ANALYSIS_SUFFIX}.csv"), index=False)

if __name__ == "__main__":
    compute_mitf_regulon()
