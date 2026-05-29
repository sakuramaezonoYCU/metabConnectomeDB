import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX
import os
import pandas as pd

def compute_stat3_network():
    """
    Computes the network overlap between the strictly conserved pan-cancer genes
    and the STAT3 transcriptional targets from the ChEA 2022 database (U87 cells).
    """
    # File paths
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    chea_path = os.path.join(base_dir, 'input', 'databases', 'chea', 'ChEA_2022.txt')
    pan_cancer_genes_path = os.path.join(base_dir, 'output', 'pan_cancer_meta_results', f'pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv')
    
    output_dir = os.path.join(base_dir, 'output', f"deepdive_conserved_metabGeneSig{ANALYSIS_SUFFIX}", 'stat3_network')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'stat3_u87_targets_strictly_conserved.csv')
    
    print(f"Reading pan-cancer genes from: {pan_cancer_genes_path}")
    df_genes = pd.read_csv(pan_cancer_genes_path)
    pan_cancer_genes = set(df_genes['Strictly_Conserved_Gene'].dropna().str.upper())
    print(f"Loaded {len(pan_cancer_genes)} pan-cancer genes.")
    
    print(f"Reading ChEA 2022 database from: {chea_path}")
    stat3_targets = set()
    
    with open(chea_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if not parts:
                continue
            
            clean_parts = [p for p in parts if p.strip()]
            if not clean_parts:
                continue
            
            sig_name = clean_parts[0]
            if sig_name.startswith('STAT3 ') and ' U87 ' in sig_name:
                print(f"Found STAT3 U87 signature: {sig_name}")
                targets = clean_parts[1:]
                stat3_targets.update([t.upper() for t in targets])
                
    print(f"Loaded {len(stat3_targets)} STAT3 targets from U87 cells.")
    
    # Intersect
    overlapping_genes = pan_cancer_genes.intersection(stat3_targets)
    print(f"Found {len(overlapping_genes)} overlapping genes: {sorted(list(overlapping_genes))}")
    
    # Save Network Edges
    results = []
    for gene in sorted(list(overlapping_genes)):
        results.append({
            'Source': 'STAT3',
            'Target': gene,
            'Interaction': 'Transcriptional Regulation (ChIP-Seq)',
            'Cell_Line': 'U87',
            'Database': 'ChEA_2022'
        })
        
    df_network = pd.DataFrame(results)
    df_network.to_csv(output_path, index=False)
    print(f"Saved STAT3 network to {output_path}")

if __name__ == '__main__':
    compute_stat3_network()
