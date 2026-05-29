import os
import json
import requests
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import pandas as pd
from typing import Dict, List, Set, Tuple

def load_json_cache(filepath: str) -> dict:
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            return json.load(f)
    return {}

def save_json_cache(filepath: str, data: dict):
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)

def get_kegg_bulk_mapping(filepath: str, delimiter='\t') -> dict:
    """Reads a local KEGG bulk mapping file and returns a dictionary {source: [targets]}."""
    print(f"    Reading {filepath}...")
    if not os.path.exists(filepath):
        print(f"      Warning: Local file not found: {filepath}")
        return {}
        
    mapping = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split(delimiter)
            if len(parts) == 2:
                src, tgt = parts[0], parts[1]
                if src not in mapping:
                    mapping[src] = []
                mapping[src].append(tgt)
    return mapping

def fetch_kegg_reactions(reaction_ids: List[str]) -> dict:
    """Fetches KEGG reactions in batches of 10 and parses the EQUATION field."""
    results = {}
    
    for i in range(0, len(reaction_ids), 10):
        batch = reaction_ids[i:i+10]
        query = "+".join(batch)
        print(f"    Fetching batch {i+1} to {min(i+10, len(reaction_ids))}...")
        url = f"https://rest.kegg.jp/get/{query}"
        
        try:
            response = requests.get(url, timeout=30, verify=False)
        except Exception as e:
            print(f"      Warning: Batch {query} failed: {e}")
            continue
            
        if response.status_code != 200:
            print(f"      Warning: Batch {query} failed.")
            continue
            
        current_entry = None
        for line in response.text.split('\n'):
            if line.startswith("ENTRY"):
                parts = line.split()
                if len(parts) >= 2:
                    current_entry = "rn:" + parts[1]
                    results[current_entry] = {"substrates": [], "products": []}
            elif line.startswith("EQUATION") and current_entry:
                equation = line.replace("EQUATION", "").strip()
                # Parse equation e.g. "C05167 + C00001 <=> C00161 + C00014"
                if "<=>" in equation:
                    left, right = equation.split("<=>")
                elif "=>" in equation:
                    left, right = equation.split("=>")
                else:
                    continue
                
                # Extract C-numbers
                subs = [word for word in left.split() if word.startswith("C") and len(word) == 6 and word[1:].isdigit()]
                prods = [word for word in right.split() if word.startswith("C") and len(word) == 6 and word[1:].isdigit()]
                
                # Prefix with cpd:
                results[current_entry]["substrates"] = [f"cpd:{s}" for s in subs]
                results[current_entry]["products"] = [f"cpd:{p}" for p in prods]
                
    return results

def run_kegg_enrichment():
    print("================================================================================")
    print("🧪 Starting KEGG Enzyme Product/Substrate Enrichment Pipeline...")
    print("================================================================================")
    
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_dir = os.path.join(base_dir, 'input')
    db_dir = os.path.join(input_dir, 'databases')
    output_dir = os.path.join(base_dir, 'output')
    
    # 1. Load Caches
    print("\\n📊 Step 1: Loading required caches...")
    uniprot_cache = load_json_cache(os.path.join(input_dir, 'uniprot_annotations_cache.json'))
    hmdb_to_chebi_cache = load_json_cache(os.path.join(input_dir, 'hmdb_to_chebi_cache.json'))
    
    # Reverse mapping: ChEBI -> HMDB
    chebi_to_hmdb = {}
    for hmdb, chebi in hmdb_to_chebi_cache.items():
        if chebi:
            if chebi not in chebi_to_hmdb:
                chebi_to_hmdb[chebi] = []
            chebi_to_hmdb[chebi].append(hmdb)
    
    # 2. Load KEGG Bulk Mappings from local files
    print("\\n🔗 Step 2: Loading KEGG bulk mappings...")
    # UniProt -> KEGG Gene (Target: KEGG Gene, Source: up:...)
    kegg_uniprot_hsa = get_kegg_bulk_mapping(os.path.join(db_dir, "kegg_hsa_uniprot.txt"))
    kegg_uniprot_mmu = get_kegg_bulk_mapping(os.path.join(db_dir, "kegg_mmu_uniprot.txt"))
    uniprot_to_kegg = {**kegg_uniprot_hsa, **kegg_uniprot_mmu}
    
    # KEGG Gene -> EC
    kegg_gene_to_ec_hsa = get_kegg_bulk_mapping(os.path.join(db_dir, "kegg_hsa_ec.txt"))
    kegg_gene_to_ec_mmu = get_kegg_bulk_mapping(os.path.join(db_dir, "kegg_mmu_ec.txt"))
    kegg_gene_to_ec = {**kegg_gene_to_ec_hsa, **kegg_gene_to_ec_mmu}
    
    # EC -> Reaction
    ec_to_rn = get_kegg_bulk_mapping(os.path.join(db_dir, "kegg_ec_rn.txt"))
    
    # KEGG Compound -> ChEBI
    # Source: cpd:..., Target: chebi:...
    cpd_to_chebi = get_kegg_bulk_mapping(os.path.join(db_dir, "kegg_cpd_chebi.txt"))
    
    # 3. Identify all required KEGG reactions
    print("\\n🧬 Step 3: Mapping Target Symbols to KEGG Reactions...")
    symbol_to_rn = {}
    all_needed_rns = set()
    
    for symbol, uniprot_data in uniprot_cache.items():
        accession = uniprot_data.get('accession')
        if not accession: continue
            
        up_key = f"up:{accession}"
        if up_key in uniprot_to_kegg:
            symbol_to_rn[symbol] = set()
            kegg_genes = uniprot_to_kegg[up_key]
            
            for kgene in kegg_genes:
                if kgene in kegg_gene_to_ec:
                    ecs = kegg_gene_to_ec[kgene]
                    for ec in ecs:
                        if ec in ec_to_rn:
                            rns = ec_to_rn[ec]
                            symbol_to_rn[symbol].update(rns)
                            all_needed_rns.update(rns)
                            
    print(f"  Mapped {len(symbol_to_rn)} target symbols to {len(all_needed_rns)} unique KEGG reactions.")
    
    # 4. Fetch missing KEGG Reactions
    print("\\n⚗️ Step 4: Fetching KEGG Reaction definitions...")
    reaction_cache_file = os.path.join(input_dir, 'kegg_reaction_cache.json')
    reaction_cache = load_json_cache(reaction_cache_file)
    
    missing_rns = [rn for rn in all_needed_rns if rn not in reaction_cache]
    if missing_rns:
        print(f"  Fetching {len(missing_rns)} missing KEGG reactions...")
        new_reactions = fetch_kegg_reactions(missing_rns)
        reaction_cache.update(new_reactions)
        save_json_cache(reaction_cache_file, reaction_cache)
    else:
        print("  -> All KEGG reactions already cached.")
        
    # 5. Process Target-Metabolite pairs
    print("\\n🎯 Step 5: Applying KEGG enzyme product/substrate annotations...")
    
    # Target Symbol -> Set of HMDB Substrates / Products
    target_kegg_map = {}
    for symbol, rns in symbol_to_rn.items():
        subs_hmdb = set()
        prods_hmdb = set()
        for rn in rns:
            if rn in reaction_cache:
                for cpd in reaction_cache[rn].get("substrates", []):
                    if cpd in cpd_to_chebi:
                        for chebi in cpd_to_chebi[cpd]:
                            # cpd_to_chebi gives 'chebi:1234'; hmdb_to_chebi_cache uses bare number '1234'
                            # Try bare number first, then 'CHEBI:1234' as fallback
                            bare_id = chebi.replace("chebi:", "").replace("CHEBI:", "")
                            for chebi_key in [bare_id, f"CHEBI:{bare_id}", chebi]:
                                if chebi_key in chebi_to_hmdb:
                                    subs_hmdb.update(chebi_to_hmdb[chebi_key])
                                    break
                for cpd in reaction_cache[rn].get("products", []):
                    if cpd in cpd_to_chebi:
                        for chebi in cpd_to_chebi[cpd]:
                            bare_id = chebi.replace("chebi:", "").replace("CHEBI:", "")
                            for chebi_key in [bare_id, f"CHEBI:{bare_id}", chebi]:
                                if chebi_key in chebi_to_hmdb:
                                    prods_hmdb.update(chebi_to_hmdb[chebi_key])
                                    break
        target_kegg_map[symbol] = {
            "substrates": subs_hmdb,
            "products": prods_hmdb
        }
        
    master_files = [
        os.path.join(db_dir, 'human_database_merge_unique_metab_target_pairs.csv'),
        os.path.join(db_dir, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv'),
        os.path.join(db_dir, 'mouse_database_merge_unique_metab_target_pairs.csv'),
        os.path.join(db_dir, 'mouse_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv'),
        os.path.join(output_dir, 'human_database_merge_unique_metab_target_pairs.csv'),
        os.path.join(output_dir, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv'),
        os.path.join(output_dir, 'mouse_database_merge_unique_metab_target_pairs.csv'),
        os.path.join(output_dir, 'mouse_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')
    ]
    
    for fpath in master_files:
        if not os.path.exists(fpath): continue
        print(f"\\n  Processing: '{os.path.basename(fpath)}'...")
        
        df = pd.read_csv(fpath)
        # The target column is named 'Target' (not 'Target_Symbol')
        target_col = None
        for candidate in ['Target', 'Target_Gene', 'Target_Symbol']:
            if candidate in df.columns:
                target_col = candidate
                break
        if target_col is None or 'HMDB_ID' not in df.columns:
            print(f"    Skipping: missing Target or HMDB_ID column (found: {list(df.columns[:10])})")
            continue
            
        annotated = 0
        def get_kegg_annotation(row, _target_col=target_col):
            nonlocal annotated
            raw_symbol = row[_target_col]
            hmdb = row['HMDB_ID']
            if pd.isna(raw_symbol) or pd.isna(hmdb):
                return None
            # Handle comma-separated multi-gene entries like 'ICAM1,HAVCR1'
            symbols = [s.strip() for s in str(raw_symbol).split(',')]
            is_sub_any = False
            is_prod_any = False
            for symbol in symbols:
                if symbol in target_kegg_map:
                    is_sub_any = is_sub_any or (hmdb in target_kegg_map[symbol]["substrates"])
                    is_prod_any = is_prod_any or (hmdb in target_kegg_map[symbol]["products"])
            
            if is_sub_any and is_prod_any:
                annotated += 1
                return "substrate; product"
            elif is_sub_any:
                annotated += 1
                return "substrate"
            elif is_prod_any:
                annotated += 1
                return "product"
            return None
            
        df['KEGG_enzyme product/substrate'] = df.apply(get_kegg_annotation, axis=1)
        
        pct = (annotated / len(df)) * 100 if len(df) > 0 else 0
        print(f"    -> KEGG annotations: {annotated}/{len(df)} pairs ({pct:.2f}%)")
        
        df.to_csv(fpath, index=False)
        
    print("\\n================================================================================")
    print("🎉 KEGG enrichment pipeline complete!")
    print("================================================================================")

if __name__ == "__main__":
    run_kegg_enrichment()
