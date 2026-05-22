import os
import re
import pandas as pd

def split_genes(target_str):
    if pd.isna(target_str) or not str(target_str).strip():
        return []
    parts = re.split(r'[;|,|/|\|]+', str(target_str))
    symbols = []
    for p in parts:
        p_clean = p.strip().upper()
        if p_clean and p_clean not in ["NAN", "NONE", "NULL"]:
            symbols.append(p_clean)
    return symbols

# 1. Load enzyme mapping databases
def load_enzyme_mappings():
    mappings = {}  # (organism, hmdb_id, gene) -> set of roles
    mappings_by_name = {}  # (organism, metab_name, gene) -> set of roles
    
    # Files to parse
    # Organism: human
    # Cellinker2 human (no header)
    cellinker_human = "input/databases/Cellinker2/Homo sapiens enzyme.txt"
    if os.path.exists(cellinker_human):
        print("Parsing Cellinker2 human...")
        with open(cellinker_human, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.strip("\n").split("\t")
                if len(parts) >= 9:
                    hmdb_id = parts[0].strip().upper()
                    metab_name = parts[1].strip().upper()
                    gene = parts[5].strip().upper()
                    role = parts[8].strip().lower()
                    if hmdb_id and gene and role in ["s", "p", "s+p"]:
                        k = ("human", hmdb_id, gene)
                        if k not in mappings:
                            mappings[k] = set()
                        mappings[k].add(role)
                    if metab_name and gene and role in ["s", "p", "s+p"]:
                        k_name = ("human", metab_name, gene)
                        if k_name not in mappings_by_name:
                            mappings_by_name[k_name] = set()
                        mappings_by_name[k_name].add(role)
                        
    # MRCLinkDB human (has header)
    mrclink_human = "input/databases/MRCLinkDB/Homo sapiens enzyme.txt"
    if os.path.exists(mrclink_human):
        print("Parsing MRCLinkDB human...")
        df = pd.read_csv(mrclink_human, sep="\t", encoding="utf-8", encoding_errors="replace")
        for _, row in df.iterrows():
            hmdb_id = str(row.get("HMDB_ID", "")).strip().upper()
            metab_name = str(row.get("METABOLITE_NAME", "")).strip().upper()
            gene = str(row.get("GENE_NAME", "")).strip().upper()
            role = str(row.get("enzyme product/substrate", "")).strip().lower()
            if hmdb_id and gene and role in ["s", "p", "s+p"]:
                k = ("human", hmdb_id, gene)
                if k not in mappings:
                    mappings[k] = set()
                mappings[k].add(role)
            if metab_name and gene and role in ["s", "p", "s+p"]:
                k_name = ("human", metab_name, gene)
                if k_name not in mappings_by_name:
                    mappings_by_name[k_name] = set()
                mappings_by_name[k_name].add(role)

    # Cellinker2 mouse (has header)
    cellinker_mouse = "input/databases/Cellinker2/Mus musculus enzyme.txt"
    if os.path.exists(cellinker_mouse):
        print("Parsing Cellinker2 mouse...")
        df = pd.read_csv(cellinker_mouse, sep="\t", encoding="utf-8", encoding_errors="replace")
        for _, row in df.iterrows():
            hmdb_id = str(row.get("HMDB_ID", "")).strip().upper()
            metab_name = str(row.get("METABOLITE_NAME", "")).strip().upper()
            gene = str(row.get("Mouse_gene symbol", "")).strip().upper()
            role = str(row.get("enzyme product/substrate", "")).strip().lower()
            if hmdb_id and gene and role in ["s", "p", "s+p"]:
                k = ("mouse", hmdb_id, gene)
                if k not in mappings:
                    mappings[k] = set()
                mappings[k].add(role)
            if metab_name and gene and role in ["s", "p", "s+p"]:
                k_name = ("mouse", metab_name, gene)
                if k_name not in mappings_by_name:
                    mappings_by_name[k_name] = set()
                mappings_by_name[k_name].add(role)

    # MRCLinkDB mouse (has header)
    mrclink_mouse = "input/databases/MRCLinkDB/Mus musculus enzyme.txt"
    if os.path.exists(mrclink_mouse):
        print("Parsing MRCLinkDB mouse...")
        df = pd.read_csv(mrclink_mouse, sep="\t", encoding="utf-8", encoding_errors="replace")
        for _, row in df.iterrows():
            hmdb_id = str(row.get("HMDB_ID", "")).strip().upper()
            metab_name = str(row.get("METABOLITE_NAME", "")).strip().upper()
            gene = str(row.get("Mouse_gene symbol", "")).strip().upper()
            role = str(row.get("enzyme product/substrate", "")).strip().lower()
            if hmdb_id and gene and role in ["s", "p", "s+p"]:
                k = ("mouse", hmdb_id, gene)
                if k not in mappings:
                    mappings[k] = set()
                mappings[k].add(role)
            if metab_name and gene and role in ["s", "p", "s+p"]:
                k_name = ("mouse", metab_name, gene)
                if k_name not in mappings_by_name:
                    mappings_by_name[k_name] = set()
                mappings_by_name[k_name].add(role)
                
    return mappings, mappings_by_name

# Consolidated mapping logic
def resolve_role(org, hmdb_id, metab_name, synonyms_str, gene, mappings, mappings_by_name):
    roles = set()
    
    # Try exact HMDB ID match
    if pd.notna(hmdb_id):
        h_clean = str(hmdb_id).strip().upper()
        k = (org, h_clean, gene)
        if k in mappings:
            for r in mappings[k]:
                if r == "s+p":
                    roles.add("s")
                    roles.add("p")
                else:
                    roles.add(r)
                    
    # Try Metabolite Name match
    if pd.notna(metab_name):
        n_clean = str(metab_name).strip().upper()
        k_name = (org, n_clean, gene)
        if k_name in mappings_by_name:
            for r in mappings_by_name[k_name]:
                if r == "s+p":
                    roles.add("s")
                    roles.add("p")
                else:
                    roles.add(r)
                    
    # Try Synonyms match
    if pd.notna(synonyms_str):
        for syn in re.split(r'[;,\/|\|]+', str(synonyms_str)):
            syn_clean = syn.strip().upper()
            if syn_clean:
                k_syn = (org, syn_clean, gene)
                if k_syn in mappings_by_name:
                    for r in mappings_by_name[k_syn]:
                        if r == "s+p":
                            roles.add("s")
                            roles.add("p")
                        else:
                            roles.add(r)
                            
    if "s" in roles and "p" in roles:
        return "s+p"
    elif "s" in roles:
        return "s"
    elif "p" in roles:
        return "p"
    return None

# Load mappings
mappings, mappings_by_name = load_enzyme_mappings()
print(f"Loaded {len(mappings)} mappings by HMDB ID and {len(mappings_by_name)} mappings by name.")

# Test on a master file
test_file = "input/databases/human_database_merge_unique_metab_target_pairs.csv"
if os.path.exists(test_file):
    df = pd.read_csv(test_file)
    print(f"Loaded master file with {len(df)} rows.")
    
    matched_count = 0
    original_non_null = df["enzyme product/substrate"].notna().sum()
    print(f"Original non-null count: {original_non_null}")
    
    for idx, row in df.iterrows():
        hmdb_id = row.get("HMDB_ID")
        metab_name = row.get("Metabolite_Name")
        synonyms = row.get("Synonyms")
        target_val = row.get("Target")
        
        genes = split_genes(target_val)
        row_roles = set()
        for g in genes:
            role = resolve_role("human", hmdb_id, metab_name, synonyms, g, mappings, mappings_by_name)
            if role:
                row_roles.add(role)
                
        if row_roles:
            matched_count += 1
            
    print(f"Enriched matching count using raw enzyme files: {matched_count} / {len(df)} rows ({matched_count/len(df)*100:.2f}%)")
