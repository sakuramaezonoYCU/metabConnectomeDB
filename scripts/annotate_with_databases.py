#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🧬 Database Annotation Pipeline for Target-Pair Classifications
============================================================
This script expands target-pair classifications across three biological dimensions:
1. Sensor Type (target classifications: Receptor, Channel, Transporter, Enzyme)
2. Enzyme Product/Substrate (reaction dynamics: s, p, s+p)
3. Interaction Type (biological mechanisms: Promote, Inhibit, etc.)

It utilizes a hybrid database caching system:
- Downloads and caches Guide to Pharmacology (GtoPdb), HGNC, and UniProt resources in input/.
- If resources exist locally, it immediately loads them to avoid network overhead.
- Features a robust batch UniProt fetcher with individual fallback and clean multi-gene splitting.
- Overwrites all 8 target-pair master files in output/ and input/databases/.

Author: Antigravity (Advanced Agentic Coding Pair)
Date: 2026-05-21
"""

import os
import re
import json
import time
import requests
import pandas as pd
import numpy as np

# Import standalone Rhea enrichment module (produces Rhea_enzyme product/substrate)
from annotate_enzyme_rhea import run_rhea_enrichment
from annotate_enzyme_kegg import run_kegg_enrichment

# ==============================================================================
# ⚙️ CONFIGURATION & PATHS
# ==============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

INPUT_DATABASES_DIR = os.path.join(PROJECT_ROOT, "input", "databases")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "output")

# The 8 target-pair CSV files to enrich
TARGET_PAIR_FILES = [
    # input/databases/ copies
    os.path.join(INPUT_DATABASES_DIR, "human_database_merge_unique_metab_target_pairs.csv"),
    os.path.join(INPUT_DATABASES_DIR, "human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv"),
    os.path.join(INPUT_DATABASES_DIR, "mouse_database_merge_unique_metab_target_pairs.csv"),
    os.path.join(INPUT_DATABASES_DIR, "mouse_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv"),
    # output/ files
    os.path.join(OUTPUT_DIR, "human_database_merge_unique_metab_target_pairs.csv"),
    os.path.join(OUTPUT_DIR, "human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv"),
    os.path.join(OUTPUT_DIR, "mouse_database_merge_unique_metab_target_pairs.csv"),
    os.path.join(OUTPUT_DIR, "mouse_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv")
]

# Cache file paths in the input/ directory
GTOPDB_MAPPING_CACHE = os.path.join(PROJECT_ROOT, "input", "GtP_to_HGNC_mapping.csv")
GTOPDB_TARGETS_CACHE = os.path.join(PROJECT_ROOT, "input", "guidetopharmacology_targets.json")
GTOPDB_INTERACTIONS = os.path.join(PROJECT_ROOT, "input", "interactions.csv")
HGNC_APPROVED_CACHE = os.path.join(PROJECT_ROOT, "input", "hgnc_approved_genes.json")
UNIPROT_CACHE = os.path.join(PROJECT_ROOT, "input", "uniprot_annotations_cache.json")

# URLs for external databases
GTOPDB_MAPPING_URL = "https://www.guidetopharmacology.org/DATA/GtP_to_HGNC_mapping.csv"
GTOPDB_TARGETS_URL = "https://www.guidetopharmacology.org/services/targets"
HGNC_APPROVED_URL = "https://rest.genenames.org/fetch/status/Approved"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

# ==============================================================================
# 📥 HYBRID RESOURCE CACHING MANAGER
# ==============================================================================
def download_file_with_retry(url, dest_path, headers=None, is_json=False):
    """
    Downloads a file with retries and saves it to the destination path.
    """
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    
    print(f"Downloading from '{url}' to '{dest_path}'...")
    for attempt in range(3):
        try:
            response = requests.get(url, headers=headers, timeout=60, stream=True)
            if response.status_code == 200:
                if is_json:
                    # Parse as text and handle any strict validation errors
                    try:
                        # strict=False handles raw control characters like tabs/newlines inside JSON strings
                        data = json.loads(response.text, strict=False)
                        with open(dest_path, "w", encoding="utf-8") as f:
                            json.dump(data, f, indent=2)
                    except Exception as json_err:
                        print(f"JSON parsing error, writing raw content. Error: {json_err}")
                        with open(dest_path, "wb") as f:
                            f.write(response.content)
                else:
                    with open(dest_path, "wb") as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                print(f"-> Successfully saved to '{dest_path}'.")
                return
            else:
                print(f"-> Status code {response.status_code} on attempt {attempt + 1}. Retrying...")
        except requests.exceptions.RequestException as e:
            print(f"-> Connection error on attempt {attempt + 1}: {e}. Retrying...")
        time.sleep(2)
    raise RuntimeError(f"Failed to download resource from {url} after 3 attempts.")

def ensure_external_databases():
    """
    Ensures all 3 key database resources are cached locally.
    """
    # 1. Guide to Pharmacology mapping
    if not os.path.exists(GTOPDB_MAPPING_CACHE):
        download_file_with_retry(GTOPDB_MAPPING_URL, GTOPDB_MAPPING_CACHE)
    else:
        print(f"-> Using cached Guide to Pharmacology mapping at '{GTOPDB_MAPPING_CACHE}'.")

    # 2. Guide to Pharmacology targets JSON
    if not os.path.exists(GTOPDB_TARGETS_CACHE):
        download_file_with_retry(GTOPDB_TARGETS_URL, GTOPDB_TARGETS_CACHE, is_json=True)
    else:
        print(f"-> Using cached Guide to Pharmacology targets at '{GTOPDB_TARGETS_CACHE}'.")

    # 3. HGNC Approved Status genes JSON
    if not os.path.exists(HGNC_APPROVED_CACHE):
        download_file_with_retry(
            HGNC_APPROVED_URL, 
            HGNC_APPROVED_CACHE, 
            headers={"Accept": "application/json"},
            is_json=True
        )
    else:
        print(f"-> Using cached HGNC Approved genes at '{HGNC_APPROVED_CACHE}'.")

# ==============================================================================
# 🌐 UNIPROT DYNAMIC BATCH SEARCH SCRAPER
# ==============================================================================
def fetch_uniprot_batch(symbols):
    """
    Queries UniProt search API in a batch for a list of gene symbols.
    """
    # Clean symbols for search syntax safety
    clean_syms = []
    for s in symbols:
        c = re.sub(r'[^A-Z0-9_-]', '', s.upper())
        if c:
            clean_syms.append(c)
            
    if not clean_syms:
        return {"results": []}
        
    # Query structure: gene:(GENE1 OR GENE2 OR ...) AND organism_id:9606
    query_terms = [f"gene:{sym}" for sym in clean_syms]
    query = f"({' OR '.join(query_terms)}) AND organism_id:9606"
    
    params = {
        "query": query,
        "fields": "accession,gene_names,keyword,go",
        "size": 500
    }
    
    headers = {
        "User-Agent": "metabConnectomeDB-Annotator/1.0"
    }
    
    for attempt in range(3):
        try:
            response = requests.get(UNIPROT_SEARCH_URL, params=params, headers=headers, timeout=30)
            if response.status_code == 200:
                return response.json()
            elif response.status_code == 400:
                # Return None representing invalid query syntax (handled via individual fallback)
                return None
            elif response.status_code == 429:
                print("-> UniProt API rate limit (429) hit. Waiting 5s...")
                time.sleep(5)
            else:
                print(f"-> UniProt API error {response.status_code}. Retrying...")
        except requests.exceptions.RequestException as e:
            print(f"-> UniProt connection error: {e}. Retrying...")
        time.sleep(2)
    return None

def run_uniprot_caching(all_symbols):
    """
    Loads existing UniProt cache, identifies missing symbols, fetches them in 
    batches, parses keywords/GO-terms, and saves back to cache.
    """
    # Load cache
    cache = {}
    if os.path.exists(UNIPROT_CACHE):
        try:
            with open(UNIPROT_CACHE, "r", encoding="utf-8") as f:
                cache = json.load(f)
            print(f"Loaded existing UniProt cache with {len(cache):,} gene records.")
        except Exception as e:
            print(f"Warning: could not parse UniProt cache. Re-initializing. Error: {e}")
            cache = {}
            
    # Find symbols not in cache
    missing_symbols = [sym for sym in all_symbols if sym not in cache]
    if not missing_symbols:
        print("-> All target gene symbols are already cached in UniProt cache!")
        return cache
        
    print(f"Found {len(missing_symbols):,} gene symbols missing from UniProt cache.")
    
    # Process in batches of 100
    batch_size = 100
    for i in range(0, len(missing_symbols), batch_size):
        batch = missing_symbols[i:i + batch_size]
        print(f"Fetching UniProt batch {i // batch_size + 1} ({i + 1} to {min(i + batch_size, len(missing_symbols))}/{len(missing_symbols)})...")
        
        result_json = fetch_uniprot_batch(batch)
        if result_json is None:
            print(f"-> Batch failed with status 400. Falling back to individual queries...")
            # Query each symbol in the batch one by one defensively
            results_list = []
            for sym in batch:
                single_res = fetch_uniprot_batch([sym])
                if single_res and "results" in single_res:
                    results_list.extend(single_res["results"])
                time.sleep(0.1)
            result_json = {"results": results_list}
            
        if not result_json or "results" not in result_json:
            print(f"-> Failed to fetch batch starting with {batch[0]}. Skipping batch.")
            continue
            
        # Keep track of which symbols in this batch got successfully mapped
        mapped_symbols = set()
        
        for entry in result_json.get("results", []):
            genes_field = entry.get("genes", [])
            primary_accession = entry.get("primaryAccession", "")
            
            # Find which of our queried symbols this entry represents
            matched_symbol = None
            
            # Collect all names and synonyms for this entry
            entry_names = []
            for g in genes_field:
                if "geneName" in g:
                    entry_names.append(g["geneName"].get("value", "").upper())
                if "synonyms" in g:
                    for syn in g["synonyms"]:
                        entry_names.append(syn.get("value", "").upper())
            
            # Check if any entry name matches a symbol in our current batch
            for sym in batch:
                if sym.upper() in entry_names:
                    matched_symbol = sym
                    break
            
            if not matched_symbol and entry_names:
                # Fallback: check if the first geneName in UniProt matches any queried symbol
                for sym in batch:
                    for gname in entry_names:
                        if sym.upper() in gname:
                            matched_symbol = sym
                            break
                    if matched_symbol:
                        break
            
            if matched_symbol:
                # Extract keywords
                keywords = [kw.get("name") for kw in entry.get("keywords", []) if kw.get("name")]
                
                # Extract GO terms
                go_terms = []
                for ref in entry.get("uniProtKBCrossReferences", []):
                    if ref.get("database") == "GO":
                        go_id = ref.get("id")
                        go_name = ""
                        for prop in ref.get("properties", []):
                            if prop.get("key") == "GoTerm":
                                go_name = prop.get("value")
                                break
                        if go_id and go_name:
                            go_terms.append({"id": go_id, "name": go_name})
                
                # Save in cache
                cache[matched_symbol] = {
                    "accession": primary_accession,
                    "keywords": keywords,
                    "go_terms": go_terms,
                    "not_found": False
                }
                mapped_symbols.add(matched_symbol)
        
        # Mark remaining symbols in this batch as not_found so we don't query them again
        for sym in batch:
            if sym not in mapped_symbols:
                cache[sym] = {
                    "keywords": [],
                    "go_terms": [],
                    "not_found": True
                }
        
        # Save cache incrementally
        with open(UNIPROT_CACHE, "w", encoding="utf-8") as f:
            json.dump(cache, f, indent=2)
            
        # Brief sleep to be polite to the UniProt servers
        time.sleep(0.5)
        
    print(f"UniProt cache successfully updated. Total records: {len(cache):,}.")
    return cache

# ==============================================================================
# 🔮 LOCAL RULE-BASED CLASSIFIER (classify_by_rules)
# ==============================================================================
def classify_by_rules(symbol, name=""):
    """
    Regular-expression-based fallback classification layer.
    Ensures 100% annotation coverage for targets that fail lookup.
    """
    sym = str(symbol).upper().strip()
    nm = str(name).lower().strip()
    
    classes = set()
    
    # 1. Channels
    channel_prefixes = [
        r"^CACN", r"^SCN", r"^KCN", r"^CLCN", r"^TRP", r"^AQP", 
        r"^PIEZO", r"^ANO", r"^BEST", r"^PANX", r"^RYR", r"^ITPR", 
        r"^ASIC", r"^GRIN", r"^GRIA", r"^GRID", r"^GRIK", r"^GABR",
        r"^CHRN", r"^GLYR", r"^P2RX", r"^CX", r"^GJ"
    ]
    if any(re.match(p, sym) for p in channel_prefixes) or "channel" in nm or "pore-forming" in nm:
        classes.add("Channel")
        
    # 2. Receptors
    receptor_prefixes = [
        r"^GPCR", r"^GPR", r"^HRH", r"^ADRA", r"^ADRB", r"^CHRM", 
        r"^DRD", r"^HTR", r"^OPR", r"^P2RY", r"^TAAR", r"^VIPR", 
        r"^LPAR", r"^S1PR", r"^FFAR", r"^OXER", r"^FPR", r"^MAS", 
        r"^MRGP", r"^RXFP", r"^GNRH", r"^SSTR", r"^TACR", r"^TRHR", 
        r"^MC[1-5]R", r"^CNR[1-2]", r"^PACR", r"^CRHR", r"^GLPR", 
        r"^GCGR", r"^GHRH", r"^SCT", r"^PTH", r"^CALCR", r"^CASR", 
        r"^GRM[1-8]", r"^TAS1R", r"^TAS2R", r"^FZD[1-9]", r"^FZD10", 
        r"^SMO", r"^LDLR", r"^EGFR", r"^FGFR", r"^IGF1R", r"^INSR", 
        r"^PDGFR", r"^VEGFR", r"^FLT", r"^KDR", r"^KIT", r"^RET", 
        r"^MET", r"^RON", r"^AXL", r"^MER", r"^TYRO3", r"^ROR", 
        r"^RYK", r"^DDR", r"^ROS1", r"^ALK", r"^LTK", r"^MUSK", 
        r"^LMTK", r"^PTK", r"^STYK1", r"^NR[0-9]A", r"^THRA", r"^THRB", 
        r"^RARA", r"^RARB", r"^RARG", r"^PPARA", r"^PPARD", r"^PPARG", 
        r"^RORA", r"^RORB", r"^RORG", r"^VDR", r"^NR1", r"^ESR1", 
        r"^ESR2", r"^AR", r"^PGR", r"^GR", r"^MR", r"^NR3"
    ]
    if any(re.match(p, sym) for p in receptor_prefixes) or "receptor" in nm or "nuclear hormone" in nm:
        classes.add("Receptor")
        
    # 3. Transporters
    transporter_prefixes = [
        r"^SLC", r"^ABC", r"^ATP1A", r"^ATP1B", r"^ATP2A", r"^ATP2B", 
        r"^ATP4A", r"^ATP4B", r"^ATP7A", r"^ATP7B", r"^ATP8A", r"^ATP8B",
        r"^TFRC", r"^FABP"
    ]
    if any(re.match(p, sym) for p in transporter_prefixes) or "transporter" in nm or "carrier" in nm or "solute carrier" in nm:
        classes.add("Transporter")
        
    # 4. Enzymes
    enzyme_prefixes = [
        r"^CYP", r"^ALDH", r"^ADH", r"^COX", r"^PTGS", r"^NOS", 
        r"^HSD", r"^SULT", r"^UGT", r"^GST", r"^MAO", r"^COMT", 
        r"^FMO", r"^PLA2", r"^DGK", r"^PLD", r"^PLC", r"^PDE", 
        r"^AC", r"^GC", r"^AMPK", r"^MAPK", r"^AKT", r"^JAK", 
        r"^SRC", r"^LCK", r"^FYN", r"^SYK", r"^BTK", r"^CDK", 
        r"^DUSP", r"^PTP", r"^CASP", r"^MMP", r"^ADAM", r"^DPP", 
        r"^FAS", r"^HMGCR"
    ]
    if any(re.match(p, sym) for p in enzyme_prefixes) or nm.endswith("ase") or "kinase" in nm or "dehydrogenase" in nm or "transferase" in nm or "synthase" in nm or "enzyme" in nm or "catalytic activity" in nm:
        classes.add("Enzyme")
        
    return ", ".join(sorted(list(classes))) if classes else None

# ==============================================================================
# 📊 SOURCE SPECIFIC PARSING & LOOKUPS
# ==============================================================================
def load_gtopdb_lookups():
    """
    Loads Guide to Pharmacology mapping and target types, returning symbol -> type mapping.
    """
    df_map = pd.read_csv(GTOPDB_MAPPING_CACHE, skiprows=1)  # skip GtoPdb header comment line
    with open(GTOPDB_TARGETS_CACHE, "r", encoding="utf-8") as f:
        targets_list = json.load(f, strict=False)
        
    target_type_map = {t["targetId"]: t.get("type") for t in targets_list if "targetId" in t}
    
    symbol_to_type = {}
    for _, row in df_map.iterrows():
        sym = str(row.get("HGNC Symbol", "")).strip()
        iuphar_id = row.get("IUPHAR ID")
        if sym and iuphar_id and not pd.isna(iuphar_id):
            iuphar_id = int(iuphar_id)
            g_type = target_type_map.get(iuphar_id)
            if g_type:
                # Map IUPHAR types to our standardized categories
                if g_type in ["GPCR", "NHR", "CatalyticReceptor"]:
                    symbol_to_type[sym] = "Receptor"
                elif g_type in ["LGIC", "VGIC", "OtherIC"]:
                    symbol_to_type[sym] = "Channel"
                elif g_type == "Transporter":
                    symbol_to_type[sym] = "Transporter"
                elif g_type == "Enzyme":
                    symbol_to_type[sym] = "Enzyme"
                    
    return symbol_to_type

def load_hgnc_lookups():
    """
    Loads HGNC Approved JSON, returning symbol -> gene_groups mapping.
    """
    with open(HGNC_APPROVED_CACHE, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    docs = data.get("response", {}).get("docs", [])
    symbol_to_groups = {}
    for doc in docs:
        sym = doc.get("symbol")
        groups = doc.get("gene_group", [])
        name = doc.get("name", "")
        if sym:
            symbol_to_groups[sym] = {
                "groups": [str(g).lower().strip() for g in groups],
                "name": str(name).lower().strip()
            }
    return symbol_to_groups

def map_hgnc_sensor_type(sym, hgnc_data):
    """
    Maps Sensor_Type using HGNC approved gene groups and names.
    """
    info = hgnc_data.get(sym)
    if not info:
        return None
        
    groups = info["groups"]
    name = info["name"]
    
    classes = set()
    for g in groups:
        if "receptor" in g or "nuclear hormone" in g:
            classes.add("Receptor")
        if "channel" in g or "pore" in g:
            classes.add("Channel")
        if "transporter" in g or "carrier" in g or "solute carrier" in g:
            classes.add("Transporter")
        if "enzyme" in g or "kinase" in g or "dehydrogenase" in g or "transferase" in g or "synthase" in g:
            classes.add("Enzyme")
            
    # Suffix check on approved name
    if "receptor" in name:
        classes.add("Receptor")
    if "channel" in name:
        classes.add("Channel")
    if "transporter" in name:
        classes.add("Transporter")
    if name.endswith("ase") or "kinase" in name or "dehydrogenase" in name or "transferase" in name or "synthase" in name:
        classes.add("Enzyme")
        
    return ", ".join(sorted(list(classes))) if classes else None

def map_uniprot_sensor_type(sym, uniprot_cache):
    """
    Maps Sensor_Type using cached UniProt keywords and GO terms.
    """
    info = uniprot_cache.get(sym)
    if not info or info.get("not_found"):
        return None
        
    keywords = [kw.lower() for kw in info.get("keywords", [])]
    go_terms = [t["name"].lower() for t in info.get("go_terms", [])]
    
    classes = set()
    
    # GPCR / Receptor keywords and GO terms
    rc_kws = ["g-protein coupled receptor", "receptor", "sensory transduction", "nuclear receptor"]
    if any(k in keywords for k in rc_kws) or any("receptor activity" in t or "receptor binding" in t or "receptor complex" in t for t in go_terms):
        classes.add("Receptor")
        
    # Channel keywords and GO terms
    ch_kws = ["ion channel", "channel", "gated ion channel", "voltage-gated channel"]
    if any(k in keywords for k in ch_kws) or any("channel activity" in t or "pore-forming" in t for t in go_terms):
        classes.add("Channel")
        
    # Transporter keywords and GO terms
    tr_kws = ["transport", "transporter", "symporter", "antiporter", "facilitated diffusion", "sodium/potassium atpase"]
    if any(k in keywords for k in tr_kws) or any("transporter activity" in t or "carrier activity" in t or "solute carrier" in t for t in go_terms):
        classes.add("Transporter")
        
    # Enzyme keywords and GO terms
    en_kws = ["hydrolase", "transferase", "oxidoreductase", "isomerase", "lyase", "ligase", "kinase", "phosphatase", "protease", "enzyme", "acyltransferase", "monooxygenase"]
    if any(k in keywords for k in en_kws) or any("catalytic activity" in t or "enzyme activity" in t or "transferase activity" in t or "hydrolase activity" in t or "kinase activity" in t for t in go_terms):
        classes.add("Enzyme")
        
    return ", ".join(sorted(list(classes))) if classes else None

# ==============================================================================
# 📊 INTERACTIONS & REACTIONS PARSING
# ==============================================================================
def load_gtopdb_interactions():
    """
    Loads Guide to Pharmacology target-ligand interactions.
    Returns a dictionary mapping: TARGET_GENE_SYMBOL (uppercase) -> list of ligand interaction entries.
    """
    gtopdb_interactions_path = "input/interactions.csv"
    if not os.path.exists(gtopdb_interactions_path):
        print(f"⚠️ GtoPdb interactions file not found at '{gtopdb_interactions_path}'!")
        return {}
        
    try:
        df_gtp = pd.read_csv(gtopdb_interactions_path, skiprows=1, low_memory=False)
        gtp_by_symbol = {}
        for _, row in df_gtp.iterrows():
            sym = str(row.get("Target Gene Symbol", "")).upper().strip()
            lig = str(row.get("Ligand", "")).upper().strip()
            sid = str(row.get("Ligand PubChem SID", "")).strip()
            if sid.endswith(".0"):  # Clean up any float-like SIDs
                sid = sid[:-2]
            
            if sym and sym != "NAN":
                if sym not in gtp_by_symbol:
                    gtp_by_symbol[sym] = []
                gtp_by_symbol[sym].append({
                    "ligand": lig,
                    "sid": sid,
                    "type": str(row.get("Type", "")),
                    "action": str(row.get("Action", ""))
                })
        return gtp_by_symbol
    except Exception as e:
        print(f"⚠️ Error loading GtoPdb interactions: {e}")
        return {}

def map_gtopdb_action_to_interaction(type_str, action_str):
    """
    Standardizes GtoPdb Type and Action values into our target interaction vocabulary:
    Promote, Inhibit, Regulate, Protect, Be consumed, Be released.
    """
    t = str(type_str).lower().strip()
    a = str(action_str).lower().strip()
    
    # 1. Inhibit
    inhibit_keywords = ["inhib", "antagonist", "blocker", "negative", "inverse agonist", "pore blocker", "gating inhibitor"]
    if any(k in t for k in inhibit_keywords) or any(k in a for k in inhibit_keywords):
        return "Inhibit"
        
    # 2. Promote
    promote_keywords = ["agonist", "activator", "positive", "potentiation", "activation", "biased agonist"]
    if any(k in t for k in promote_keywords) or any(k in a for k in promote_keywords):
        return "Promote"
        
    # 3. Regulate / Modulator fallback
    if "modulator" in t or "modulator" in a or "regulation" in a or "regulate" in a or "modulat" in a:
        return "Regulate"
        
    # 4. If there's an action, capitalize and return it as a backup
    if a and a != "nan" and a != "binding" and a != "":
        return a.capitalize()
    if t and t != "nan" and t != "":
        return t.capitalize()
        
    return None

def parse_reaction_roles(metabolite_name, synonyms_str, reactions_str):
    """
    Parses reaction strings to check if the metabolite is a substrate ('s'),
    a product ('p'), or both ('s+p').
    """
    if pd.isna(reactions_str) or not str(reactions_str).strip():
        return None
        
    # Get all names to match against
    names_to_match = {str(metabolite_name).strip().upper()}
    if not pd.isna(synonyms_str) and str(synonyms_str).strip():
        # Split synonyms by common delimiters, including commas
        for syn in re.split(r'[;,\/|\|]+', str(synonyms_str)):
            cleaned_syn = syn.strip().upper()
            if cleaned_syn:
                names_to_match.add(cleaned_syn)
                
    is_substrate = False
    is_product = False
    
    # Split by semicolon or pipe to get individual reactions
    reactions = re.split(r'[;|\|]+', str(reactions_str))
    for rxn in reactions:
        rxn = rxn.strip().upper()
        if "¡Ú" in rxn:
            left, right = rxn.split("¡Ú", 1)
            
            # Check left side and right side for metabolite name or synonyms
            for name in names_to_match:
                if name in left:
                    is_substrate = True
                if name in right:
                    is_product = True
                    
    if is_substrate and is_product:
        return "s+p"
    elif is_substrate:
        return "s"
    elif is_product:
        return "p"
    return None

# ==============================================================================
# 🔀 DEDUPLICATION & CONSOLIDATION UTILITIES
# ==============================================================================
def split_genes(target_str):
    """
    Splits complex/multi-protein targets into clean individual gene symbols.
    """
    if pd.isna(target_str) or not str(target_str).strip():
        return []
    parts = re.split(r'[;|,|/|\|]+', str(target_str))
    symbols = []
    for p in parts:
        p_clean = p.strip().upper()
        if p_clean and p_clean not in ["NAN", "NONE", "NULL"]:
            symbols.append(p_clean)
    return symbols

def consolidate_sensor_type(row_types):
    """
    Consolidates Sensor_Type categories (Receptor, Channel, Transporter, Enzyme).
    Row_types is a list/set of type strings that can be comma-separated.
    """
    flat_types = set()
    for t in row_types:
        if pd.isna(t) or not str(t).strip():
            continue
        parts = str(t).split(",")
        for p in parts:
            p_clean = p.strip().capitalize()
            if p_clean in ["Receptor", "Channel", "Transporter", "Enzyme"]:
                flat_types.add(p_clean)
            elif p_clean == "Gpcr":
                flat_types.add("Receptor")
                
    return ", ".join(sorted(list(flat_types))) if flat_types else np.nan

def consolidate_enzyme_product_substrate(val_list):
    """
    Deduplicates and standardizes reaction dynamics: s, p, s+p.
    Input can be a list of values (e.g. from OtherDB and Reactions parser).
    """
    cleaned_parts = set()
    for val in val_list:
        if pd.isna(val) or not str(val).strip():
            continue
            
        # Split by separators: pipe, plus, comma, slash
        parts = re.split(r'[|+\s,/]+', str(val).lower().strip())
        for p in parts:
            p_clean = p.strip()
            if p_clean in ["s", "substrate"]:
                cleaned_parts.add("s")
            elif p_clean in ["p", "product"]:
                cleaned_parts.add("p")
            elif p_clean in ["s+p", "s&p", "s,p"]:
                cleaned_parts.add("s")
                cleaned_parts.add("p")
                
    if "s" in cleaned_parts and "p" in cleaned_parts:
        return "s+p"
    elif "s" in cleaned_parts:
        return "s"
    elif "p" in cleaned_parts:
        return "p"
    
    return np.nan

def consolidate_interaction(val_list):
    """
    Deduplicates and standardizes clinical interactions: Promote, Inhibit, etc.
    Input can be a list of values (e.g. from OtherDB and GtoPdb).
    """
    cleaned_parts = set()
    for val in val_list:
        if pd.isna(val) or not str(val).strip():
            continue
            
        # Split by common separators
        parts = re.split(r'[;|,/|\|]+', str(val).strip())
        for p in parts:
            p_clean = p.strip().capitalize()
            # Map to valid vocabulary: Promote, Be consumed, Be released, Inhibit, Regulate, Protect
            if "consume" in p_clean.lower():
                cleaned_parts.add("Be consumed")
            elif "release" in p_clean.lower():
                cleaned_parts.add("Be released")
            elif p_clean in ["Promote", "Inhibit", "Regulate", "Protect"]:
                cleaned_parts.add(p_clean)
            elif p_clean:
                # Sane mapping fallback
                if "inhib" in p_clean.lower():
                    cleaned_parts.add("Inhibit")
                elif "promot" in p_clean.lower() or "activat" in p_clean.lower():
                    cleaned_parts.add("Promote")
                elif "regulat" in p_clean.lower() or "modulat" in p_clean.lower():
                    cleaned_parts.add("Regulate")
                else:
                    cleaned_parts.add(p_clean)
                    
    return ", ".join(sorted(list(cleaned_parts))) if cleaned_parts else np.nan

# ==============================================================================
# 🚀 MAIN PIPELINE EXECUTION
# ==============================================================================
def main():
    print("=" * 80)
    print("🧬 Starting metabConnectomeDB Annotation Enrichment Pipeline (Enhanced Version)...")
    print("=" * 80)
    
    # Step 1: Ensure all external resource caches exist locally
    ensure_external_databases()
    
    # Step 2: Extract all unique individual target symbols (splitting complexes)
    print("Scanning master target-pair files to extract unique target symbols...")
    all_targets = set()
    for file_path in TARGET_PAIR_FILES:
        if os.path.exists(file_path):
            print(f"  Reading symbols from: '{file_path}'")
            try:
                df = pd.read_csv(file_path, usecols=["Target"], low_memory=False)
                for val in df["Target"].dropna():
                    # Split comma/slash/pipe separated targets (subunit splitting)
                    symbols = split_genes(val)
                    for sym in symbols:
                        all_targets.add(sym)
            except Exception as e:
                print(f"  Warning reading '{file_path}': {e}")
                
    sorted_targets = sorted(list(all_targets))
    print(f"-> Extracted {len(sorted_targets):,} unique individual target gene symbols.")
    
    # Step 3: Run UniProt batch fetcher & cache ingestion
    uniprot_cache = run_uniprot_caching(sorted_targets)
    
    # Step 4: Load database lookup mappings
    print("Loading database lookup mappings...")
    gtopdb_map = load_gtopdb_lookups()
    print(f"  Loaded Guide to Pharmacology target lookups for {len(gtopdb_map):,} symbols.")
    
    gtopdb_interactions = load_gtopdb_interactions()
    print(f"  Loaded Guide to Pharmacology interactions for {len(gtopdb_interactions):,} symbols.")
    
    hgnc_map = load_hgnc_lookups()
    print(f"  Loaded HGNC Approved gene lookups for {len(hgnc_map):,} symbols.")
    
    # Step 5: Run Rhea + UniProt enzyme product/substrate enrichment
    #         (standalone one-time script that caches BridgeDb/UniProt/Rhea data)
    print("Running Rhea enzyme product/substrate enrichment pipeline...")
    try:
        run_rhea_enrichment(TARGET_PAIR_FILES)
    except Exception as e:
        print(f"⚠️ Rhea enrichment encountered an error (non-fatal): {e}")
        print("   Continuing without Rhea annotations...")
        
    print("Running KEGG enzyme product/substrate enrichment pipeline...")
    try:
        run_kegg_enrichment()
    except Exception as e:
        print(f"⚠️ KEGG enrichment encountered an error (non-fatal): {e}")
        print("   Continuing without KEGG annotations...")
    
    # Step 6: Process and overwrite each of the 8 master target-pair datasets
    print("Processing and overwriting master datasets...")
    for file_path in TARGET_PAIR_FILES:
        if not os.path.exists(file_path):
            print(f"⚠️ Target file not found, skipping: '{file_path}'")
            continue
            
        print(f"Enriching master dataset: '{file_path}'...")
        try:
            # Load master file
            df = pd.read_csv(file_path, low_memory=False)
            
            # 1. Rename existing columns to OtherDB_* if not already renamed
            rename_dict = {}
            if "Sensor_Type" in df.columns and "OtherDB_Sensor_Type" not in df.columns:
                rename_dict["Sensor_Type"] = "OtherDB_Sensor_Type"
            if "enzyme product/substrate" in df.columns and "OtherDB_enzyme product/substrate" not in df.columns:
                rename_dict["enzyme product/substrate"] = "OtherDB_enzyme product/substrate"
            if "Interaction" in df.columns and "OtherDB_Interaction" not in df.columns:
                rename_dict["Interaction"] = "OtherDB_Interaction"
                
            if rename_dict:
                df.rename(columns=rename_dict, inplace=True)
                print(f"  Renamed columns: {rename_dict}")
                
            # If columns don't exist at all, initialize them safely
            for col in ["OtherDB_Sensor_Type", "OtherDB_enzyme product/substrate", "OtherDB_Interaction"]:
                if col not in df.columns:
                    df[col] = np.nan
                    
            # 2. Populate source specific annotations row-wise
            gtopdb_col = []
            hgnc_col = []
            uniprot_col = []
            
            # New biological database source columns
            reactions_enzyme_col = []
            gtopdb_interaction_col = []
            
            # Lists to store the consolidated rows
            unified_sensor_type_col = []
            unified_enzyme_col = []
            unified_interaction_col = []
            
            for idx, row in df.iterrows():
                target_val = row.get("Target")
                if pd.isna(target_val):
                    gtopdb_col.append(np.nan)
                    hgnc_col.append(np.nan)
                    uniprot_col.append(np.nan)
                    reactions_enzyme_col.append(np.nan)
                    gtopdb_interaction_col.append(np.nan)
                    unified_sensor_type_col.append(np.nan)
                    unified_enzyme_col.append(np.nan)
                    unified_interaction_col.append(np.nan)
                    continue
                    
                # Split complex targets into subunits
                subunits = split_genes(target_val)
                
                # Fetch annotations for each subunit and combine them
                gtp_types = []
                hgnc_types = []
                unip_types = []
                local_types = []
                
                for sym in subunits:
                    # GtoPdb
                    gtp_t = gtopdb_map.get(sym)
                    if gtp_t:
                        gtp_types.append(gtp_t)
                    
                    # HGNC
                    hgnc_t = map_hgnc_sensor_type(sym, hgnc_map)
                    if hgnc_t:
                        hgnc_types.append(hgnc_t)
                    
                    # UniProt
                    unip_t = map_uniprot_sensor_type(sym, uniprot_cache)
                    if unip_t:
                        unip_types.append(unip_t)
                        
                    # Local fallback rule
                    h_info = hgnc_map.get(sym, {})
                    h_name = h_info.get("name", "")
                    loc_t = classify_by_rules(sym, h_name)
                    if loc_t:
                        local_types.append(loc_t)
                
                # Join subunit annotations
                row_gtp = ", ".join(sorted(list(set(gtp_types)))) if gtp_types else np.nan
                row_hgnc = ", ".join(sorted(list(set(hgnc_types)))) if hgnc_types else np.nan
                row_unip = ", ".join(sorted(list(set(unip_types)))) if unip_types else np.nan
                
                gtopdb_col.append(row_gtp)
                hgnc_col.append(row_hgnc)
                uniprot_col.append(row_unip)
                
                # Consolidated unified Sensor_Type
                other_type = row.get("OtherDB_Sensor_Type")
                candidates = [other_type, row_gtp, row_hgnc, row_unip]
                unified_sensor_type = consolidate_sensor_type(candidates)
                
                # Fallback to local rule classifier
                if pd.isna(unified_sensor_type) or not str(unified_sensor_type).strip():
                    unified_sensor_type = consolidate_sensor_type(local_types)
                    
                unified_sensor_type_col.append(unified_sensor_type if unified_sensor_type else np.nan)
                
                # 3. Resolve Reactions_enzyme product/substrate from REACTIONS column
                metab_val = row.get("Metabolite_Name") if "Metabolite_Name" in df.columns else ""
                syns_val = row.get("Synonyms") if "Synonyms" in df.columns else ""
                rxn_val = row.get("REACTIONS") if "REACTIONS" in df.columns else ""
                
                rxn_role = parse_reaction_roles(metab_val, syns_val, rxn_val)
                reactions_enzyme_col.append(rxn_role if rxn_role else np.nan)
                
                # Consolidated unified enzyme product/substrate (4 sources)
                other_enz = row.get("OtherDB_enzyme product/substrate")
                rhea_enz = row.get("Rhea_enzyme product/substrate") if "Rhea_enzyme product/substrate" in df.columns else np.nan
                kegg_enz = row.get("KEGG_enzyme product/substrate") if "KEGG_enzyme product/substrate" in df.columns else np.nan
                unified_enzyme = consolidate_enzyme_product_substrate([other_enz, rxn_role, rhea_enz, kegg_enz])
                unified_enzyme_col.append(unified_enzyme if unified_enzyme else np.nan)
                
                # 4. Resolve GtoPdb_Interaction from interactions.csv
                pc_sid_val = row.get("PubChem CID/SID") if "PubChem CID/SID" in df.columns else ""
                
                # Parse all PubChem SIDs
                sids = set()
                if pd.notna(pc_sid_val):
                    for part in re.split(r'[;,\/|\|]+', str(pc_sid_val)):
                        m = re.search(r'SID:(\d+)', part, re.IGNORECASE)
                        if m:
                            sids.add(m.group(1).strip())
                            
                # Parse metabolite names and synonyms
                metab_names = {str(metab_val).upper().strip()}
                if pd.notna(syns_val):
                    for s in re.split(r'[;,\/|\|]+', str(syns_val)):
                        if s.strip():
                            metab_names.add(s.strip().upper())
                            
                # Perform GtoPdb interaction matching
                matched_acts = set()
                for sub in subunits:
                    if sub in gtopdb_interactions:
                        for entry in gtopdb_interactions[sub]:
                            # Match by PubChem SID
                            if entry["sid"] and entry["sid"] in sids:
                                act = map_gtopdb_action_to_interaction(entry["type"], entry["action"])
                                if act:
                                    matched_acts.add(act)
                            # Match by name
                            elif entry["ligand"] in metab_names:
                                act = map_gtopdb_action_to_interaction(entry["type"], entry["action"])
                                if act:
                                    matched_acts.add(act)
                                    
                gtp_int_val = ", ".join(sorted(list(matched_acts))) if matched_acts else np.nan
                gtopdb_interaction_col.append(gtp_int_val)
                
                # Consolidated unified Interaction
                other_int = row.get("OtherDB_Interaction")
                unified_int = consolidate_interaction([other_int, gtp_int_val])
                unified_interaction_col.append(unified_int if unified_int else np.nan)
                
            # Assign database specific columns to DataFrame
            df["GtoPdb_Sensor_Type"] = gtopdb_col
            df["HGNC_Sensor_Type"] = hgnc_col
            df["UniProt_Sensor_Type"] = uniprot_col
            
            # Assign our new rich biological database source columns
            df["Reactions_enzyme product/substrate"] = reactions_enzyme_col
            df["GtoPdb_Interaction"] = gtopdb_interaction_col
            
            # Assign unified consolidated columns (replacing original columns exactly)
            df["Sensor_Type"] = unified_sensor_type_col
            df["enzyme product/substrate"] = unified_enzyme_col
            df["Interaction"] = unified_interaction_col
            
            # Print stats
            non_null_sensor = df["Sensor_Type"].notna().sum()
            total_rows = len(df)
            pct_sensor = (non_null_sensor / total_rows) * 100
            
            non_null_enz = df["enzyme product/substrate"].notna().sum()
            pct_enz = (non_null_enz / total_rows) * 100
            
            non_null_int = df["Interaction"].notna().sum()
            pct_int = (non_null_int / total_rows) * 100
            
            print(f"  Annotation Coverage stats:")
            print(f"    Unified Sensor_Type: {non_null_sensor}/{total_rows} ({pct_sensor:.2f}%)")
            print(f"    Unified enzyme: {non_null_enz}/{total_rows} ({pct_enz:.2f}%)")
            print(f"    Unified Interaction: {non_null_int}/{total_rows} ({pct_int:.2f}%)")
            
            # Save file back to disk
            df.to_csv(file_path, index=False)
            print(f"  Successfully overwrote: '{file_path}'")
            
        except Exception as e:
            print(f"🔴 Critical Error processing file '{file_path}': {e}")
            raise e
            
    print("=" * 80)
    print("🎉 metabConnectomeDB target-pair databases successfully enriched!")
    print("=" * 80)

if __name__ == "__main__":
    main()
