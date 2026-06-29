#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🧪 Rhea + UniProt Enzyme Product/Substrate Enrichment Pipeline
===============================================================
Standalone enrichment script that annotates metabolite-target pairs with
enzyme product/substrate roles (substrate 's', product 'p', or both 's+p')
using external biochemical reaction databases:

    1. BridgeDb API   → maps HMDB IDs to ChEBI IDs
    2. UniProt API     → fetches CATALYTIC ACTIVITY comments (Rhea + ChEBI participants)
    3. Rhea SPARQL     → resolves reaction directionality (left=substrate, right=product)

PIPELINE CLASSIFICATION:
    ▸ ONE-TIME enrichment script: Run once to build local caches, then re-run
      only when new metabolites or targets are added to the dataset.
    ▸ Caches are incremental: re-running skips already-cached entries.
    ▸ Called from annotate_with_databases.py during the main annotation pipeline.

CACHES PRODUCED (all under input/):
    - input/hmdb_to_chebi_cache.json     ← HMDB ID → ChEBI ID mapping
    - input/rhea_catalytic_cache.json    ← UniProt accession → catalytic activities + Rhea participants
    - input/rhea_reaction_cache.json     ← Rhea reaction ID → {substrates: [ChEBI], products: [ChEBI]}

Author: Antigravity (Advanced Agentic Coding Pair)
Date: 2026-05-22
"""

import os
import re
import json
import time
import requests
import pandas as pd
import numpy as np

# Load config to get TARGET_PAIR_FILES
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "input", "pipeline.config.json"), "r") as f:
    PIPELINE_CONFIG = json.load(f)

TARGET_PAIR_BASENAMES = PIPELINE_CONFIG.get("ANNOTATION_DATABASES", {}).get("TARGET_PAIR_FILES", [])

try:
    import requests
except ImportError:
    import urllib.request
    import urllib.parse
    requests = None  # Fallback to urllib

# ==============================================================================
# ⚙️ CONFIGURATION
# ==============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

# Cache paths
HMDB_CHEBI_CACHE = os.path.join(PROJECT_ROOT, "input", "hmdb_to_chebi_cache.json")
RHEA_CATALYTIC_CACHE = os.path.join(PROJECT_ROOT, "input", "rhea_catalytic_cache.json")
RHEA_REACTION_CACHE = os.path.join(PROJECT_ROOT, "input", "rhea_reaction_cache.json")
UNIPROT_CACHE = os.path.join(PROJECT_ROOT, "input", "uniprot_annotations_cache.json")

# API endpoints
BRIDGEDB_BASE = "https://webservice.bridgedb.org"
UNIPROT_API_BASE = "https://rest.uniprot.org/uniprotkb"
RHEA_SPARQL_ENDPOINT = "https://sparql.rhea-db.org/sparql"

# Species configuration
import json
with open(os.path.join(PROJECT_ROOT, "input", "pipeline.config.json"), "r") as __f:
    __cfg = json.load(__f)
SPECIES_CONFIG = __cfg.get("ANNOTATION_DATABASES", {}).get("SPECIES_CONFIG", {})

# Rate limiting
BRIDGEDB_DELAY = 0.2   # seconds between BridgeDb requests
UNIPROT_DELAY = 0.5    # seconds between UniProt requests
RHEA_DELAY = 0.3       # seconds between Rhea SPARQL queries


# ==============================================================================
# 🔧 HTTP HELPERS
# ==============================================================================
def _http_get_json(url, headers=None, timeout=30):
    """GET request returning parsed JSON, with retry logic."""
    default_headers = {"User-Agent": "metabConnectomeDB-RheaEnrichment/1.0"}
    if headers:
        default_headers.update(headers)
    
    for attempt in range(3):
        try:
            if requests:
                resp = requests.get(url, headers=default_headers, timeout=timeout)
                if resp.status_code == 200:
                    return resp.json()
                elif resp.status_code == 429:
                    print(f"  -> Rate limited (429). Waiting 5s...")
                    time.sleep(5)
                elif resp.status_code == 404:
                    return None
                else:
                    print(f"  -> HTTP {resp.status_code} on attempt {attempt + 1}")
            else:
                req = urllib.request.Request(url, headers=default_headers)
                resp = urllib.request.urlopen(req, timeout=timeout)
                return json.loads(resp.read().decode())
        except Exception as e:
            print(f"  -> Error on attempt {attempt + 1}: {e}")
            raise
        time.sleep(2)
    return None


def _http_get_text(url, headers=None, timeout=30):
    """GET request returning raw text, with retry logic."""
    default_headers = {"User-Agent": "metabConnectomeDB-RheaEnrichment/1.0"}
    if headers:
        default_headers.update(headers)
    
    for attempt in range(3):
        try:
            if requests:
                resp = requests.get(url, headers=default_headers, timeout=timeout)
                if resp.status_code == 200:
                    return resp.text
                elif resp.status_code == 404:
                    return None
            else:
                req = urllib.request.Request(url, headers=default_headers)
                resp = urllib.request.urlopen(req, timeout=timeout)
                return resp.read().decode()
        except Exception as e:
            print(f"  -> Error on attempt {attempt + 1}: {e}")
            raise
        time.sleep(2)
    return None


def _load_json_cache(path):
    """Loads a JSON cache file if it exists, returns empty dict otherwise."""
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            print(f"  Warning: could not parse cache '{path}': {e}")
            raise
    return {}


def _save_json_cache(cache, path):
    """Saves a dictionary to a JSON cache file."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(cache, f, indent=2, ensure_ascii=False)


# ==============================================================================
# 🔗 STEP 1: HMDB → ChEBI MAPPING (via BridgeDb)
# ==============================================================================
def build_hmdb_chebi_mapping(hmdb_ids, species="human"):
    """
    Maps HMDB IDs to ChEBI IDs using the BridgeDb web service.
    
    Args:
        hmdb_ids: List of HMDB identifiers (e.g., ['HMDB0000898', ...])
        species: 'human' or 'mouse' (BridgeDb species context)
    
    Returns:
        Dictionary {HMDB_ID: ChEBI_ID} (ChEBI ID without 'CHEBI:' prefix)
    """
    cache = _load_json_cache(HMDB_CHEBI_CACHE)
    
    # Find HMDB IDs not yet cached
    missing = [h for h in hmdb_ids if h not in cache]
    if not missing:
        print(f"  -> All {len(hmdb_ids)} HMDB IDs already cached in ChEBI mapping.")
        return cache
    
    species_label = SPECIES_CONFIG.get(species, SPECIES_CONFIG["human"])["bridgedb_species"]
    print(f"  Fetching ChEBI mappings for {len(missing)} HMDB IDs via BridgeDb ({species_label})...")
    
    for i, hmdb_id in enumerate(missing):
        if (i + 1) % 25 == 0 or (i + 1) == len(missing):
            print(f"    [{i + 1}/{len(missing)}] Processing {hmdb_id}...")
        
        url = f"{BRIDGEDB_BASE}/{species_label}/xrefs/Ch/{hmdb_id}"
        text = _http_get_text(url, timeout=15)
        
        if text:
            chebi_id = None
            for line in text.strip().split("\n"):
                parts = line.strip().split("\t")
                if len(parts) >= 2 and parts[1] == "ChEBI":
                    # Extract numeric ChEBI ID
                    raw_id = parts[0].strip()
                    # Handle both "29009" and "CHEBI:29009" formats
                    chebi_id = raw_id.replace("CHEBI:", "")
                    break
            
            if chebi_id:
                cache[hmdb_id] = chebi_id
            else:
                cache[hmdb_id] = None  # Mark as queried but no ChEBI found
        else:
            cache[hmdb_id] = None
        
        time.sleep(BRIDGEDB_DELAY)
    
    _save_json_cache(cache, HMDB_CHEBI_CACHE)
    
    found = sum(1 for v in cache.values() if v is not None)
    print(f"  -> ChEBI mapping complete: {found}/{len(cache)} HMDB IDs have ChEBI IDs.")
    return cache


# ==============================================================================
# 🧬 STEP 2: UNIPROT CATALYTIC ACTIVITY ENRICHMENT
# ==============================================================================
def fetch_uniprot_catalytic_activities(accessions_by_symbol):
    """
    Fetches CATALYTIC ACTIVITY comments from UniProt for each accession,
    extracting Rhea reaction IDs and ChEBI participant IDs.
    
    Args:
        accessions_by_symbol: dict {gene_symbol: uniprot_accession}
    
    Returns:
        Dictionary {gene_symbol: { ... }}
    """
    import concurrent.futures
    import threading
    
    cache = _load_json_cache(RHEA_CATALYTIC_CACHE)
    
    # Find symbols not yet cached
    missing = {sym: acc for sym, acc in accessions_by_symbol.items() if sym not in cache}
    if not missing:
        print(f"  -> All {len(accessions_by_symbol)} symbols already cached in catalytic activity cache.")
        return cache
    
    print(f"  Fetching catalytic activities for {len(missing)} UniProt accessions in parallel...")
    
    lock = threading.Lock()
    processed_count = [0]
    total_missing = len(missing)
    
    def fetch_one(sym, acc):
        url = f"{UNIPROT_API_BASE}/{acc}?format=json"
        data = _http_get_json(url, headers={"Accept": "application/json"}, timeout=20)
        
        result = None
        if data:
            activities = []
            for comment in data.get("comments", []):
                if comment.get("commentType") == "CATALYTIC ACTIVITY":
                    rxn = comment.get("reaction", {})
                    rhea_ids = []
                    chebi_ids = []
                    for xref in rxn.get("reactionCrossReferences", []):
                        db = xref.get("database", "")
                        xid = xref.get("id", "")
                        if db == "Rhea" and not xid.startswith("RHEA-COMP"):
                            rhea_ids.append(xid.replace("RHEA:", ""))
                        elif db == "ChEBI":
                            chebi_ids.append(xid.replace("CHEBI:", ""))
                    
                    activities.append({
                        "reaction_name": rxn.get("name", ""),
                        "ec_number": rxn.get("ecNumber", ""),
                        "rhea_ids": rhea_ids,
                        "chebi_ids": chebi_ids,
                    })
            
            result = {
                "accession": acc,
                "catalytic_activities": activities,
            }
        else:
            result = {
                "accession": acc,
                "catalytic_activities": [],
                "fetch_failed": True,
            }
            
        with lock:
            cache[sym] = result
            processed_count[0] += 1
            i = processed_count[0]
            if i % 50 == 0 or i == total_missing:
                print(f"    [{i}/{total_missing}] Fetched {sym} ({acc})...")
            if i % 50 == 0 or i == total_missing:
                _save_json_cache(cache, RHEA_CATALYTIC_CACHE)
                
        time.sleep(0.1) # Small delay to be polite to API
        return result

    # Run in parallel with 10 workers
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        futures = []
        for sym, acc in missing.items():
            futures.append(executor.submit(fetch_one, sym, acc))
        concurrent.futures.wait(futures)
    
    _save_json_cache(cache, RHEA_CATALYTIC_CACHE)
    
    has_activities = sum(1 for v in cache.values() if v.get("catalytic_activities"))
    print(f"  -> Catalytic activity cache complete: {has_activities}/{len(cache)} symbols have reactions.")
    return cache


# ==============================================================================
# ⚗️ STEP 3: RHEA DIRECTIONAL RESOLUTION (substrate vs product)
# ==============================================================================
def resolve_rhea_directions(rhea_ids):
    """
    For each Rhea reaction ID, queries the Rhea SPARQL endpoint to determine
    which ChEBI compounds are on the left side (substrates) vs right side (products).
    
    Args:
        rhea_ids: List of Rhea master reaction IDs (numeric strings)
    
    Returns:
        Dictionary {rhea_id: {
            'substrates': [chebi_id, ...],  (left-side participants)
            'products': [chebi_id, ...],    (right-side participants)
        }}
    """
    cache = _load_json_cache(RHEA_REACTION_CACHE)
    
    missing = [rid for rid in rhea_ids if rid not in cache]
    if not missing:
        print(f"  -> All {len(rhea_ids)} Rhea reactions already cached.")
        return cache
    
    print(f"  Resolving directionality for {len(missing)} Rhea reactions via SPARQL...")
    
    # Process in batches of 50 for efficiency
    batch_size = 50
    for batch_start in range(0, len(missing), batch_size):
        batch = missing[batch_start:batch_start + batch_size]
        batch_end = min(batch_start + batch_size, len(missing))
        print(f"    [{batch_start + 1}-{batch_end}/{len(missing)}] Querying batch...")
        
        # Build SPARQL VALUES clause for batch
        values_str = " ".join(f'"RHEA:{rid}"' for rid in batch)
        
        query = f"""
PREFIX rh: <http://rdf.rhea-db.org/>
SELECT ?accession ?side ?chebi WHERE {{
  VALUES ?accession {{ {values_str} }}
  ?reaction rh:accession ?accession .
  ?reaction rh:side ?side .
  ?side rh:contains ?participant .
  ?participant rh:compound ?compound .
  ?compound rh:chebi ?chebi .
}}
"""
        
        try:
            if requests:
                resp = requests.get(
                    RHEA_SPARQL_ENDPOINT,
                    params={"query": query, "format": "json"},
                    headers={
                        "Accept": "application/sparql-results+json",
                        "User-Agent": "metabConnectomeDB-RheaEnrichment/1.0",
                    },
                    timeout=60,
                )
                if resp.status_code == 200:
                    result = resp.json()
                else:
                    print(f"    -> SPARQL error {resp.status_code}, falling back to individual queries...")
                    result = None
            else:
                params = urllib.parse.urlencode({"query": query, "format": "json"})
                req = urllib.request.Request(
                    f"{RHEA_SPARQL_ENDPOINT}?{params}",
                    headers={
                        "Accept": "application/sparql-results+json",
                        "User-Agent": "metabConnectomeDB-RheaEnrichment/1.0",
                    },
                )
                resp = urllib.request.urlopen(req, timeout=60)
                result = json.loads(resp.read().decode())
        except Exception as e:
            print(f"    -> SPARQL batch error: {e}. Falling back to individual queries...")
            result = None
            raise
        
        if result:
            # Parse batch results
            for binding in result.get("results", {}).get("bindings", []):
                accession = binding["accession"]["value"]  # e.g., "RHEA:25286"
                rid = accession.replace("RHEA:", "")
                side_uri = binding["side"]["value"]         # e.g., "...25286_L" or "...25286_R"
                chebi_uri = binding["chebi"]["value"]       # e.g., "http://purl.obolibrary.org/obo/CHEBI_15377"
                
                # Extract ChEBI numeric ID
                chebi_match = re.search(r"CHEBI[_:](\d+)", chebi_uri)
                if not chebi_match:
                    continue
                chebi_id = chebi_match.group(1)
                
                # Determine side: _L = substrates (left), _R = products (right)
                if rid not in cache:
                    cache[rid] = {"substrates": [], "products": []}
                
                if f"_{rid}_L" in side_uri or side_uri.endswith(f"{rid}_L"):
                    if chebi_id not in cache[rid]["substrates"]:
                        cache[rid]["substrates"].append(chebi_id)
                elif f"_{rid}_R" in side_uri or side_uri.endswith(f"{rid}_R"):
                    if chebi_id not in cache[rid]["products"]:
                        cache[rid]["products"].append(chebi_id)
        else:
            # Fallback: query each reaction individually
            for rid in batch:
                if rid in cache:
                    continue
                _resolve_single_rhea(rid, cache)
                time.sleep(RHEA_DELAY)
        
        # Mark any batch IDs that returned no results as empty
        for rid in batch:
            if rid not in cache:
                cache[rid] = {"substrates": [], "products": []}
        
        # Incremental save
        _save_json_cache(cache, RHEA_REACTION_CACHE)
        time.sleep(RHEA_DELAY)
    
    has_data = sum(1 for v in cache.values() if v.get("substrates") or v.get("products"))
    print(f"  -> Rhea direction cache complete: {has_data}/{len(cache)} reactions have substrate/product data.")
    return cache


def _resolve_single_rhea(rid, cache):
    """Fallback: resolve a single Rhea reaction via SPARQL."""
    query = f"""
PREFIX rh: <http://rdf.rhea-db.org/>
SELECT ?side ?chebi WHERE {{
  ?reaction rh:accession "RHEA:{rid}" .
  ?reaction rh:side ?side .
  ?side rh:contains ?participant .
  ?participant rh:compound ?compound .
  ?compound rh:chebi ?chebi .
}}
"""
    try:
        if requests:
            resp = requests.get(
                RHEA_SPARQL_ENDPOINT,
                params={"query": query, "format": "json"},
                headers={
                    "Accept": "application/sparql-results+json",
                    "User-Agent": "metabConnectomeDB-RheaEnrichment/1.0",
                },
                timeout=30,
            )
            result = resp.json() if resp.status_code == 200 else None
        else:
            params = urllib.parse.urlencode({"query": query, "format": "json"})
            req = urllib.request.Request(
                f"{RHEA_SPARQL_ENDPOINT}?{params}",
                headers={
                    "Accept": "application/sparql-results+json",
                    "User-Agent": "metabConnectomeDB-RheaEnrichment/1.0",
                },
            )
            resp = urllib.request.urlopen(req, timeout=30)
            result = json.loads(resp.read().decode())
    except Exception:
        result = None
        raise
    
    cache[rid] = {"substrates": [], "products": []}
    if result:
        for binding in result.get("results", {}).get("bindings", []):
            side_uri = binding["side"]["value"]
            chebi_uri = binding["chebi"]["value"]
            chebi_match = re.search(r"CHEBI[_:](\d+)", chebi_uri)
            if not chebi_match:
                continue
            chebi_id = chebi_match.group(1)
            
            if f"_{rid}_L" in side_uri or side_uri.endswith(f"{rid}_L"):
                if chebi_id not in cache[rid]["substrates"]:
                    cache[rid]["substrates"].append(chebi_id)
            elif f"_{rid}_R" in side_uri or side_uri.endswith(f"{rid}_R"):
                if chebi_id not in cache[rid]["products"]:
                    cache[rid]["products"].append(chebi_id)


# ==============================================================================
# 🎯 STEP 4: ANNOTATION MATCHING
# ==============================================================================
def determine_enzyme_role(target_symbol, hmdb_id, hmdb_chebi_map, catalytic_cache, rhea_cache):
    """
    For a given (target, metabolite) pair, determines whether the metabolite
    is a substrate, product, or both in the target enzyme's known reactions.
    
    Returns: 's', 'p', 's+p', or None
    """
    # Get metabolite's ChEBI ID
    chebi_id = hmdb_chebi_map.get(hmdb_id)
    if not chebi_id:
        return None
    
    # Get target's catalytic activities
    cat_data = catalytic_cache.get(target_symbol)
    if not cat_data or not cat_data.get("catalytic_activities"):
        return None
    
    is_substrate = False
    is_product = False
    
    for activity in cat_data["catalytic_activities"]:
        # Quick check: is our ChEBI even in this reaction's participants?
        if chebi_id not in activity.get("chebi_ids", []):
            continue
        
        # Use Rhea directional data for precise substrate/product assignment
        for rhea_id in activity.get("rhea_ids", []):
            rxn_data = rhea_cache.get(rhea_id)
            if not rxn_data:
                continue
            
            if chebi_id in rxn_data.get("substrates", []):
                is_substrate = True
            if chebi_id in rxn_data.get("products", []):
                is_product = True
    
    if is_substrate and is_product:
        return "s+p"
    elif is_substrate:
        return "s"
    elif is_product:
        return "p"
    return None


# ==============================================================================
# 🚀 MAIN PIPELINE
# ==============================================================================
def run_rhea_enrichment(target_pair_files=None):
    """
    Main entry point. Runs the full Rhea enrichment pipeline on all
    target-pair CSVs, adding a 'Rhea_enzyme product/substrate' column.
    
    Can be called standalone or imported from annotate_with_databases.py.
    
    Args:
        target_pair_files: List of CSV file paths to enrich. If None, uses defaults.
    
    Returns:
        Dictionary of caches: {hmdb_chebi, catalytic, rhea_reactions}
    """
    print("=" * 80)
    print("🧪 Starting Rhea + UniProt Enzyme Product/Substrate Enrichment Pipeline...")
    print("=" * 80)
    
    if target_pair_files is None:
        input_db_dir = os.path.join(PROJECT_ROOT, "input", "databases")
        output_dir = os.path.join(PROJECT_ROOT, "output")
        target_pair_files = []
        for _f in TARGET_PAIR_BASENAMES:
            target_pair_files.append(os.path.join(input_db_dir, _f))
            target_pair_files.append(os.path.join(output_dir, _f))
    
    # ── Collect all unique HMDB IDs and target symbols across all files ──
    print("\n📊 Step 0: Scanning all target-pair files for unique identifiers...")
    all_hmdb_ids = set()
    all_targets = set()
    existing_files = []
    
    for fp in target_pair_files:
        if not os.path.exists(fp):
            continue
        existing_files.append(fp)
        try:
            df = pd.read_csv(fp, usecols=["HMDB_ID", "Target"], low_memory=False)
            all_hmdb_ids.update(df["HMDB_ID"].dropna().unique())
            all_targets.update(df["Target"].dropna().unique())
        except Exception as e:
            print(f"  Warning: could not scan '{fp}': {e}")
            raise e
    
    print(f"  Found {len(all_hmdb_ids)} unique HMDB IDs and {len(all_targets)} unique targets across {len(existing_files)} files.")
    
    # ── Step 1: HMDB → ChEBI mapping ──
    print("\n🔗 Step 1: Building HMDB → ChEBI mapping via BridgeDb...")
    hmdb_chebi_map = build_hmdb_chebi_mapping(sorted(all_hmdb_ids), species="human")
    
    # ── Step 2: UniProt catalytic activity enrichment ──
    print("\n🧬 Step 2: Fetching UniProt catalytic activities...")
    
    # Load existing UniProt annotations cache to get accession mappings
    uniprot_cache = _load_json_cache(UNIPROT_CACHE)
    
    # Build symbol → accession map for targets that have UniProt entries
    accessions_map = {}
    for sym in sorted(all_targets):
        entry = uniprot_cache.get(sym, {})
        acc = entry.get("accession")
        if acc and not entry.get("not_found"):
            accessions_map[sym] = acc
    
    print(f"  {len(accessions_map)}/{len(all_targets)} targets have UniProt accessions cached.")
    catalytic_cache = fetch_uniprot_catalytic_activities(accessions_map)
    
    # ── Step 3: Collect all unique Rhea IDs and resolve directions ──
    print("\n⚗️ Step 3: Resolving Rhea reaction directionality...")
    all_rhea_ids = set()
    for sym_data in catalytic_cache.values():
        for activity in sym_data.get("catalytic_activities", []):
            all_rhea_ids.update(activity.get("rhea_ids", []))
    
    print(f"  Found {len(all_rhea_ids)} unique Rhea reaction IDs to resolve.")
    rhea_cache = resolve_rhea_directions(sorted(all_rhea_ids))
    
    # ── Step 4: Apply annotations to each file ──
    print("\n🎯 Step 4: Applying Rhea enzyme product/substrate annotations...")
    
    for fp in existing_files:
        # Detect species from filename
        fname = os.path.basename(fp).lower()
        species = "mouse" if "mouse" in fname else "human"
        species_label = SPECIES_CONFIG[species]["label"]
        
        print(f"\n  Processing: '{os.path.basename(fp)}' ({species_label})...")
        try:
            df = pd.read_csv(fp, low_memory=False)
            
            rhea_col = []
            matches_found = 0
            
            for _, row in df.iterrows():
                target = row.get("Target")
                hmdb_id = row.get("HMDB_ID")
                
                if pd.isna(target) or pd.isna(hmdb_id):
                    rhea_col.append(np.nan)
                    continue
                
                role = determine_enzyme_role(
                    str(target).strip(),
                    str(hmdb_id).strip(),
                    hmdb_chebi_map,
                    catalytic_cache,
                    rhea_cache,
                )
                
                if role:
                    rhea_col.append(role)
                    matches_found += 1
                else:
                    rhea_col.append(np.nan)
            
            df["Rhea_enzyme product/substrate"] = rhea_col
            
            # Save back
            df.to_csv(fp, index=False)
            
            total = len(df)
            pct = (matches_found / total) * 100 if total > 0 else 0
            print(f"    -> Rhea annotations: {matches_found}/{total} pairs ({pct:.2f}%)")
            
        except Exception as e:
            print(f"    ❌ Error processing '{fp}': {e}")
            raise
    
    print("\n" + "=" * 80)
    print("🎉 Rhea enrichment pipeline complete!")
    print("=" * 80)
    
    return {
        "hmdb_chebi": hmdb_chebi_map,
        "catalytic": catalytic_cache,
        "rhea_reactions": rhea_cache,
    }


# ==============================================================================
# 📦 HELPER: Get Rhea annotations for use by annotate_with_databases.py
# ==============================================================================
def get_rhea_role(target_symbol, hmdb_id):
    """
    Convenience function for annotate_with_databases.py to look up a single
    (target, metabolite) pair's enzyme role from the cached data.
    
    Returns 's', 'p', 's+p', or None.
    Requires caches to be populated first (via run_rhea_enrichment).
    """
    hmdb_chebi = _load_json_cache(HMDB_CHEBI_CACHE)
    catalytic = _load_json_cache(RHEA_CATALYTIC_CACHE)
    rhea_rxns = _load_json_cache(RHEA_REACTION_CACHE)
    
    return determine_enzyme_role(target_symbol, hmdb_id, hmdb_chebi, catalytic, rhea_rxns)


if __name__ == "__main__":
    run_rhea_enrichment()
