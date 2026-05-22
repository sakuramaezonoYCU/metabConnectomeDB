#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
📥 External Database Resource Downloader & Caching Manager
==========================================================
This script fetches all required external biological datasets from GtoPdb, 
HGNC, and other resources and caches them locally in the `input/` folder.
This ensures subsequent runs of the annotation pipeline are entirely offline.

Author: Antigravity
Date: 2026-05-21
"""

import os
import json
import time
import requests

# URLs for external databases
RESOURCES = {
    "GTOPDB_MAPPING": {
        "url": "https://www.guidetopharmacology.org/DATA/GtP_to_HGNC_mapping.csv",
        "dest": "input/GtP_to_HGNC_mapping.csv",
        "is_json": False
    },
    "GTOPDB_TARGETS": {
        "url": "https://www.guidetopharmacology.org/services/targets",
        "dest": "input/guidetopharmacology_targets.json",
        "is_json": True
    },
    "GTOPDB_INTERACTIONS": {
        "url": "https://www.guidetopharmacology.org/DATA/interactions.csv",
        "dest": "input/interactions.csv",
        "is_json": False
    },
    "HGNC_APPROVED": {
        "url": "https://rest.genenames.org/fetch/status/Approved",
        "dest": "input/hgnc_approved_genes.json",
        "headers": {"Accept": "application/json"},
        "is_json": True
    }
}

def download_file(url, dest_path, headers=None, is_json=False):
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
                    try:
                        # strict=False handles control characters inside JSON strings
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
                return True
            else:
                print(f"-> Status code {response.status_code} on attempt {attempt + 1}. Retrying...")
        except requests.exceptions.RequestException as e:
            print(f"-> Connection error on attempt {attempt + 1}: {e}. Retrying...")
        time.sleep(2)
    
    print(f"❌ Failed to download resource from {url} after 3 attempts.")
    return False

def main():
    print("=" * 80)
    print("📥 Starting External Database Resource Downloader...")
    print("=" * 80)
    
    success_count = 0
    total_resources = len(RESOURCES)
    
    for key, info in RESOURCES.items():
        url = info["url"]
        dest = info["dest"]
        is_json = info["is_json"]
        headers = info.get("headers")
        
        # Check if already cached
        if os.path.exists(dest):
            print(f"-> Resource '{key}' is already cached at '{dest}' (Size: {os.path.getsize(dest) / 1024 / 1024:.2f} MB).")
            success_count += 1
            continue
            
        success = download_file(url, dest, headers=headers, is_json=is_json)
        if success:
            success_count += 1
            
    print("=" * 80)
    print(f"🎉 Downloader finished. {success_count}/{total_resources} resources are cached locally in 'input/'.")
    print("=" * 80)

if __name__ == "__main__":
    main()
