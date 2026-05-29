#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
📚 PubMed Metadata Scraper for Unique Metabolite-Target Pairs
============================================================
This script extracts all unique PubMed IDs (PMIDs) from the unique metabolite-target
pair database, web-scrapes their metadata (Title, Abstract, Journal, Year, Authors)
from the official NCBI Entrez Utilities (E-Fetch) API, and saves the results in an
incremental fashion to prevent redundant API calls.

Requirements:
  - requests, pandas, numpy
  - [Optional] input/ncbi_api_key.txt containing NCBI API Key for higher rate limits (10 req/s vs 3 req/s)

Author: Antigravity (Advanced Agentic Coding Pair)
Date: 2026-05-21
"""

import os
import re
import xml.etree.ElementTree as ET
import time
import requests
import pandas as pd
import numpy as np

# ==============================================================================
# ⚙️ CONFIGURATION & PARAMETERS
# ==============================================================================
INPUT_DATABASE = "output/human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv"
OUTPUT_FILE = "input/pubmed_results.csv"
API_KEY_FILE = "input/ncbi_api_key.txt"
BATCH_SIZE = 50  # PubMed E-Fetch supports batch querying

# Flagged outputs for incomplete data
MISSING_TITLE_FILE = "input/missing_PMID_title_info.csv"
MISSING_ABSTRACT_FILE = "input/missing_PMID_abstract_info.csv"

# ==============================================================================
# 🔍 EXTRACT UNIQUE PMIDs FROM DATABASE
# ==============================================================================
def extract_unique_db_pmids(csv_path):
    """
    Extracts, cleans, and deduplicates all unique PMIDs across PMID, Evidence,
    and Text_Evidence columns in the database.
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Database file not found at: {csv_path}")
        
    print(f"Reading database to extract PMIDs: '{csv_path}'...")
    df = pd.read_csv(csv_path, low_memory=False)
    all_pmids = set()
    
    # 1. Extract from PMID and Evidence columns
    for col in ['PMID', 'Evidence']:
        if col in df.columns:
            for val in df[col].dropna():
                val_str = str(val).strip()
                # Split by common separators: semicolon, comma, pipe, space, slash
                for chunk in re.split(r'[;|,|\s|/|\|]+', val_str):
                    chunk = chunk.strip()
                    # Clean decimals from float conversions (e.g., 12345.0 -> 12345)
                    if chunk.endswith('.0'):
                        chunk = chunk[:-2]
                    if chunk.isdigit() and len(chunk) >= 5:
                        all_pmids.add(chunk)
                    elif chunk and chunk.lower() not in ['nan', 'none', 'null']:
                        # Regex match for digits inside a string
                        clean_num = re.findall(r'\b\d{5,}\b', chunk)
                        for num in clean_num:
                            all_pmids.add(num)
                            
    # 2. Extract from Text_Evidence column (typically text paragraphs with inline citations)
    if 'Text_Evidence' in df.columns:
        for val in df['Text_Evidence'].dropna():
            val_str = str(val).strip()
            matches = re.findall(r'(?i)(?:pmid|pubmed)\s*:?\s*(\d+)', val_str)
            for m in matches:
                if len(m) >= 5:
                    all_pmids.add(m)
                    
    sorted_pmids = sorted(list(all_pmids), key=lambda x: int(x) if x.isdigit() else x)
    print(f"-> Extracted {len(sorted_pmids):,} unique PMIDs from database.")
    return sorted_pmids

# ==============================================================================
# 🌐 WEB-SCRAPE METADATA FROM NCBI PUBMED
# ==============================================================================
def fetch_pubmed_metadata(pmid_list):
    """
    Fetches Title, Abstract, Journal, Year, and Authors for a list of PMIDs
    via NCBI's E-Fetch API. Supports incremental fetching and rate limiting.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
    # 1. Load existing results if file exists for incremental fetching
    if os.path.exists(OUTPUT_FILE):
        print(f"Loading existing scraped data from '{OUTPUT_FILE}'...")
        existing_df = pd.read_csv(OUTPUT_FILE)
        existing_df["PMID"] = existing_df["PMID"].astype(str)
        scraped_pmids = set(existing_df["PMID"].unique())
        print(f"-> Found {len(scraped_pmids):,} already-scraped PMIDs.")
    else:
        existing_df = pd.DataFrame(columns=["PMID", "Title", "Abstract", "Journal", "Year", "Author"])
        scraped_pmids = set()
        
    # 2. Filter list to scrape only new PMIDs
    to_scrape = [p for p in pmid_list if str(p) not in scraped_pmids]
    if not to_scrape:
        print("🎉 All PMIDs are already scraped and up-to-date! No web-scraping needed.")
        return existing_df
        
    print(f"🚀 Preparing to scrape {len(to_scrape):,} new PMIDs from PubMed...")
    
    # 3. Read NCBI API key if available
    ncbi_api_key = None
    if os.path.exists(API_KEY_FILE):
        try:
            with open(API_KEY_FILE, "r") as file:
                ncbi_api_key = file.read().strip()
            print("🔑 NCBI API Key loaded successfully. Higher rate limits enabled (10 req/s).")
        except Exception as e:
            print(f"⚠️ Warning loading API Key: {e}. Defaulting to standard rate limit (3 req/s).")
    else:
        print("💡 NCBI API Key not found. Using standard rate limits (3 req/s).")
        
    # 4. Batch fetching with rate limiting
    new_results = []
    pmid_batches = [to_scrape[i:i + BATCH_SIZE] for i in range(0, len(to_scrape), BATCH_SIZE)]
    
    # Set request delay based on whether we have an API key
    delay = 0.12 if ncbi_api_key else 0.35
    
    for batch_idx, batch in enumerate(pmid_batches):
        print(f"   -> Fetching batch {batch_idx + 1}/{len(pmid_batches)} ({len(batch)} PMIDs)...")
        
        params = {
            "db": "pubmed",
            "id": ",".join(batch),
            "retmode": "xml"
        }
        if ncbi_api_key:
            params["api_key"] = ncbi_api_key
            
        try:
            response = requests.get(base_url, params=params, timeout=20)
            response.raise_for_status()
            
            # Parse XML Response
            root = ET.fromstring(response.text)
            for article in root.findall(".//PubmedArticle"):
                pmid_node = article.find(".//MedlineCitation/PMID")
                if pmid_node is None or not pmid_node.text:
                    continue
                pmid = pmid_node.text.strip()
                
                # Title
                title = article.findtext(".//ArticleTitle", default="N/A").strip()
                
                # Abstract (can contain nested tags or multiple AbstractText sections)
                abstract_texts = []
                for abs_text in article.findall(".//AbstractText"):
                    if abs_text.text:
                        abstract_texts.append(abs_text.text.strip())
                abstract = " ".join(abstract_texts) if abstract_texts else "No abstract available"
                
                # Journal Title
                journal_tag = article.find(".//Journal/Title")
                journal = journal_tag.text.strip() if (journal_tag is not None and journal_tag.text) else "N/A"
                # Strip parentheses details if present
                if journal != "N/A" and " (" in journal:
                    journal = journal.split(" (", 1)[0]
                    
                # Publication Year
                pub_date = article.find(".//PubDate")
                year = "N/A"
                if pub_date is not None:
                    year_tag = pub_date.find("Year")
                    if year_tag is not None and year_tag.text:
                        year = year_tag.text.strip()
                    else:
                        medline_date_tag = pub_date.find("MedlineDate")
                        if medline_date_tag is not None and medline_date_tag.text:
                            year = medline_date_tag.text.split()[0].strip()
                            
                # Authors (format: "Last Initials, Last Initials")
                authors = []
                author_list = article.find(".//AuthorList")
                if author_list is not None:
                    for author in author_list.findall("Author"):
                        last_name = author.findtext("LastName")
                        fore_name = author.findtext("ForeName")
                        initials = author.findtext("Initials")
                        collective = author.findtext("CollectiveName")
                        
                        if last_name and initials:
                            authors.append(f"{last_name} {initials}")
                        elif last_name and fore_name:
                            authors.append(f"{last_name} {fore_name}")
                        elif collective:
                            authors.append(collective)
                        elif last_name:
                            authors.append(last_name)
                author_str = ", ".join(authors) if authors else "N/A"
                
                new_results.append({
                    "PMID": pmid,
                    "Title": title,
                    "Abstract": abstract,
                    "Journal": journal,
                    "Year": year,
                    "Author": author_str
                })
                
        except Exception as e:
            print(f"❌ Error fetching batch starting with {batch[0]}: {e}")
            
        # Respect NCBI API rate limit compliance
        time.sleep(delay)
        
    # 5. Merge and Save updated CSV
    if new_results:
        new_df = pd.DataFrame(new_results)
        updated_df = pd.concat([existing_df, new_df], ignore_index=True)
        
        # Post-process, drop duplicates and order by numeric PMID
        updated_df["PMID"] = updated_df["PMID"].astype(str)
        updated_df = updated_df.dropna(subset=["PMID"])
        
        # Sort and clean duplicates
        updated_df = (
            updated_df.assign(pmid_num=pd.to_numeric(updated_df["PMID"], errors='coerce'))
            .sort_values(by="pmid_num")
            .drop(columns=["pmid_num"])
            .drop_duplicates(subset="PMID", keep="last")
        )
        
        # Ensure parent folder exists
        os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
        updated_df.to_csv(OUTPUT_FILE, index=False)
        print(f"🎉 Incremental update complete! Scraped results saved to '{OUTPUT_FILE}' ({len(updated_df):,} total PMIDs).")
        return updated_df
    else:
        print("⚠️ No new metadata could be fetched.")
        return existing_df

# ==============================================================================
# 🔍 FLAG AND EXPORT INCOMPLETE INFORMATION FOR MANUAL VERIFICATION
# ==============================================================================
def verify_scraped_data(df):
    """
    Checks titles and abstracts for completeness (e.g. ending punctuation)
    and saves flagged PMIDs to specific files for manual review.
    """
    print("\n🔍 Conducting completeness checks on scraped data...")
    
    # 1. Incomplete Title Check (No ending period or question mark)
    df["Title"] = df["Title"].fillna("N/A")
    missing_title_mask = df["Title"].isna() | ~df["Title"].str.endswith((".", "?")).fillna(False)
    incomplete_titles = df[missing_title_mask]
    
    print(f"   -> Found {len(incomplete_titles)} titles that might be incomplete or need manual check.")
    incomplete_titles[["PMID", "Title", "Journal", "Year"]].to_csv(MISSING_TITLE_FILE, index=False)
    print(f"      Saved flagged titles to '{MISSING_TITLE_FILE}'")
    
    # 2. Incomplete Abstract Check (No standard ending)
    df["Abstract"] = df["Abstract"].fillna("")
    endings = ('.', '?', 'available', ')', 'Review', '...')
    missing_abstract_mask = ~df["Abstract"].str.endswith(endings)
    incomplete_abstracts = df[missing_abstract_mask]
    
    print(f"   -> Found {len(incomplete_abstracts)} abstracts that might be incomplete.")
    incomplete_abstracts[["PMID", "Title", "Abstract"]].to_csv(MISSING_ABSTRACT_FILE, index=False)
    print(f"      Saved flagged abstracts to '{MISSING_ABSTRACT_FILE}'")

# ==============================================================================
# 🚀 MAIN SCRIPT EXECUTION
# ==============================================================================
def main():
    try:
        # Extract unique PMIDs from database
        pmid_list = extract_unique_db_pmids(INPUT_DATABASE)
        
        # Webscrape metadata
        scraped_df = fetch_pubmed_metadata(pmid_list)
        
        # Perform validation checks
        verify_scraped_data(scraped_df)
        
    except Exception as e:
        print(f"❌ Execution failed: {e}")

if __name__ == "__main__":
    main()
