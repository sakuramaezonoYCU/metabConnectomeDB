import pandas as pd
import os
import sys
import json
import re
import argparse
import numpy as np
import time
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup

if 'scripts' not in sys.path and '.' not in sys.path:
    sys.path.append('scripts')
    sys.path.append('.')

# Load explicit parameters from config
try:
    with open('input/pipeline.config.json', 'r') as cf:
        pipeline_config = json.load(cf)
except FileNotFoundError:
    print("Error: input/pipeline.config.json not found.")
    sys.exit(1)

# Strict config extraction (will hard crash via KeyError if missing)
_p45 = pipeline_config["PHASE_4_5_META_VALIDATION"]
SKEW_THRESHOLD = _p45["SKEW_THRESHOLD"]
SUBCLONE_SD_MULTIPLIER = _p45["SUBCLONE_SD_MULTIPLIER"]

try:
    from pan_cancer_config import ANALYSIS_SUFFIX, CANCER_CAP
    CANCERS = list(CANCER_CAP.keys())
except Exception as e:
    raise ImportError(f"Failed to load dynamic parameters from pan_cancer_config. Hardcoding parameters is strictly prohibited. Error: {e}")

MD_OUT_PATH = 'output/AI_summary_and_insights.md'
PAPERS_DIR = 'papers'
NCBI_API_KEY_FILE = 'input/ncbi_api_key.txt'

def df_to_markdown(df):
    header = '| ' + ' | '.join(df.columns) + ' |'
    separator = '| ' + ' | '.join(['---'] * len(df.columns)) + ' |'
    rows = []
    for _, row in df.iterrows():
        formatted_row = []
        for x in row.values:
            if isinstance(x, (int, np.integer)):
                formatted_row.append(f"{x:,}")
            elif isinstance(x, (float, np.floating)) and not pd.isna(x):
                if x.is_integer():
                    formatted_row.append(f"{int(x):,}")
                else:
                    formatted_row.append(str(x))
            else:
                val = str(x)
                if '|' in val:
                    val = re.sub(r'\s*\|+\s*', ', ', val)
                formatted_row.append(val)
        rows.append('| ' + ' | '.join(formatted_row) + ' |')
    return '\n'.join([header, separator] + rows)

def scrape_notebook_output(html_path, extractors):
    if not os.path.exists(html_path):
        return {k: f"(Not yet generated: {os.path.basename(html_path)})" for k in extractors}
    with open(html_path, 'r', encoding='utf-8') as f:
        soup = BeautifulSoup(f, 'html.parser')
        
    text = soup.get_text()
    
    results = {}
    for k, pattern in extractors.items():
        match = re.search(pattern, text)
        if match:
            val = match.group(1).strip()
            clean_val = val.replace(',', '')
            if clean_val.isdigit():
                val = f"{int(clean_val):,}"
            results[k] = val
        else:
            results[k] = f"(Could not find metric for pattern: {pattern[:30]}...)"
    return results

CONTEXT_FILE = 'output/.cumulative_ai_context.json'
cumulative_ai_context = []
if os.path.exists(CONTEXT_FILE):
    try:
        with open(CONTEXT_FILE, 'r') as f:
            cumulative_ai_context = json.load(f)
    except Exception as e:
        print(f"Warning: Failed to load existing context: {e}")

# ==============================================================================
# 📄 GEMINI FILE API: Upload papers/ PDFs once and cache handles
# ==============================================================================
_uploaded_paper_handles = None  # Global cache

def _get_gemini_client():
    """Returns a configured Gemini client. Hard-crashes if secret is missing."""
    secret_path = os.path.join('input', '.geminiSecret')
    if not os.path.exists(secret_path):
        return None
    with open(secret_path, 'r') as f:
        api_key = f.read().strip()
    if not api_key or "PASTE_YOUR_GEMINI_API_KEY_HERE" in api_key:
        return None
    from google import genai
    return genai.Client(api_key=api_key)

def upload_papers_once(client):
    """Uploads all PDFs from papers/ directory to Gemini File API. Cached globally."""
    global _uploaded_paper_handles
    if _uploaded_paper_handles is not None:
        return _uploaded_paper_handles

    _uploaded_paper_handles = []
    if not os.path.isdir(PAPERS_DIR):
        print(f"    [i] No '{PAPERS_DIR}/' directory found. Skipping literature upload.")
        return _uploaded_paper_handles

    pdf_files = sorted([f for f in os.listdir(PAPERS_DIR) if f.endswith('.pdf')])
    if not pdf_files:
        print(f"    [i] No PDFs found in '{PAPERS_DIR}/'. Skipping literature upload.")
        return _uploaded_paper_handles

    print(f"    [📄] Uploading {len(pdf_files)} reference PDFs to Gemini File API...")
    for pdf_name in pdf_files:
        pdf_path = os.path.join(PAPERS_DIR, pdf_name)
        try:
            handle = client.files.upload(file=pdf_path)
            _uploaded_paper_handles.append(handle)
            print(f"        ✓ Uploaded: {pdf_name}")
            time.sleep(4)  # Rate limit (15 RPM free tier)
        except Exception as e:
            print(f"        ✗ Failed to upload {pdf_name}: {e}")
    return _uploaded_paper_handles

def _parse_html_to_text(html_path):
    """
    Parses an HTML file locally using BeautifulSoup to extract only text and tables.
    Strips all scripts, styles, and embedded images (which consume millions of tokens).
    """
    try:
        from bs4 import BeautifulSoup
        with open(html_path, 'r', encoding='utf-8') as f:
            html_content = f.read()
        soup = BeautifulSoup(html_content, 'html.parser')
        # Remove non-text elements
        for script in soup(["script", "style", "img", "svg"]):
            script.extract()
        # Extract text with line breaks
        text = soup.get_text(separator='\n', strip=True)
        return text
    except Exception as e:
        print(f"        ✗ Failed to parse HTML {os.path.basename(html_path)}: {e}")
        return ""


# ==============================================================================
# 🔬 PMID VERIFICATION: NCBI E-Fetch + Semantic Validation
# ==============================================================================

def _load_ncbi_api_key():
    """Returns the NCBI API key if available, else None."""
    if os.path.exists(NCBI_API_KEY_FILE):
        try:
            with open(NCBI_API_KEY_FILE, 'r') as f:
                key = f.read().strip()
            return key if key else None
        except Exception:
            return None
    return None

def _extract_pmid_claims(text):
    """
    Extracts (pmid, surrounding_claim_sentence) pairs from AI-generated text.
    Returns a list of dicts: [{'pmid': '12345678', 'claim': 'The sentence containing the PMID.', 'match': 'PMID:12345678'}]
    """
    results = []
    # Find all PMID references in various formats
    pmid_pattern = re.compile(r'(?:PMID\s*:?\s*(\d{5,})|PubMed\s*(?:ID)?\s*:?\s*(\d{5,}))', re.IGNORECASE)
    
    # Split text into sentences for claim extraction
    sentences = re.split(r'(?<=[.!?])\s+', text)
    
    for sentence in sentences:
        for match in pmid_pattern.finditer(sentence):
            pmid = match.group(1) or match.group(2)
            results.append({
                'pmid': pmid,
                'claim': sentence.strip(),
                'match': match.group(0)
            })
    return results

def _fetch_pubmed_abstracts(pmid_list, ncbi_api_key=None):
    """
    Batch-fetches Title and Abstract from NCBI E-Fetch for a list of PMIDs.
    Returns a dict: {pmid: {'title': ..., 'abstract': ...}} or {pmid: None} if not found.
    """
    import requests
    
    if not pmid_list:
        return {}
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    results = {}
    
    # Batch up to 50 at a time
    batch_size = 50
    for i in range(0, len(pmid_list), batch_size):
        batch = pmid_list[i:i + batch_size]
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
            root = ET.fromstring(response.text)
            
            found_pmids = set()
            for article in root.findall(".//PubmedArticle"):
                pmid_node = article.find(".//MedlineCitation/PMID")
                if pmid_node is None or not pmid_node.text:
                    continue
                pmid = pmid_node.text.strip()
                found_pmids.add(pmid)
                
                title = article.findtext(".//ArticleTitle", default="N/A").strip()
                abstract_texts = []
                for abs_text in article.findall(".//AbstractText"):
                    if abs_text.text:
                        abstract_texts.append(abs_text.text.strip())
                abstract = " ".join(abstract_texts) if abstract_texts else "No abstract available"
                
                results[pmid] = {'title': title, 'abstract': abstract}
            
            # Mark PMIDs not found as None
            for pmid in batch:
                if pmid not in found_pmids:
                    results[pmid] = None
                    
        except Exception as e:
            print(f"        [!] NCBI E-Fetch error: {e}")
            for pmid in batch:
                if pmid not in results:
                    results[pmid] = None
        
        # Respect NCBI rate limits
        delay = 0.12 if ncbi_api_key else 0.35
        time.sleep(delay)
    
    return results

def _batch_semantic_verify_pmids(client, items_to_verify):
    """
    Uses a single secondary Gemini call to batch verify multiple PMIDs.
    items_to_verify is a list of dicts: [{'pmid': '...', 'claim': '...', 'title': '...', 'abstract': '...'}]
    Returns a dict mapping PMID to boolean: {'12345678': True, '87654321': False}
    """
    if not items_to_verify:
        return {}
        
    from google.genai import types
    import json
    
    prompt = "You are a scientific citation auditor. Your ONLY job is to determine if cited papers actually support the claims made about them.\n\n"
    prompt += "Below is a list of claims and the corresponding cited paper's title and abstract.\n"
    prompt += "For each item, answer YES if the abstract provides evidence that supports or is directly relevant to the claim, or NO if it does not.\n\n"
    
    for i, item in enumerate(items_to_verify):
        prompt += f"--- ITEM {i+1} ---\n"
        prompt += f"PMID: {item['pmid']}\n"
        prompt += f"CLAIM: \"{item['claim']}\"\n"
        prompt += f"PAPER TITLE: {item['title']}\n"
        prompt += f"PAPER ABSTRACT: {item['abstract']}\n\n"
        
    prompt += "Output your response as a valid JSON object where the keys are the PMIDs (as strings) and the values are boolean (true for YES, false for NO). "
    prompt += "Do not include markdown blocks or any other text. Output strict JSON only, e.g. {\"12345678\": true, \"87654321\": false}."

    try:
        response = client.models.generate_content(
            model='gemini-2.5-flash',
            contents=prompt,
            config=types.GenerateContentConfig(
                temperature=0.0,
                response_mime_type="application/json",
            )
        )
        
        # Parse the JSON response
        text = response.text.strip()
        if text.startswith("```json"):
            text = text[7:]
        if text.endswith("```"):
            text = text[:-3]
        
        results = json.loads(text.strip())
        
        # Ensure all PMIDs are strings in the results
        return {str(k): bool(v) for k, v in results.items()}
    except Exception as e:
        print(f"        [!] Batch semantic verification failed: {e}")
        # On failure, conservatively reject all
        return {str(item['pmid']): False for item in items_to_verify}

def verify_and_format_pmids(text, client):
    """
    Post-processes AI-generated text to:
    1. Extract all cited PMIDs
    2. Verify they exist via NCBI E-Fetch
    3. Semantically verify the claim matches the paper's abstract (in one batch call)
    4. Format valid PMIDs as clickable links with titles
    5. Flag invalid/hallucinated PMIDs
    
    Returns the cleaned text with verified links and a verification summary.
    """
    claims = _extract_pmid_claims(text)
    if not claims:
        return text, "No PMIDs cited."
    
    # Deduplicate PMIDs for batch fetching
    unique_pmids = list(set(c['pmid'] for c in claims))
    print(f"    [🔬] Verifying {len(unique_pmids)} unique PMIDs against NCBI PubMed...")
    
    ncbi_key = _load_ncbi_api_key()
    abstracts = _fetch_pubmed_abstracts(unique_pmids, ncbi_key)
    
    # Prepare batch for semantic verification
    items_to_verify = []
    replacements = {}
    verified_count = 0
    failed_count = 0
    verification_log = []
    
    for claim_info in claims:
        pmid = claim_info['pmid']
        claim_text = claim_info['claim']
        original_match = claim_info['match']
        
        paper_data = abstracts.get(pmid)
        
        if paper_data is None:
            # PMID does not exist in PubMed at all
            replacements[original_match] = f"~~{original_match}~~ [⚠️ PMID Not Found in PubMed]"
            failed_count += 1
            verification_log.append(f"  ✗ {original_match}: Does not exist in PubMed database")
        else:
            items_to_verify.append({
                'pmid': pmid,
                'claim': claim_text,
                'title': paper_data['title'],
                'abstract': paper_data['abstract'],
                'original_match': original_match
            })
            
    # Run batch semantic verification
    if items_to_verify:
        print(f"    [🔬] Running batch semantic verification for {len(items_to_verify)} citations...")
        time.sleep(2)  # Short delay before hitting Gemini again
        verification_results = _batch_semantic_verify_pmids(client, items_to_verify)
        
        for item in items_to_verify:
            pmid = item['pmid']
            original_match = item['original_match']
            title = item['title']
            is_valid = verification_results.get(pmid, False)
            
            if is_valid:
                clean_title = title.rstrip('.')
                link = f"[PMID:{pmid} - {clean_title}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)"
                replacements[original_match] = link
                verified_count += 1
                verification_log.append(f"  ✓ PMID:{pmid} — \"{clean_title}\" — Claim verified")
            else:
                replacements[original_match] = f"~~{original_match}~~ [⚠️ Citation does not support claim — removed]"
                failed_count += 1
                verification_log.append(f"  ✗ PMID:{pmid} — \"{title}\" — Claim NOT supported by abstract")
    
    # Apply replacements to text
    processed_text = text
    for original, replacement in replacements.items():
        processed_text = processed_text.replace(original, replacement)
    
    summary = f"PMID Verification: {verified_count} verified, {failed_count} removed/flagged out of {len(claims)} citations."
    print(f"    [🔬] {summary}")
    for log_line in verification_log:
        print(f"        {log_line}")
    
    return processed_text, summary


# ==============================================================================
# 🤖 GEMINI INTERPRETATION: Deep Research with File API + PMID Verification
# ==============================================================================

def ask_gemini_interpretation(markdown_text, phase, html_paths=None):
    """
    Queries Gemini for a scientific interpretation of the raw data.
    
    Enhancements over original:
    - Uploads PDFs from papers/ directory for deep biological cross-referencing
    - Uploads relevant HTML reports for the current phase
    - Post-processes response to verify all cited PMIDs via NCBI E-Fetch
    - Semantically validates that cited papers actually support the claims made
    """
    global cumulative_ai_context
    
    client = _get_gemini_client()
    if client is None:
        return "\n> [!WARNING]\n> **AI Interpretation Skipped**: `input/.geminiSecret` not found or invalid.\n"
         
    try:
        from google.genai import types
        
        # Upload reference papers (cached after first call)
        paper_handles = upload_papers_once(client)
        
        # Parse phase-specific HTML reports locally into text
        parsed_html_texts = ""
        if html_paths:
            existing_paths = [p for p in html_paths if os.path.exists(p)]
            if existing_paths:
                print(f"    [📊] Parsing {len(existing_paths)} HTML reports locally to extract text and tables...")
                for path in existing_paths:
                    parsed_text = _parse_html_to_text(path)
                    if parsed_text:
                        parsed_html_texts += f"\n\n--- EXTRACTED DATA FROM HTML REPORT: {os.path.basename(path)} ---\n"
                        # Limit to first 25000 chars per file just in case it's still large
                        parsed_html_texts += parsed_text[:25000] + "\n..."
        
        # Build cumulative context string (condensed to avoid token bloat)
        context_str = ""
        if cumulative_ai_context:
            context_str = "PREVIOUS PHASES CONTEXT & FINDINGS:\n"
            for past_phase in cumulative_ai_context:
                context_str += f"--- Phase {past_phase['phase']} ---\n"
                context_str += f"Raw Data Summary:\n{past_phase['data'][:500]}...\n"
                context_str += f"AI Interpretation:\n{past_phase['interpretation'][:1000]}...\n\n"
        
        # Build the file context description
        file_context_note = ""
        if paper_handles:
            file_context_note += f"\nYou have been provided {len(paper_handles)} reference PDF papers (via File API). "
            file_context_note += "Cross-reference your interpretation against findings in these papers. "
            file_context_note += "When citing findings from these papers, use their actual PMIDs.\n"
        if parsed_html_texts:
            file_context_note += f"\nBelow, you are provided with text extracted from {len(existing_paths)} Jupyter notebook HTML reports. "
            file_context_note += "These contain the FULL analytical outputs including plots, statistical tests, and intermediate results. "
            file_context_note += "Reference specific results, figures, and statistical outputs from these reports in your interpretation.\n"
        
        prompt = f"""You are an expert computational biologist analyzing single-cell metabolism and pan-cancer data.
Below is the raw markdown data summary from Phase {phase} of our metabConnectomeDB pipeline.

{file_context_note}

{context_str}

YOUR TASK:
Provide a deeply scientific, data-driven interpretation of the results below. You MUST:
1. Summarize the key quantitative findings (do NOT just copy tables — synthesize patterns and highlight the most significant results).
2. Explain the biological significance of these findings, cross-referencing with the provided PDF literature and HTML reports where applicable.
3. When citing external literature to support biological claims, format PMIDs exactly as: PMID:12345678
   - Only cite PMIDs you are confident are real and relevant.
   - Every PMID you cite will be programmatically verified against PubMed. Hallucinated or irrelevant citations will be automatically removed.

SCIENTIFIC INTEGRITY POLICY (ABSOLUTE):
- DO NOT fabricate, guess, or mock any data, metrics, or biological mechanisms.
- If the data is sparse or inconclusive, state that explicitly.
- Do NOT claim causation from correlational data.
- Only reference biological mechanisms that are directly supported by the data tables, the provided HTML reports, or the provided PDF literature.

RESPONSE FORMAT (use these exact headers):
### 1. NOVEL FINDINGS
[Synthesize the most important findings. Reference specific metrics from the data. Cross-reference with the provided literature PDFs and HTML reports.]

### 2. PROPOSED RESEARCH QUESTIONS
[List 2-3 deep biological questions raised by these specific results.]

### 3. SUGGESTED NEXT STEPS
[Actionable computational or experimental next steps grounded in the data.]

Data to interpret:
{markdown_text}

{parsed_html_texts}
"""
        
        # We only pass paper_handles into the contents array now, not HTMLs
        contents = paper_handles + [prompt]
        
        attempt = 1
        max_attempts = 5
        
        while attempt < max_attempts:
            try:
                response = client.models.generate_content(
                    model='gemini-2.5-flash',
                    contents=contents,
                    config=types.GenerateContentConfig(
                        safety_settings=[
                            types.SafetySetting(category=types.HarmCategory.HARM_CATEGORY_DANGEROUS_CONTENT, threshold=types.HarmBlockThreshold.BLOCK_NONE),
                            types.SafetySetting(category=types.HarmCategory.HARM_CATEGORY_HATE_SPEECH, threshold=types.HarmBlockThreshold.BLOCK_NONE),
                            types.SafetySetting(category=types.HarmCategory.HARM_CATEGORY_HARASSMENT, threshold=types.HarmBlockThreshold.BLOCK_NONE),
                            types.SafetySetting(category=types.HarmCategory.HARM_CATEGORY_SEXUALLY_EXPLICIT, threshold=types.HarmBlockThreshold.BLOCK_NONE),
                        ]
                    )
                )
                break  # Success!
            except Exception as e:
                error_str = str(e)
                attempt += 1
                if any(code in error_str for code in ['403', '429', '503', '500', '502', '504']):
                    backoff_time = min(15 * (2 ** (attempt - 1)), 300)
                    print(f"    [!] Gemini API error ({e}). Retrying in {backoff_time}s (Attempt {attempt}/{max_attempts})...")
                    time.sleep(backoff_time)
                else:
                    print(f"    [!] Gemini API Error ({e}). Retrying in 30s (Attempt {attempt}/{max_attempts})...")
                    time.sleep(30)
        
        if attempt == max_attempts:
            return f"\n> [!WARNING]\n> **AI Interpretation Failed**: Reached max attempts ({max_attempts}) due to API errors.\n"
        
        raw_response = response.text
        
        # === PMID Verification Pass ===
        print(f"    [🔬] Running PMID verification pass on Phase {phase} response...")
        verified_response, verification_summary = verify_and_format_pmids(raw_response, client)
        
        # Save context for future phases (use raw response for context, verified for output)
        cumulative_ai_context.append({
            'phase': phase,
            'data': markdown_text,
            'interpretation': raw_response
        })
        try:
            with open(CONTEXT_FILE, 'w') as f:
                json.dump(cumulative_ai_context, f)
        except Exception as e:
            print(f"Warning: Failed to save context: {e}")
        
        # Prefix every line with '> ' for markdown alert block
        formatted_text = "> " + verified_response.replace('\n', '\n> ')
        verification_note = f"\n> \n> ---\n> *{verification_summary}*"
        return f"\n> [!NOTE]\n> **Data-Driven AI Interpretation (with PMID Verification)**\n{formatted_text}{verification_note}\n"
        
    except ImportError:
        return "\n> [!WARNING]\n> **AI Interpretation Skipped**: `google-genai` library not installed. Please run `pip install google-genai`.\n"
    except Exception as e:
        return f"\n> [!WARNING]\n> **AI Interpretation Failed**: Error calling Gemini API: {e}\n"


def ask_gemini_batch_interpretation(phase_data_map):
    """
    Sends ALL phase data in a SINGLE API request and returns per-phase interpretations.
    
    phase_data_map: dict of {phase_num: (markdown_content, html_paths_list)}
    Returns: dict of {phase_num: formatted_interpretation_string}
    
    This approach uses exactly 1 API request instead of 6, completely bypassing
    the gemini-2.5-flash 20-request/day free tier limit.
    """
    client = _get_gemini_client()
    if client is None:
        return {p: "\n> [!WARNING]\n> **AI Interpretation Skipped**: `input/.geminiSecret` not found or invalid.\n" for p in phase_data_map}
    
    try:
        from google.genai import types
        
        # Upload reference papers (cached after first call)
        paper_handles = upload_papers_once(client)
        
        # Build the mega-prompt with all phase data
        all_phases_text = ""
        for phase_num in sorted(phase_data_map.keys()):
            content, html_paths = phase_data_map[phase_num]
            all_phases_text += f"\n\n{'='*80}\n"
            all_phases_text += f"PHASE {phase_num} DATA\n"
            all_phases_text += f"{'='*80}\n\n"
            all_phases_text += content
            
            # Parse and append HTML text for this phase
            if html_paths:
                existing = [p for p in html_paths if os.path.exists(p)]
                if existing:
                    print(f"    [📊] Parsing {len(existing)} HTML reports for Phase {phase_num}...")
                    for path in existing:
                        parsed_text = _parse_html_to_text(path)
                        if parsed_text:
                            all_phases_text += f"\n--- EXTRACTED DATA FROM HTML REPORT: {os.path.basename(path)} ---\n"
                            all_phases_text += parsed_text[:25000] + "\n"
        
        phase_numbers = sorted(phase_data_map.keys())
        
        file_context_note = ""
        if paper_handles:
            file_context_note += f"\nYou have been provided {len(paper_handles)} reference PDF papers (via File API). "
            file_context_note += "Cross-reference your interpretation against findings in these papers. "
            file_context_note += "When citing findings from these papers, use their actual PMIDs.\n"
        
        prompt = f"""You are an expert computational biologist analyzing single-cell metabolism and pan-cancer data.
Below is the COMPLETE raw data from ALL {len(phase_numbers)} phases of our metabConnectomeDB pipeline.

{file_context_note}

YOUR TASK:
For EACH phase, provide a deeply scientific, data-driven interpretation. You MUST:
1. Synthesize the key quantitative findings for that phase (do NOT just copy tables — highlight the most significant results).
2. Explain the biological significance, cross-referencing with the provided PDF literature where applicable.
3. As you interpret later phases, actively reference and build upon findings from earlier phases (cumulative insight).
4. When citing external literature, format PMIDs exactly as: PMID:12345678
   - Only cite PMIDs you are confident are real and relevant.
   - Every PMID you cite will be programmatically verified against PubMed.

SCIENTIFIC INTEGRITY POLICY (ABSOLUTE):
- DO NOT fabricate, guess, or mock any data, metrics, or biological mechanisms.
- If the data is sparse or inconclusive, state that explicitly.
- Do NOT claim causation from correlational data.
- Only reference biological mechanisms that are directly supported by the data tables, the provided HTML reports, or the provided PDF literature.

RESPONSE FORMAT (CRITICAL — you MUST use these exact delimiters for each phase):

For each phase, output:

=== PHASE N INTERPRETATION ===
### 1. NOVEL FINDINGS
[Synthesize the most important findings for this phase. Reference specific metrics. Cross-reference with literature.]

### 2. PROPOSED RESEARCH QUESTIONS
[List 2-3 deep biological questions raised by this phase's results.]

### 3. SUGGESTED NEXT STEPS
[Actionable computational or experimental next steps grounded in the data.]
=== END PHASE N ===

Replace N with the actual phase number ({', '.join(str(p) for p in phase_numbers)}).
You MUST produce one delimited block for EACH phase listed above.

{all_phases_text}
"""
        
        contents = paper_handles + [prompt]
        
        print(f"    [🤖] Sending single batch request for {len(phase_numbers)} phases...")
        attempt = 1
        max_attempts = 5
        
        while attempt < max_attempts:
            try:
                response = client.models.generate_content(
                    model='gemini-2.5-flash',
                    contents=contents,
                    config=types.GenerateContentConfig(
                        safety_settings=[
                            types.SafetySetting(category=types.HarmCategory.HARM_CATEGORY_DANGEROUS_CONTENT, threshold=types.HarmBlockThreshold.BLOCK_NONE),
                            types.SafetySetting(category=types.HarmCategory.HARM_CATEGORY_HATE_SPEECH, threshold=types.HarmBlockThreshold.BLOCK_NONE),
                            types.SafetySetting(category=types.HarmCategory.HARM_CATEGORY_HARASSMENT, threshold=types.HarmBlockThreshold.BLOCK_NONE),
                            types.SafetySetting(category=types.HarmCategory.HARM_CATEGORY_SEXUALLY_EXPLICIT, threshold=types.HarmBlockThreshold.BLOCK_NONE),
                        ]
                    )
                )
                break
            except Exception as e:
                error_str = str(e)
                attempt += 1
                if any(code in error_str for code in ['403', '429', '503', '500', '502', '504']):
                    backoff_time = min(15 * (2 ** (attempt - 1)), 300)
                    print(f"    [!] Gemini API error ({e}). Retrying in {backoff_time}s (Attempt {attempt}/{max_attempts})...")
                    time.sleep(backoff_time)
                else:
                    print(f"    [!] Gemini API Error ({e}). Retrying in 30s (Attempt {attempt}/{max_attempts})...")
                    time.sleep(30)
        
        if attempt == max_attempts:
            return {p: "\n> [!WARNING]\n> **AI Interpretation Failed**: Reached max attempts due to API errors.\n" for p in phase_data_map}
        
        raw_response = response.text
        print(f"    [✓] Received batch response ({len(raw_response)} chars)")
        
        # === PMID Verification Pass (single pass on entire response) ===
        print(f"    [🔬] Running PMID verification pass on batch response...")
        verified_response, verification_summary = verify_and_format_pmids(raw_response, client)
        
        # Split the response into per-phase blocks
        results = {}
        for phase_num in phase_numbers:
            start_marker = f"=== PHASE {phase_num} INTERPRETATION ==="
            end_marker = f"=== END PHASE {phase_num} ==="
            
            start_idx = verified_response.find(start_marker)
            end_idx = verified_response.find(end_marker)
            
            if start_idx != -1 and end_idx != -1:
                phase_text = verified_response[start_idx + len(start_marker):end_idx].strip()
            elif start_idx != -1:
                # End marker missing — grab until next phase or end
                next_start = verified_response.find(f"=== PHASE", start_idx + len(start_marker))
                if next_start != -1:
                    phase_text = verified_response[start_idx + len(start_marker):next_start].strip()
                else:
                    phase_text = verified_response[start_idx + len(start_marker):].strip()
            else:
                phase_text = f"*Phase {phase_num} interpretation not found in batch response.*"
            
            # Format as markdown blockquote
            formatted_text = "> " + phase_text.replace('\n', '\n> ')
            verification_note = f"\n> \n> ---\n> *{verification_summary}*"
            results[phase_num] = f"\n> [!NOTE]\n> **Data-Driven AI Interpretation (with PMID Verification)**\n{formatted_text}{verification_note}\n"
        
        return results
        
    except ImportError:
        return {p: "\n> [!WARNING]\n> **AI Interpretation Skipped**: `google-genai` library not installed.\n" for p in phase_data_map}
    except Exception as e:
        return {p: f"\n> [!WARNING]\n> **AI Interpretation Failed**: Error calling Gemini API: {e}\n" for p in phase_data_map}


def append_to_md(content):
    if not os.path.exists(MD_OUT_PATH):
        with open(MD_OUT_PATH, 'w') as f:
            f.write('# AI Summary and Insights (Pipeline Execution)\n\n')
            f.write('This document provides a high-level summary of the latest `metabConnectomeDB` data integration and analytical results.\n\n')
    
    with open(MD_OUT_PATH, 'a') as f:
        f.write(content + '\n')

def build_phase_1():
    print("Building Phase 1 Summary...")
    content = '## Phase 1: Database Exploration and Reporting\n\n'
    content += '> [!NOTE]\n'
    content += '> Metabolites are mapped to the official HMDB names while genes are mapped to the official HGNC gene symbols.\n\n'
    
    # ---------------- metab_targetPair_analysis ----------------
    content += '### 1.1 Core Interaction Network Metrics\n'
    content += 'Source: `input/metab_target_database.csv`\n'
    content += 'Script: `scripts/metab_targetPair_analysis.ipynb`\n'
    content += 'Scraper: `scripts/tmp_build_md.py`\n'
    content += 'Output: `output/metab_targetPair_analysis_full_report.html`\n\n'
    
    s1_metrics = scrape_notebook_output(
        'output/metab_targetPair_analysis_full_report.html',
        {
            'total_pairs': r"Total unique interaction pairs in database:\s*([\d,]+)",
            'total_pmids': r"Total unique PMIDs cited:\s*([\d,]+)",
            'unannotated': r"Number of unannotated values:\s*([\d,]+ \(\s*[\d\.]+\s*%\))",
            'cancer_annotated': r"🎉 Successfully saved ([\d,]+) cancer-annotated interaction pairs",
            'hc_pairs': r"High-confidence CCC-ready pairs:\s*([\d,]+)",
            'sensor_unified': r"\*\s*Consolidated Unified \(with regex fallbacks\):\s*([\d,]+ / [\d,]+ pairs \(\s*[\d.]+\s*%\))",
            'sensor_receptors': r"\*\s*Receptor\s*:\s*([\d,]+ occurrences \(\s*[\d.]+\s*%\))",
            'sensor_enzymes': r"\*\s*Enzyme\s*:\s*([\d,]+ occurrences \(\s*[\d.]+\s*%\))"
        }
    )
    
    content += f"- **Total Unique Interaction Pairs:** {s1_metrics['total_pairs']}\n"
    content += f"- **High-Confidence CCC-Ready Pairs (3+ DBs):** {s1_metrics['hc_pairs']}\n"
    content += f"- **Cancer-Annotated Pairs:** {s1_metrics['cancer_annotated']}\n"
    content += f"- **Literature Evidence:** {s1_metrics['total_pmids']} unique PMIDs successfully matched.\n"
    content += f"- **Sensor Annotation Coverage:** {s1_metrics['sensor_unified']}\n"
    content += f"  - Receptors: {s1_metrics['sensor_receptors']}\n"
    content += f"  - Enzymes: {s1_metrics['sensor_enzymes']}\n\n"
    
    # ---------------- unique_metab_data_exploration ----------------
    content += '### 1.2 Metabolite Biochemical Catalog\n'
    content += 'Source: `input/unique_metab_database.csv`\n'
    content += 'Script: `scripts/unique_metab_data_exploration.ipynb`\n'
    content += 'Scraper: `scripts/tmp_build_md.py`\n'
    content += 'Output: `output/unique_metab_data_exploration_full_report.html`\n\n'
    
    s2_metrics = scrape_notebook_output(
        'output/unique_metab_data_exploration_full_report.html',
        {
            'unique_metabs': r"Dataset loaded:\s*(\d+)\s*unique HMDB metabolites",
            'tier_1': r"\*\s*Tier 1 \(High.*?\):\s*(\d+ metabolites \(\s*[\d.]+\s*%\))",
            'tier_2': r"\*\s*Tier 2 \(Medium.*?\):\s*(\d+ metabolites \(\s*[\d.]+\s*%\))",
            'tier_3': r"\*\s*Tier 3 \(Low.*?\):\s*(\d+ metabolites \(\s*[\d.]+\s*%\))",
            'small_mol': r"([0-9,]+/[0-9,]+ \(\s*[\d\.]+\s*%\)) metabolites are < 500 Da",
            'median': r"Median mass:\s*([0-9\.]+ Da)",
            'range': r"Range:\s*([0-9\.]+ \S+ [0-9\.]+ Da)",
            'sulfur': r"S:\s*(\d+ metabolites \(\s*[\d.]+\s*%\))",
            'nitrogen': r"N:\s*(\d+ metabolites \(\s*[\d.]+\s*%\))"
        }
    )
    
    content += f"- **Unique Metabolites:** {s2_metrics['unique_metabs']}\n"
    content += f"- **Confidence Tier Breakdown:**\n"
    content += f"  - Tier 1 (High): {s2_metrics['tier_1']}\n"
    content += f"  - Tier 2 (Medium): {s2_metrics['tier_2']}\n"
    content += f"  - Tier 3 (Low): {s2_metrics['tier_3']}\n"
    content += f"- **Biophysics:** {s2_metrics['small_mol']} | Median Mass: {s2_metrics['median']} | Range: {s2_metrics['range']}\n"
    content += f"- **Elemental Composition:** Nitrogen-containing: {s2_metrics['nitrogen']} | Sulfur-containing: {s2_metrics['sulfur']}\n\n"

    print("Querying Gemini API for Phase 1 Interpretation...")
    phase1_htmls = [
        'output/metab_targetPair_analysis_full_report.html',
        'output/unique_metab_data_exploration_full_report.html',
    ]
    return content, phase1_htmls


def build_phase_2():
    print("Building Phase 2 & 3 Summary (Per-Cancer Integration)...")
    content = '## Phase 2 & 3: Single-Cell Transcriptome Integration\n\n'
    
    for c in CANCERS:
        content += f"### {c.capitalize()} Single-Cell Profiling\n"
        
        cap = CANCER_CAP[c]
        html_integration = f"output/{c}_results/cancer_cellxgene_integration_{cap}.html"
        
        import glob
        
        # Integration Metrics
        s_int = scrape_notebook_output(
            html_integration,
            {
                'cells_dl': r"Retrieved metadata for ([\d,]+) cells",
                'genes_dl': r"Downloading expression matrix.*for (\d+) cells and (\d+) genes"
            }
        )
        
        # Find DE HTML
        html_prim_met_matches = glob.glob(f"output/{c}_results/primary_vs_metastasis_*_{cap}.html")
        report_path = html_prim_met_matches[0] if html_prim_met_matches else f"output/{c}_results/primary_vs_metastasis_MISSING.html"
        
        content += f'Source: `input/pipeline.config.json` (CELLxGENE Data)\n'
        content += f'Script: `scripts/run_cancer_pipeline.py` -> `scripts/{c}_primary_vs_metastasis.ipynb`\n'
        content += f'Scraper: `scripts/tmp_build_md.py`\n'
        content += f'Output: `{report_path}`\n\n'
        
        # DE Metrics
        s_de = scrape_notebook_output(
            report_path,
            {
                'cells_prim': r"Primary cells:\s*([\d,]+)",
                'cells_met': r"Metastatic cells:\s*([\d,]+)",
                'total_analyzed': r"Total Metabolic Target Genes analyzed:\s*([\d,]+)",
                'up_met': r"Up in Metastasis\s+([\d,]+)"
            }
        )
        
        content += f"- **Total Cells Analyzed:** {s_int['cells_dl']}\n"
        content += f"- **Primary/Metastatic Split:** {s_de['cells_prim']} Primary vs {s_de['cells_met']} Metastatic\n"
        content += f"- **Differential Expression:** Of {s_de['total_analyzed']} metabolic targets tested, **{s_de['up_met']}** were significantly upregulated in metastasis.\n\n"

        # Dynamically inject the new CSVs from the pipeline patches
        orphan_csv = f"output/{c}_results/immune_evasion_orphan_metabolic_candidates.csv"
        if os.path.exists(orphan_csv):
            # Strict read. If file exists but is empty/malformed, this will intentionally crash
            df_orphan = pd.read_csv(orphan_csv)
            if not df_orphan.empty:
                content += f"#### Orphan Metabolic Immune Evasion Candidates\n\n"
                content += df_to_markdown(df_orphan.head(5)) + "\n\n"

        upset_csv = f"output/{c}_results/primary_vs_metastatic_convergence_upset_{c}.csv"
        if os.path.exists(upset_csv):
            # Strict read. If file exists but is empty/malformed, this will intentionally crash
            df_upset = pd.read_csv(upset_csv)
            if not df_upset.empty:
                content += f"#### Metastatic Convergence (Niche Target Up-Regulation)\n\n"
                content += df_to_markdown(df_upset) + "\n\n"

    print("Querying Gemini API for Phase 2 Interpretation...")
    import glob
    phase2_htmls = []
    for c in CANCERS:
        cap = CANCER_CAP[c]
        phase2_htmls += glob.glob(f"output/{c}_results/primary_vs_metastasis_*_{cap}.html")
        phase2_htmls += glob.glob(f"output/{c}_results/orphan_immune_*_{cap}.html")
    return content, phase2_htmls

def build_phase_3():
    # Because Phase 2 & 3 are conceptually merged in the output as requested by the user,
    # Phase 3 simply outputs the summary tables from `generate_ai_summary_tables.py`
    print("Building Phase 3 Summary (Pan-Cancer Tables)...")
    content = '## Phase 3: Dataset Summary Tables\n\n'

    content += 'Source: `output/pan_cancer_meta_results/*.csv` and `output/ai_summary_tables/*.csv`\n'
    content += 'Script: `scripts/generate_ai_summary_tables.py`\n'
    content += 'Scraper: `scripts/tmp_build_md.py`\n'
    content += 'Output: Markdown tables generated below\n\n'

    content += '### 3.1 Cell Type Composition\n\n'
    counts_csv = f'output/pan_cancer_meta_results/cell_type_counts{ANALYSIS_SUFFIX}.csv'
    if os.path.exists(counts_csv):
        df_counts = pd.read_csv(counts_csv)
        if 'Dataset' in df_counts.columns:
            df_counts.rename(columns={'Dataset': 'Cancer'}, inplace=True)
        for col in df_counts.columns:
            if col != 'Cancer':
                df_counts[col] = df_counts[col].apply(lambda x: f"{int(x):,}" if pd.notnull(x) else x)
        content += df_to_markdown(df_counts) + '\n\n'
    else:
        content += f"*(Pending: {counts_csv})*\n\n"

    content += '### 3.1.5 Disease Counts (From HTML Scrape)\n\n'
    disease_csv = f'output/pan_cancer_meta_results/disease_counts_pan_cancer{ANALYSIS_SUFFIX}.csv'
    if os.path.exists(disease_csv):
        df_disease = pd.read_csv(disease_csv)
        
        if 'Cell_Count' in df_disease.columns:
            # Remove string artifacts like '...'
            df_disease = df_disease[~df_disease['Cell_Count'].astype(str).str.contains(r'^\.+$')]
            # Convert to numeric, coercing errors to NaN
            df_disease['Cell_Count'] = pd.to_numeric(df_disease['Cell_Count'], errors='coerce')
            # Drop NaNs and 0s
            df_disease = df_disease[df_disease['Cell_Count'] > 0]
            # Format nicely
            df_disease['Cell_Count'] = df_disease['Cell_Count'].astype(int).apply(lambda x: f"{x:,}")
            
        content += df_to_markdown(df_disease) + '\n\n'
    else:
        content += f"*(Pending: {disease_csv})*\n\n"

    content += '### 3.2 Cell-Cell Communication (CCC) Potential\n\n'
    df_ccc_path = f'output/ai_summary_tables/ccc_potential_summary{ANALYSIS_SUFFIX}.csv'
    if os.path.exists(df_ccc_path):
        df = pd.read_csv(df_ccc_path)
        content += df_to_markdown(df) + '\n\n'
    else:
        content += f"*(Pending: {df_ccc_path})*\n\n"

    content += '### 3.3 Immune Evasion and CCC Quantification (LIANA+ analysis)\n\n'
    quant_path = f'output/pan_cancer_meta_results/immune_evasion_ccc_quantification{ANALYSIS_SUFFIX}.csv'
    if os.path.exists(quant_path):
        df_quant = pd.read_csv(quant_path)
        content += df_to_markdown(df_quant) + '\n\n'
    else:
        content += f"*(Pending: {quant_path})*\n\n"

    print("Querying Gemini API for Phase 3 Interpretation...")
    phase3_htmls = [
        f'output/pan_cancer_meta_results/pan_cancer_meta_analysis_report.html',
    ]
    return content, phase3_htmls

def build_phase_4():
    print("Building Phase 4 Summary...")
    content = '## Phase 4: Pan-Cancer Meta-Analysis & Signatures\n\n'

    content += '### 4.1 Strict and Relaxed Intersection Metrics\n'
    content += 'Source: `output/*_results/metab_targetPair_analysis_DE_full.csv`\n'
    content += 'Script: `scripts/pan_cancer_meta_analysis.ipynb`\n'
    content += 'Scraper: `scripts/tmp_build_md.py`\n'
    content += 'Output: `output/pan_cancer_meta_results/pan_cancer_meta_analysis_report.html`\n\n'
    
    html_meta = f"output/pan_cancer_meta_results/pan_cancer_meta_analysis_report.html"
    s_meta = scrape_notebook_output(
        html_meta,
        {
            'strict_up': r"(\d+) genes commonly upregulated across all \d+ cancers",
            'relaxed_combos': r"(combinations of \d+ out of \d+ cancers)"
        }
    )
    from pan_cancer_config import CANCERS_TO_RUN as CANCERS
    num_cancers = len(CANCERS)

    content += f"- **Strictly Conserved (All {num_cancers} Cancers):** {s_meta['strict_up']}\n"
    content += f"- **Relaxed Threshold:** Checked for {s_meta['relaxed_combos']}\n\n"
    content += '### 4.2 Conserved Metastatic Gene Signatures\n\n'
    content += '> [!IMPORTANT]\n'
    content += f'> **Methodology Note (MAX CANCER - 1 Rule):** If the strictly conserved signature across all {num_cancers} cancers yields exactly 0 genes, the pipeline automatically falls back to the union of all {num_cancers-1}-cancer combinations. This ensures a robust, viable meta-signature is still evaluated downstream, as per the established project methodology.\n\n'
    meta_results_dir = 'output/pan_cancer_meta_results'
    if os.path.exists(meta_results_dir):
        sig_files = [f for f in os.listdir(meta_results_dir) if f.startswith('pan_cancer_signature_') and f.endswith(f'{ANALYSIS_SUFFIX}.csv')]
        if not sig_files:
            content += '*No pan-cancer signature files found.*\n\n'
        else:
            for sig_file in sorted(sig_files):
                sig_path = os.path.join(meta_results_dir, sig_file)
                df_sig = pd.read_csv(sig_path)
                
                # IMPORTANT: Strict adherence to data columns (Fixing the Target -> Gene bug)
                if 'Gene' not in df_sig.columns and 'Target' not in df_sig.columns and 'Strictly_Conserved_Gene' not in df_sig.columns:
                    raise KeyError(f"CRITICAL ERROR: Expected 'Gene', 'Target', or 'Strictly_Conserved_Gene' column missing from {sig_file}. Halting to prevent data falsification.")
                    
                if 'Gene' in df_sig.columns:
                    gene_col = 'Gene'
                elif 'Strictly_Conserved_Gene' in df_sig.columns:
                    gene_col = 'Strictly_Conserved_Gene'
                else:
                    gene_col = 'Target'
                
                cancers_in_sig = sig_file.replace('pan_cancer_signature_', '').replace(f'{ANALYSIS_SUFFIX}.csv', '').split('_')
                num_genes = len(df_sig)
                
                content += f"Source: `output/*_results/metab_targetPair_analysis_DE_full.csv`\n"
                content += f"Script: `scripts/pan_cancer_meta_analysis.ipynb`\n"
                content += f"Scraper: `scripts/tmp_build_md.py`\n"
                content += f"Output: `{sig_path}`\n\n"
                
                content += f'**Conserved across {len(cancers_in_sig)} cancers ({", ".join(cancers_in_sig)}):** '
                content += f'{num_genes} target genes\n'
                
                if num_genes > 0:
                    genes = df_sig[gene_col].tolist()
                    content += f'**Genes**: {", ".join(genes)}\n\n'
                else:
                    content += '**Genes**: None\n\n'
    else:
        content += '*pan_cancer_meta_results directory not found.*\n\n'

    content += '### 4.3 Metastatic Enrichment Summary\n\n'
    df_met_path = f'output/ai_summary_tables/metastatic_enrichment_summary{ANALYSIS_SUFFIX}.csv'
    if os.path.exists(df_met_path):
        content += df_to_markdown(pd.read_csv(df_met_path)) + '\n\n'

    content += '### 4.4 Immune Evasion Summary\n\n'
    df_imm_path = f'output/ai_summary_tables/immune_evasion_summary{ANALYSIS_SUFFIX}.csv'
    if os.path.exists(df_imm_path):
        content += df_to_markdown(pd.read_csv(df_imm_path)) + '\n\n'

    content += '### 4.5 Cancer-Specific Unique Signatures\n\n'
    df_uniq_path = f'output/ai_summary_tables/cancer_specific_unique_signatures{ANALYSIS_SUFFIX}.csv'
    if os.path.exists(df_uniq_path):
        content += df_to_markdown(pd.read_csv(df_uniq_path)) + '\n\n'

    content += '### 4.6 Pan-Cancer Conserved CCC Links\n\n'
    ccc_links_path = f'output/pan_cancer_meta_results/pan_cancer_conserved_ccc_links{ANALYSIS_SUFFIX}.csv'
    if os.path.exists(ccc_links_path):
        content += df_to_markdown(pd.read_csv(ccc_links_path).head(10)) + '\n\n'

    content += '### 4.7 Conserved Gene Functional Annotations\n\n'
    gene_annot_path = f'output/ai_summary_tables/conserved_gene_directed_signature_annotation{ANALYSIS_SUFFIX}.csv'
    if os.path.exists(gene_annot_path):
        content += df_to_markdown(pd.read_csv(gene_annot_path)) + '\n\n'

    print("Querying Gemini API for Phase 4 Interpretation...")
    import glob
    phase4_htmls = [
        'output/pan_cancer_meta_results/pan_cancer_meta_analysis_report.html',
    ] + glob.glob(f'output/pan_cancer_meta_results/predictive_signature_biomarker_*.html')
    return content, phase4_htmls

def build_phase_5():
    print("Building Phase 5 Summary...")
    content = '## Phase 5: Druggability Axis Analysis\n\n'

    content += 'Source: `output/pan_cancer_meta_results/conserved_target_genes.csv`\n'
    content += 'Script: `scripts/druggability_axis_analysis.ipynb`\n'
    content += 'Scraper: `scripts/tmp_build_md.py`\n'
    content += f'Output: `output/druggability/druggability_axis_analysis{ANALYSIS_SUFFIX}.html`\n\n'

    html_drug = f"output/druggability/druggability_axis_analysis{ANALYSIS_SUFFIX}.html"
    
    content += '### 5.1 Pipeline Output Metrics\n\n'
    s_drug = scrape_notebook_output(
        html_drug,
        {
            'dg_genes': r"\[DGIdb\] Querying interactions for (\d+) genes",
            'dg_hits': r"\[DGIdb\] Found (\d+) interactions",
            'ot_hits': r"\[OpenTargets\] Found (\d+) known drug indications",
            'trac_genes': r"\[Tractability\] Found tractability data for (\d+) genes",
            'depmap_genes': r"\[DepMap\].*Found columns for:\s*\[(.*?)\]",
            'strict_interactions': r"Total drug interactions for strictly conserved genes:\s*(\d+)",
            'broad_interactions': r"Total drug interactions for broadly conserved genes:\s*(\d+)"
        }
    )
    
    content += f"- **DGIdb Queries:** Checked {s_drug['dg_genes']} genes, found {s_drug['dg_hits']} interactions.\n"
    content += f"- **OpenTargets Pipeline:** Identified {s_drug['ot_hits']} known drug indications.\n"
    content += f"- **Tractability Assessment:** Data found for {s_drug['trac_genes']} targets.\n"
    content += f"- **DepMap Co-dependency:** Computed synergy matrix for genes: [{s_drug['depmap_genes']}].\n"
    content += f"- **Cross-Pipeline Totals:**\n"
    content += f"  - Strictly Conserved Targets: {s_drug['strict_interactions']} total drug interactions.\n"
    content += f"  - Broadly Conserved Targets: {s_drug['broad_interactions']} total drug interactions.\n\n"

    content += '### 5.2 Druggable Targets Summary (CSV)\n\n'
    drug_csv = f'output/druggability/druggability_axis{ANALYSIS_SUFFIX}_drug_targets.csv'
    if os.path.exists(drug_csv):
        df_drugs = pd.read_csv(drug_csv)
        if not df_drugs.empty:
            required_cols = ['Target_Gene', 'Drug_Name', 'Interaction_Type']
            missing_cols = [c for c in required_cols if c not in df_drugs.columns]
            if missing_cols:
                raise KeyError(f"CRITICAL ERROR: Expected columns {missing_cols} missing from {drug_csv}. Halting to prevent falsification.")
                
            total_targets = df_drugs['Target_Gene'].nunique()
            total_drugs = df_drugs['Drug_Name'].nunique()
            
            content += f"- **Total Unique Druggable Targets:** {total_targets}\n"
            content += f"- **Total Unique Drugs Identified:** {total_drugs}\n\n"
            
            top_targets = df_drugs['Target_Gene'].value_counts().head(5)
            content += "**Top Targets by Number of Drugs:**\n"
            for t, c in top_targets.items():
                content += f"- {t} ({c} drugs)\n"
            content += "\n"
            
            top_interactions = df_drugs['Interaction_Type'].value_counts().head(5)
            content += "**Primary Interaction Types:**\n"
            for t, c in top_interactions.items():
                content += f"- {t} ({c} instances)\n"
            content += "\n"
        else:
            content += "*(No druggable targets found in dataset)*\n\n"
    else:
        content += f"*(Pending: {drug_csv})*\n\n"

    print("Querying Gemini API for Phase 5 Interpretation...")
    import glob
    phase5_htmls = glob.glob('output/druggability/druggability_axis_*.html')
    return content, phase5_htmls


def build_phase_6():
    import pandas as pd
    print("Building Phase 6 Summary...")
    content = '## Phase 6: Gene Signature Validation\n\n'

    content += 'Source: `output/pan_cancer_meta_results/pan_cancer_signature_*`\n'
    content += 'Script: `scripts/run_validation_phase.py` (which orchestrates TCGA, MassSpec, Spatial, and Single-Cell predictive validation)\n'
    content += 'Scraper: `scripts/tmp_build_md.py`\n'
    content += f'Output: `Multiple dynamic output directories`\n\n'

    content += '### 6.1 Clinical Prognostic Power (TCGA)\n\n'
    content += '**Methodology:** Each signature was evaluated against TCGA survival cohorts via Cox proportional hazard regression. '
    content += 'The table below reads the actual `true_signature_metrics.csv` files, counts the number of cohorts with P < 0.05, '
    content += 'calculates the median Hazard Ratio across cohorts, and ranks the top 5 signatures by prognostic significance. '
    content += 'HR > 1 = high expression associates with worse survival; HR < 1 = protective.\n\n'
    tcga_dir = 'output/tcga_validation'
    if os.path.exists(tcga_dir):
        sig_dirs = [d for d in os.listdir(tcga_dir) if os.path.isdir(os.path.join(tcga_dir, d))]
        tcga_summary = []
        for sig_name in sorted(sig_dirs):
            csv_path = os.path.join(tcga_dir, sig_name, 'true_signature_metrics.csv')
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                sig_count = len(df[df['P_Value'] < 0.05])
                median_hr = df['Hazard_Ratio'].median()
                tcga_summary.append({'Signature': sig_name, 'Significant_Cohorts (P<0.05)': sig_count, 'Median_Hazard_Ratio': median_hr})
        if tcga_summary:
            tcga_df = pd.DataFrame(tcga_summary).sort_values(by=['Significant_Cohorts (P<0.05)', 'Median_Hazard_Ratio'], ascending=[False, False])
            content += f"Evaluated {len(tcga_summary)} signatures. **Top 5 Prognostic Signatures:**\n\n"
            content += df_to_markdown(tcga_df.head(5)) + '\n\n'
        else:
            content += '*(No valid TCGA metrics found)*\n\n'
    else:
        content += '*(TCGA validation results pending)*\n\n'

    content += '### 6.2 Pre-Metastatic Single-Cell Subclone Identification\n\n'
    content += '**Methodology:** Pre-computed `pre_metastatic_subclone_summary` percentages are extracted, converted to numeric values, '
    content += 'and sorted descending to highlight the subclones with the highest pre-metastatic fraction. '
    content += 'Higher % indicates a larger fraction of cells in the primary tumor already expressing the metastatic gene signature.\n\n'
    subclone_csv = f'output/pan_cancer_meta_results/pre_metastatic_subclone_summary{ANALYSIS_SUFFIX}.csv'
    if os.path.exists(subclone_csv):
        df_subclone = pd.read_csv(subclone_csv)
        if 'Pre-Metastatic Subclone (%)' in df_subclone.columns:
            # Clean and sort
            df_subclone['Numeric_Pct'] = df_subclone['Pre-Metastatic Subclone (%)'].astype(str).str.replace(r'[~%]', '', regex=True)
            df_subclone['Numeric_Pct'] = pd.to_numeric(df_subclone['Numeric_Pct'], errors='coerce')
            df_subclone = df_subclone.sort_values(by='Numeric_Pct', ascending=False).drop(columns=['Numeric_Pct'])
            content += f"Evaluated {len(df_subclone)} subclone profiles. **Top 5 Pre-Metastatic Subclones:**\n\n"
            content += df_to_markdown(df_subclone.head(5)) + '\n\n'
        else:
            content += df_to_markdown(df_subclone.head(5)) + '\n\n'
    else:
        content += '*(Single-cell subclone validation results pending)*\n\n'

    content += '### 6.3 Spatial Transcriptomics Validation (Visium)\n\n'
    content += '**Methodology:** Each signature\'s spatial coherence is assessed using Visium 10x spatial transcriptomics data. '
    content += 'Moran\'s I (spatial autocorrelation) and its P-value are computed per tissue section. '
    content += 'The table counts samples with significant spatial clustering (Moran\'s I > 0.1 and P < 0.05) and reports the average Moran\'s I. '
    content += '**Interpretation:** Moran\'s I ranges from -1 to 1. Values near 0 indicate random distribution. '
    content += 'Values > 0.1 = moderate clustering; > 0.3 = strong clustering; > 0.5 = very strong spatial co-localization of the signature genes within the tissue architecture.\n\n'
    spatial_dir = 'output/spatial_verification'
    if os.path.exists(spatial_dir):
        sig_dirs = [d for d in os.listdir(spatial_dir) if os.path.isdir(os.path.join(spatial_dir, d))]
        spatial_summary = []
        for sig_name in sorted(sig_dirs):
            csv_path = os.path.join(spatial_dir, sig_name, 'visium_spatial_clustering_summary.csv')
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                if 'Morans_Pval' in df.columns and 'Morans_I' in df.columns:
                    sig_count = len(df[(df['Morans_Pval'] < 0.05) & (df['Morans_I'] > 0.1)])
                    avg_moran = df['Morans_I'].mean()
                    spatial_summary.append({'Signature': sig_name, 'Significant_Spatial_Samples': sig_count, 'Average_Morans_I': avg_moran})
        if spatial_summary:
            spatial_df = pd.DataFrame(spatial_summary).sort_values(by=['Significant_Spatial_Samples', 'Average_Morans_I'], ascending=[False, False])
            content += f"Evaluated {len(spatial_summary)} signatures. **Top 5 Spatially Conserved Signatures:**\n\n"
            content += df_to_markdown(spatial_df.head(5)) + '\n\n'
        else:
             content += '*(No valid Spatial metrics found)*\n\n'
    else:
        content += '*(Spatial validation results pending)*\n\n'

    content += '### 6.4 MassSpec Metabolomics Validation\n\n'
    content += '**Methodology:** Each signature\'s genes are cross-referenced against clinical MassSpec metabolomics data from `per_gene_metabolite_profile.csv`. '
    content += 'The table counts the total genes in the signature, how many have at least one metabolite directly detected in MassSpec clinical data (Detected_in_MassSpec > 0), '
    content += 'and computes the detection rate (%). Higher detection rate = stronger orthogonal metabolomic evidence for the signature\'s biological activity.\n\n'
    massspec_dir = 'output/massspec_metabolomics'
    if os.path.exists(massspec_dir):
        sig_dirs = [d for d in os.listdir(massspec_dir) if os.path.isdir(os.path.join(massspec_dir, d))]
        ms_summary = []
        for sig_name in sorted(sig_dirs):
            csv_path = os.path.join(massspec_dir, sig_name, 'per_gene_metabolite_profile.csv')
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                total_genes = len(df)
                detected = len(df[df['Detected_in_MassSpec'] > 0])
                pct = (detected / total_genes * 100) if total_genes > 0 else 0
                ms_summary.append({'Signature': sig_name, 'Total_Genes': total_genes, 'Genes_with_MassSpec_Hits': detected, 'Detection_Rate_Pct': pct})
        if ms_summary:
            ms_df = pd.DataFrame(ms_summary).sort_values(by='Detection_Rate_Pct', ascending=False)
            content += f"Evaluated {len(ms_summary)} signatures. **Top 5 Signatures by MassSpec Validation Rate:**\n\n"
            content += df_to_markdown(ms_df.head(5)) + '\n\n'
        else:
             content += '*(No valid MassSpec metrics found)*\n\n'
    else:
        content += '*(MassSpec validation results pending)*\n\n'

    content += '### 6.5 Advanced Downstream Meta-Analyses\n\n'
    import glob
    druggability_html = "output/druggability/druggability_axis_*.html"
    druggability_matches = glob.glob(druggability_html)
    druggability_path = druggability_matches[0] if druggability_matches else "output/druggability/druggability_axis_analysis_report.html"

    advanced_reports = {
        'Druggability Axis Analysis': druggability_path,
        'Ovarian Serotonin Immune Evasion': 'output/ovarian_serotonin_immune_evasion_report.html',
        'Oxygen Tension Analysis': 'output/oxygen_tension_analysis_report.html',
        'MITF Regulon Expansion': 'output/mitf_regulon_expansion_report.html',
        'Serotonin Axis Spatial Mapping': 'output/serotonin_axis_spatial_mapping_report.html',
        'Deep-Dive Conserved Metab Gene Sig': 'output/deepdive_conserved_metabGeneSig/deepdive_conserved_metabGeneSig_report.html',
        'CAMP Pan-Cancer Integration': 'output/camp_pancancer_integration_report.html',
        'Master Regulator Analysis': 'output/master_regulator_analysis_report.html'
    }
    
    for title, path in advanced_reports.items():
        if os.path.exists(path):
            content += f"- **{title}**: Generated Successfully (`{path}`)\n"
            
            if title == 'Druggability Axis Analysis':
                drug_targets = glob.glob('output/druggability/druggable_targets_strictly_conserved_*.csv')
                if drug_targets:
                    df_drug = pd.read_csv(sorted(drug_targets)[0])
                    if not df_drug.empty and 'Gene' in df_drug.columns:
                        gene_counts = df_drug['Gene'].value_counts()
                        content += f"  - *Strictly Conserved Druggable Targets*: {len(gene_counts)} genes targeted.\n"
                        for gene, count in gene_counts.items():
                            content += f"    - **{gene}**: {count} available drugs/compounds.\n"
                
                interactions = glob.glob('output/druggability/druggability_axis_*_drug_targets.csv')
                if interactions:
                    df_int = pd.read_csv(sorted(interactions)[0])
                    total_interactions = len(df_int)
                    content += f"  - *Total Database Interactions Found*: {total_interactions}\n"

            elif title == 'Master Regulator Analysis':
                tf_csv = 'output/master_regulator_results/tf_enrichment_results.csv'
                if os.path.exists(tf_csv):
                    import pandas as pd
                    df_tf = pd.read_csv(tf_csv)
                    content += "  - *Top 5 Transcription Factors*:\n"
                    top_tfs = df_tf.head(5)
                    for _, row in top_tfs.iterrows():
                        tf = row['Transcription Factor']
                        pval = row['Adjusted P-value']
                        score = row['Combined Score']
                        genes = row['Overlapping Genes']
                        content += f"    - **{tf}**: Adj P-Value = {pval:.2e}, Score = {score:.1f}. Regulates: {genes}\n"

            elif title == 'MITF Regulon Expansion':
                mitf_csvs = glob.glob('output/mitf_regulon/mitf_metabolic_regulon_pairs*.csv')
                if mitf_csvs:
                    import pandas as pd
                    df_mitf = pd.read_csv(sorted(mitf_csvs)[0])
                    if 'Metabolite_Name' in df_mitf.columns:
                        top_metabs = df_mitf['Metabolite_Name'].value_counts().head(5)
                        content += f"  - *Top 5 Metabolites Regulated by MITF Network*:\n"
                        for metab, count in top_metabs.items():
                            targets = df_mitf[df_mitf['Metabolite_Name'] == metab]['Target'].unique()
                            content += f"    - **{metab}**: Regulated via {len(targets)} enzymes (e.g., {', '.join(targets[:3])})\n"
                    
            elif title == 'Oxygen Tension Analysis':
                oxy_csvs = glob.glob('output/oxygen_tension/oxygen_tension_correlation_results.csv')
                if oxy_csvs:
                    import pandas as pd
                    df_oxy = pd.read_csv(sorted(oxy_csvs)[0])
                    content += f"  - *Metabolic Adaptation Profiles*:\n"
                    for _, row in df_oxy.iterrows():
                        content += f"    - **{row['Cancer']}**: OXPHOS/Glycolysis Ratio = {row['OXPHOS_Glycolysis_Ratio']:.3f}, HIF1 LFC = {row['Mean_HIF1_LFC']:.3f}\n"

            elif title == 'CAMP Pan-Cancer Integration':
                camp_csvs = glob.glob('output/camp_integration/metabolite_immune_covariation_*.csv')
                if camp_csvs:
                    import pandas as pd
                    df_camp = pd.read_csv(sorted(camp_csvs)[0])
                    if 'Unnamed: 0' in df_camp.columns:
                        df_camp = df_camp.set_index('Unnamed: 0')
                    corrs = df_camp.unstack().reset_index()
                    corrs.columns = ['Immune_Cell', 'Metabolite', 'Correlation']
                    
                    pos_corrs = corrs.sort_values('Correlation', ascending=False).head(3)
                    neg_corrs = corrs.sort_values('Correlation', ascending=True).head(3)
                    
                    content += f"  - *Top Positive (Immune-Activating) Covariations*:\n"
                    for _, row in pos_corrs.iterrows():
                        content += f"    - **{row['Metabolite']}** with **{row['Immune_Cell']}**: R = {row['Correlation']:.3f}\n"
                    
                    content += f"  - *Top Negative (Immune-Suppressing) Covariations*:\n"
                    for _, row in neg_corrs.iterrows():
                        content += f"    - **{row['Metabolite']}** with **{row['Immune_Cell']}**: R = {row['Correlation']:.3f}\n"

            elif title == 'Serotonin Axis Spatial Mapping':
                spat_csvs = glob.glob('output/serotonin_axis_spatial_mapping/primary_vs_metastatic_immune_evasion_summary.csv')
                if spat_csvs:
                    import pandas as pd
                    df_spat = pd.read_csv(sorted(spat_csvs)[0])
                    content += f"  - *Spatial TME Shifts (Primary → Metastatic)*:\n"
                    prim = df_spat[df_spat['Niche'] == 'Primary'].iloc[0] if not df_spat[df_spat['Niche'] == 'Primary'].empty else None
                    met = df_spat[df_spat['Niche'] == 'Metastatic'].iloc[0] if not df_spat[df_spat['Niche'] == 'Metastatic'].empty else None
                    
                    if prim is not None and met is not None:
                        content += f"    - **Macrophage Count**: {prim['Macrophage_Count']} → {met['Macrophage_Count']}\n"
                        content += f"    - **T/NK Cell Count**: {prim['T_NK_Count']} → {met['T_NK_Count']}\n"
                        content += f"    - **Exhaustion Score**: {prim['Mean_Exhaustion_Score']:.3f} → {met['Mean_Exhaustion_Score']:.3f}\n"
                        content += f"    - **Treg Score**: {prim['Mean_Treg_Score']:.3f} → {met['Mean_Treg_Score']:.3f}\n"

            elif title == 'Deep-Dive Conserved Metab Gene Sig':
                tcga_csvs = glob.glob('output/tcga_validation/*/true_signature_metrics.csv')
                if tcga_csvs:
                    import pandas as pd
                    df_tcga = pd.read_csv(sorted(tcga_csvs)[0])
                    content += f"  - *TCGA Survival Validation*:\n"
                    if 'P_Value' in df_tcga.columns:
                        top_tcga = df_tcga.sort_values('P_Value').head(3)
                        for _, row in top_tcga.iterrows():
                            # The signature metrics has "TCGA_Cohort", "Hazard_Ratio", "P_Value", "N_Samples"
                            # I'll just skip 'Significance_Level' since I didn't see it in the subset above
                            content += f"    - **{row['TCGA_Cohort']}**: HR = {row['Hazard_Ratio']:.3f}, P = {row['P_Value']:.2e}\n"
        else:
            content += f"- **{title}**: *(Pending: {path})*\n"
    content += '\n'

    # --- Specific Extraction for ML Prognostic Classifier ---
    import glob
    content += "### 6.6 ML Prognostic Classifiers (Clinical OS)\n\n"
    content += '**Methodology:** Cox Proportional Hazards models were trained on each gene signature against clinical Overall Survival data. '
    content += 'The optimal L1/L2 penalizer was selected via 5-fold cross-validation and applied to the final model to prevent overfitting. '
    content += '**Interpretation of Optimal Penalizer:** A low penalizer (e.g., 0.01 - 0.1) means the gene signature is robust and predictive without much shrinkage. '
    content += 'A high penalizer (e.g., 0.5 - 1.0) means the model required heavy regularization, suggesting the genes are either highly correlated with each other or weakly predictive on their own. '
    content += 'The table reads the actual `ml_metrics.csv` files from each cohort, sorts all evaluated signatures by Test C-Index, '
    content += 'and displays only the Top 3 best-performing signatures. '
    content += 'C-Index > 0.5 = better than random; > 0.6 = clinically meaningful; > 0.7 = strong discriminative power.\n\n'
    
    metric_files = glob.glob('output/ml_prognostic_results/*/*/ml_metrics.csv')
    if not metric_files:
        content += "*(No ML metrics CSV files found. Please ensure notebook execution finished.)*\n\n"
    else:
        for m_file in sorted(metric_files):
            cohort = os.path.basename(os.path.dirname(m_file))
            dataset_source = os.path.basename(os.path.dirname(os.path.dirname(m_file)))
            content += f"#### Cohort: {dataset_source.upper()} - {cohort.upper()}\n"
            content += f"*(Full metrics available in: `{m_file}`)*\n\n"
            
            # Strict read. If file exists but is empty/malformed, this will intentionally crash
            df = pd.read_csv(m_file)
            total_sigs = len(df)
            if total_sigs == 0:
                content += "*(No signatures were successfully evaluated)*\n\n"
                continue
            
            # Sort by Test C-Index descending and grab Top 3
            df_sorted = df.sort_values(by="Test_C_Index", ascending=False).head(3)
            
            # Calculate Median
            median_test = df["Test_C_Index"].median()
            median_cv = df["CV_C_Index"].median()
            
            content += f"**Aggregate Performance ({total_sigs} signatures evaluated):**\n"
            content += f"- Median CV C-Index: {median_cv:.3f}\n"
            content += f"- Median Test C-Index: {median_test:.3f}\n\n"
            
            content += "**Top 3 Performing Signatures (by Test C-Index):**\n"
            content += df_to_markdown(df_sorted) + "\n\n"

    # --- Specific Extraction for Ovarian Serotonin Immune Evasion ---
    ov_html = 'output/ovarian_serotonin_immune_evasion_report.html'
    content += "### 6.7 Ovarian Serotonin Immune Evasion\n\n"
    content += 'Source: `output/ovarian_results/` and Serotonin/TAM annotations\n'
    content += 'Script: `scripts/ovarian_serotonin_immune_evasion.ipynb`\n'
    content += 'Scraper: `scripts/tmp_build_md.py`\n'
    content += f'Output: `{ov_html}`\n\n'
    if os.path.exists(ov_html):
        s_ov = scrape_notebook_output(ov_html, {
            'total_cells': r"Total cells loaded:\s*([\d,]+)",
            'h5ad_shape': r"AnnData object with n_obs × n_vars = ([\d,]+ × [\d,]+)"
        })
        content += f"- **Total Spatial/scRNA Cells Loaded:** {s_ov['total_cells']}\n"
        content += f"- **Raw Dataset Shape:** {s_ov['h5ad_shape']}\n\n"
    else:
        content += f"*(Pending HTML: {ov_html})*\n\n"

    content += "**Primary vs Metastatic Immune Evasion Niche Summary:**\n"
    macs_csv = 'output/serotonin_axis_spatial_mapping/plot_data_macs.csv'
    tnk_csv = 'output/serotonin_axis_spatial_mapping/plot_data_tnk.csv'
    
    if os.path.exists(macs_csv) and os.path.exists(tnk_csv):
        df_macs = pd.read_csv(macs_csv)
        df_tnk = pd.read_csv(tnk_csv)
        
        try:
            plot_config = pipeline_config.get("SEROTONIN_AXIS", {}).get("PLOTTING_CONFIG", [])
            if plot_config:
                for p_cfg in plot_config:
                    score_col = p_cfg['score_col']
                    title = p_cfg['title']
                    target = p_cfg['cell_target']
                    
                    df = df_macs if target == 'df_macs' else df_tnk
                    if score_col in df.columns:
                        prim_mean = df[df['Niche'] == 'Primary'][score_col].mean()
                        met_mean = df[df['Niche'] == 'Metastatic'][score_col].mean()
                        content += f"- **{title}:** {prim_mean:.3f} (Primary) vs {met_mean:.3f} (Metastatic).\n"
                content += "\n"
            else:
                content += "*(No PLOTTING_CONFIG found in pipeline config)*\n\n"
        except Exception as e:
            content += f"*(Error computing summary statistics: {e})*\n\n"
    else:
        content += "*(Pending: plot_data_macs.csv and plot_data_tnk.csv)*\n\n"

    content += "**Visium Spatial Pearson Correlation (HTR7 vs Immune Suppression):**\n"
    v_csv = 'output/serotonin_axis_spatial_mapping/visium_immune_evasion_summary.csv'
    if os.path.exists(v_csv):
        v_df = pd.read_csv(v_csv)
        avg_r = v_df['Pearson_R'].mean()
        avg_p = v_df['Pearson_P'].mean()
        num_samples = len(v_df)
        content += f"- **Aggregated over {num_samples} Spatial Samples:** Mean Pearson R = {avg_r:.3f}, Mean P-value = {avg_p:.2e}\n\n"
    else:
        content += "*(Pending: visium_immune_evasion_summary.csv)*\n\n"

    print("Querying Gemini API for Phase 6 Interpretation...")
    import glob
    phase6_htmls = [
        'output/visium_spatial_validation_report.html',
        'output/ovarian_serotonin_immune_evasion_report.html',
        'output/oxygen_tension_analysis_report.html',
        'output/mitf_regulon_expansion_report.html',
        'output/serotonin_axis_spatial_mapping_report.html',
        'output/master_regulator_analysis_report.html',
        'output/camp_pancancer_integration_report.html',
        'output/deepdive_conserved_metabGeneSig/deepdive_conserved_metabGeneSig_report.html',
    ] + glob.glob('output/*_ml_prognostic_classifier_report.html')
    return content, phase6_htmls

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Modular AI Summary & Insight Builder")
    parser.add_argument('--phase', type=str, choices=['1', '2', '3', '4', '5', '6', 'all'],
                        help="Which phase results to build and interpret (1, 2, 3, 4, 5, 6, or all)")
    parser.add_argument('--reset', action='store_true', help="Delete existing MD file to start fresh")
    parser.add_argument('--check-failed', action='store_true', help="Check existing MD file for failed AI API phases and rerun only them")
    args = parser.parse_args()

    if not args.phase and not args.check_failed:
        parser.error("You must specify either --phase or --check-failed")

    phases_to_run = []

    if args.check_failed:
        if os.path.exists(MD_OUT_PATH):
            with open(MD_OUT_PATH, 'r') as f:
                content = f.read()
            for p in range(1, 7):
                if f"## Phase {p}:" in content:
                    # Check if AI failed in that block
                    # A naive check: if "**AI Interpretation Failed**" or missing entirely
                    phase_block_start = content.find(f"## Phase {p}:")
                    next_phase_start = content.find(f"## Phase {p+1}:", phase_block_start)
                    if next_phase_start == -1:
                        block = content[phase_block_start:]
                    else:
                        block = content[phase_block_start:next_phase_start]
                    
                    if "**AI Interpretation Failed**" in block:
                        print(f"Phase {p} AI interpretation failed previously. Re-queueing...")
                        phases_to_run.append(str(p))
                else:
                    print(f"Phase {p} missing entirely. Re-queueing...")
                    phases_to_run.append(str(p))
            
            if not phases_to_run:
                print("No failed AI phases detected. All good!")
                sys.exit(0)
        else:
            print(f"{MD_OUT_PATH} not found. Running all phases.")
            phases_to_run = ['all']
    else:
        if args.reset or args.phase == 'all':
            if os.path.exists(MD_OUT_PATH):
                os.remove(MD_OUT_PATH)
            if os.path.exists(CONTEXT_FILE):
                os.remove(CONTEXT_FILE)
        phases_to_run = ['all'] if args.phase == 'all' else [args.phase]

    phase_data = {}
    if '1' in phases_to_run or 'all' in phases_to_run:
        phase_data[1] = build_phase_1()
    if '2' in phases_to_run or 'all' in phases_to_run:
        phase_data[2] = build_phase_2()
    if '3' in phases_to_run or 'all' in phases_to_run:
        phase_data[3] = build_phase_3()
    if '4' in phases_to_run or 'all' in phases_to_run:
        phase_data[4] = build_phase_4()
    if '5' in phases_to_run or 'all' in phases_to_run:
        phase_data[5] = build_phase_5()
    if '6' in phases_to_run or 'all' in phases_to_run:
        phase_data[6] = build_phase_6()
    
    if not phase_data:
        print("No phase data collected. Exiting.")
        sys.exit(1)
    
    # === SINGLE BATCH API CALL ===
    print(f"\n{'='*60}")
    print(f"Sending ALL {len(phase_data)} phases to Gemini in a SINGLE API request...")
    print(f"{'='*60}")
    
    interpretations = ask_gemini_batch_interpretation(
        {p: (content, htmls) for p, (content, htmls) in phase_data.items()}
    )
    
    # === STITCH interpretations back into phase content and write MD ===
    for phase_num in sorted(phase_data.keys()):
        content, _ = phase_data[phase_num]
        interp = interpretations.get(phase_num, "\n> [!WARNING]\n> **AI Interpretation Missing**\n")
        content += interp + '\n---\n'
        append_to_md(content)
    
    print(f"\n✓ AI Summary written to {MD_OUT_PATH}")
