import os
import glob
import re
import pandas as pd
from bs4 import BeautifulSoup
from pan_cancer_config import CANCERS_TO_RUN, ANALYSIS_SUFFIX, get_liana_csv_paths

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(os.path.dirname(BASE_DIR), 'output', 'pan_cancer_meta_results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

results = []

def find_file(folder, suffix):
    search_path = os.path.join(os.path.dirname(BASE_DIR), 'output', folder, f"*{suffix}*")
    matches = glob.glob(search_path)
    if matches:
        return matches[0]
    return None

def parse_disease_counts(html_file, cancer, cap):
    """Parses the disease_counts table from the integration HTML."""
    if not os.path.exists(html_file):
        return pd.DataFrame(columns=[])
    
    try:
        tables = pd.read_html(html_file)
        for df in tables:
            if isinstance(df.columns, pd.MultiIndex):
                df.columns = [col[-1] if isinstance(col, tuple) else col for col in df.columns]
                
            cols_lower = [str(c).lower() for c in df.columns]
            if 'disease' in df.columns or 'disease' in cols_lower:
                clean_df = df.copy()
                clean_df.columns = ['Disease', 'Cell_Count'] if len(df.columns) == 2 else df.columns
                clean_df['Cancer'] = cancer.capitalize()
                return clean_df
    except Exception as e:
        print(f"[{cancer.capitalize()}] Warning: Could not parse disease table from {html_file}: {e}")
        return pd.DataFrame(columns=[])
        raise

def scrape_integration_html(html_file, ccc_file):
    """Scrapes cancer_cellxgene_integration HTML for missing metrics."""
    metrics = {
        'Total_Inferred_CCC_Links': '0',
        'Total_Cells': '0',
        'Cell_Types': '0',
        'Genes_gt_10pct_Detection': '0'
    }
    
    if ccc_file and os.path.exists(ccc_file):
        try:
            df = pd.read_csv(ccc_file)
            metrics['Total_Inferred_CCC_Links'] = str(len(df))
        except Exception as e:
            pass
            raise
            
    if not os.path.exists(html_file):
        return metrics

    with open(html_file, 'r', encoding='utf-8') as f:
        text = BeautifulSoup(f, 'html.parser').get_text()

    if metrics['Total_Inferred_CCC_Links'] == '0':
        m_ccc = re.search(r"LIANA\+ found ([\d,]+) interactions", text)
        if not m_ccc:
            m_ccc = re.search(r"Successfully inferred ([\d,]+) metabolic cell-cell communication links", text)
        if m_ccc: metrics['Total_Inferred_CCC_Links'] = m_ccc.group(1).replace(',', '')

    m_cells = re.search(r"Total cells:\s*([\d,]+)", text)
    if m_cells: metrics['Total_Cells'] = m_cells.group(1).replace(',', '')

    m_ctypes = re.search(r"Cell types:\s*([\d,]+)", text)
    if m_ctypes: metrics['Cell_Types'] = m_ctypes.group(1).replace(',', '')

    m_det = re.search(r"Genes with >10% detection:\s*([\d,]+)", text)
    if m_det: metrics['Genes_gt_10pct_Detection'] = m_det.group(1).replace(',', '')

    return metrics

def scrape_pvm_html(html_file, cancer):
    """Scrapes primary_vs_metastasis HTML for missing metrics and uses CSVs for LIANA."""
    metrics = {
        'Primary_LIANA_Targets': '0',
        'Meta_LIANA_Targets': '0',
        'Pan_Metastatic_Conserved_Targets': '0',
        'Site_Specific_Upregulated_Summary': 'None'
    }
    
    # Use CSV files to get exact LIANA target counts
    try:
        paths = get_liana_csv_paths(cancer)
        
        if paths['primary'] and os.path.exists(paths['primary']):
            df_prim = pd.read_csv(paths['primary'])
            prim_targs = set(df_prim['receptor_complex'].dropna().unique()).union(set(df_prim['ligand_complex'].dropna().unique()))
            metrics['Primary_LIANA_Targets'] = str(len(prim_targs))
            
        if paths['meta'] and os.path.exists(paths['meta']):
            df_meta = pd.read_csv(paths['meta'])
            meta_targs = set(df_meta['receptor_complex'].dropna().unique()).union(set(df_meta['ligand_complex'].dropna().unique()))
            metrics['Meta_LIANA_Targets'] = str(len(meta_targs))
    except Exception as e:
        print(f"[{cancer.capitalize()}] Warning: Could not read LIANA CSVs: {e}")
        raise
        
    if not os.path.exists(html_file):
        return metrics

    with open(html_file, 'r', encoding='utf-8') as f:
        text = BeautifulSoup(f, 'html.parser').get_text()

    if metrics['Primary_LIANA_Targets'] == '0':
        m_prim_liana = re.search(r"Primary LIANA results:\s*([\d,]+)", text)
        if m_prim_liana: metrics['Primary_LIANA_Targets'] = m_prim_liana.group(1).replace(',', '')

    if metrics['Meta_LIANA_Targets'] == '0':
        m_meta_liana = re.search(r"Meta LIANA results:\s*([\d,]+)", text)
        if m_meta_liana: metrics['Meta_LIANA_Targets'] = m_meta_liana.group(1).replace(',', '')

    m_cons = re.search(r"Pan-Metastatic Conserved Targets \((\d+) genes\)", text)
    if m_cons: metrics['Pan_Metastatic_Conserved_Targets'] = m_cons.group(1).replace(',', '')

    site_specific = []
    for match in re.finditer(r"-\s+(.*?):\s+(\d+)\s+upregulated metabolic targets", text):
        site = match.group(1).strip()
        count = match.group(2).strip()
        site_specific.append(f"{site}:{count}")
    
    if site_specific:
        metrics['Site_Specific_Upregulated_Summary'] = " | ".join(site_specific)

    return metrics

def scrape_omi_html(html_file):
    """Scrapes orphan_immune HTML for missing metrics."""
    metrics = {
        'Highly_Enriched_Immune_Targets': '0',
        'Overlapping_Orphan_Interactions': '0'
    }
    if not os.path.exists(html_file):
        return metrics

    with open(html_file, 'r', encoding='utf-8') as f:
        text = BeautifulSoup(f, 'html.parser').get_text()

    m_enrich = re.search(r"Found ([\d,]+) unique highly enriched immune target genes", text)
    if m_enrich: metrics['Highly_Enriched_Immune_Targets'] = m_enrich.group(1).replace(',', '')

    m_over = re.search(r"Found ([\d,]+) interactions overlapping across multiple tumor sites", text)
    if m_over: metrics['Overlapping_Orphan_Interactions'] = m_over.group(1).replace(',', '')

    return metrics

def quantify_immune_evasion():
    print("Quantifying Immune Evasion and CCC Microenvironments (including exhaustive HTML scraping)...")
    
    results_parent = os.path.dirname(BASE_DIR)
    pan_cancer_stats = []
    all_disease_counts = []
    
    for cancer in CANCERS_TO_RUN:
        prefix = cancer.lower()
        if prefix == 'ovary':
            prefix = 'ovarian'
        cancer_dir = os.path.join(results_parent, 'output', f"{prefix}_results")
        
        ccc_file = find_file(f"{prefix}_results", 'cellxgene_communication_potential.csv')
        int_html = glob.glob(os.path.join(cancer_dir, "cancer_cellxgene_integration*.html"))
        int_html = int_html[0] if int_html else ""
        int_metrics = scrape_integration_html(int_html, ccc_file)
        
        disease_df = parse_disease_counts(int_html, cancer, "")
        if not disease_df.empty:
            all_disease_counts.append(disease_df)
            
        pvm_html = glob.glob(os.path.join(cancer_dir, "primary_vs_metastasis_*.html"))
        pvm_html = pvm_html[0] if pvm_html else ""
        pvm_metrics = scrape_pvm_html(pvm_html, cancer)
        
        omi_html = glob.glob(os.path.join(cancer_dir, "orphan_immune_*.html"))
        omi_html = omi_html[0] if omi_html else ""
        omi_metrics = scrape_omi_html(omi_html)

        stats = {
            'Cancer': cancer.capitalize(),
            **int_metrics,
            **pvm_metrics,
            **omi_metrics
        }
        pan_cancer_stats.append(stats)
        
        print(f"[{cancer.capitalize()}] CCC Links: {int_metrics['Total_Inferred_CCC_Links']} | "
              f"Primary LIANA Targets: {pvm_metrics['Primary_LIANA_Targets']} | "
              f"Meta LIANA Targets: {pvm_metrics['Meta_LIANA_Targets']} | "
              f"Conserved Targets: {pvm_metrics['Pan_Metastatic_Conserved_Targets']} | "
              f"Overlapping Orphan: {omi_metrics['Overlapping_Orphan_Interactions']}")

    out_csv = os.path.join(OUTPUT_DIR, f'immune_evasion_ccc_quantification{ANALYSIS_SUFFIX}.csv')
    df_out = pd.DataFrame(pan_cancer_stats)
    df_out.to_csv(out_csv, index=False)
    
    if all_disease_counts:
        combined_disease = pd.concat(all_disease_counts, ignore_index=True)
        combined_disease.to_csv(os.path.join(OUTPUT_DIR, f'disease_counts_pan_cancer{ANALYSIS_SUFFIX}.csv'), index=False)
            
    print(f"Saved comprehensive HTML metrics quantification to {out_csv}")

if __name__ == "__main__":
    quantify_immune_evasion()
