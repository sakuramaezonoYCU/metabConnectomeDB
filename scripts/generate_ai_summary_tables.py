"""
Purpose: Aggregates all the extracted outputs into formatted CSVs in output/ai_summary_tables/ for the AI insights generation.
"""
import os
import sys
import pandas as pd
import numpy as np

# Try to import ANALYSIS_SUFFIX from config
try:
    from pan_cancer_config import ANALYSIS_SUFFIX
except ImportError:
    # If run from outside scripts dir, try adding scripts dir to path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.append(script_dir)
    try:
        from pan_cancer_config import ANALYSIS_SUFFIX
    except ImportError:
        ANALYSIS_SUFFIX = ''

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(os.path.dirname(BASE_DIR), 'output')
META_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'pan_cancer_meta_results')
SUMMARY_TABLES_DIR = os.path.join(OUTPUT_DIR, 'ai_summary_tables')

os.makedirs(SUMMARY_TABLES_DIR, exist_ok=True)

# Helper config for the 5 cancers
# NOTE: The size suffix might differ depending on the dataset used.
# Checking for 100k or 500k files.
def get_cancer_prefix(cancer_name):
    from pan_cancer_config import normalize_cancer_name
    return normalize_cancer_name(cancer_name)

def find_file(prefix, folder, pattern, cancer_name):
    from pan_cancer_config import CANCER_CAP, normalize_cancer_name
    normalized_cancer = normalize_cancer_name(cancer_name)
    cap = CANCER_CAP[normalized_cancer]
    folder_path = os.path.join(OUTPUT_DIR, folder)
    if not os.path.exists(folder_path):
        return None
    for f in os.listdir(folder_path):
        if pattern in f and f.endswith('.csv'):
            return os.path.join(folder_path, f)
    return None

def generate_metastatic_enrichment_table():
    from pan_cancer_config import CANCERS_TO_RUN
    cancers = [c.capitalize() for c in CANCERS_TO_RUN]
    data = []
    
    for cancer in cancers:
        prefix = get_cancer_prefix(cancer)
        folder = f"{prefix}_results"
        
        file_path = find_file(prefix, folder, 'DE_metabolic_targets', cancer)
        if not file_path:
            continue
            
        df = pd.read_csv(file_path)
        # Count Up, Down, Not Sig
        if 'Regulation' in df.columns:
            up = len(df[df['Regulation'] == 'Up in Metastasis'])
            down = len(df[df['Regulation'] == 'Up in Primary'])
            ns = len(df[df['Regulation'] == 'Not Significant'])
        else:
            # Fallback calculation if column missing
            sig_thresh = 0.05
            lfc_thresh = 0.5
            if 'pvals_adj' in df.columns and 'logfoldchanges' in df.columns:
                up = len(df[(df['pvals_adj'] < sig_thresh) & (df['logfoldchanges'] > lfc_thresh)])
                down = len(df[(df['pvals_adj'] < sig_thresh) & (df['logfoldchanges'] < -lfc_thresh)])
                ns = len(df) - up - down
            else:
                up, down, ns = 0, 0, 0
                
        ratio = up / max(down, 1) # Avoid division by zero
        
        # Fetch oxygen tension data
        from pan_cancer_config import CANCER_PO2_CSV_MAPPING, _project_root, normalize_cancer_name
        csv_path = os.path.join(_project_root(), "input", "pO2_guide_24588669.csv")
        o2_tumour = "N/A"
        o2_normal = "N/A"
        try:
            po2_df = pd.read_csv(csv_path)
            csv_name = CANCER_PO2_CSV_MAPPING.get(normalize_cancer_name(cancer), "")
            if csv_name and not po2_df.empty:
                row = po2_df[po2_df['Tumour type'] == csv_name]
                if not row.empty:
                    o2_tumour = row['Median % oxygen_tumour'].values[0]
                    o2_normal = row['Median % oxygen_normal'].values[0]
                    if pd.isna(o2_normal) and pd.notna(o2_tumour):
                        o2_normal = 6.0
        except Exception:
            pass
            raise
        
        data.append({
            'Cancer': cancer,
            'Up in Metastasis': up,
            'Up in Primary': down,
            'Not Significant': ns,
            'Metastasis/Primary Ratio': f"{ratio:.2f}x",
            'Tumor pO2 (%)': o2_tumour,
            'Normal pO2 (%)': o2_normal
        })
        
    if data:
        summary_df = pd.DataFrame(data)
        out_path = os.path.join(SUMMARY_TABLES_DIR, f'metastatic_enrichment_summary{ANALYSIS_SUFFIX}.csv')
        summary_df.to_csv(out_path, index=False)
        print(f"Generated: {out_path}")

def generate_immune_evasion_table():
    from pan_cancer_config import CANCERS_TO_RUN
    cancers = [c.capitalize() for c in CANCERS_TO_RUN]
    data = []
    
    for cancer in cancers:
        prefix = get_cancer_prefix(cancer)
        folder = f"{prefix}_results"
        
        file_path = find_file(prefix, folder, 'immune_evasion_orphan_metabolic_candidates', cancer)
        if not file_path:
            continue
            
        df = pd.read_csv(file_path)
        total = len(df)
        if total > 0 and 'gene' in df.columns:
            unique_targets = df['gene'].nunique()
        elif total > 0 and 'Target_Gene' in df.columns:
            unique_targets = df['Target_Gene'].nunique()
        elif total > 0 and 'Target' in df.columns:
            unique_targets = df['Target'].nunique()
        else:
            unique_targets = 0
            
        data.append({
            'Cancer': cancer,
            'Total Candidates': total,
            'Unique Targets': unique_targets
        })
        
    if data:
        summary_df = pd.DataFrame(data)
        out_path = os.path.join(SUMMARY_TABLES_DIR, f'immune_evasion_summary{ANALYSIS_SUFFIX}.csv')
        summary_df.to_csv(out_path, index=False)
        print(f"Generated: {out_path}")

def generate_ccc_potential_table():
    from pan_cancer_config import CANCERS_TO_RUN
    cancers = [c.capitalize() for c in CANCERS_TO_RUN]
    data = []
    
    for cancer in cancers:
        prefix = get_cancer_prefix(cancer)
        folder = f"{prefix}_results"
        
        file_path = find_file(prefix, folder, 'cellxgene_communication_potential', cancer)
        if not file_path:
            continue
            
        df = pd.read_csv(file_path)
        total_targets = len(df)
        if 'Number_of_Expressing_Cell_Types' in df.columns:
            max_ct = df['Number_of_Expressing_Cell_Types'].max()
            mean_ct = df['Number_of_Expressing_Cell_Types'].mean()
        elif 'N_Cell_Types' in df.columns:
            max_ct = df['N_Cell_Types'].max()
            mean_ct = df['N_Cell_Types'].mean()
        else:
            max_ct, mean_ct = 0, 0
            
        data.append({
            'Cancer': cancer,
            'Metabolic Target Genes': total_targets,
            'Max Cell Types/Target': max_ct,
            'Mean Cell Types/Target': f"{mean_ct:.1f}"
        })
        
    if data:
        summary_df = pd.DataFrame(data)
        out_path = os.path.join(SUMMARY_TABLES_DIR, f'ccc_potential_summary{ANALYSIS_SUFFIX}.csv')
        summary_df.to_csv(out_path, index=False)
        print(f"Generated: {out_path}")
        
def generate_annotated_signature():
    """Generates the signature annotation from data, completely removing AI hallucinated biological roles."""
    
    # 1. Load the 21 conserved genes
    conserved_csv = os.path.join(META_RESULTS_DIR, f'pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv')
    if not os.path.exists(conserved_csv):
        print(f"File not found: {conserved_csv}")
        return
        
    sig_df = pd.read_csv(conserved_csv)
    genes = sig_df['Strictly_Conserved_Gene'].tolist()
    import subprocess
    
    # The signature genes are dynamic per run, so we always re-fetch the latest 
    # OpenTargets and UniProt annotations for the current exact signature.
    print(f"Triggering scripts/fetch_uniprot_roles.py for {len(genes)} dynamic genes...")
    subprocess.run(["python", "scripts/fetch_uniprot_roles.py"])
    
    print(f"Triggering scripts/fetch_opentargets.py for {len(genes)} dynamic genes...")
    subprocess.run(["python", "scripts/fetch_opentargets.py"])

    # Load UniProt biological roles if available
    uniprot_roles = {}
    uniprot_csv = os.path.join(OUTPUT_DIR, f'uniprot_biological_roles{ANALYSIS_SUFFIX}.csv')
    if os.path.exists(uniprot_csv):
        u_df = pd.read_csv(uniprot_csv)
        for _, row in u_df.iterrows():
            uniprot_roles[row['Gene']] = row['UniProt_Biological_Role']
            
    # Load OpenTargets diseases if available
    opentargets_diseases = {}
    ot_csv = os.path.join(OUTPUT_DIR, f'opentargets_diseases{ANALYSIS_SUFFIX}.csv')
    if os.path.exists(ot_csv):
        ot_df = pd.read_csv(ot_csv)
        for _, row in ot_df.iterrows():
            opentargets_diseases[row['Gene']] = row['OpenTargets_Diseases']
            
    # 2. Load the primary database to extract real metabolic associations and literature PMIDs
    db_path = os.path.join(OUTPUT_DIR, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')
    if not os.path.exists(db_path):
        print(f"Database not found: {db_path}")
        return
        
    db_df = pd.read_csv(db_path, low_memory=False)
    
    # The database maps multiple targets per row, so we need to explode or search carefully
    # Target column can be comma/semicolon separated.
    db_df['Target'] = db_df['Target'].astype(str).str.split(r'[,;]')
    db_df = db_df.explode('Target')
    db_df['Target'] = db_df['Target'].str.strip()
    
    # Filter the database to just our 21 genes
    db_filtered = db_df[db_df['Target'].isin(genes)]
    
    results = []
    for gene in genes:
        role_text = uniprot_roles.get(gene, "Requires literature review")
        disease_text = opentargets_diseases.get(gene, "Not fetched")
        gepia_link = f"http://gepia.cancer-pku.cn/detail.php?gene={gene}"
        
        gene_data = db_filtered[db_filtered['Target'] == gene]
        if gene_data.empty:
            results.append({
                "Gene": gene,
                "Key Metabolite(s)": "Unknown in DB",
                "Source_Database": "",
                "Sensor_Type": "",
                "Rhea_Reaction": "",
                "PMID(s)": "",
                "Biological Role": role_text,
                "Top OpenTargets Diseases": disease_text,
                "MRCLinkDB_Disease": "",
                "GEPIA Link": gepia_link
            })
            continue
            
        # Get unique metabolites
        metabolites = gene_data['Metabolite_Name'].dropna().unique()
        met_str = ", ".join(metabolites)
        
        # Get Sensor Type
        sensor_types = gene_data['Sensor_Type'].dropna().unique()
        flat_types = set()
        for st in sensor_types:
            flat_types.update([s.strip() for s in str(st).split(',') if s.strip() != 'nan'])
        sensor_str = ", ".join(sorted(list(flat_types)))
        
        # Get unique rhea reactions from dynamically enriched columns
        if 'Rhea_enzyme product/substrate' in gene_data.columns:
            rhea = gene_data['Rhea_enzyme product/substrate'].dropna().unique()
        elif 'enzyme product/substrate' in gene_data.columns:
            rhea = gene_data['enzyme product/substrate'].dropna().unique()
        elif 'rhea_reaction' in gene_data.columns:
            rhea = gene_data['rhea_reaction'].dropna().unique()
        else:
            rhea = []
        rhea_str = "; ".join([str(x) for x in rhea if pd.notna(x) and str(x).strip() != ''])
        
        # Get unique PMIDs
        pmids = gene_data['PMID'].dropna().unique()
        pmids_str = ";".join([str(int(x)) for x in pmids if str(x).replace('.0','').isdigit()])
        
        # Get Diseases
        if 'Disease' in gene_data.columns:
            db_diseases = gene_data['Disease'].dropna().unique()
            db_diseases_str = ", ".join(sorted([str(d).strip() for d in db_diseases if str(d).strip() and str(d).lower() != 'nan']))
        else:
            db_diseases_str = ""
            
        # Get Source Database
        if 'database' in gene_data.columns:
            source_dbs = gene_data['database'].dropna().unique()
            source_db_str = ", ".join(sorted([str(d).strip() for d in source_dbs if str(d).strip() and str(d).lower() != 'nan']))
        else:
            source_db_str = ""
        
        results.append({
            "Gene": gene,
            "Key Metabolite(s)": met_str.lower(),
            "Source_Database": source_db_str,
            "Sensor_Type": sensor_str,
            "Rhea_Reaction": rhea_str,
            "PMID(s)": pmids_str,
            "Biological Role": role_text,
            "Top OpenTargets Diseases": disease_text,
            "MRCLinkDB_Disease": db_diseases_str,
            "GEPIA Link": gepia_link
        })
        
    df = pd.DataFrame(results)
    out_path = os.path.join(SUMMARY_TABLES_DIR, f'conserved_gene_directed_signature_annotation{ANALYSIS_SUFFIX}.csv')
    df.to_csv(out_path, index=False)
    print(f"Generated data-driven annotation: {out_path}")

def generate_dataset_overview(suffix, out_name):
    counts_csv = os.path.join(META_RESULTS_DIR, f'cell_type_counts{suffix}.csv')
    if not os.path.exists(counts_csv):
        print(f"Cell counts file not found: {counts_csv}")
        return
        
    df = pd.read_csv(counts_csv)
    df.rename(columns={'Dataset': 'Cancer'}, inplace=True)
    
    # Format numbers with commas
    for col in df.columns:
        if col != 'Cancer':
            df[col] = df[col].apply(lambda x: f"{int(x):,}" if pd.notnull(x) else x)
            
    out_path = os.path.join(SUMMARY_TABLES_DIR, f'{out_name}.csv')
    df.to_csv(out_path, index=False)
    print(f"Generated data-driven overview: {out_path}")

def generate_subclone_summary(suffix, out_name):
    from pan_cancer_config import CANCER_CAP
    cancers = list(CANCER_CAP.keys())
    
    data = []
    for c in cancers:
        file_path = os.path.join(META_RESULTS_DIR, f"{c}_primary_signature_scores{suffix}.csv")
        if not os.path.exists(file_path):
            continue
        df = pd.read_csv(file_path)
        count = len(df)
        if count == 0:
            continue
            
        score_cols = [col for col in df.columns if col.startswith('Metastatic_Signature_Score')]
        if not score_cols:
            continue
            
        for score_col in score_cols:
            scores = df[score_col]
            skew = scores.skew()
            
            if skew > 0.5:
                dist = "Right-skewed"
            elif skew < -0.5:
                dist = "Left-skewed"
            else:
                dist = "Symmetric"
                
            # Mathematical extraction of highly metastatic subclone
            # Defined as cells with scores > Mean + 1 Standard Deviation
            mean_score = scores.mean()
            std_score = scores.std()
            subclone_pct = (scores > (mean_score + std_score)).mean() * 100
            
            sig_name = score_col.replace('Metastatic_Signature_Score_', '')
            if sig_name == 'Metastatic_Signature_Score':
                sig_name = 'Conserved Pan-Cancer'
                
            data.append({
                'Cancer': f"**{c.capitalize()}**",
                'Signature': sig_name,
                'Primary Cells Scored': f"{count:,}",
                'Score Distribution': dist,
                'Pre-Metastatic Subclone (%)': f"{subclone_pct:.1f}% (> +1 SD)"
            })
    if data:
        summary_df = pd.DataFrame(data)
        out_path = os.path.join(SUMMARY_TABLES_DIR, f'{out_name}.csv')
        summary_df.to_csv(out_path, index=False)
        print(f"Generated: {out_path}")

def generate_unique_signatures():
    # Dynamically read the upset plot intersection data to get UNIQUE genes per cancer
    upset_csv = os.path.join(META_RESULTS_DIR, f'upset_plot_data{ANALYSIS_SUFFIX}.csv')
    if not os.path.exists(upset_csv):
        print(f"Upset plot data not found: {upset_csv}")
        return
        
    upset_df = pd.read_csv(upset_csv)
    
    # Load primary metabolite database to get REAL metabolite associations
    db_path = os.path.join(OUTPUT_DIR, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')
    if not os.path.exists(db_path):
        print(f"Database not found: {db_path}")
        return
    db_df = pd.read_csv(db_path, low_memory=False)
    db_df['Target'] = db_df['Target'].astype(str).str.split(r'[,;]')
    db_df = db_df.explode('Target')
    db_df['Target'] = db_df['Target'].str.strip()
    
    # Find unique genes (appearing in exactly 1 cancer)
    gene_counts = upset_df.groupby('Up_Regulated_Gene')['Cancer_Type'].count()
    unique_genes = gene_counts[gene_counts == 1].index
    
    # Filter upset_df to only include the unique genes
    unique_df = upset_df[upset_df['Up_Regulated_Gene'].isin(unique_genes)]
    
    cancers = upset_df['Cancer_Type'].unique()
    
    results = []
    for cancer in cancers:
        cancer_unique_genes = unique_df[unique_df['Cancer_Type'] == cancer]['Up_Regulated_Gene'].tolist()
        
        # Take first 5 for the table
        top_genes = cancer_unique_genes[:5]
        
        # Find metabolites for these top genes from the primary database
        top_metabs = set()
        for gene in top_genes:
            gene_data = db_df[db_df['Target'] == gene]
            if not gene_data.empty:
                for m in gene_data['Metabolite_Name'].dropna().unique():
                    top_metabs.add(m)
                    
        results.append({
            "Cancer": cancer,
            "Unique Metastatic Targets": len(cancer_unique_genes),
            "Top 5 Unique Genes": ", ".join(top_genes),
            "Associated Metabolites (Database)": ", ".join(list(top_metabs)[:10]) # Show up to 10
        })
        
    df = pd.DataFrame(results)
    out_path = os.path.join(SUMMARY_TABLES_DIR, f'cancer_specific_unique_signatures{ANALYSIS_SUFFIX}.csv')
    df.to_csv(out_path, index=False)
    print(f"Generated data-driven unique signatures: {out_path}")

def main():
    print("Generating AI Summary Tables from raw data and databases...")
    generate_metastatic_enrichment_table()
    generate_immune_evasion_table()
    generate_ccc_potential_table()
    generate_annotated_signature()
    generate_unique_signatures()
    
    # Generate dynamically based on the current run's configuration
    generate_dataset_overview(ANALYSIS_SUFFIX, f'dataset_overview{ANALYSIS_SUFFIX}')
    generate_subclone_summary(ANALYSIS_SUFFIX, f'subclone_summary{ANALYSIS_SUFFIX}')
    
    print(f"\nAll summary CSVs have been saved to: {SUMMARY_TABLES_DIR}")

if __name__ == '__main__':
    main()
