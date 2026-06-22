import pandas as pd

def generate_v1_metrics():
    print("Generating Version 1 Metrics...")
    
    # 1. Total Pairs and Directionality
    df = pd.read_csv('output/human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv', low_memory=False)
    total_pairs = len(df)
    
    has_rhea = df['Rhea_enzyme product/substrate'].notna()
    has_kegg = df['KEGG_enzyme product/substrate'].notna()
    has_direction = has_rhea | has_kegg
    directionality_gap = 100 - (has_direction.sum() / total_pairs * 100)
    
    # 2. B-cell targets from Pan-Cancer 100k
    de_file = 'output/ovarian_results/ovary_abdomen_omentum_uterus_100k_whole_transcriptome_2025-11-08_DE_genome_wide.csv'
    df_de = pd.read_csv(de_file)
    b_cell = df_de[df_de['group'] == 'B cell'].sort_values('scores', ascending=False)
    top_15_b_cell = b_cell.head(15)['names'].tolist()
    
    # 3. Biophysical Classes
    mw = df['MONO_MASS'].astype(float)
    under_300 = mw[mw < 300].count()
    mid_300_750 = mw[(mw >= 300) & (mw <= 750)].count()
    over_750 = mw[mw > 750].count()
    
    with open('output/version_1_metrics.txt', 'w') as f:
        f.write("=== Version 1 Verified Metrics ===\n\n")
        f.write(f"1. Total Unique Interaction Pairs: {total_pairs}\n")
        f.write(f"2. Literature Evidence (from metab_targetPair_analysis_full_report.html): 524 pairs (6.10%)\n")
        f.write(f"3. Directionality Gap (Missing Rhea/KEGG equations): {directionality_gap:.2f}%\n")
        f.write(f"4. Top 15 Highly Enriched Targets in B Cells (Ovarian 100k subset):\n")
        f.write(f"   {', '.join(top_15_b_cell)}\n")
        f.write(f"5. Biophysical Classifications based on Molecular Weight:\n")
        f.write(f"   <300 Da (Paracrine/Soluble): {under_300} pairs\n")
        f.write(f"   300-750 Da (GPCR/Hormonal): {mid_300_750} pairs\n")
        f.write(f"   >750 Da (Juxtacrine/Vesicular): {over_750} pairs\n")
        
    print("Metrics saved to output/version_1_metrics.txt")

if __name__ == '__main__':
    generate_v1_metrics()
