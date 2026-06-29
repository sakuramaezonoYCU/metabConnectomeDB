import pandas as pd
import os
import matplotlib.pyplot as plt
from upsetplot import from_contents, plot
import networkx as nx

from pan_cancer_config import CANCER_CAP, ANALYSIS_SUFFIX, get_de_csv_path

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
META_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'pan_cancer_meta_results')
os.makedirs(META_RESULTS_DIR, exist_ok=True)

from pan_cancer_config import CANCERS_TO_RUN as CANCERS

def generate_upset_plot():
    contents = {}
    pan_cancer_genes = None
    upset_data_rows = []
    
    for cancer in CANCERS:
        res_file = get_de_csv_path(cancer)
        if os.path.exists(res_file):
            df = pd.read_csv(res_file)
            up_genes = df[df['Significance'] == 'Up in Metastasis']['names'].tolist()
            contents[cancer.capitalize()] = up_genes
            
            for gene in up_genes:
                upset_data_rows.append({'Cancer_Type': cancer.capitalize(), 'Up_Regulated_Gene': gene})
            
            if pan_cancer_genes is None:
                pan_cancer_genes = set(up_genes)
            else:
                pan_cancer_genes = pan_cancer_genes.intersection(up_genes)

    if not contents:
        print("No DE data found for UpSet plot.")
        return []

    # Export CSV
    df_upset = pd.DataFrame(upset_data_rows)
    csv_path = os.path.join(META_RESULTS_DIR, f'upset_plot_data{ANALYSIS_SUFFIX}.csv')
    df_upset.to_csv(csv_path, index=False)
    print(f"Saved UpSet plot data to {csv_path}")

    # Generate UpSet plot
    upset_data = from_contents(contents)
    fig = plt.figure(figsize=(10, 6))
    plot(upset_data, fig=fig, sort_by='degree', sort_categories_by=None)
    plt.title(f'Overlap of Up-Regulated Metastatic Metabolic Genes Across {len(CANCERS)} Cancers')
    
    plot_path = os.path.join(META_RESULTS_DIR, f'upset_plot{ANALYSIS_SUFFIX}.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved UpSet plot to {plot_path}")
    
    from itertools import combinations

    # Dynamically find conserved signatures, relaxing the cancer overlap threshold if <= 1 gene is found
    target_overlap = len(CANCERS)
    final_genes = []
    
    while target_overlap > 0:
        combo_genes = set()
        
        # Save each combination and aggregate their union
        for combo in combinations([c.capitalize() for c in CANCERS], target_overlap):
            if all(c in contents for c in combo):
                c_genes = set(contents[combo[0]])
                for c in combo[1:]:
                    c_genes = c_genes.intersection(set(contents[c]))
                
                if c_genes:
                    combo_name = "_".join(combo)
                    df_combo = pd.DataFrame({'Strictly_Conserved_Gene': list(c_genes)})
                    combo_csv = os.path.join(META_RESULTS_DIR, f'pan_cancer_signature_{combo_name}{ANALYSIS_SUFFIX}.csv')
                    df_combo.to_csv(combo_csv, index=False)
                    print(f"Saved {target_overlap}-cancer signature {combo_name} with {len(c_genes)} genes.")
                    combo_genes.update(c_genes)
                    
        final_genes = list(combo_genes)
        
        if len(final_genes) > 1:
            if target_overlap == len(CANCERS):
                print(f"Saved strict {target_overlap}-cancer signature ({len(final_genes)} genes) as the primary conserved signature.")
            else:
                print(f"Saved the union of {target_overlap}-cancer combinations ({len(final_genes)} genes) as the primary conserved signature (Relaxation Rule applied).")
            break
        else:
            if target_overlap == len(CANCERS):
                print(f"Note: There are only {len(final_genes)} genes commonly upregulated across all {len(CANCERS)} cancers, which is insufficient for meaningful network analysis.")
            else:
                print(f"Note: There are only {len(final_genes)} genes conserved across {target_overlap} cancers.")
            print(f"Relaxing threshold to {target_overlap - 1} cancers...")
            target_overlap -= 1
        
    df_pan = pd.DataFrame({'Strictly_Conserved_Gene': final_genes})
    pan_csv = os.path.join(META_RESULTS_DIR, f'pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv')
    df_pan.to_csv(pan_csv, index=False)
    
    # Save an explicitly named alias expected by the downstream markdown reports
    conserved_csv = os.path.join(META_RESULTS_DIR, 'conserved_target_genes.csv')
    df_pan.to_csv(conserved_csv, index=False)
    
    return final_genes

def generate_network_plot(pan_cancer_genes):
    if not pan_cancer_genes:
        raise ValueError("No pan-cancer genes identified to build network for. UpSet plot generation likely failed.")
        
    pairs_file = os.path.join(OUTPUT_DIR, 'human_database_merge_unique_metab_target_pairs.csv')
        
    if not os.path.exists(pairs_file):
        raise FileNotFoundError(f"Could not find pairs file to build network. Checked: {pairs_file}")

    df_pairs = pd.read_csv(pairs_file, low_memory=False)
    
    target_col = 'Target'
    metab_col = 'Metabolite_Name'
    
    # Ensure targets are exploded if they are piped (e.g. 'A | B | C')
    if df_pairs[target_col].astype(str).str.contains(r' \| ').any():
        df_pairs = df_pairs.assign(**{target_col: df_pairs[target_col].astype(str).str.split(r' \| ')}).explode(target_col)

    df_net = df_pairs[df_pairs[target_col].isin(pan_cancer_genes)].copy()
    
    # Export Network CSV
    csv_path = os.path.join(META_RESULTS_DIR, f'metabolite_target_network_edges{ANALYSIS_SUFFIX}.csv')
    df_net[[metab_col, target_col]].drop_duplicates().to_csv(csv_path, index=False)
    print(f"Saved Network edge data to {csv_path}")
    
    G = nx.Graph()
    metabolites = set()
    targets = set()
    
    for _, row in df_net.iterrows():
        m = row[metab_col]
        t = row[target_col]
        G.add_node(m, bipartite=0)
        G.add_node(t, bipartite=1)
        G.add_edge(m, t)
        metabolites.add(m)
        targets.add(t)
        
    plt.figure(figsize=(14, 10))
    pos = nx.spring_layout(G, k=0.5, iterations=50)
    
    nx.draw_networkx_nodes(G, pos, nodelist=list(metabolites), node_color='skyblue', node_size=300, alpha=0.8, label='Metabolite')
    nx.draw_networkx_nodes(G, pos, nodelist=list(targets), node_color='lightcoral', node_size=500, alpha=0.9, label='Gene Target')
    nx.draw_networkx_edges(G, pos, alpha=0.4, edge_color='gray')
    nx.draw_networkx_labels(G, pos, font_size=8, font_family="sans-serif")
    
    plt.title(f'Pan-Cancer {len(pan_cancer_genes)}-Gene Conserved Target-Metabolite Network', size=16)
    plt.axis('off')
    plt.legend(scatterpoints=1, loc='upper left')
    
    plot_path = os.path.join(META_RESULTS_DIR, f'metabolite_target_network{ANALYSIS_SUFFIX}.png')
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved Network plot to {plot_path}")

if __name__ == '__main__':
    pan_genes = generate_upset_plot()
    generate_network_plot(pan_genes)
