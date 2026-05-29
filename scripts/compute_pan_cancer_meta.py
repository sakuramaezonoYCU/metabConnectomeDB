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

CANCERS = ['breast', 'colorectal', 'lung', 'melanoma', 'ovarian']

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
    plt.title('Overlap of Up-Regulated Metastatic Metabolic Genes Across 5 Cancers')
    
    plot_path = os.path.join(META_RESULTS_DIR, f'upset_plot{ANALYSIS_SUFFIX}.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved UpSet plot to {plot_path}")
    
    # Save the strictly conserved genes to a separate CSV
    pan_cancer_genes_list = list(pan_cancer_genes)
    df_pan = pd.DataFrame({'Strictly_Conserved_Gene': pan_cancer_genes_list})
    pan_csv = os.path.join(META_RESULTS_DIR, f'pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv')
    df_pan.to_csv(pan_csv, index=False)
    
    return pan_cancer_genes_list

def generate_network_plot(pan_cancer_genes):
    if not pan_cancer_genes:
        print("No pan-cancer genes to build network for.")
        return
        
    pairs_file = os.path.join(OUTPUT_DIR, 'human_metab_target_pairs_cancer_annotated.csv')
    if not os.path.exists(pairs_file):
        pairs_file = os.path.join(OUTPUT_DIR, 'human_database_merge_unique_metab_target_pairs.csv')
        
    if not os.path.exists(pairs_file):
        print(f"Could not find pairs file to build network.")
        return

    df_pairs = pd.read_csv(pairs_file)
    
    target_col = 'Target' if 'Target' in df_pairs.columns else 'Target_Gene'
    metab_col = 'Metabolite' if 'Metabolite' in df_pairs.columns else 'Metabolite_Name'
    
    if target_col not in df_pairs.columns or metab_col not in df_pairs.columns:
        target_col = df_pairs.columns[1]
        metab_col = df_pairs.columns[0]

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
