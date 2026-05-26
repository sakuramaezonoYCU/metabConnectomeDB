import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import os

# Set beautiful styling for publication-ready figures
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Inter', 'Helvetica', 'DejaVu Sans'],
    'axes.edgecolor': '#cccccc',
    'axes.linewidth': 0.8,
    'grid.color': '#eeeeee',
    'grid.linestyle': '--',
    'figure.titlesize': 20,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
})

# Load the candidate CSV
workspace_dir = '.'
output_dir = os.path.join(workspace_dir, 'output')
csv_path = os.path.join(output_dir, 'immune_evasion_orphan_metabolic_candidates.csv')

print(f"Loading candidate data from {csv_path}...")
df = pd.read_csv(csv_path)

print(f"Total rows loaded: {len(df)}")

# Create a clean Label for the interaction
df['Interaction'] = df['Metabolite_Name'] + ' → ' + df['Target']

# ----------------------------------------------------
# 1. ANALYSIS: Group by Interaction and Cell Type
# ----------------------------------------------------
# Grouping by both Interaction and the Immune Cell Type (group) to see robustness
grouped = df.groupby(['Interaction', 'Metabolite_Name', 'Target', 'group']).agg(
    cancers_count=('Cancer', 'nunique'),
    mean_logfc=('logfoldchanges', 'mean'),
    max_logfc=('logfoldchanges', 'max'),
    cancers_list=('Cancer', lambda x: ", ".join(sorted(x.unique())))
).reset_index()

# Sort by consistency and strength
grouped = grouped.sort_values(by=['cancers_count', 'mean_logfc'], ascending=False)

# Get overall interaction stats (ignoring cell type) for general ranking
overall_grouped = df.groupby(['Interaction', 'Metabolite_Name', 'Target']).agg(
    cancers_count=('Cancer', 'nunique'),
    mean_logfc=('logfoldchanges', 'mean'),
    max_logfc=('logfoldchanges', 'max')
).reset_index().sort_values(by=['cancers_count', 'mean_logfc'], ascending=False)

print("\nTop 10 overall most robust orphan metabolic immune evasion interactions:")
print(overall_grouped.head(10)[['Interaction', 'cancers_count', 'mean_logfc']])

# Save these summarized stats to outputs for user review
overall_grouped.to_csv(os.path.join(output_dir, 'robust_orphan_interaction_rankings.csv'), index=False)

# ----------------------------------------------------
# VISUALIZATION 1: Cross-Cancer Consistency Bar Plot
# ----------------------------------------------------
plt.figure(figsize=(12, 8))
top_overall = overall_grouped.head(20)

# Create a horizontal bar plot colored by average logFC
norm = plt.Normalize(top_overall['mean_logfc'].min(), top_overall['mean_logfc'].max())
sm = plt.cm.ScalarMappable(cmap="magma", norm=norm)
sm.set_array([])

bars = plt.barh(
    top_overall['Interaction'][::-1], 
    top_overall['cancers_count'][::-1],
    color=[plt.cm.magma(norm(val)) for val in top_overall['mean_logfc'][::-1]],
    edgecolor='none',
    height=0.7
)

cbar = plt.colorbar(sm, ax=plt.gca(), orientation='vertical')
cbar.set_label('Mean Enrichment ($Log_2$ Fold Change)', rotation=270, labelpad=20)

plt.xlabel('Cross-Cancer Consistency (Number of Significant Cancers)', fontsize=14, fontweight='bold')
plt.ylabel('Orphan Metabolic Interaction Pair', fontsize=14, fontweight='bold')
plt.title('Top 20 Most Robust Orphan Metabolic Immune Evasion Checkpoints\n(Present Across the Most Tumor Types)', fontsize=16, fontweight='bold', pad=15)
plt.xlim(0, top_overall['cancers_count'].max() + 0.5)
plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
plt.tight_layout()

vis1_path = os.path.join(output_dir, 'robust_orphan_interactions_bar.png')
plt.savefig(vis1_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved robust bar plot to {vis1_path}")

# ----------------------------------------------------
# VISUALIZATION 2: Cell-Type Specificity Dot Plot (Bubble Plot)
# ----------------------------------------------------
# Filter grouped dataset to only include the top 20 most robust interactions
top_interactions = overall_grouped.head(20)['Interaction'].tolist()
filtered_grouped = grouped[grouped['Interaction'].isin(top_interactions)]

plt.figure(figsize=(13, 10))

# Scatter (Dot Plot)
scatter = sns.scatterplot(
    data=filtered_grouped,
    x='group',
    y='Interaction',
    size='cancers_count',
    hue='mean_logfc',
    palette='magma',
    sizes=(80, 500),
    alpha=0.85,
    edgecolor='black',
    linewidth=0.5
)

# Customizing legends
handles, labels = scatter.get_legend_handles_labels()

# Separate hue and size legend items
legend1 = plt.legend(
    handles[1:5], labels[1:5], 
    title='Mean $Log_2$ FC', 
    bbox_to_anchor=(1.02, 1), 
    loc='upper left',
    frameon=True,
    facecolor='white'
)
plt.gca().add_artist(legend1)

# Finding where sizes start in legend handles (it depends on how seaborn structures it)
# We can just let seaborn handle it or make a cleaner legend
plt.legend(
    handles[6:], labels[6:], 
    title='Cancer Count', 
    bbox_to_anchor=(1.02, 0.6), 
    loc='upper left',
    frameon=True,
    facecolor='white'
)

plt.xlabel('Immune Cell Population', fontsize=14, fontweight='bold', labelpad=10)
plt.ylabel('Orphan Metabolic Interaction Pair', fontsize=14, fontweight='bold', labelpad=10)
plt.title('Immune Cell Specificity of Robust Orphan Metabolic Interactions', fontsize=16, fontweight='bold', pad=15)
plt.grid(True, which='both', linestyle=':', alpha=0.6)
plt.tight_layout()

vis2_path = os.path.join(output_dir, 'immune_cell_specificity_dotplot.png')
plt.savefig(vis2_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved dot plot to {vis2_path}")

# ----------------------------------------------------
# VISUALIZATION 3: Bipartite Connectome Network
# ----------------------------------------------------
# Focus on top 10 interactions to keep the network highly readable and elegant
top_10_interactions = overall_grouped.head(10)

G = nx.Graph()

# Add nodes with bipartite attribute
metabolites = top_10_interactions['Metabolite_Name'].unique().tolist()
receptors = top_10_interactions['Target'].unique().tolist()

for m in metabolites:
    G.add_node(m, bipartite=0, label=m)
for r in receptors:
    G.add_node(r, bipartite=1, label=r)

# Add edges with weights and colors
for _, row in top_10_interactions.iterrows():
    G.add_edge(
        row['Metabolite_Name'], 
        row['Target'], 
        weight=row['mean_logfc'],
        consistency=row['cancers_count']
    )

plt.figure(figsize=(14, 8))

# Positional layout for bipartite graph
pos = {}
pos.update((node, (1, index)) for index, node in enumerate(metabolites))
pos.update((node, (2, index * len(metabolites) / len(receptors))) for index, node in enumerate(receptors))

# Get edge weights and consistency for drawing
edges = G.edges()
weights = [G[u][v]['weight'] * 1.5 for u, v in edges]
consistencies = [G[u][v]['consistency'] for u, v in edges]

# Draw nodes
nx.draw_networkx_nodes(
    G, pos, 
    nodelist=metabolites, 
    node_color='#3f2d54', 
    node_size=1000, 
    alpha=0.9
)
nx.draw_networkx_nodes(
    G, pos, 
    nodelist=receptors, 
    node_color='#c0392b', 
    node_size=1000, 
    alpha=0.9
)

# Draw edges with color mapping based on consistency
edge_cmap = plt.cm.plasma
edge_colors = consistencies

draw_edges = nx.draw_networkx_edges(
    G, pos, 
    width=weights, 
    edge_color=edge_colors, 
    edge_cmap=edge_cmap,
    edge_vmin=min(consistencies),
    edge_vmax=max(consistencies),
    alpha=0.7
)

# Colorbar for edge consistency
cbar = plt.colorbar(draw_edges, ax=plt.gca(), orientation='vertical', pad=0.08)
cbar.set_label('Cross-Cancer Consistency (Cancers Count)', rotation=270, labelpad=20)

# Draw node labels with offset for readability
metab_labels = {n: n for n in metabolites}
rec_labels = {n: n for n in receptors}

# Shift labels slightly to avoid overlap with nodes
pos_metab_labels = {k: (v[0] - 0.05, v[1]) for k, v in pos.items() if k in metabolites}
pos_rec_labels = {k: (v[0] + 0.05, v[1]) for k, v in pos.items() if k in receptors}

nx.draw_networkx_labels(G, pos_metab_labels, labels=metab_labels, font_size=11, font_weight='bold', horizontalalignment='right')
nx.draw_networkx_labels(G, pos_rec_labels, labels=rec_labels, font_size=11, font_weight='bold', horizontalalignment='left')

plt.xlim(0.7, 2.3)
plt.axis('off')
plt.title('Top 10 Orphan Metabolic-Immune Connectome Network\n(Edge Thickness: Enrichment Strength, Edge Color: Consistency across cancers)', fontsize=16, fontweight='bold', pad=15)
plt.tight_layout()

vis3_path = os.path.join(output_dir, 'orphan_metabolic_immune_connectome_network.png')
plt.savefig(vis3_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved connectome network to {vis3_path}")
print("🎉 Success: All three publications-quality figures have been successfully generated!")
