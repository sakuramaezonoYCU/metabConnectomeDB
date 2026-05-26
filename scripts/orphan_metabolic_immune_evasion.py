#!/usr/bin/env python
# coding: utf-8

# # Metabolic "Orphan" Predicting Immune Evasion
# 
# **Research Question 1:** Can metabolic "orphan" interactions predict immune evasion across tumor types?
# 
# This notebook systematically cross-references our computationally predicted, literature-sparse ("Tier 2/3") metabolic target pairs against immune cell populations (such as B-cells, macrophages, and dendritic cells) in our cancer datasets. 
# 
# By identifying uncharacterized metabolic ligands that are heavily upregulated in these immune populations, we can potentially discover novel, unpatented immune-checkpoints driven by the metabolic microenvironment.

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import sys

# Define workspace
workspace_dir = '../'
output_dir = os.path.join(workspace_dir, 'output')
sys.path.append(output_dir)

# We use the parse_md_tables script to extract previously computed DE tables from HTML reports
from parse_md_tables import TableExtractor

print("Loading comprehensive target pairs...")
pairs_df = pd.read_csv(os.path.join(workspace_dir, 'output', 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv'))

# Explode Target column just in case it has commas
pairs_df['Target'] = pairs_df['Target'].astype(str).str.split(r'[,;]')
pairs_df = pairs_df.explode('Target')
pairs_df['Target'] = pairs_df['Target'].str.strip()

print(f"Total unique pairs loaded: {len(pairs_df)}")


# ## 1. Isolate the "Orphan" (Tier 2 & 3 Pairs)
# 
# We define "Orphan" pairs as those that are supported by 2-3 databases but lack heavy, redundant validation (which would make them Tier 1).

# In[2]:


def calc_tier(db_str):
    if pd.isna(db_str): return 'Tier 3 (Low)'
    count = len(str(db_str).split(','))
    if count >= 4: return 'Tier 1 (High)'
    elif count >= 2: return 'Tier 2 (Medium)'
    else: return 'Tier 3 (Low)'

pairs_df['Pair_Confidence_Tier'] = pairs_df['database'].apply(calc_tier)

orphan_metabolic_pairs = pairs_df[pairs_df['Pair_Confidence_Tier'].isin(['Tier 2 (Medium)', 'Tier 3 (Low)'])]
print(f"Orphan pairs identified: {len(orphan_metabolic_pairs)}")


# ## 2. Extract Highly Enriched Immune Targets
# 
# We parse the existing `cancer_*.html` reports to extract the genome-wide Differential Expression (DE) tables. We specifically filter for targets highly upregulated (`Log2FC > 1.0`, `padj < 0.05`) in immune clusters like B cells, Macrophages, and Dendritic cells.

# In[3]:


import glob

de_files = glob.glob(os.path.join(output_dir, '*_DE_genome_wide.csv'))
immune_enriched_genes = []

print("Parsing DE tables from CSVs...")
for f in de_files:
    # Get cancer name from the filename: "breast_50150_DE_genome_wide.csv" -> "breast-cancer"
    basename = os.path.basename(f)
    cancer_name = basename.replace('cancer_', '').split('_')[0]

    try:
        de_table = pd.read_csv(f)

        de_table['logfoldchanges'] = pd.to_numeric(de_table['logfoldchanges'], errors='coerce')
        de_table['pvals_adj'] = pd.to_numeric(de_table['pvals_adj'], errors='coerce')

        # Filter for immune cells
        if 'group' in de_table.columns:
            immune_mask = de_table['group'].str.contains('B cell|Macrophage|mononuclear phagocyte|dendritic cell|T cell', case=False, na=False)
            sig_mask = (de_table['logfoldchanges'] > 1.0) & (de_table['pvals_adj'] < 0.05)
            filtered = de_table[immune_mask & sig_mask].copy()
            filtered['Cancer'] = cancer_name

            immune_enriched_genes.append(filtered)

    except Exception as e:
        print(f"Error parsing {f}: {e}")

if len(immune_enriched_genes) > 0:
    immune_de_df = pd.concat(immune_enriched_genes, ignore_index=True)
    print(f"Found {len(immune_de_df['names'].unique())} unique highly enriched immune target genes across all cancers.")
else:
    immune_de_df = pd.DataFrame()
    print("No immune enriched targets found (or no CSVs found).")


# ## 3. Map Orphan Metabolites to Immune Targets
# 
# We cross-reference our immune DE genes with the Orphan pairs, and visualize the top hits.

# In[4]:


# Inner merge to find the intersection
orphan_metabolic_immune_df = pd.merge(
    immune_de_df, 
    orphan_metabolic_pairs, 
    left_on='names', 
    right_on='Target', 
    how='inner'
)

print(f"Identified {len(orphan_metabolic_immune_df)} critical immune evasion orphan interactions.")

import networkx as nx

# Grouping by both Interaction and Immune Cell Type to calculate consistency across tumor types
df = orphan_metabolic_immune_df.copy()
df['Interaction'] = df['Metabolite_Name'] + ' → ' + df['Target']

grouped = df.groupby(['Interaction', 'Metabolite_Name', 'Target', 'group']).agg(
    cancers_count=('Cancer', 'nunique'),
    mean_logfc=('logfoldchanges', 'mean'),
    max_logfc=('logfoldchanges', 'max'),
    cancers_list=('Cancer', lambda x: ", ".join(sorted(x.unique())))
).reset_index()

grouped = grouped.sort_values(by=['cancers_count', 'mean_logfc'], ascending=False)

overall_grouped = df.groupby(['Interaction', 'Metabolite_Name', 'Target']).agg(
    cancers_count=('Cancer', 'nunique'),
    mean_logfc=('logfoldchanges', 'mean'),
    max_logfc=('logfoldchanges', 'max')
).reset_index().sort_values(by=['cancers_count', 'mean_logfc'], ascending=False)

overall_grouped.to_csv(os.path.join(output_dir, 'robust_orphan_interaction_rankings.csv'), index=False)

# ----------------------------------------------------
# Figure 1: Cross-Cancer Consistency Bar Plot
# ----------------------------------------------------
plt.figure(figsize=(7.5, 4))
top_overall = overall_grouped.head(15)

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

cbar = plt.colorbar(sm, ax=plt.gca(), orientation='vertical', pad=0.02, shrink=0.8)
cbar.set_label('Mean Enrichment ($Log_2$ Fold Change)', rotation=270, labelpad=12, fontsize=8.5)
cbar.ax.tick_params(labelsize=7.5)

plt.xlabel('Cross-Cancer Consistency (Number of Significant Cancers)', fontsize=9, fontweight='bold')
plt.ylabel('Orphan Metabolic Interaction Pair', fontsize=9, fontweight='bold')
plt.title('Top 15 Most Robust Orphan Metabolic Immune Evasion Checkpoints\n(Present Across the Most Tumor Types)', fontsize=10, fontweight='bold', pad=8)
plt.xlim(0, top_overall['cancers_count'].max() + 0.5)
plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
plt.gca().tick_params(axis='both', which='major', labelsize=8)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'robust_orphan_interactions_bar.png'), dpi=300, bbox_inches='tight')
plt.show()

# ----------------------------------------------------
# Figure 2: Cell-Type Specificity Dot Plot (Bubble Plot)
# ----------------------------------------------------
top_interactions = overall_grouped.head(15)['Interaction'].tolist()
filtered_grouped = grouped[grouped['Interaction'].isin(top_interactions)]

fig, ax = plt.subplots(figsize=(8, 4.8))
norm = plt.Normalize(filtered_grouped['mean_logfc'].min(), filtered_grouped['mean_logfc'].max())
mapper = plt.cm.ScalarMappable(norm=norm, cmap='magma')
mapper.set_array([])

scatter = sns.scatterplot(
    data=filtered_grouped,
    x='group',
    y='Interaction',
    size='cancers_count',
    hue='mean_logfc',
    hue_norm=norm,
    palette='magma',
    sizes=(40, 240),
    alpha=0.85,
    edgecolor='black',
    linewidth=0.5,
    ax=ax,
    legend='brief'
)

handles, labels = ax.get_legend_handles_labels()
size_handles = []
size_labels = []
is_size_section = False
for h, l in zip(handles, labels):
    if l == 'cancers_count':
        is_size_section = True
        continue
    if l == 'mean_logfc':
        is_size_section = False
        continue
    if is_size_section:
        size_handles.append(h)
        size_labels.append(l)

if not size_handles:
    size_handles = handles[6:]
    size_labels = labels[6:]

ax.legend(
    size_handles, size_labels,
    title='Cancer Consistency',
    bbox_to_anchor=(1.02, 1.0),
    loc='upper left',
    frameon=True,
    facecolor='white',
    fontsize=7.5,
    title_fontsize=8
)

cbar = fig.colorbar(mapper, ax=ax, orientation='vertical', pad=0.02, shrink=0.5, anchor=(0.0, 0.0))
cbar.set_label('Mean Enrichment ($Log_2$ Fold Change)', rotation=270, labelpad=12, fontsize=8.5)
cbar.ax.tick_params(labelsize=7.5)

plt.xlabel('Immune Cell Population', fontsize=9, fontweight='bold', labelpad=3)
plt.ylabel('Orphan Metabolic Interaction Pair', fontsize=9, fontweight='bold', labelpad=3)
plt.title('Immune Cell Specificity of Robust Orphan Metabolic Interactions', fontsize=10, fontweight='bold', pad=8)
plt.grid(True, which='both', linestyle=':', alpha=0.5)
ax.tick_params(axis='both', which='major', labelsize=8)
plt.tight_layout(rect=[0, 0, 0.82, 1])
plt.savefig(os.path.join(output_dir, 'immune_cell_specificity_dotplot.png'), dpi=300, bbox_inches='tight')
plt.show()

# ----------------------------------------------------
# Figure 3: Metabolic-Immune Connectome Network Graph
# ----------------------------------------------------
top_10_interactions = overall_grouped.head(10)
G = nx.Graph()

metabolites = top_10_interactions['Metabolite_Name'].unique().tolist()
receptors = top_10_interactions['Target'].unique().tolist()

for m in metabolites: G.add_node(m, bipartite=0, label=m)
for r in receptors: G.add_node(r, bipartite=1, label=r)

for _, row in top_10_interactions.iterrows():
    G.add_edge(row['Metabolite_Name'], row['Target'], weight=row['mean_logfc'], consistency=row['cancers_count'])

plt.figure(figsize=(8, 4.2))
pos = {}
pos.update((node, (1, index)) for index, node in enumerate(metabolites))
pos.update((node, (2, index * len(metabolites) / len(receptors))) for index, node in enumerate(receptors))

edges = G.edges()
weights = [G[u][v]['weight'] * 1.0 for u, v in edges]
consistencies = [G[u][v]['consistency'] for u, v in edges]

nx.draw_networkx_nodes(G, pos, nodelist=metabolites, node_color='#3f2d54', node_size=400, alpha=0.9)
nx.draw_networkx_nodes(G, pos, nodelist=receptors, node_color='#c0392b', node_size=400, alpha=0.9)

draw_edges = nx.draw_networkx_edges(
    G, pos, width=weights, edge_color=consistencies, edge_cmap=plt.cm.plasma,
    edge_vmin=min(consistencies), edge_vmax=max(consistencies), alpha=0.7
)

cbar = plt.colorbar(draw_edges, ax=plt.gca(), orientation='vertical', pad=0.02, shrink=0.8)
cbar.set_label('Cross-Cancer Consistency (Cancers Count)', rotation=270, labelpad=12, fontsize=8.5)
cbar.ax.tick_params(labelsize=7.5)

metab_labels = {n: n for n in metabolites}
rec_labels = {n: n for n in receptors}
pos_metab_labels = {k: (v[0] - 0.08, v[1]) for k, v in pos.items() if k in metabolites}
pos_rec_labels = {k: (v[0] + 0.08, v[1]) for k, v in pos.items() if k in receptors}

nx.draw_networkx_labels(G, pos_metab_labels, labels=metab_labels, font_size=7.5, font_weight='bold', horizontalalignment='right')
nx.draw_networkx_labels(G, pos_rec_labels, labels=rec_labels, font_size=7.5, font_weight='bold', horizontalalignment='left')

plt.xlim(0.3, 2.7)
plt.margins(y=0.1)
plt.axis('off')
plt.title('Top 10 Orphan Metabolic-Immune Connectome Network\n(Edge Thickness: Enrichment Strength, Edge Color: Consistency)', fontsize=10, fontweight='bold', pad=8)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'orphan_metabolic_immune_connectome_network.png'), dpi=300, bbox_inches='tight')
plt.show()


# In[5]:


# Save the comprehensive list for manual review
output_csv = os.path.join(output_dir, 'immune_evasion_orphan_metabolic_candidates.csv')
orphan_metabolic_immune_df.to_csv(output_csv, index=False)
print(f"Saved highly confident orphan candidates to {output_csv}")

# Display top 10
orphan_metabolic_immune_df[['Cancer', 'group', 'names', 'logfoldchanges', 'Metabolite_Name', 'Pair_Confidence_Tier']].sort_values(by='logfoldchanges', ascending=False).head(10)


# ## Automated HTML Export
# 
# Generate a styled HTML report of this notebook for sharing and viewing.

# In[12]:


import subprocess
import sys
import os

notebook_filename = 'orphan_metabolic_immune_evasion.ipynb'

# remove extension for nbconvert
output_base = output_csv.replace(".csv", "")

print(f"Executing full notebook HTML export for '{notebook_filename}'...")

# Get Jupyter binary
jupyter_bin = os.path.join(os.path.dirname(sys.executable), 'jupyter')
if not os.path.exists(jupyter_bin):
    jupyter_bin = 'jupyter'

cmd_html = [
    jupyter_bin,
    "nbconvert",
    "--to", "html",
    notebook_filename,
    "--output", output_base
]

res_html = subprocess.run(cmd_html, capture_output=True, text=True)

if res_html.returncode == 0:
    print("🎉 SUCCESS: Notebook successfully exported as a styled HTML report!")
    print(f"   -> Saved to: '{output_base}.html'")
else:
    print("❌ HTML export failed.")
    print(res_html.stderr)


# In[ ]:




