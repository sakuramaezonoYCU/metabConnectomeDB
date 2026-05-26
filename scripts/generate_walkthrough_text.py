import base64
import os

workspace_dir = '.'
output_dir = os.path.join(workspace_dir, 'output')

def get_base64(filepath):
    with open(filepath, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')

print("Converting figures to base64...")
b64_fig1 = get_base64(os.path.join(output_dir, 'robust_orphan_interactions_bar.png'))
b64_fig2 = get_base64(os.path.join(output_dir, 'immune_cell_specificity_dotplot.png'))
b64_fig3 = get_base64(os.path.join(output_dir, 'orphan_metabolic_immune_connectome_network.png'))

print("Generating walkthrough text...")
walkthrough_content = f"""# Walkthrough - Advanced Visualizations and Consistency Analysis for Metabolic Immune Evasion

I have successfully updated **Section 3: Map Orphan Metabolites to Immune Targets** of your Jupyter notebook [orphan_metabolic_immune_evasion.ipynb](file:///Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/orphan_metabolic_immune_evasion.ipynb). The original, cluttered scatter plot has been replaced with a robust cross-cancer consistency grouping analysis and three publications-grade, highly-informative figures.

The continuous variables (like Enrichment $Log_2$ Fold Change) are now modeled with continuous color ranges (colorbars), while discrete variables (like Cancer consistency count) are modeled using separate categorical Legends.

---

## 1. Top Robust Checkpoints Discovered

Our grouping analysis evaluated the 934 candidates across B-cells, macrophages, and dendritic cells, ranking them by **Consistency** (the number of different cancers they are active in) and **Average Enrichment strength**.

The top 10 most consistent metabolic-immune checkpoints are:
1. **nitric oxide → HLA-DRA** (Significant across 6 cancers, mean $Log_2$ FC = 5.25)
2. **guanosine triphosphate → CD74** (Significant across 6 cancers, mean $Log_2$ FC = 5.07)
3. **n-acetylglucosamine → CD74** (Significant across 6 cancers, mean $Log_2$ FC = 5.07)
4. **acetaldehyde → HLA-DRB1** (Significant across 6 cancers, mean $Log_2$ FC = 4.14)
5. **adenosine triphosphate → HLA-DRB1** (Significant across 6 cancers, mean $Log_2$ FC = 4.14)
6. **d-galactose → HLA-DRB1** (Significant across 6 cancers, mean $Log_2$ FC = 4.14)
7. **d-glucose → HLA-DRB1** (Significant across 6 cancers, mean $Log_2$ FC = 4.14)
8. **d-mannose → HLA-DRB1** (Significant across 6 cancers, mean $Log_2$ FC = 4.14)
9. **formate → HLA-DRB1** (Significant across 6 cancers, mean $Log_2$ FC = 4.14)
10. **fumarate → HLA-DRB1** (Significant across 6 cancers, mean $Log_2$ FC = 4.14)

---

## 2. Publications-Quality Visualizations

### Figure 1: Cross-Cancer Consistency Bar Plot (The "Robustness" Chart)
This figure ranks the top 15 checkpoints by their occurrence across different tumor types (cancers). The continuous color mapping (colorbar) shows the average enrichment strength across these cancers.

![Robust Orphan Checkpoints Bar Plot](data:image/png;base64,{b64_fig1})

---

### Figure 2: Cell-Type Specificity Dot Plot (The "Checkpoint Profile" Bubble Chart)
This bubble plot shows where these robust checkpoints are active across immune cell populations (B-cells, macrophages, dendritic cells).
* **Dot Size**: Discrete Legend representing the number of cancer datasets in which the interaction is active.
* **Dot Color**: Continuous colorbar representing the average enrichment strength ($Log_2$ Fold Change).

![Immune Cell Specificity Dot Plot](data:image/png;base64,{b64_fig2})

---

### Figure 3: Bipartite Connectome Network Graph
This systems-biology network connects the top 10 metabolic ligands (left) to their target receptors and immune populations (right).
* **Edge Thickness**: Direct mapping to the average enrichment fold change.
* **Edge Color**: Continuous colorbar showing the cross-cancer consistency count.

![Metabolic-Immune Connectome Network Graph](data:image/png;base64,{b64_fig3})

---

## 3. Styled Interactive HTML Report

I have compiled the updated notebook cells and created a premium styled HTML candidate report in [immune_evasion_orphan_metabolic_candidates.html](file:///Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/output/immune_evasion_orphan_metabolic_candidates.html) (currently active in your editor!).

This report features:
- **Interactive Tabs** to instantly switch between the three high-resolution figures.
- **Glassmorphism Design & Premium Typography** using Inter and Helvetica.
- **Interactive Table Badge Styling** mapping Tier 2 and Tier 3 candidate pairs automatically.
"""

output_walkthrough_path = os.path.join(output_dir, 'walkthrough_content.md')
print(f"Writing walkthrough content to {output_walkthrough_path}...")
with open(output_walkthrough_path, 'w', encoding='utf-8') as f:
    f.write(walkthrough_content)

print("🎉 Walkthrough content successfully written!")
