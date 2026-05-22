#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
🧬 Surgical Jupyter Notebook Annotation Analysis Cell Updater
============================================================
This script programmatically modifies `scripts/metab_targetPair_analysis.ipynb`
to update the markdown and code cells for Sections 5.1, 5.2, and 5.3.
It replaces them with enriched details about GtoPdb, HGNC, UniProt, and OtherDB coverage,
the Regex Rules Classifier guide, and advanced multi-panel plots.

Author: Antigravity (Advanced Agentic Coding Pair)
Date: 2026-05-21
"""

import json
import os

def edit_notebook():
    notebook_path = 'scripts/metab_targetPair_analysis.ipynb'
    if not os.path.exists(notebook_path):
        print(f"❌ Error: Notebook not found at {notebook_path}")
        return
        
    print(f"📖 Reading notebook: {notebook_path}")
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = json.load(f)
        
    modified_count = 0
    
    # Define the new cell contents
    
    cell_15_markdown = """---

## 5. Interaction Biology

### 5.1. Sensor Type Distribution

#### Purpose
To categorize the functional roles of target proteins (Receptors, Channels, Transporters, Enzymes) mediating metabolic communication. Metabolites act through distinct biological mechanisms: triggering intracellular signaling (Receptors), gating selective transport (Channels), internalizing cargo (Transporters), or undergoing chemical processing (Enzymes).

#### Database Origins & Multi-Database Mapping
To achieve high-fidelity and comprehensive functional profiling, we split and enriched target classifications across multiple database layers:
1. **`OtherDB_Sensor_Type`**: The original classifications integrated from the merged databases of **MEBOCOST** (metabolite-sensor maps), **CellPhoneDB v5** (ligand-receptor metadata), and **Celllinker2** (transporter-mediated pairings).
2. **`GtoPdb_Sensor_Type`**: Classifications retrieved from the IUPHAR/BPS **Guide to Pharmacology** mapping, mapping gene symbols to established pharmacological target classes.
3. **`HGNC_Sensor_Type`**: Gene group tags and standardized functional definitions retrieved from the **HGNC (HUGO Gene Nomenclature Committee)** Approved database.
4. **`UniProt_Sensor_Type`**: Classifications extracted dynamically from **UniProt KB** keywords and Gene Ontology (GO) terms using a high-throughput API batch query pipeline (enriched with 2,000+ cached gene symbols).
5. **`Sensor_Type` (Unified)**: The consolidated, deduplicated functional annotation that merges all the above sources and resolves multi-valued annotations row-wise (e.g., when a target acts as both a receptor and enzyme).

#### 🔮 Local Rule Classifier (classify_by_rules) Guide
To ensure **100% annotation completeness** for novel, rare, or unmapped targets, a regular-expression-based heuristic classifier serves as the final fallback layer in the pipeline. These rules recognize gene family prefixes and standard suffixes to assign correct functional roles:
- **Channels**: Matches prefixes such as `CACN*`, `SCN*`, `KCN*`, `CLCN*`, `TRP*`, `AQP*`, `ANO*`, `BEST*`, `PANX*`, `RYR*`, `ITPR*`, `ASIC*`, `GRIN*`, `GRIA*`, `GABR*`, `CHRN*`, `GLYR*`, `P2RX*`, `CX*`, `GJ*` or descriptors containing `"channel"` or `"pore"`.
- **Receptors**: Matches GPCR and nuclear receptor families including `GPR*`, `HRH*`, `ADR*`, `DRD*`, `HTR*`, `OPR*`, `P2RY*`, `LPAR*`, `S1PR*`, `FFAR*`, `CNR*`, `FZD*`, `SMO*`, `EGFR*`, `FGFR*`, `VEGFR*`, `NR*`, `THRA*`, `PPARA*`, `ESR*`, `AR*`, `GR*`, `MR*` or names containing `"receptor"`.
- **Transporters**: Matches solute carriers and active transporters such as `SLC*`, `ABC*`, `ATP1A*`, `ATP2A*`, `ATP4A*`, `ATP7A*`, `TFRC*`, `FABP*` or names containing `"transporter"`, `"carrier"`, or `"solute carrier"`.
- **Enzymes**: Matches metabolizing and signaling enzyme families including `CYP*`, `ALDH*`, `ADH*`, `COX*`, `PTGS*`, `NOS*`, `HSD*`, `SULT*`, `UGT*`, `GST*`, `MAO*`, `COMT*`, `PLA2*`, `PDE*`, `MAPK*`, `AKT*`, `JAK*`, `CDK*`, `CASP*`, `MMP*`, `DPP*` or suffixes ending in `"-ase"` (e.g. kinases, transferases, synthases, dehydrogenases)."""

    cell_16_code = """# 1. Print Coverage Comparison
print("=== Sensor Type Annotation Coverage ===")
cols_to_check = {
    "OtherDB_Sensor_Type": "Original DBs (OtherDB)",
    "GtoPdb_Sensor_Type": "Guide to Pharmacology (GtoPdb)",
    "HGNC_Sensor_Type": "HGNC approved gene groups",
    "UniProt_Sensor_Type": "UniProt KB annotations",
    "Sensor_Type": "Consolidated Unified (with regex fallbacks)"
}

coverage_data = {}
for col, label in cols_to_check.items():
    if col in df.columns:
        cnt = df[col].notna().sum()
        pct = (cnt / len(df)) * 100
        print(f"  * {label:42}: {cnt:5,} / {len(df):5,} pairs ({pct:6.2f}%)")
        coverage_data[label] = cnt
    else:
        print(f"  * {label:42}: Column '{col}' not found!")

# 2. Visualize Coverage Comparison
plt.figure(figsize=(12, 5))
cov_series = pd.Series(coverage_data).sort_values(ascending=True)
colors_cov = plt.cm.get_cmap('plasma')(np.linspace(0.2, 0.8, len(cov_series)))
bars = plt.barh(cov_series.index, cov_series.values, color=colors_cov, edgecolor='gray', height=0.6)
plt.title('Annotation Coverage Comparison by Database Source', fontweight='bold', fontsize=15, pad=15)
plt.xlabel('Number of Annotated Pairs', fontweight='bold', fontsize=12)
plt.xlim(0, len(df) * 1.15)
for bar in bars:
    width = bar.get_width()
    plt.text(width + 100, bar.get_y() + bar.get_height()/2, f"{width:,}\\n({width/len(df)*100:.1f}%)", 
             va='center', ha='left', fontsize=10, fontweight='bold')
plt.tight_layout()
plt.show()

# 3. Explode and count multi-valued consolidated Sensor_Type values
sensor_elements = []
for val in df['Sensor_Type'].dropna():
    parts = [p.strip() for p in str(val).split(',') if p.strip()]
    sensor_elements.extend(parts)

st_counts = pd.Series(sensor_elements).value_counts()
print("\\n=== Consolidated Unified Sensor Type Distribution (exploded multi-values) ===")
for cat, count in st_counts.items():
    print(f"  * {cat:12}: {count:5,} occurrences ({count/len(df)*100:6.2f}%)")

# 4. Plot Unified Sensor Type Distribution
plt.figure(figsize=(10, 5))
colors_dist = sns.color_palette('viridis', len(st_counts))
bars_dist = plt.bar(st_counts.index, st_counts.values, color=colors_dist, edgecolor='gray', width=0.4)
plt.title(f'Consolidated Unified Sensor Type Distribution (Exploded)', fontweight='bold', fontsize=15, pad=15)
plt.ylabel('Occurrence Count', fontweight='bold', fontsize=12)
plt.ylim(0, max(st_counts.values) * 1.15)
for bar in bars_dist:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval + 100, f"{yval:,}\\n({yval/len(df)*100:.1f}%)", 
             va='bottom', ha='center', fontsize=10, fontweight='bold')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()"""

    cell_17_markdown = """### 5.2. Enzyme Product/Substrate Relationships

#### Purpose
Focuses specifically on the subset of interactions classified as "Enzymes" to determine if the target protein treats the metabolite as a **substrate** (consuming it) or a **product** (synthesizing/releasing it). This defines the localized metabolic flux in the TME.

#### Database Origins & Multi-Database Mapping
Biochemical classifications are mapped across different layers to determine how target enzymes interact with metabolic signals:
1. **`OtherDB_enzyme product/substrate`**: Original classifications retrieved from the enzyme table of **Celllinker2** and **MRCLinkDB** (`Homo sapiens enzyme.txt` and `Mus musculus enzyme.txt`), which map metabolic enzymes to their corresponding substrate (consumed) or product (synthesized/released) metabolites, as well as MEBOCOST mappings.
2. **`Rhea_enzyme product/substrate`**: External enrichment via the **Rhea** biochemical reaction database (expert-curated, non-redundant reactions used as the standard vocabulary for UniProt enzyme annotation). The enrichment pipeline uses three external APIs:
   - **BridgeDb** → maps HMDB IDs to ChEBI identifiers
   - **UniProt CATALYTIC ACTIVITY** → fetches enzyme reaction annotations with Rhea reaction IDs and ChEBI participant IDs per target gene
   - **Rhea SPARQL** → resolves reaction directionality (left side = substrates, right side = products) for each Rhea reaction
   - Matching: if the metabolite's ChEBI ID appears on the left side of a Rhea reaction catalyzed by the target enzyme → `s` (substrate); right side → `p` (product); both → `s+p`
3. **`KEGG_enzyme product/substrate`**: External enrichment via the **KEGG** REST API. Maps Target Symbols to UniProt IDs, then to KEGG Genes, then to KEGG Reactions. Parses EQUATIONs for KEGG Compounds (C-numbers) and maps back to HMDB IDs using ChEBI cross-references.
4. **`enzyme product/substrate` (Unified)**: The consolidated, row-wise standardized and deduplicated biochemical classification merging `OtherDB`, `Reactions` column parser, `Rhea`, and `KEGG`. Values are mapped to `s` (substrate), `p` (product), or `s+p` (both substrate and product roles for reactions where the target participates in multiple ways)."""

    cell_18_code = """# 1. Print Coverage Comparison for Enzyme annotations
print("=== Enzyme Product/Substrate Annotation Coverage ===")
cols_to_check_enz = {
    "OtherDB_enzyme product/substrate": "Original Database (OtherDB)",
    "Rhea_enzyme product/substrate": "Rhea (UniProt + BridgeDb + SPARQL)",
    "KEGG_enzyme product/substrate": "KEGG (REST API)",
    "enzyme product/substrate": "Consolidated Unified"
}

coverage_data_enz = {}
for col, label in cols_to_check_enz.items():
    if col in df.columns:
        cnt = df[col].notna().sum()
        pct = (cnt / len(df)) * 100
        print(f"  * {label:28}: {cnt:5,} / {len(df):5,} pairs ({pct:6.2f}%)")
        coverage_data_enz[label] = cnt
    else:
        print(f"  * {label:28}: Column '{col}' not found!")

# 2. Plot Consolidated Enzyme-Metabolite Relationship Types Distribution
eps = df['enzyme product/substrate'].value_counts(dropna=False)
nan_count = df['enzyme product/substrate'].isna().sum()
nan_percentage = df['enzyme product/substrate'].isna().mean() * 100

print(f"\\nUnified Enzyme Relationships stats:")
print(f"  * Number of unannotated values: {nan_count:,} ({nan_percentage:.2f}%)")

eps_clean = df['enzyme product/substrate'].value_counts(dropna=True)
plt.figure(figsize=(8, 4.5))
colors_eps = sns.color_palette('Set2', len(eps_clean))
bars_eps = plt.bar(eps_clean.index, eps_clean.values, color=colors_eps, edgecolor='gray', width=0.3)
plt.title(f'Enzyme-Metabolite Relationship Types (Unified)', fontweight='bold', fontsize=14, pad=15)
plt.ylabel('Occurrence Count', fontweight='bold', fontsize=11)
plt.xlabel('Relationship (s = substrate, p = product, s+p = both)', fontweight='bold', fontsize=11)
plt.ylim(0, max(eps_clean.values) * 1.15)
for bar in bars_eps:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval + 10, f"{yval:,}", 
             va='bottom', ha='center', fontsize=10, fontweight='bold')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()"""

    cell_19_markdown = """### 5.3. Interaction Type (Promote/Inhibit/Release/Consume)

#### Purpose
Beyond simple binding, this section analyzes the functional outcome of the interaction in clinical oncology and cellular signaling contexts. It defines whether the metabolite-target pairing triggers activation (promote), downregulation (inhibit), secretion (release), or uptake (consume).

#### Database Origins & Multi-Database Mapping
Interaction dynamics are parsed and harmonized to analyze biological effects in the tumor microenvironment:
1. **`OtherDB_Interaction`**: Interaction directions originally retrieved from the clinical oncology tables of **MRCLinkDB** (`Metabolite-cell interaction.txt`), representing curated experimentally validated clinical and microenvironmental effects.
2. **`GtoPdb_Interaction`**: Interaction directions derived from the **IUPHAR/BPS Guide to Pharmacology** (`input/interactions.csv`), where target–ligand interaction Type and Action fields (e.g., Agonist, Antagonist, Inhibitor) are mapped to our standardized vocabulary using `map_gtopdb_action_to_interaction()` in `annotate_with_databases.py`. Matching is performed by PubChem SID or metabolite name against each target gene symbol's known ligand interactions.
3. **`Interaction` (Unified)**: The consolidated column produced by merging `OtherDB_Interaction` and `GtoPdb_Interaction` via `consolidate_interaction()`. Values are standardized to: `Promote`, `Inhibit`, `Be released`, `Be consumed`, `Regulate`, `Protect`."""

    cell_20_code = """# 1. Print Coverage Comparison for Interaction annotations
print("=== Interaction Type Annotation Coverage ===")
cols_to_check_int = {
    "OtherDB_Interaction": "MRCLinkDB (OtherDB)",
    "GtoPdb_Interaction": "Guide to Pharmacology (GtoPdb)",
    "Interaction": "Consolidated Unified"
}

coverage_data_int = {}
for col, label in cols_to_check_int.items():
    if col in df.columns:
        cnt = df[col].notna().sum()
        pct = (cnt / len(df)) * 100
        print(f"  * {label:28}: {cnt:5,} / {len(df):5,} pairs ({pct:6.2f}%)")
        coverage_data_int[label] = cnt
    else:
        print(f"  * {label:28}: Column '{col}' not found!")

# 2. Plot Consolidated Interaction Type Distribution
inter = df['Interaction'].value_counts(dropna=False)
nan_count = df['Interaction'].isna().sum()
nan_percentage = df['Interaction'].isna().mean() * 100
print(f"\\nUnified Interaction Directionality stats:")
print(f"  * Number of unannotated values: {nan_count:,} ({nan_percentage:.2f}%)")

inter_clean = df['Interaction'].value_counts(dropna=True)
colors_int = {
    'Promote': '#e74c3c', 
    'Inhibit': '#3498db', 
    'Be released': '#2ecc71', 
    'Be consumed': '#f39c12', 
    'Regulate': '#9b59b6', 
    'Protect': '#1abc9c'
}

plt.figure(figsize=(10, 5))
# Generate horizontal bar plot
bars_int = plt.barh(inter_clean.index, inter_clean.values, 
                    color=[colors_int.get(i, '#95a5a6') for i in inter_clean.index], 
                    edgecolor='gray', height=0.6)
plt.title(f'Unified Interaction Directionality Distribution', fontweight='bold', fontsize=14, pad=15)
plt.xlabel('Count', fontweight='bold', fontsize=11)
plt.ylabel('Interaction Type', fontweight='bold', fontsize=11)
plt.xlim(0, max(inter_clean.values) * 1.15)
for bar in bars_int:
    width = bar.get_width()
    plt.text(width + 5, bar.get_y() + bar.get_height()/2, f"{width:,}", 
             va='center', ha='left', fontsize=10, fontweight='bold')
plt.grid(axis='x', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()"""

    # Helper function to split text block into line array (Jupyter format)
    def to_jupyter_lines(text):
        return [line + '\n' for line in text.split('\n')]
        
    # Replace content of cells programmatically based on index checks
    cells = nb['cells']
    
    # 5.1 Markdown
    if '### 5.1. Sensor Type Distribution' in ''.join(cells[15]['source']):
        cells[15]['source'] = to_jupyter_lines(cell_15_markdown)
        # remove trailing split newline if needed
        if cells[15]['source'][-1] == '\n':
            cells[15]['source'].pop()
        modified_count += 1
        print("✏️ Programmatically updated Cell 15 (Section 5.1 Markdown)")
    else:
        print("⚠️ Warning: Cell 15 did not match '### 5.1. Sensor Type Distribution'")
        
    # 5.1 Code
    if "Sensor_Type" in ''.join(cells[16]['source']):
        cells[16]['source'] = to_jupyter_lines(cell_16_code)
        if cells[16]['source'][-1] == '\n':
            cells[16]['source'].pop()
        modified_count += 1
        print("✏️ Programmatically updated Cell 16 (Section 5.1 Code)")
    else:
        print("⚠️ Warning: Cell 16 did not match expected code content")
        
    # 5.2 Markdown
    if '### 5.2. Enzyme Product/Substrate Relationships' in ''.join(cells[17]['source']):
        cells[17]['source'] = to_jupyter_lines(cell_17_markdown)
        if cells[17]['source'][-1] == '\n':
            cells[17]['source'].pop()
        modified_count += 1
        print("✏️ Programmatically updated Cell 17 (Section 5.2 Markdown)")
    else:
        print("⚠️ Warning: Cell 17 did not match '### 5.2. Enzyme Product/Substrate Relationships'")
        
    # 5.2 Code
    if "enzyme product/substrate" in ''.join(cells[18]['source']):
        cells[18]['source'] = to_jupyter_lines(cell_18_code)
        if cells[18]['source'][-1] == '\n':
            cells[18]['source'].pop()
        modified_count += 1
        print("✏️ Programmatically updated Cell 18 (Section 5.2 Code)")
    else:
        print("⚠️ Warning: Cell 18 did not match expected code content")
        
    # 5.3 Markdown
    if '### 5.3. Interaction Type (Promote/Inhibit/Release/Consume)' in ''.join(cells[19]['source']):
        cells[19]['source'] = to_jupyter_lines(cell_19_markdown)
        if cells[19]['source'][-1] == '\n':
            cells[19]['source'].pop()
        modified_count += 1
        print("✏️ Programmatically updated Cell 19 (Section 5.3 Markdown)")
    else:
        print("⚠️ Warning: Cell 19 did not match '### 5.3. Interaction Type (Promote/Inhibit/Release/Consume)'")
        
    # 5.3 Code
    if "Interaction" in ''.join(cells[20]['source']):
        cells[20]['source'] = to_jupyter_lines(cell_20_code)
        if cells[20]['source'][-1] == '\n':
            cells[20]['source'].pop()
        modified_count += 1
        print("✏️ Programmatically updated Cell 20 (Section 5.3 Code)")
    else:
        print("⚠️ Warning: Cell 20 did not match expected code content")

    # Save notebook back
    with open(notebook_path, 'w', encoding='utf-8') as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)
        
    print(f"\n🎉 Finished programmatically updating {modified_count} cells in {notebook_path}!")

if __name__ == '__main__':
    edit_notebook()
