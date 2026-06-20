#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
📝 Notebook Addition Script (Literature Evidence & Temporal Analysis)
===================================================================
This script programmatically inserts/updates a premium, publication-grade 
literature evidence and temporal (PMID publication year) analysis section 
into the Jupyter Notebook scripts/metab_targetPair_analysis.ipynb.

Author: Antigravity (Advanced Agentic Coding Pair)
Date: 2026-05-21
"""

import json
import os

def update_notebook():
    notebook_path = 'scripts/metab_targetPair_analysis.ipynb'
    if not os.path.exists(notebook_path):
        print(f"Error: {notebook_path} not found.")
        return

    print(f"Loading notebook '{notebook_path}'...")
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = json.load(f)

    # 1. Prepare new Markdown Cell with Temporal Analysis details
    markdown_cell = {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "---\n",
            "\n",
            "### 4.2. Literature Evidence, Temporal Dynamics & PubMed References Validation\n",
            "\n",
            "### Purpose\n",
            "While database consensus indicates standard curation agreement, empirical validation from peer-reviewed scientific literature is the ultimate standard of biological confidence. This section evaluates metabolite-target interaction pairs against empirical literature support, cleaning and parsing PubMed IDs (PMIDs) and their corresponding **publication years** to explore the historical depth and recency of research backing our database.\n",
            "\n",
            "### Analytical Pillars\n",
            "1. **Evidence Rate by Tier:** Verifying if our confidence tiers (Tiers 1, 2, 3) map to active literature. Higher tiers should exhibit high citation rates, validating the classification logic.\n",
            "2. **Research Density (Top Hubs):** Identifying the specific metabolite-target pairs with the highest PMID citation counts to reveal the central scientific communication nodes.\n",
            "3. **Temporal Dynamics (Publication Years):** Exploring the historical trajectory of research from legacy classical findings (e.g., 1980s or earlier) to cutting-edge recent discoveries (2020-2026), mapping active vs. historical interaction nodes."
        ]
    }

    # 2. Prepare new Code Cell with robust parsing, calculations, temporal stats, and Three-Panel visualization
    code_cell = {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# 📚 4.2. LITERATURE EVIDENCE (PMID/REFERENCES) & TEMPORAL ANALYSIS\n",
            "# ==============================================================================\n",
            "import re\n",
            "import os\n",
            "import numpy as np\n",
            "import pandas as pd\n",
            "import matplotlib.pyplot as plt\n",
            "import seaborn as sns\n",
            "\n",
            "# ⚙️ CONFIGURABLE PARAMETER\n",
            "TOP_N_PAIRS = 10  # Customize the number of top referenced pairs to display in Panel B\n",
            "\n",
            "print(f\"Running Literature Evidence & Temporal Analysis (Top {TOP_N_PAIRS} Pairs in Panel B)...\\n\")\n",
            "\n",
            "# 1. Load PubMed Scraped Metadata for Publication Year mapping\n",
            "pubmed_results_path = '../input/pubmed_results.csv'\n",
            "if os.path.exists(pubmed_results_path):\n",
            "    pubmed_df = pd.read_csv(pubmed_results_path, low_memory=False)\n",
            "    # Standardize PMID to string for robust mapping\n",
            "    pubmed_df['PMID'] = pubmed_df['PMID'].astype(str).str.strip()\n",
            "    \n",
            "    # Standardize Year column (extract numeric 4-digit year)\n",
            "    def clean_year(y):\n",
            "        if pd.isna(y):\n",
            "            return np.nan\n",
            "        y_str = str(y).strip()\n",
            "        matches = re.findall(r'\\b(19\\d{2}|20\\d{2})\\b', y_str)\n",
            "        if matches:\n",
            "            return int(matches[0])\n",
            "        return np.nan\n",
            "        \n",
            "    pubmed_df['Clean_Year'] = pubmed_df['Year'].apply(clean_year)\n",
            "    pmid_to_year = dict(zip(pubmed_df['PMID'], pubmed_df['Clean_Year']))\n",
            "    print(f\"\\u2705 Loaded {len(pubmed_df):,} unique PMIDs from PubMed metadata with publication years.\")\n",
            "else:\n",
            "    pmid_to_year = {}\n",
            "    print(\"\\u26a0\\ufe0f Warning: 'input/pubmed_results.csv' not found. Temporal analysis will be limited.\")\n",
            "\n",
            "# 2. Clean, Coalesce, and Standardize PMIDs Row-wise\n",
            "def clean_row_pmids(row):\n",
            "    all_pmids = set()\n",
            "    # Extract PMIDs from PMID and Evidence columns\n",
            "    for col in ['PMID', 'Evidence']:\n",
            "        if col in row and pd.notna(row[col]):\n",
            "            val_str = str(row[col]).strip()\n",
            "            # Split by common separators: semicolon, comma, pipe, space, slash\n",
            "            for chunk in re.split(r'[;|,\\s/|]+', val_str):\n",
            "                chunk = chunk.strip()\n",
            "                if chunk.endswith('.0'):\n",
            "                    chunk = chunk[:-2]\n",
            "                if chunk.isdigit() and len(chunk) >= 5:\n",
            "                    all_pmids.add(chunk)\n",
            "                elif chunk and chunk.lower() not in ['nan', 'none', 'null']:\n",
            "                    clean_num = re.findall(r'\\b\\d{5,}\\b', chunk)\n",
            "                    for num in clean_num:\n",
            "                        all_pmids.add(num)\n",
            "                        \n",
            "    # Extract PMIDs from Text_Evidence (typically paragraph text with citations)\n",
            "    if 'Text_Evidence' in row and pd.notna(row['Text_Evidence']):\n",
            "        val_str = str(row['Text_Evidence']).strip()\n",
            "        matches = re.findall(r'(?i)(?:pmid|pubmed)\\s*:?\\s*(\\d+)', val_str)\n",
            "        for m in matches:\n",
            "            if len(m) >= 5:\n",
            "                all_pmids.add(m)\n",
            "                \n",
            "    if not all_pmids:\n",
            "        return np.nan\n",
            "    return ';'.join(sorted(list(all_pmids), key=lambda x: int(x) if x.isdigit() else x))\n",
            "\n",
            "# Apply the cleaning function to the active dataframe df\n",
            "df['Cleaned_PMID'] = df.apply(clean_row_pmids, axis=1)\n",
            "df['evidence_count'] = df['Cleaned_PMID'].apply(lambda x: len([p for p in str(x).split(';') if p.strip()]) if pd.notna(x) else 0)\n",
            "df['has_evidence'] = df['evidence_count'] > 0\n",
            "\n",
            "# 3. Map PMIDs to Publication Years for Temporal Analytics\n",
            "def get_row_years(pmid_str):\n",
            "    if pd.isna(pmid_str):\n",
            "        return []\n",
            "    years = []\n",
            "    for p in str(pmid_str).split(';'):\n",
            "        p = p.strip()\n",
            "        if p in pmid_to_year and pd.notna(pmid_to_year[p]):\n",
            "            years.append(int(pmid_to_year[p]))\n",
            "    return sorted(years)\n",
            "\n",
            "df['evidence_years'] = df['Cleaned_PMID'].apply(get_row_years)\n",
            "df['earliest_year'] = df['evidence_years'].apply(lambda x: min(x) if x else np.nan)\n",
            "df['latest_year'] = df['evidence_years'].apply(lambda x: max(x) if x else np.nan)\n",
            "df['median_year'] = df['evidence_years'].apply(lambda x: np.nanmedian(x) if x else np.nan)\n",
            "\n",
            "# 4. Calculate Global Temporal Metrics\n",
            "all_pmids_set = set()\n",
            "for p in df['Cleaned_PMID'].dropna():\n",
            "    all_pmids_set.update([x.strip() for x in str(p).split(';') if x.strip()])\n",
            "num_unique_pmids = len(all_pmids_set)\n",
            "\n",
            "unique_pmid_years = []\n",
            "for p in all_pmids_set:\n",
            "    if p in pmid_to_year and pd.notna(pmid_to_year[p]):\n",
            "        unique_pmid_years.append(int(pmid_to_year[p]))\n",
            "\n",
            "total_pairs = len(df)\n",
            "with_evidence = df['has_evidence'].sum()\n",
            "pct_evidence = (with_evidence / total_pairs) * 100 if total_pairs > 0 else 0\n",
            "\n",
            "print(f\"\\u2705 Literature Analysis Complete:\")\n",
            "print(f\"   -> Total unique interaction pairs in database: {total_pairs:,}\")\n",
            "print(f\"   -> Pairs supported by literature: {with_evidence:,} ({pct_evidence:.2f}%)\")\n",
            "print(f\"   -> Total unique PMIDs cited: {num_unique_pmids:,}\")\n",
            "\n",
            "if unique_pmid_years:\n",
            "    min_yr = min(unique_pmid_years)\n",
            "    max_yr = max(unique_pmid_years)\n",
            "    med_yr = int(np.nanmedian(unique_pmid_years)) if not np.isnan(np.nanmedian(unique_pmid_years)) else np.nan\n",
            "    print(f\"   -> Temporal span of cited literature: {min_yr} to {max_yr} (Median: {med_yr})\")\n",
            "    \n",
            "    # Recency statistics (2020 onwards)\n",
            "    recent_pmids = sum(1 for y in unique_pmid_years if y >= 2020)\n",
            "    recent_pairs = sum(1 for y in df['latest_year'].dropna() if y >= 2020)\n",
            "    print(f\"   -> Modern research (2020+): {recent_pmids} unique PMIDs supporting {recent_pairs} pairs\")\n",
            "else:\n",
            "    med_yr = np.nan\n",
            "\n",
            "# 5. Create Premium Three-Panel Figure\n",
            "sns.set_theme(style=\"whitegrid\")\n",
            "fig, axes = plt.subplots(1, 3, figsize=(22, 6.5))\n",
            "\n",
            "# --- Panel A: Literature Evidence Rate by Confidence Tier ---\n",
            "tier_stats = df.groupby('Confidence_Tier')['has_evidence'].mean() * 100\n",
            "tier_stats = tier_stats.reindex(['Tier 1 (High)', 'Tier 2 (Medium)', 'Tier 3 (Low)']).fillna(0)\n",
            "\n",
            "colors_a = ['#1a73e8', '#8ab4f8', '#ffad46']  # Harmonious custom palette\n",
            "bars_a = axes[0].bar(tier_stats.index, tier_stats.values, color=colors_a, edgecolor='black', width=0.45)\n",
            "\n",
            "# Add values on top of bars\n",
            "for bar in bars_a:\n",
            "    yval = bar.get_height()\n",
            "    axes[0].text(bar.get_x() + bar.get_width()/2.0, yval + 1.5, f'{yval:.1f}%', ha='center', va='bottom', fontweight='bold', fontsize=11)\n",
            "\n",
            "axes[0].set_ylim(0, 108)\n",
            "axes[0].set_ylabel('Pairs with Literature Evidence (%)', fontweight='bold', fontsize=12)\n",
            "axes[0].set_title('Literature Evidence Rate by Confidence Tier', fontweight='bold', fontsize=13, pad=15)\n",
            "axes[0].tick_params(axis='both', which='major', labelsize=11)\n",
            "\n",
            "# --- Panel B: Top N Pairs by Reference Count ---\n",
            "top_pairs = df.sort_values(by='evidence_count', ascending=False).head(TOP_N_PAIRS).copy()\n",
            "top_pairs['Pair_Label'] = top_pairs['Metabolite_Name'] + ' \\u2194 ' + top_pairs['Target']\n",
            "\n",
            "colors_b = sns.color_palette(\"flare\", n_colors=TOP_N_PAIRS)\n",
            "bars_b = axes[1].barh(top_pairs['Pair_Label'], top_pairs['evidence_count'], color=colors_b, edgecolor='black')\n",
            "\n",
            "# Add values on right of horizontal bars\n",
            "for bar in bars_b:\n",
            "    xval = bar.get_width()\n",
            "    axes[1].text(xval + max(1, xval*0.02), bar.get_y() + bar.get_height()/2.0, f'{int(xval):,}', ha='left', va='center', fontweight='bold', fontsize=10)\n",
            "\n",
            "axes[1].invert_yaxis()  # top-down ranking\n",
            "axes[1].set_xlabel('Number of Unique PMIDs', fontweight='bold', fontsize=12)\n",
            "axes[1].set_title(f'Top {TOP_N_PAIRS} Pairs by Reference Count', fontweight='bold', fontsize=13, pad=15)\n",
            "axes[1].tick_params(axis='both', which='major', labelsize=10)\n",
            "\n",
            "# --- Panel C: Temporal Trend / Reference Year Distribution ---\n",
            "if unique_pmid_years:\n",
            "    # Premium histogram with a smooth KDE line for cited paper years\n",
            "    sns.histplot(unique_pmid_years, bins=max(10, min(30, max_yr - min_yr)), kde=True, ax=axes[2], color='#2ca02c', edgecolor='black', alpha=0.6)\n",
            "    \n",
            "    # Highlight and label the median year with a distinct indicator line\n",
            "    if pd.notna(med_yr):\n",
            "        axes[2].axvline(med_yr, color='#d62728', linestyle='--', linewidth=2, label=f'Median: {int(med_yr)}')\n",
            "        axes[2].legend(fontsize=11, loc='upper left')\n",
            "    \n",
            "    axes[2].set_xlabel('Publication Year', fontweight='bold', fontsize=12)\n",
            "    axes[2].set_ylabel('Number of Unique PMIDs', fontweight='bold', fontsize=12)\n",
            "    axes[2].set_title('Publication Year Distribution of Cited PMIDs', fontweight='bold', fontsize=13, pad=15)\n",
            "    axes[2].tick_params(axis='both', which='major', labelsize=11)\n",
            "else:\n",
            "    axes[2].text(0.5, 0.5, \"No Year Information Available\", ha='center', va='center', fontsize=14, color='gray')\n",
            "    axes[2].set_title('Publication Year Distribution', fontweight='bold', fontsize=13, pad=15)\n",
            "\n",
            "plt.tight_layout()\n",
            "plt.show()\n",
            "\n",
            "# 6. Detailed Tier breakdown with Median Publication Year\n",
            "print(\"\\n=== Literature Evidence & Temporal Statistics by Tier ===\")\n",
            "tier_counts = df.groupby('Confidence_Tier')['has_evidence'].agg(['count', 'sum'])\n",
            "for tier in ['Tier 1 (High)', 'Tier 2 (Medium)', 'Tier 3 (Low)']:\n",
            "    if tier in tier_counts.index:\n",
            "        cnt = tier_counts.loc[tier, 'count']\n",
            "        ev = tier_counts.loc[tier, 'sum']\n",
            "        pct = (ev / cnt) * 100 if cnt > 0 else 0\n",
            "        \n",
            "        # Calculate median publication year for references in this tier\n",
            "        tier_df = df[df['Confidence_Tier'] == tier]\n",
            "        tier_years = []\n",
            "        for yrs in tier_df['evidence_years'].dropna():\n",
            "            tier_years.extend(yrs)\n",
            "            \n",
            "        if tier_years:\n",
            "            tier_med = np.nanmedian(tier_years)\n",
            "            med_yr_str = f\"Median Reference Year: {int(tier_med)}\" if pd.notna(tier_med) else \"No reference years\"\n",
            "        else:\n",
            "            med_yr_str = \"No reference years\"\n",
            "        print(f\"  * {tier:18}: {ev:4,} / {cnt:5,} pairs ({pct:6.2f}%) | {med_yr_str}\")\n",
            "\n",
            "# 7. Show temporal landmark references (Oldest vs Newest in database)\n",
            "if unique_pmid_years and os.path.exists(pubmed_results_path):\n",
            "    print(\"\\n=== Landmark References in Literature Evidence ===\")\n",
            "    \n",
            "    # Sort PMIDs by year from metadata\n",
            "    valid_meta = pubmed_df.dropna(subset=['Clean_Year']).copy()\n",
            "    if not valid_meta.empty:\n",
            "        # Find row with minimum year\n",
            "        oldest_row = valid_meta.loc[valid_meta['Clean_Year'].idxmin()]\n",
            "        print(f\"  * Oldest Reference ({int(oldest_row['Clean_Year'])}): PMID {oldest_row['PMID']} - \\\"{oldest_row['Title']}\\\"\")\n",
            "        print(f\"    Journal: {oldest_row['Journal']} | Authors: {oldest_row['Author']}\\n\")\n",
            "        \n",
            "        # Find row with maximum year\n",
            "        newest_row = valid_meta.loc[valid_meta['Clean_Year'].idxmax()]\n",
            "        print(f\"  * Newest Reference ({int(newest_row['Clean_Year'])}): PMID {newest_row['PMID']} - \\\"{newest_row['Title']}\\\"\")\n",
            "        print(f\"    Journal: {newest_row['Journal']} | Authors: {newest_row['Author']}\")"
        ]
    }

    # 3. Find insertion index after Section 4.1 code cell
    target_index = -1
    for i, cell in enumerate(nb['cells']):
        if cell['cell_type'] == 'markdown':
            content = ''.join(cell['source'])
            if '### 4.1. Cross-Database Validation of Pairs' in content:
                # The code cell is at index i + 1
                # We insert the new cells after that, i.e., starting at index i + 2
                target_index = i + 2
                break

    if target_index != -1:
        # Check if 4.2 already exists to prevent duplicate insertions
        already_exists_idx = -1
        for i, cell in enumerate(nb['cells']):
            if cell['cell_type'] == 'markdown':
                content = ''.join(cell['source'])
                if '### 4.2. Literature Evidence' in content:
                    already_exists_idx = i
                    break

        if already_exists_idx == -1:
            print(f"Inserting new literature evidence cells at index {target_index}...")
            nb['cells'].insert(target_index, markdown_cell)
            nb['cells'].insert(target_index + 1, code_cell)
            print("Successfully inserted cells.")
        else:
            print(f"Literature evidence cells already exist. Overwriting cells at index {already_exists_idx} and {already_exists_idx + 1} with temporal analysis version...")
            nb['cells'][already_exists_idx] = markdown_cell
            nb['cells'][already_exists_idx + 1] = code_cell
            print("Successfully updated cells.")
    else:
        print("Warning: Could not find Section 4.1 insertion marker.")

    # 4. Find and update the Summary Statistics cell to include year details safely
    summary_index = -1
    for i, cell in enumerate(nb['cells']):
        if cell['cell_type'] == 'code':
            content = ''.join(cell['source'])
            if 'summary = {' in content and 'Unique HMDB metabolites' in content:
                summary_index = i
                break

    if summary_index != -1:
        print(f"Updating Summary Statistics scorecard cell at index {summary_index}...")
        nb['cells'][summary_index]['source'] = [
            "summary = {\n",
            "    'Total interaction pairs': f'{len(df):,}',\n",
            "    'Unique HMDB metabolites': df['HMDB_ID'].nunique(),\n",
            "    'Unique metabolite names': df['Metabolite_Name'].nunique(),\n",
            "    'Unique targets': df['Target'].nunique(),\n",
            "    'Cancer-annotated pairs': cancer_mask.sum(),\n",
            "    'CCC-ready pairs (sender+receiver)': ccc_ready.sum(),\n",
            "    'High-conf CCC-ready (3+ DB)': len(hc),\n",
            "    'With Sensor_Type annotation': df['Sensor_Type'].notna().sum(),\n",
            "    'With Disease annotation': df['Disease'].notna().sum(),\n",
            "    'With scCellFie score': df['scCellFie_value'].notna().sum(),\n",
            "    'With Evidence_Score': df['Evidence_Score'].notna().sum(),\n",
            "    'Pairs with Literature Evidence (PMID)': df['has_evidence'].sum() if 'has_evidence' in df.columns else 0,\n",
            "    'Total unique PMIDs cited': len(set(p.strip() for p in df['Cleaned_PMID'].dropna().str.split(';').explode().unique() if p.strip())) if 'Cleaned_PMID' in df.columns else 0,\n",
            "    'Median Publication Year of Citations': int(np.nanmedian([pmid_to_year[p] for p in set(p.strip() for p in df['Cleaned_PMID'].dropna().str.split(';').explode().unique() if p.strip()) if 'pmid_to_year' in globals() and p in pmid_to_year and pd.notna(pmid_to_year[p])])) if ('Cleaned_PMID' in df.columns and 'pmid_to_year' in globals() and any(p in pmid_to_year and pd.notna(pmid_to_year[p]) for p in set(p.strip() for p in df['Cleaned_PMID'].dropna().str.split(';').explode().unique() if p.strip()))) else 'N/A',\n",
            "}\n",
            "print(pd.DataFrame.from_dict(summary, orient='index', columns=['Value']))\n"
        ]
        print("Successfully updated Summary Statistics scorecard.")
    else:
        print("Warning: Could not find Summary Statistics cell.")

    # Save modified notebook
    print(f"Saving notebook back to '{notebook_path}'...")
    with open(notebook_path, 'w', encoding='utf-8') as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)
    print("Notebook saved successfully!")

if __name__ == '__main__':
    update_notebook()
