import json

with open("/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/oxygen_tension_analysis.ipynb", "r") as f:
    nb = json.load(f)

# Find the cell we want to replace
for cell in nb.get("cells", []):
    if cell["cell_type"] == "code":
        source = "".join(cell.get("source", []))
        if "plt.figure(figsize=(8, 6))" in source and "sns.scatterplot" in source:
            # Replace the source
            cell["source"] = [
                "fig, axes = plt.subplots(1, 2, figsize=(16, 6))\n",
                "\n",
                "# Clean NaN if any for both hypotheses\n",
                "df_plot_normal = df_res.dropna(subset=['O2_Tension_Normal_Pct', 'OXPHOS_Glycolysis_Ratio']).copy()\n",
                "df_plot_tumour = df_res.dropna(subset=['O2_Tension_Tumour_Pct', 'OXPHOS_Glycolysis_Ratio']).copy()\n",
                "\n",
                "# Hypotheses 1: Normal Tissue Physioxia\n",
                "sns.scatterplot(data=df_plot_normal, x='O2_Tension_Normal_Pct', y='OXPHOS_Glycolysis_Ratio', hue='Cancer', s=200, palette='Set1', ax=axes[0])\n",
                "sns.regplot(data=df_plot_normal, x='O2_Tension_Normal_Pct', y='OXPHOS_Glycolysis_Ratio', scatter=False, color='grey', line_kws={\"linestyle\":\"--\"}, ax=axes[0])\n",
                "axes[0].set_title(\"Normal Tissue Physioxia vs OXPHOS/Glycolysis\", fontsize=14, weight='bold')\n",
                "axes[0].set_xlabel(\"Normal Tissue Oxygen Tension (%)\", fontsize=12)\n",
                "axes[0].set_ylabel(\"OXPHOS / Glycolysis Ratio\", fontsize=12)\n",
                "axes[0].set_yscale(\"log\")\n",
                "axes[0].axhline(1.0, color='black', linewidth=1, linestyle=':')\n",
                "axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')\n",
                "\n",
                "if len(df_plot_normal) > 1:\n",
                "    r_n, p_n = stats.pearsonr(df_plot_normal['O2_Tension_Normal_Pct'], np.log2(df_plot_normal['OXPHOS_Glycolysis_Ratio']))\n",
                "    axes[0].annotate(f\"Pearson r: {r_n:.2f}\\np-value: {p_n:.3f}\", xy=(0.05, 0.85), xycoords='axes fraction', bbox=dict(boxstyle=\"round,pad=0.3\", fc=\"white\", ec=\"gray\", alpha=0.8))\n",
                "\n",
                "# Hypotheses 2: Tumour Tissue Hypoxia\n",
                "sns.scatterplot(data=df_plot_tumour, x='O2_Tension_Tumour_Pct', y='OXPHOS_Glycolysis_Ratio', hue='Cancer', s=200, palette='Set1', ax=axes[1])\n",
                "sns.regplot(data=df_plot_tumour, x='O2_Tension_Tumour_Pct', y='OXPHOS_Glycolysis_Ratio', scatter=False, color='grey', line_kws={\"linestyle\":\"--\"}, ax=axes[1])\n",
                "axes[1].set_title(\"Tumour Hypoxia vs OXPHOS/Glycolysis\", fontsize=14, weight='bold')\n",
                "axes[1].set_xlabel(\"Tumour Oxygen Tension (%)\", fontsize=12)\n",
                "axes[1].set_ylabel(\"OXPHOS / Glycolysis Ratio\", fontsize=12)\n",
                "axes[1].set_yscale(\"log\")\n",
                "axes[1].axhline(1.0, color='black', linewidth=1, linestyle=':')\n",
                "axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')\n",
                "\n",
                "if len(df_plot_tumour) > 1:\n",
                "    r_t, p_t = stats.pearsonr(df_plot_tumour['O2_Tension_Tumour_Pct'], np.log2(df_plot_tumour['OXPHOS_Glycolysis_Ratio']))\n",
                "    axes[1].annotate(f\"Pearson r: {r_t:.2f}\\np-value: {p_t:.3f}\", xy=(0.05, 0.85), xycoords='axes fraction', bbox=dict(boxstyle=\"round,pad=0.3\", fc=\"white\", ec=\"gray\", alpha=0.8))\n",
                "\n",
                "plt.tight_layout()\n",
                "\n",
                "# output path\n",
                "out_dir = os.path.join(os.path.abspath(os.path.join(os.getcwd(), '..')), \"output\", \"oxygen_tension\")\n",
                "os.makedirs(out_dir, exist_ok=True)\n",
                "filename = \"oxygen_tension_correlation\"+ ANALYSIS_SUFFIX + \".png\"\n",
                "output_plot = os.path.join(out_dir,filename)\n",
                "plt.savefig(output_plot, dpi=300, bbox_inches='tight')\n",
                "print(f\"Saved plot to {output_plot}\")\n",
                "plt.show()\n",
                "\n",
                "from IPython.display import display, Markdown\n",
                "md_str = \"### Oxygen Tension Data Citations\\n\"\n",
                "for idx, row in df_res.iterrows():\n",
                "    cancer = row['Cancer']\n",
                "    pmids = str(row['PMID_Reference']).split(',')\n",
                "    links = []\n",
                "    for p in pmids:\n",
                "        p = p.strip()\n",
                "        if p and p != 'nan' and p != 'None':\n",
                "            links.append(f\"[PMID {p}](https://pubmed.ncbi.nlm.nih.gov/{p}/)\")\n",
                "    if links:\n",
                "        md_str += f\"- **{cancer}**: {', '.join(links)}\\n\"\n",
                "display(Markdown(md_str))"
            ]

with open("/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/oxygen_tension_analysis.ipynb", "w") as f:
    json.dump(nb, f, indent=1)

print("Notebook patched successfully!")
