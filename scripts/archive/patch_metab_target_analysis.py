import json
import os

notebook_path = "metab_targetPair_analysis.ipynb"

with open(notebook_path, "r") as f:
    nb = json.load(f)

for cell in nb.get("cells", []):
    if cell.get("cell_type") == "markdown":
        new_source = []
        for line in cell.get("source", []):
            if "Local Rule Classifier (classify_by_rules)" in line:
                break
            if "To ensure **100% annotation completeness**" in line:
                break
            if line.startswith("- **Channels**: Matches prefixes"):
                break
            if line.startswith("- **Receptors**: Matches GPCR"):
                break
            if line.startswith("- **Transporters**: Matches solute"):
                break
            if line.startswith("- **Enzymes**: Matches metabolizing"):
                break
            new_source.append(line)
        # Strip the last empty lines if any
        while new_source and new_source[-1] == "\n":
            new_source.pop()
        cell["source"] = new_source

    if cell.get("cell_type") == "code":
        source_text = "".join(cell.get("source", []))
        if "Annotation Coverage Comparison" in source_text and "plt.barh" in source_text:
            new_code = [
                "# 2. Visualize Coverage Comparison as a Pie Chart\n",
                "plt.figure(figsize=(8, 8))\n",
                "consolidated_label = \"Consolidated Unified (with regex fallbacks)\"\n",
                "annotated_count = coverage_data.get(consolidated_label, 0)\n",
                "unannotated_count = len(df) - annotated_count\n",
                "labels = [f\"Annotated ({annotated_count:,})\", f\"Unannotated ({unannotated_count:,})\"]\n",
                "sizes = [annotated_count, unannotated_count]\n",
                "colors = ['#ff9999', '#66b3ff']\n",
                "explode = (0.1, 0)\n",
                "plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False, startangle=140, textprops={'fontsize': 12, 'fontweight': 'bold'})\n",
                "plt.title('Annotation Coverage (Consolidated Unified)', fontweight='bold', fontsize=15, pad=15)\n",
                "plt.tight_layout()\n",
                "plt.show()\n"
            ]
            
            # Find where to replace
            lines = cell.get("source", [])
            out_lines = []
            skip = False
            for line in lines:
                if "# 2. Visualize Coverage Comparison" in line:
                    skip = True
                    out_lines.extend(new_code)
                elif "# 3. Explode and count multi-valued" in line:
                    skip = False
                    out_lines.append("\n")
                    out_lines.append(line)
                elif not skip:
                    out_lines.append(line)
            cell["source"] = out_lines

with open(notebook_path, "w") as f:
    json.dump(nb, f, indent=1)
    f.write("\n")
