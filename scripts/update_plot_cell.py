#!/usr/bin/env python3
"""
Update the orphan_metabolic_immune_evasion.ipynb plot cell.

- Replaces the seaborn scatterplot that used `size='logfoldchanges'`
  with a plain matplotlib scatter that maps `logfoldchanges` → color
  via the viridis colormap.
- Adds a colorbar labeled “Log2 Fold Change”.
- Preserves any existing imports and data preparation steps.
"""

import json
import pathlib
import sys

NOTEBOOK_PATH = pathlib.Path(
    "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/"
    "Personal/05_Python_repositories/metabConnectomeDB/scripts/"
    "orphan_metabolic_immune_evasion.ipynb"
)

# ----------------------------------------------------------------------
# Helper to locate the target cell (the one that creates the figure)
# ----------------------------------------------------------------------
def find_plot_cell(nb_json: dict) -> int | None:
    """
    Returns the index of the code cell that contains the original
    `plt.figure(figsize=(14, 10))` line.  If the notebook structure changes,
    the search will still succeed as long as that line is present.
    """
    for idx, cell in enumerate(nb_json.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        if "plt.figure(figsize=(14, 10))" in source:
            return idx
    return None


# ----------------------------------------------------------------------
# New plotting block (multi‑line string, exactly as you’d type in a notebook)
# ----------------------------------------------------------------------
NEW_PLOT_BLOCK = """\
plt.figure(figsize=(14, 10))
# Use a continuous viridis colormap to represent log2 fold changes
sc = plt.scatter(
    x=top_hits['Cell_Type'] if 'Cell_Type' in top_hits.columns else top_hits['group'],
    y=top_hits['Interaction'],
    c=top_hits['logfoldchanges'],
    cmap='viridis',
    s=200,          # fixed marker size for visual clarity
    alpha=0.8,
    edgecolors='w',
)
plt.title("Top Orphan Metabolic Interactions Enriched in Immune Cells", fontsize=16)
plt.xlabel("Immune Cell Population")
plt.ylabel("Metabolite → Receptor Pair")
cbar = plt.colorbar(sc)
cbar.set_label('Log2 Fold Change')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
"""

# ----------------------------------------------------------------------
def main() -> None:
    if not NOTEBOOK_PATH.is_file():
        sys.exit(f"❌ Notebook not found at {NOTEBOOK_PATH}")

    # Load the notebook JSON
    with NOTEBOOK_PATH.open("r", encoding="utf-8") as f:
        nb = json.load(f)

    cell_idx = find_plot_cell(nb)
    if cell_idx is None:
        sys.exit("❌ Could not locate the plotting cell – check that the notebook "
                 "contains `plt.figure(figsize=(14, 10))`.")

    # Replace the entire source of the cell with the new block
    nb["cells"][cell_idx]["source"] = [line + "\n" for line in NEW_PLOT_BLOCK.splitlines()]

    # Write back (overwrites the original notebook)
    with NOTEBOOK_PATH.open("w", encoding="utf-8") as f:
        json.dump(nb, f, indent=1)

    print(f"✅ Updated cell #{cell_idx} in {NOTEBOOK_PATH}")

if __name__ == "__main__":
    main()

