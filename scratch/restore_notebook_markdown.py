import json

with open('scripts/metab_targetPair_analysis.ipynb', 'r') as f:
    nb = json.load(f)

markdown_to_add = [
    "\n",
    "#### Local Rule Classifier (`classify_by_rules`) Guide\n",
    "To ensure maximum annotation coverage when external databases lack functional assignments, the pipeline implements a robust structural naming heuristic based on standard HGNC nomenclatures:\n",
    "- **Channels**: Matches prefixes like `CACN`, `SCN`, `KCN`, `TRP`, `AQP`, `PIEZO`, `GABR`, `CHRN`, `GLYR`, etc., or names containing 'channel' or 'pore-forming'.\n",
    "- **Receptors**: Matches GPCR and Nuclear Receptor families (`HRH`, `ADRA`, `HTR`, `OPR`, `MC1R`, `TAS1R`, `EGFR`, `IGF1R`, `NR1`, `ESR1`, etc.), or names containing 'receptor'.\n",
    "- **Transporters**: Matches carrier families like `SLC` (solute carriers), `ABC` (ATP-binding cassettes), `ATP`, `TFRC`, `FABP`, or names containing 'transporter' or 'carrier'.\n",
    "- **Enzymes**: Matches catalytic families (`CYP`, `ALDH`, `COX`, `NOS`, `MAPK`, `CDK`, `CASP`, etc.) and suffixes/keywords like `ase`, `kinase`, `synthase`, `dehydrogenase`.\n"
]

for cell in nb['cells']:
    if cell['cell_type'] == 'markdown':
        src = ''.join(cell['source'])
        if "Classifications extracted dynamically from **UniProt KB**" in src:
            # Check if not already added
            if "#### Local Rule Classifier" not in src:
                cell['source'].extend(markdown_to_add)
                print("Injected Local Rule Classifier explanation.")
                break

with open('scripts/metab_targetPair_analysis.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)
    f.write('\n')
