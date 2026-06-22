import json

nb_path = 'unique_metab_data_exploration.ipynb'
with open(nb_path, 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        new_source = []
        for line in cell['source']:
            # Replace Communication Modes with High-Level Mass Bins
            line = line.replace('Communication Modes (The Peaks)', 'High-Level Mass Bins')
            line = line.replace('High-Level Communication Modes (Horizontal Bar)', 'High-Level Mass Bins (Horizontal Bar)')
            line = line.replace("labels_comm = ['Paracrine', 'GPCR/Hormonal', 'Juxtacrine/Vesicular']", 
                                "labels_comm = ['<300 Da', '300-750 Da', '>750 Da']")
            
            if 'Paracrine / Soluble Transmitters (<300 Da)' in line:
                line = "    '<300 Da',\n"
            elif 'GPCR & Hormonal / Bioactive Lipids (300-750 Da)' in line:
                line = "    '300-750 Da',\n"
            elif 'Juxtacrine & Vesicular / Structural Anchors (>750 Da)' in line:
                line = "    '>750 Da'\n"
            
            line = line.replace("ax.set_title('Metabolite Distribution by Downstream Communication Mode'", 
                                "ax.set_title('Metabolite Distribution by High-Level Mass Bins'")
            
            new_source.append(line)
        cell['source'] = new_source

with open(nb_path, 'w') as f:
    json.dump(nb, f, indent=1)

print("Notebook patched successfully!")
