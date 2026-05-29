import json
import os

def update_notebook(path):
    print(f"Processing {path}...")
    with open(path, "r", encoding="utf-8") as f:
        nb = json.load(f)
    
    modified = False
    for cell in nb.get("cells", []):
        if cell.get("cell_type") == "code":
            source_lines = cell.get("source", [])
            new_lines = []
            cell_modified = False
            for line in source_lines:
                # Replacement 1: Dataset loaded in unique_metab_data_exploration.ipynb
                if "Dataset loaded: {len(df):,} unique metabolites" in line:
                    line = line.replace(
                        "print(f'Dataset loaded: {len(df):,} unique metabolites, {df.shape[1]} columns')",
                        "print(f'Dataset loaded: {df[\"HMDB_ID\"].nunique():,} unique HMDB metabolites ({df[\"Metabolite_Name\"].nunique()} names, {len(df):,} total rows), {df.shape[1]} columns')"
                    )
                    cell_modified = True
                
                # Replacement 2: Summary table in unique_metab_data_exploration.ipynb
                elif "'Total unique metabolites': len(df)," in line:
                    line = line.replace(
                        "'Total unique metabolites': len(df),",
                        "'Total unique metabolites (HMDB IDs)': df['HMDB_ID'].nunique(),\n    'Total unique metabolite names': df['Metabolite_Name'].nunique(),\n    'Total database entries': len(df),"
                    )
                    cell_modified = True
                
                # Replacement 3: Loaded count in metab_targetPair_analysis.ipynb
                elif "print(f'Loaded {len(df):,} interaction pairs | {df[\"Metabolite_Name\"].nunique()} metabolites | {df[\"Target\"].nunique()} targets | {df.shape[1]} columns')" in line:
                    line = line.replace(
                        "print(f'Loaded {len(df):,} interaction pairs | {df[\"Metabolite_Name\"].nunique()} metabolites | {df[\"Target\"].nunique()} targets | {df.shape[1]} columns')",
                        "print(f'Loaded {len(df):,} interaction pairs | {df[\"HMDB_ID\"].nunique()} unique HMDB metabolites ({df[\"Metabolite_Name\"].nunique()} names) | {df[\"Target\"].nunique()} targets | {df.shape[1]} columns')"
                    )
                    cell_modified = True
                
                # Replacement 4: Remaining counts in metab_targetPair_analysis.ipynb
                elif 'print(f"   -> Unique metabolites remaining: {df[\'Metabolite_Name\'].nunique():,}")' in line:
                    line = line.replace(
                        'print(f"   -> Unique metabolites remaining: {df[\'Metabolite_Name\'].nunique():,}")',
                        'print(f"   -> Unique HMDB metabolites remaining: {df[\'HMDB_ID\'].nunique():,} ({df[\'Metabolite_Name\'].nunique():,} names)")'
                    )
                    cell_modified = True
                
                # Replacement 5: Export counts dict in metab_targetPair_analysis.ipynb
                elif "'Unique metabolites': df['Metabolite_Name'].nunique()," in line:
                    line = line.replace(
                        "'Unique metabolites': df['Metabolite_Name'].nunique(),",
                        "'Unique HMDB metabolites': df['HMDB_ID'].nunique(),\n        'Unique metabolite names': df['Metabolite_Name'].nunique(),"
                    )
                    cell_modified = True
                
                new_lines.append(line)
            
            if cell_modified:
                cell["source"] = new_lines
                modified = True
                
    if modified:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(nb, f, indent=1, ensure_ascii=False)
        print(f"✅ Successfully updated {path}")
    else:
        print(f"ℹ️ No replacements made in {path}")

if __name__ == "__main__":
    update_notebook("scripts/unique_metab_data_exploration.ipynb")
    update_notebook("scripts/metab_targetPair_analysis.ipynb")
