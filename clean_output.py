import os
import shutil
import glob

output_dir = "output"

# 1. Clean up duplicate HTMLs in breast_results
breast_results = os.path.join(output_dir, "breast_results")
if os.path.exists(breast_results):
    dup_html = os.path.join(breast_results, "primary_vs_metastasis_DE_metabolic_targets.html")
    if os.path.exists(dup_html):
        print(f"Removing duplicate HTML: {dup_html}")
        os.remove(dup_html)

# 2. Move files from output/ to respective folders and remove prefixes
for f in os.listdir(output_dir):
    if not os.path.isfile(os.path.join(output_dir, f)):
        continue
        
    if f.startswith("cancer_"):
        # Example: cancer_breast-cancer_breast_100k...
        parts = f.split("_", 2)
        if len(parts) >= 3:
            cancer_slug = parts[1] # e.g. breast-cancer
            rest = parts[2]
            
            # Map cancer_slug to short name for results folder
            short_name = cancer_slug.split("-")[0] # breast
            target_dir = os.path.join(output_dir, f"{short_name}_results")
            
            os.makedirs(target_dir, exist_ok=True)
            
            src = os.path.join(output_dir, f)
            dst = os.path.join(target_dir, rest)
            print(f"Moving {src} -> {dst}")
            shutil.move(src, dst)

# 3. Remove old breast_cancer dir if it's empty or just has checkpoints
breast_cancer = os.path.join(output_dir, "breast_cancer")
if os.path.exists(breast_cancer):
    print(f"Removing old directory: {breast_cancer}")
    shutil.rmtree(breast_cancer)

print("Cleanup complete.")
