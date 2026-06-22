import os
import glob
import shutil

cancer_to_primary = {
    'lung cancer, lung adenocarcinoma': ['lung'],
    'breast cancer': ['breast', 'mammary gland'],
    'colorectal cancer, colorectal carcinoma || metastatic malignant neoplasm': ['colon', 'large intestine'],
    'melanoma, metastatic melanoma': ['skin of body'],
    'ovarian cancer, malignant ovarian serous tumor': ['ovary']
}

output_dir = "output"

for disease_name in cancer_to_primary.keys():
    cancer_name_safe = str(disease_name).replace(" ", "_").replace(",", "").replace("/", "_") + "_results"
    target_dir = os.path.join(output_dir, cancer_name_safe)
    os.makedirs(target_dir, exist_ok=True)
    
    # Match by the primary tissue prefixes
    primary_tissues = cancer_to_primary[disease_name]
    tissue_slug = "_".join(t.replace(" ", "-") for t in primary_tissues)
    
    # Find h5ad files in output/ that start with tissue_slug
    for h5ad_file in glob.glob(os.path.join(output_dir, f"{tissue_slug}*.h5ad")):
        target_path = os.path.join(target_dir, os.path.basename(h5ad_file))
        print(f"Moving {os.path.basename(h5ad_file)} -> {cancer_name_safe}/")
        shutil.move(h5ad_file, target_path)

print("Done moving h5ad files.")
