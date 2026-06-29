import os

def main():
    print("=" * 60)
    print("KIDNEY RCC SPATIAL DATA DOWNLOAD")
    print("=" * 60)
    print("Please manually download the spatial dataset from Kalogirou 2025:")
    print("Zenodo DOI: 10.5281/zenodo.16833780")
    print("URL: https://zenodo.org/records/16833780")
    print()
    
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    target_dir = os.path.join(base_dir, 'input', 'spatial', 'Zenodo_16833780_RCC')
    
    print(f"Once downloaded, please extract the data files (e.g., .rds, .h5ad, or matrix files) into:")
    print(f" -> {target_dir}")
    print()
    
    os.makedirs(target_dir, exist_ok=True)
    
    print("If the data is provided as a Seurat object (.rds), you may need to use R to convert it to a scanpy-compatible .h5ad file.")
    print("Example R command:")
    print("  library(Seurat)")
    print("  library(SeuratDisk)")
    print("  seu <- readRDS('data.rds')")
    print("  SaveH5Seurat(seu, filename = 'data.h5Seurat')")
    print("  Convert('data.h5Seurat', dest = 'h5ad')")
    print("=" * 60)

if __name__ == "__main__":
    main()
