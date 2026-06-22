import json

def main():
    nb_path = "scripts/pan_cancer_meta_analysis.ipynb"
    with open(nb_path, "r") as f:
        nb = json.load(f)
        
    for i, c in enumerate(nb["cells"]):
        if i < 15:
            continue
        print(f"CELL {i} ({c['cell_type']})")
        src = "".join(c.get("source", []))
        print(src[:200] + ("..." if len(src) > 200 else ""))
        print("---")

if __name__ == "__main__":
    main()
