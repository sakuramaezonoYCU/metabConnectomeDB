import pandas as pd
import cellxgene_census
import requests
import sys

# Parameters
CENSUS_VERSION = "2023-12-15"
ORGANISM = "Homo sapiens"
TARGET_PAIRS_PATH = "output/human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv"

print("1. Loading target genes from target pairs CSV...")
try:
    target_df = pd.read_csv(TARGET_PAIRS_PATH)
    target_genes = sorted(target_df['Target'].dropna().unique().tolist())
    target_500 = target_genes[:500]
    print(f"Total target genes in DB: {len(target_genes)}")
    print(f"Checking the first 500 genes (which are queried in Section 4)...")
except Exception as e:
    print(f"Error loading target pairs CSV: {e}")
    sys.exit(1)

print("\n2. Connecting to CellxGene Census to get all human gene feature names...")
try:
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        var_df = cellxgene_census.get_var(
            census,
            organism=ORGANISM,
            column_names=["soma_joinid", "feature_name"]
        )
    census_genes = set(var_df["feature_name"].tolist())
    print(f"Loaded {len(census_genes):,} unique genes from CellxGene Census.")
except Exception as e:
    print(f"Error connecting to CellxGene Census: {e}")
    sys.exit(1)

# Find unmapped genes
mapped_500 = [g for g in target_500 if g in census_genes]
unmapped_500 = [g for g in target_500 if g not in census_genes]

print(f"\nResults for first 500 genes:")
print(f"→ Mapped: {len(mapped_500)} / 500")
print(f"→ Unmapped: {len(unmapped_500)} / 500")

if len(unmapped_500) == 0:
    print("All 500 genes are already mapped!")
    sys.exit(0)

print(f"\n3. Querying MyGene.info to find current approved symbols for the {len(unmapped_500)} unmapped genes...")
try:
    # Batch query MyGene.info
    url = "https://mygene.info/v3/query"
    payload = {
        "q": ",".join(unmapped_500),
        "scopes": "symbol,alias",
        "fields": "symbol,name",
        "species": "human",
        "size": 1
    }
    response = requests.post(url, data=payload)
    response.raise_for_status()
    results = response.json()
    
    # Parse results
    mapping_dict = {}
    for item in results:
        query_term = item.get("query")
        if not query_term:
            continue
        
        if "notfound" in item:
            mapping_dict[query_term] = {
                "approved_symbol": None,
                "name": None,
                "status": "Not Found on MyGene.info"
            }
        else:
            approved = item.get("symbol")
            name = item.get("name")
            mapping_dict[query_term] = {
                "approved_symbol": approved,
                "name": name,
                "status": "Found"
            }
            
except Exception as e:
    print(f"Error querying MyGene.info: {e}")
    sys.exit(1)

# Analyze potential replacements
report_data = []
for gene in unmapped_500:
    mapping = mapping_dict.get(gene, {"approved_symbol": None, "name": None, "status": "Not Found"})
    approved = mapping["approved_symbol"]
    
    in_census = "No"
    action = "Keep (External/Non-coding or custom name)"
    
    if approved:
        if approved in census_genes:
            in_census = "Yes (Found!)"
            action = f"Rename '{gene}' → '{approved}'"
        elif approved == gene:
            in_census = "No"
            action = "Keep (Valid gene name, but not in Census scRNAseq expression)"
        else:
            in_census = "No"
            action = f"Approved name is '{approved}' but also missing in Census"
    else:
        action = "Unrecognized symbol; verify source"
        
    report_data.append({
        "Original_Symbol": gene,
        "Approved_Symbol": approved or "N/A",
        "In_CellxGene": in_census,
        "Action": action
    })

report_df = pd.DataFrame(report_data)

# Print a nice markdown table of results
print("\n### Unmapped Genes Analysis Table:\n")
print(report_df.to_markdown(index=False))

# Also calculate stats
can_be_mapped = report_df[report_df["In_CellxGene"] == "Yes (Found!)"]
print(f"\nStats:")
print(f"→ Unmapped genes: {len(unmapped_500)}")
print(f"→ Can be successfully mapped by updating to latest HGNC approved symbol: {len(can_be_mapped)}")
if len(can_be_mapped) > 0:
    print(f"→ Examples of successful mapping:")
    for idx, row in can_be_mapped.head(5).iterrows():
        print(f"  * {row['Original_Symbol']} → {row['Approved_Symbol']}")
        
# Save report to csv for later use
report_df.to_csv("output/unmapped_genes_census_analysis.csv", index=False)
print("\nSaved full mapping report to output/unmapped_genes_census_analysis.csv")
