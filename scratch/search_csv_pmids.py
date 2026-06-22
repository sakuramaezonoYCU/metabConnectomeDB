import csv

csv_file = 'input/pubmed_results.csv'
pathways = {
    'IDO1/Kynurenine': ['ido1', 'kynurenine', 'ahr'],
    'xCT/Glutamate': ['xct', 'slc7a11', 'glutamate', 'ferroptosis'],
    'CD73/Adenosine': ['cd73', 'nt5e', 'adenosine', 'cd39'],
    'COX-2/PGE2': ['cox-2', 'pge2', 'ptgs2'],
    'SPHK1/S1P': ['sphk1', 's1p', 'sphingosine']
}

results = {k: [] for k in pathways.keys()}

with open(csv_file, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        text = (row.get('Title', '') + " " + row.get('Abstract', '')).lower()
        pmid = row.get('PMID', '')
        title = row.get('Title', '')
        
        for pway, keywords in pathways.items():
            # Check if any keyword matches
            if any(kw in text for kw in keywords):
                results[pway].append((pmid, title))

for pway, hits in results.items():
    print(f"\n--- {pway} ---")
    for hit in hits[:5]:
        print(f"PMID: {hit[0]} | Title: {hit[1]}")
