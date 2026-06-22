import urllib.request
import urllib.parse
import json

queries = [
    "IDO1 AND kynurenine AND AHR AND tumor",
    "SLC7A11 AND ferroptosis AND cancer",
]

for q in queries:
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={urllib.parse.quote(q)}&retmode=json&retmax=3"
    with urllib.request.urlopen(url) as response:
        data = json.loads(response.read().decode())
    idlist = data['esearchresult']['idlist']
    
    if not idlist:
        print(f"No results for {q}")
        continue
        
    summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={','.join(idlist)}&retmode=json"
    with urllib.request.urlopen(summary_url) as summary_response:
        summary_data = json.loads(summary_response.read().decode())
        
    print(f"\n--- Results for: {q} ---")
    for pid in idlist:
        item = summary_data['result'][pid]
        title = item.get('title', 'No Title')
        first_author = item.get('authors', [{'name': 'Unknown'}])[0]['name']
        print(f"PMID: {pid} | {first_author} et al. | {title}")
