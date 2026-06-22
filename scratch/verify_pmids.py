import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

pmids = ['29551271', '34322129', '36496662', '30872535', '26880461', '29758241', '35121582', '36776289', '22298596']
url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={','.join(pmids)}"

try:
    with urllib.request.urlopen(url) as response:
        xml_data = response.read()
        
    root = ET.fromstring(xml_data)
    for docsum in root.findall('DocSum'):
        pmid = docsum.find('Id').text
        title = docsum.find(".//Item[@Name='Title']").text
        first_author = docsum.find(".//Item[@Name='AuthorList']/Item")
        author = first_author.text if first_author is not None else "Unknown"
        print(f"PMID: {pmid} | Author: {author} | Title: {title}")
except Exception as e:
    print(f"Error querying PubMed: {e}")
