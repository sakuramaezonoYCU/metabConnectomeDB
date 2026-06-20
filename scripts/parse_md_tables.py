import os
import glob
from html.parser import HTMLParser

class TableExtractor(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_table = False
        self.in_tr = False
        self.in_td_or_th = False
        self.current_cell = []
        self.current_row = []
        self.tables = []

    def handle_starttag(self, tag, attrs):
        if tag == 'table':
            self.in_table = True
            self.tables.append([])
        elif tag == 'tr' and self.in_table:
            self.in_tr = True
            self.current_row = []
        elif tag in ['td', 'th'] and self.in_tr:
            self.in_td_or_th = True
            self.current_cell = []

    def handle_endtag(self, tag):
        if tag == 'table':
            self.in_table = False
        elif tag == 'tr' and self.in_table:
            self.in_tr = False
            self.tables[-1].append(self.current_row)
        elif tag in ['td', 'th'] and self.in_tr:
            self.in_td_or_th = False
            self.current_row.append(' '.join(self.current_cell).strip().replace('\n', ' '))

    def handle_data(self, data):
        if self.in_td_or_th:
            stripped = data.strip()
            if stripped:
                self.current_cell.append(stripped)

def format_markdown_table(table):
    if not table: return ""
    md = []
    # Header
    md.append("| " + " | ".join(table[0]) + " |")
    # Separator
    md.append("|" + "|".join(["---"] * len(table[0])) + "|")
    # Rows
    for row in table[1:]:
        # pad row if shorter than header
        padded_row = row + [""] * (len(table[0]) - len(row))
        md.append("| " + " | ".join(padded_row) + " |")
    return "\n".join(md)

workspace_dir = '/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB'
output_dir = os.path.join(workspace_dir, 'output')
html_files = glob.glob(os.path.join(output_dir, '*.html'))

out_file = os.path.join(output_dir, 'parsed_md_tables.txt')

with open(out_file, 'w', encoding='utf-8') as out_f:
    for f in html_files:
        if not os.path.exists(f):
            continue
            
        out_f.write(f"========== TABLES FOR {os.path.basename(f)} ==========\n")
        try:
            with open(f, 'r', encoding='utf-8') as inf:
                content = inf.read()
            extractor = TableExtractor()
            extractor.feed(content)
            
            for i, tbl in enumerate(extractor.tables):
                if len(tbl) > 1: # Ignore empty tables or single row tables which might be layout
                    out_f.write(f"--- Table {i+1} ---\n")
                    out_f.write(format_markdown_table(tbl))
                    out_f.write("\n\n")
        except Exception as e:
            out_f.write(f"Error parsing: {e}\n\n")

print(f"Done parsing tables to {out_file}.")
