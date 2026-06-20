import pandas as pd
import os
import glob
import sys

workspace_dir = '/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB'
output_dir = os.path.join(workspace_dir, 'output')

html_files = glob.glob(os.path.join(output_dir, '*.html'))

out_file = os.path.join(output_dir, 'parsed_tables.txt')

with open(out_file, 'w') as out_f:
    for f in html_files:
        if not os.path.exists(f):
            continue
            
        out_f.write(f"========== TABLES FOR {os.path.basename(f)} ==========\n")
        try:
            # Match all tables using pandas
            tables = pd.read_html(f)
            if not tables:
                out_f.write("No tables found.\n\n")
                continue
            for i, tbl in enumerate(tables):
                out_f.write(f"--- Table {i+1} ---\n")
                out_f.write(tbl.to_string(index=False) + "\n\n")
        except Exception as e:
            out_f.write(f"Error parsing tables: {e}\n\n")

print(f"Done parsing tables to {out_file}.")
