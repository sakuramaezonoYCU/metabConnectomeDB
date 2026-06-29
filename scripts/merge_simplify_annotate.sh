#!/bin/bash

cd "$(dirname "$0")"

python fetch_kegg_pathways.py
python merge_dbs.py
python generate_final_outputs.py
python annotate_with_hmdb.py
python annotate_with_databases.py
