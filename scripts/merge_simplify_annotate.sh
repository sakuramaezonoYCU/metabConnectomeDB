#!/bin/bash

cd "$(dirname "$0")"

python merge_dbs_claude.py
python generate_final_outputs.py
python annotate_with_hmdb.py
