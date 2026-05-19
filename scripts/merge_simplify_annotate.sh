#!/bin/bash

cd "$(dirname "$0")"

../venv/bin/python merge_dbs_claude.py
../venv/bin/python generate_final_outputs.py
../venv/bin/python annotate_with_hmdb.py
../venv/bin/python generate_column_definitions.py
