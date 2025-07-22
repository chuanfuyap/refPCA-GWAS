#!/bin/bash
# Usage:
# bash run_meta_regintest.sh 

# NOTE: BEFORE RUNNNG THIS PLEASE MAKE SURE ALL THE COLUMNS ARE THE SAME FOR ALL SUMSTATS, NAMELY IT HAS 'VARIANT' COLUMN AND ALL THE VARIANT IDs MATCH ACROSS THE SUMSTATS FILES. 
set -e

# EDIT THIS TO THE CORRESPONDING PATHS
FILELIST="./sumstats.txt"       # Required: file containing list of summary stat files, one per line
OUTPUT="test_meta"         # Optional: output file name (default: meta_results.tsv)

# Activate your python environment if needed
# source /path/to/your/env/bin/activate


echo "Running meta-analysis for regintest summary statistics..."
echo "File list: $FILELIST"
echo "Output file: $OUTPUT"

python synregslopes.py \
    --input-filelist "$FILELIST" \
    --output-file "$OUTPUT"

echo "Meta-analysis complete: $OUTPUT"