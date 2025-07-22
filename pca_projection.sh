#!/bin/bash
# Usage: bash pca_projection.sh 

# Exit immediately if a command exits with a non-zero status.
set -e

# EDIT THIS TO THE CORRESPONDING PATHS
PLINK_PREFIX="$1"         # e.g. ukbb_genotype (without .bed/.bim/.fam)
GNOMAD_LOADINGS_HT="$2"   # e.g. /path/to/gnomad_loadings.ht
OUTPUT_PREFIX="$3"        # e.g. output/ukbb_pca

# Activate your hail/python environment if needed
# source /path/to/your/env/bin/activate

echo "Running UKBB PCA projection..."
echo "PLINK prefix: $PLINK_PREFIX"
echo "gnomAD loadings HT: $GNOMAD_LOADINGS_HT"
echo "Output prefix: $OUTPUT_PREFIX"

python project_gnomad.py \
    --input "$PLINK_PREFIX" \
    --reference "$GNOMAD_LOADINGS_HT" \
    --output "$OUTPUT_PREFIX"

echo "Projection completed. Result: ${OUTPUT_PREFIX}_pca.tsv"