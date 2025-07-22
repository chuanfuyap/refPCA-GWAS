#!/bin/bash

set -e
# Activate your hail/python environment if needed
# source /path/to/your/env/bin/activate

# EDIT THIS TO THE CORRESPONDING PATHS
CHR=20
NUM_CORES=8
CHUNKSIZE=225
SUBCHUNKSIZE=4
PHENO=""
COVAR=""
RAWFILE=""
VARIANT_FILE=""
LOCO=""
OUTPUT=""
model_type="logistic"

python ./regintest.py \
  --pheno-path $PHENO \
  --covar-path $COVAR \
  --rawfile-path $RAWFILE \
  --variant-list-path $VARIANT_FILE \
  --loco-path $LOCO \
  --output-path $OUTPUT \
  --model-type $model_type \
  --chr $CHR \
  --num-cores $NUM_CORES \
  --chunksize $CHUNKSIZE \
  --subchunksize $SUBCHUNKSIZE
