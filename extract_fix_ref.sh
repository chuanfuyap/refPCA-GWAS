#!/bin/bash

## extract intersect between query genotype with reference pruned genotype and update reference allele

# EDIT THIS TO THE CORRESPONDING PATHS
PLINK="${LOCATION TO PLINK}"
DATA="${PLINK FILE LOCATION}"
VARIANTS="./extract_gnomad_overlap.txt"
OUTPUT="${OUTPUT DIRECTORY LOCATION AND NAME}"
REFERENCE="./gnomad_ref_allele.txt"

$PLINK \
  --bfile $DATA \
  --extract $VARIANTS \
  --ref-allele $REFERENCE \
  --make-bed \
  --threads ${OPTIONAL THREADS} \
  --out $OUTPUT