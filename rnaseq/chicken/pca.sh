#!/bin/bash
set -euo pipefail

Rscript pca.R \
  --combined=/home/johan/johan/output/chicken/combined.HTseq_report \
  --metadata=/home/johan/pipeline/rnaseq/chicken/sample.csv \
  --output_base=/home/johan/johan/output/chicken \
  --min_count=10 \
  --color_by=tissue --shape_by=egg_production --label_by=sample --ellipse_by=tissue
