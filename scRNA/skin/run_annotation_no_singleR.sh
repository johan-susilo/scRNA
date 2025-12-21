#!/bin/bash
# Cell Annotation Without SingleR (due to BiocFileCache compatibility issue)
# This script uses CelliD + scCATCH + Classical Markers for annotation

set -e

cd /home/johan/pipeline/scRNA/skin

echo "================================================"
echo "Cell Annotation Pipeline (Without SingleR)"
echo "================================================"
echo ""
echo "Methods to be used:"
echo "  1. CelliD (PanglaoDB signatures)"
echo "  2. scCATCH (skin-specific markers)"
echo "  3. Classical marker plots"
echo "  4. Consensus from methods 1-2"
echo ""

RDS_FILE="/home/johan/output/skin_pmh/TN.combined_dim30.rds"
OUTPUT_DIR="/home/johan/output/skin_pmh/annotations"

# Step 1: Read RDS
echo "Step 1/5: Reading RDS file..."
Rscript cell_annotation.R \
  -r "${RDS_FILE}" \
  -s read_rds \
  -o "${OUTPUT_DIR}"

# Step 2: CelliD
echo ""
echo "Step 2/5: Running CelliD annotation..."
Rscript cell_annotation.R \
  -s celliD \
  -o "${OUTPUT_DIR}"

# Step 3: scCATCH
echo ""
echo "Step 3/5: Running scCATCH annotation..."
Rscript cell_annotation.R \
  -s scCATCH \
  -o "${OUTPUT_DIR}" \
  --tissue skin

# Step 4: Classical Markers
echo ""
echo "Step 4/5: Generating classical marker plots..."
Rscript cell_annotation.R \
  -s markers \
  -o "${OUTPUT_DIR}"

# Step 5: Consensus
echo ""
echo "Step 5/5: Generating consensus annotation..."
Rscript cell_annotation.R \
  -s consensus \
  -o "${OUTPUT_DIR}"

echo ""
echo "================================================"
echo "✓ Annotation Complete!"
echo "================================================"
echo ""
echo "Results:"
echo "  - CelliD: ${OUTPUT_DIR}/celliD/CelliD_PanglaoDB_summary.tsv"
echo "  - scCATCH: ${OUTPUT_DIR}/scCATCH/scCATCH_summary.tsv"
echo "  - Markers: ${OUTPUT_DIR}/markers/*.pdf"
echo "  - Consensus: ${OUTPUT_DIR}/consensus/consensus_annotation.tsv"
echo ""
echo "View consensus results:"
echo "  cat ${OUTPUT_DIR}/consensus/consensus_annotation.tsv"
echo ""
