#!/bin/bash
###############################################################################
# Master Script: Complete Autism Trio Pipeline with UCSC hg19
#
# Purpose: Orchestrate the complete pipeline from re-alignment to variant calling
#
# This script will:
#   1. Wait for BWA index to complete
#   2. Re-align all trio samples to UCSC hg19
#   3. Run variant calling and de novo detection
#
###############################################################################

set -euo pipefail

PIPELINE_DIR="/home/johan/pipeline/varcall"
REF_DIR="/home/johan/johan/johan/johan/reference/human"
REFERENCE="${REF_DIR}/hg19_ucsc.fa"

echo "=========================================="
echo "Complete Autism Trio Pipeline"
echo "=========================================="
echo "Date: $(date)"
echo ""

###############################################################################
# STEP 1: Verify BWA Index is Complete
###############################################################################

echo "=== STEP 1: Checking BWA index ==="
echo ""

while [ ! -f "${REFERENCE}.bwt" ]; do
    echo "Waiting for BWA index to complete..."
    echo "Current time: $(date)"
    sleep 60
done

echo "✓ BWA index complete"
echo ""

# Verify all index files exist
REQUIRED_FILES=(
    "${REFERENCE}"
    "${REFERENCE}.fai"
    "${REFERENCE}.bwt"
    "${REFERENCE}.pac"
    "${REFERENCE}.ann"
    "${REFERENCE}.amb"
    "${REFERENCE}.sa"
    "${REF_DIR}/hg19_ucsc.dict"
)

echo "Verifying all required files..."
for FILE in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "${FILE}" ]; then
        echo "ERROR: Required file not found: ${FILE}"
        exit 1
    fi
    echo "  ✓ $(basename ${FILE})"
done

echo ""
echo "✓ All reference files ready"
echo ""

###############################################################################
# STEP 2: Re-align Trio to UCSC hg19
###############################################################################

echo "=========================================="
echo "STEP 2: Re-aligning Trio to UCSC hg19"
echo "=========================================="
echo ""

bash "${PIPELINE_DIR}/realign_with_ucsc_hg19.sh"

if [ $? -ne 0 ]; then
    echo "ERROR: Re-alignment failed!"
    exit 1
fi

echo ""
echo "✓ Re-alignment complete"
echo ""

###############################################################################
# STEP 3: Variant Calling and De Novo Detection
###############################################################################

echo "=========================================="
echo "STEP 3: Variant Calling"
echo "=========================================="
echo ""

bash "${PIPELINE_DIR}/varcall_trio_ucsc.sh"

if [ $? -ne 0 ]; then
    echo "ERROR: Variant calling failed!"
    exit 1
fi

echo ""
echo "✓ Variant calling complete"
echo ""

###############################################################################
# Final Summary
###############################################################################

echo ""
echo "=========================================="
echo "PIPELINE COMPLETE!"
echo "=========================================="
echo ""
echo "All steps completed successfully:"
echo "  1. ✓ Reference preparation (UCSC hg19 with chr naming)"
echo "  2. ✓ Re-alignment of trio samples"
echo "  3. ✓ Variant calling and de novo detection"
echo ""
echo "Results location: /home/johan/output/autism/variants/"
echo ""
echo "Key output files:"
echo "  - De novo variants: denovo_high_confidence.vcf.gz"
echo "  - Summary report: denovo_analysis_report.txt"
echo "  - Variant table: denovo_variants_table.tsv"
echo ""
echo "Completed: $(date)"
echo ""
