#!/bin/bash
###############################################################################
# Master Script: Complete De Novo Variant Analysis for Autism Trio
#
# This script orchestrates the complete pipeline:
#   1. Download FASTQ files from ENA (if needed)
#   2. Preprocessing (alignment, deduplication, BQSR)
#   3. Variant calling and de novo detection
#
# Reference: Based on methodology from autism exome sequencing studies
# Coverage requirement: 20x minimum in all trio members
#
# Usage:
#   bash run_denovo_analysis.sh [--skip-download] [--skip-preprocess]
###############################################################################

set -euo pipefail

# Parse arguments
SKIP_DOWNLOAD=false
SKIP_PREPROCESS=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-download)
            SKIP_DOWNLOAD=true
            shift
            ;;
        --skip-preprocess)
            SKIP_PREPROCESS=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--skip-download] [--skip-preprocess]"
            exit 1
            ;;
    esac
done

echo "=========================================="
echo "De Novo Variant Analysis Pipeline"
echo "=========================================="
echo "Date: $(date)"
echo ""
echo "This pipeline will:"
echo "  1. Download FASTQ files from ENA"
echo "  2. Align reads to hg19 reference"
echo "  3. Perform quality control and preprocessing"
echo "  4. Call variants using GATK"
echo "  5. Identify de novo variants"
echo "  6. Generate analysis reports"
echo ""

###############################################################################
# STEP 1: Download FASTQ files from ENA
###############################################################################

if [ "${SKIP_DOWNLOAD}" = false ]; then
    echo "=========================================="
    echo "STEP 1: Downloading FASTQ files from ENA"
    echo "=========================================="
    echo ""

    bash /home/johan/pipeline/varcall/convert_sra_fixed.sh

    echo ""
    echo "✓ FASTQ download complete"
    echo ""
else
    echo "Skipping FASTQ download (--skip-download specified)"
    echo ""
fi

###############################################################################
# STEP 2: Preprocessing (Alignment, Deduplication, BQSR)
###############################################################################

if [ "${SKIP_PREPROCESS}" = false ]; then
    echo "=========================================="
    echo "STEP 2: Preprocessing"
    echo "=========================================="
    echo ""
    echo "This step will:"
    echo "  - Align reads with BWA-MEM"
    echo "  - Mark duplicates with GATK"
    echo "  - Perform base quality recalibration (BQSR)"
    echo ""
    echo "Expected time: 4-6 hours (with 60 threads)"
    echo ""

    read -p "Start preprocessing? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        bash /home/johan/pipeline/varcall/preprocess_trio_optimized.sh
        echo ""
        echo "✓ Preprocessing complete"
        echo ""
    else
        echo "Preprocessing skipped by user"
        exit 0
    fi
else
    echo "Skipping preprocessing (--skip-preprocess specified)"
    echo ""
fi

###############################################################################
# STEP 3: Variant Calling and De Novo Detection
###############################################################################

echo "=========================================="
echo "STEP 3: Variant Calling & De Novo Detection"
echo "=========================================="
echo ""
echo "This step will:"
echo "  - Call variants with GATK HaplotypeCaller"
echo "  - Perform joint genotyping"
echo "  - Apply quality filters"
echo "  - Identify de novo variants"
echo "  - Annotate variants"
echo "  - Generate analysis reports"
echo ""
echo "Expected time: 2-4 hours"
echo ""

read -p "Start variant calling and de novo detection? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    bash /home/johan/pipeline/varcall/complete_trio_pipeline.sh
    echo ""
    echo "✓ Variant calling and de novo detection complete"
    echo ""
else
    echo "Variant calling skipped by user"
    exit 0
fi

###############################################################################
# Final Summary
###############################################################################

echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo ""
echo "Results are in: ${HOME}/output/autism/variants/"
echo ""
echo "Key files:"
echo "  - De novo variants: denovo_high_confidence.vcf.gz"
echo "  - Analysis report: denovo_analysis_report.txt"
echo "  - Variant table: denovo_variants_table.tsv"
echo ""
echo "Next steps:"
echo "  1. Review the analysis report"
echo "  2. Prioritize variants by gene function"
echo "  3. Validate high-priority variants"
echo ""
echo "Completed: $(date)"
