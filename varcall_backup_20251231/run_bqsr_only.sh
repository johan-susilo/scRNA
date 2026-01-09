#!/bin/bash
###############################################################################
# Run BQSR Only - For resuming pipeline after dedup stage
###############################################################################

set -euo pipefail

OUTPUT_DIR="/home/johan/output/autism"
REF_DIR="/home/johan/johan/reference/human"
REFERENCE="${REF_DIR}/hg19.fa"
DBSNP="${REF_DIR}/dbsnp_138.hg19.vcf.gz"
MILLS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
THREADS=60

echo "=========================================="
echo "BQSR for Autism Trio"
echo "=========================================="
date
echo ""

# Activate conda environment
source "${HOME}/tool/anaconda3/etc/profile.d/conda.sh"
conda activate varcall

cd "${OUTPUT_DIR}"

###############################################################################
# Function to run BQSR on one sample
###############################################################################

run_bqsr() {
    local LABEL=$1
    local BAM_IN="${OUTPUT_DIR}/${LABEL}_dedup.bam"
    local RECAL_TABLE="${OUTPUT_DIR}/${LABEL}_recal_data.table"
    local BAM_OUT="${OUTPUT_DIR}/${LABEL}_final.bam"

    if [ -f "${BAM_OUT}" ] && [ -f "${BAM_OUT}.bai" ]; then
        echo "✓ ${LABEL} final BAM already exists, skipping"
        return 0
    fi

    echo ""
    echo "=== Processing ${LABEL} ==="
    echo ""

    # Step 1: BaseRecalibrator
    echo "Step 1/3: BaseRecalibrator..."
    gatk BaseRecalibrator \
        --input "${BAM_IN}" \
        --reference "${REFERENCE}" \
        --known-sites "${DBSNP}" \
        --known-sites "${MILLS}" \
        --output "${RECAL_TABLE}"

    echo "✓ Recalibration table created"

    # Step 2: ApplyBQSR
    echo "Step 2/3: ApplyBQSR..."
    gatk ApplyBQSR \
        --input "${BAM_IN}" \
        --reference "${REFERENCE}" \
        --bqsr-recal-file "${RECAL_TABLE}" \
        --output "${BAM_OUT}"

    echo "✓ BQSR applied"

    # Step 3: Index final BAM
    echo "Step 3/3: Indexing final BAM..."
    samtools index -@ ${THREADS} "${BAM_OUT}"

    echo "✓ ${LABEL} BQSR complete"
    ls -lh "${BAM_OUT}" "${BAM_OUT}.bai"
}

###############################################################################
# Run BQSR for all trio members
###############################################################################

run_bqsr "proband"
run_bqsr "mother"
run_bqsr "father"

echo ""
echo "=========================================="
echo "BQSR Complete for All Samples"
echo "=========================================="
date
echo ""
echo "Final BAM files:"
ls -lh "${OUTPUT_DIR}"/*_final.bam
echo ""
echo "Next step: Run variant calling"
echo "  bash /home/johan/pipeline/varcall/complete_autism_analysis.sh --skip-download --skip-preprocess"
echo ""
