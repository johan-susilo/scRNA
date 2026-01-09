#!/bin/bash
###############################################################################
# Fix Missing Reference Index Files for GATK
#
# This script creates the missing index files required for BQSR and variant
# calling with GATK.
#
# Missing files:
#   - hg19.dict (sequence dictionary)
#   - dbsnp_138.hg19.vcf.gz.tbi (VCF index)
#   - Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz.tbi (VCF index)
###############################################################################

set -euo pipefail

REF_DIR="/home/johan/johan/reference/human"
REFERENCE="${REF_DIR}/hg19.fa"
DBSNP="${REF_DIR}/dbsnp_138.hg19.vcf.gz"
MILLS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"

echo "=========================================="
echo "Fixing Reference Index Files for GATK"
echo "=========================================="
echo "Date: $(date)"
echo ""

# Activate conda environment
echo "Activating varcall environment..."
source "${HOME}/tool/anaconda3/etc/profile.d/conda.sh"
conda activate varcall
echo ""

###############################################################################
# 1. Create Sequence Dictionary
###############################################################################

echo "=== Step 1: Creating sequence dictionary ==="
echo "File: ${REF_DIR}/hg19.dict"
echo ""

if [ ! -f "${REF_DIR}/hg19.dict" ]; then
    echo "Creating hg19.dict..."

    gatk CreateSequenceDictionary \
        -R "${REFERENCE}" \
        -O "${REF_DIR}/hg19.dict"

    if [ -f "${REF_DIR}/hg19.dict" ]; then
        echo "✓ hg19.dict created successfully"
        ls -lh "${REF_DIR}/hg19.dict"
    else
        echo "✗ Failed to create hg19.dict"
        exit 1
    fi
else
    echo "✓ hg19.dict already exists"
    ls -lh "${REF_DIR}/hg19.dict"
fi

echo ""

###############################################################################
# 2. Index dbSNP VCF
###############################################################################

echo "=== Step 2: Indexing dbSNP VCF ==="
echo "File: ${DBSNP}"
echo ""

if [ ! -f "${DBSNP}.tbi" ]; then
    echo "Checking compression format..."

    # Check if file is bgzipped (required for tabix)
    if ! bgzip -t "${DBSNP}" 2>/dev/null; then
        echo "File is not bgzipped. Recompressing with bgzip..."
        echo "This may take several minutes for large VCF files..."

        # Decompress and recompress with bgzip
        gunzip -c "${DBSNP}" | bgzip -c > "${DBSNP}.tmp"
        mv "${DBSNP}.tmp" "${DBSNP}"

        echo "✓ Recompression complete"
    else
        echo "✓ File is already bgzipped"
    fi

    echo "Creating index for dbSNP..."
    tabix -p vcf "${DBSNP}"

    if [ -f "${DBSNP}.tbi" ]; then
        echo "✓ dbSNP index created successfully"
        ls -lh "${DBSNP}.tbi"
    else
        echo "✗ Failed to create dbSNP index"
        exit 1
    fi
else
    echo "✓ dbSNP index already exists"
    ls -lh "${DBSNP}.tbi"
fi

echo ""

###############################################################################
# 3. Index Mills Indels VCF
###############################################################################

echo "=== Step 3: Indexing Mills indels VCF ==="
echo "File: ${MILLS}"
echo ""

if [ ! -f "${MILLS}.tbi" ]; then
    echo "Checking compression format..."

    # Check if file is bgzipped (required for tabix)
    if ! bgzip -t "${MILLS}" 2>/dev/null; then
        echo "File is not bgzipped. Recompressing with bgzip..."

        # Decompress and recompress with bgzip
        gunzip -c "${MILLS}" | bgzip -c > "${MILLS}.tmp"
        mv "${MILLS}.tmp" "${MILLS}"

        echo "✓ Recompression complete"
    else
        echo "✓ File is already bgzipped"
    fi

    echo "Creating index for Mills indels..."
    tabix -p vcf "${MILLS}"

    if [ -f "${MILLS}.tbi" ]; then
        echo "✓ Mills indels index created successfully"
        ls -lh "${MILLS}.tbi"
    else
        echo "✗ Failed to create Mills indels index"
        exit 1
    fi
else
    echo "✓ Mills indels index already exists"
    ls -lh "${MILLS}.tbi"
fi

echo ""

###############################################################################
# Verification
###############################################################################

echo "=========================================="
echo "Verification"
echo "=========================================="
echo ""

ALL_PRESENT=true

# Check sequence dictionary
if [ -f "${REF_DIR}/hg19.dict" ]; then
    echo "✓ hg19.dict"
else
    echo "✗ hg19.dict - MISSING"
    ALL_PRESENT=false
fi

# Check dbSNP index
if [ -f "${DBSNP}.tbi" ]; then
    echo "✓ dbsnp_138.hg19.vcf.gz.tbi"
else
    echo "✗ dbsnp_138.hg19.vcf.gz.tbi - MISSING"
    ALL_PRESENT=false
fi

# Check Mills index
if [ -f "${MILLS}.tbi" ]; then
    echo "✓ Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz.tbi"
else
    echo "✗ Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz.tbi - MISSING"
    ALL_PRESENT=false
fi

echo ""

if [ "${ALL_PRESENT}" = true ]; then
    echo "=========================================="
    echo "SUCCESS! All reference indices created"
    echo "=========================================="
    echo ""
    echo "You can now run the variant calling pipeline:"
    echo ""
    echo "  bash /home/johan/pipeline/varcall/complete_autism_analysis.sh \\"
    echo "      --skip-download \\"
    echo "      --skip-preprocess"
    echo ""
    echo "This will:"
    echo "  1. Skip to BQSR (using existing dedup BAMs)"
    echo "  2. Generate final BAMs"
    echo "  3. Call variants with GATK"
    echo "  4. Detect de novo variants"
    echo "  5. Detect rare inherited variants"
    echo "  6. Generate comprehensive reports"
    echo ""
    echo "Expected runtime: 4-6 hours"
    echo ""
else
    echo "=========================================="
    echo "ERROR: Some indices failed to create"
    echo "=========================================="
    echo ""
    echo "Please check the error messages above"
    exit 1
fi

echo "Completed: $(date)"
