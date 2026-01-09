#!/bin/bash
###############################################################################
# FASTQ Download Script - ENA Method (RECOMMENDED)
# Downloads pre-converted FASTQ files from European Nucleotide Archive
# This avoids SRA Toolkit reference dependency issues entirely
#
# IMPORTANT: For aligned SRA files, use this script instead of SRA Toolkit
# See README_ENA_DOWNLOAD.md for details
###############################################################################

set -euo pipefail

OUTPUT_DIR="${HOME}/output/autism"

# Trio sample IDs
PROBAND="SRR8697636"
MOTHER="SRR8697627"
FATHER="SRR8697645"

# ENA FTP base URL
ENA_FTP="ftp://ftp.sra.ebi.ac.uk/vol1/fastq"

echo "=========================================="
echo "Downloading FASTQ from ENA"
echo "=========================================="
echo "Date: $(date)"
echo "Output: ${OUTPUT_DIR}"
echo ""
echo "NOTE: This script downloads from ENA instead of converting from SRA"
echo "      ENA provides pre-converted FASTQ files, avoiding reference issues"
echo ""

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Function to download FASTQ from ENA
download_from_ena() {
    local SRA=$1
    local LABEL=$2

    # Construct ENA path: /vol1/fastq/SRR###/00#/SRR#######/
    local DIR1="${SRA:0:6}"  # SRR869
    local DIR2="${SRA: -1}"   # Last digit
    # Pad to 3 digits: 6 -> 006, 7 -> 007, 5 -> 005
    local DIR2_PADDED=$(printf "%03d" ${DIR2})
    local ENA_PATH="${ENA_FTP}/${DIR1}/${DIR2_PADDED}/${SRA}"

    # Check if files already exist and verify
    if [ -f "${OUTPUT_DIR}/${LABEL}_1.fastq.gz" ] && [ -f "${OUTPUT_DIR}/${LABEL}_2.fastq.gz" ]; then
        echo "✓ ${LABEL} FASTQ already exists, verifying..."

        # Quick validation - check if files are readable
        if zcat "${LABEL}_1.fastq.gz" | head -4 > /dev/null 2>&1; then
            echo "  Files validated successfully"
            ls -lh "${LABEL}"_*.fastq.gz
            return 0
        else
            echo "  Files corrupted, re-downloading..."
            rm -f "${LABEL}"_*.fastq.gz
        fi
    fi

    echo "Downloading ${SRA} (${LABEL}) from ENA..."
    echo "  URL: ${ENA_PATH}"

    # Download R1 with resume support
    echo "  Downloading Read 1..."
    wget -c "${ENA_PATH}/${SRA}_1.fastq.gz" -O "${LABEL}_1.fastq.gz.tmp"
    mv "${LABEL}_1.fastq.gz.tmp" "${LABEL}_1.fastq.gz"

    # Download R2 with resume support
    echo "  Downloading Read 2..."
    wget -c "${ENA_PATH}/${SRA}_2.fastq.gz" -O "${LABEL}_2.fastq.gz.tmp"
    mv "${LABEL}_2.fastq.gz.tmp" "${LABEL}_2.fastq.gz"

    echo "✓ ${LABEL} download complete"
    ls -lh "${OUTPUT_DIR}/${LABEL}"_*.fastq.gz
}

# Download all three samples
download_from_ena "${PROBAND}" "proband"
download_from_ena "${MOTHER}" "mother"
download_from_ena "${FATHER}" "father"

echo ""
echo "=========================================="
echo "Download Complete!"
echo "=========================================="
echo ""
ls -lh "${OUTPUT_DIR}"/*.fastq.gz
echo ""

# Verify read counts
echo "=========================================="
echo "Verifying Read Counts"
echo "=========================================="
for sample in proband mother father; do
    if [ -f "${sample}_1.fastq.gz" ]; then
        echo -n "${sample}: "
        zcat ${sample}_1.fastq.gz | wc -l | awk '{printf "%'"'"'d read pairs\n", $1/4}'
    fi
done
echo ""

echo "Next: Run the preprocessing/alignment pipeline"
echo "Note: These files are ready for alignment with BWA-MEM or your chosen aligner"
