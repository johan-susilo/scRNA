#!/bin/bash
###############################################################################
# Robust Trio Download Script
# Handles common SRA download issues with multiple fallback methods
###############################################################################

set -euo pipefail

# Configuration
OUTPUT_DIR="${HOME}/output/autism"
THREADS=60

# Trio sample IDs
PROBAND="SRR8697636"
MOTHER="SRR8697627"
FATHER="SRR8697645"

echo "=========================================="
echo "Robust Trio Data Download"
echo "=========================================="
echo "Date: $(date)"
echo ""

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

###############################################################################
# Method 1: Direct fasterq-dump (bypasses prefetch)
###############################################################################

download_direct() {
    local SRA_ID=$1
    local LABEL=$2

    echo ""
    echo "=== Downloading ${LABEL} (${SRA_ID}) ==="

    # Check if already downloaded
    if [ -f "${OUTPUT_DIR}/${LABEL}_1.fastq.gz" ] && [ -f "${OUTPUT_DIR}/${LABEL}_2.fastq.gz" ]; then
        echo "✓ ${LABEL} FASTQ files already exist, skipping"
        return 0
    fi

    echo "Method 1: Direct download with fasterq-dump..."

    # Try direct fasterq-dump (downloads and converts in one step)
    if fasterq-dump \
        --split-files \
        --threads ${THREADS} \
        --progress \
        --outdir "${OUTPUT_DIR}" \
        --temp "${OUTPUT_DIR}/tmp_${SRA_ID}" \
        ${SRA_ID} 2>&1 | tee "${OUTPUT_DIR}/${SRA_ID}_download.log"; then

        echo "✓ Direct download successful"

        # Rename files
        if [ -f "${OUTPUT_DIR}/${SRA_ID}_1.fastq" ]; then
            mv "${OUTPUT_DIR}/${SRA_ID}_1.fastq" "${OUTPUT_DIR}/${LABEL}_1.fastq"
            mv "${OUTPUT_DIR}/${SRA_ID}_2.fastq" "${OUTPUT_DIR}/${LABEL}_2.fastq"

            # Compress
            echo "Compressing with pigz..."
            pigz -p ${THREADS} "${OUTPUT_DIR}/${LABEL}_1.fastq" &
            pigz -p ${THREADS} "${OUTPUT_DIR}/${LABEL}_2.fastq" &
            wait

            echo "✓ ${LABEL} complete!"
            return 0
        fi
    fi

    echo "Method 1 failed, trying Method 2..."

    ###########################################################################
    # Method 2: Using SRA Toolkit with vdb-config fix
    ###########################################################################

    echo "Method 2: Configuring SRA Toolkit..."

    # Configure SRA Toolkit to use cloud access
    vdb-config --restore-defaults
    vdb-config --set /repository/user/main/public/root="${OUTPUT_DIR}/sra_cache"
    vdb-config --prefetch-to-cwd

    # Try prefetch with smaller cache
    mkdir -p "${OUTPUT_DIR}/${SRA_ID}"

    if prefetch \
        --output-directory "${OUTPUT_DIR}/${SRA_ID}" \
        --max-size 50G \
        ${SRA_ID} 2>&1 | tee -a "${OUTPUT_DIR}/${SRA_ID}_download.log"; then

        echo "✓ Prefetch successful"

        # Convert to FASTQ
        fasterq-dump \
            --split-files \
            --threads ${THREADS} \
            --outdir "${OUTPUT_DIR}" \
            "${OUTPUT_DIR}/${SRA_ID}/${SRA_ID}.sra"

        # Rename and compress
        mv "${OUTPUT_DIR}/${SRA_ID}_1.fastq" "${OUTPUT_DIR}/${LABEL}_1.fastq"
        mv "${OUTPUT_DIR}/${SRA_ID}_2.fastq" "${OUTPUT_DIR}/${LABEL}_2.fastq"

        pigz -p ${THREADS} "${OUTPUT_DIR}/${LABEL}_1.fastq" &
        pigz -p ${THREADS} "${OUTPUT_DIR}/${LABEL}_2.fastq" &
        wait

        echo "✓ ${LABEL} complete!"
        return 0
    fi

    echo "Method 2 failed, trying Method 3..."

    ###########################################################################
    # Method 3: Download from ENA (European Nucleotide Archive)
    ###########################################################################

    echo "Method 3: Downloading from ENA mirror..."

    # ENA FTP URLs (usually more reliable than NCBI)
    ENA_FTP="ftp://ftp.sra.ebi.ac.uk/vol1/fastq"

    # Parse SRA ID for ENA path
    # SRR8697636 -> SRR869/006/SRR8697636
    PREFIX="${SRA_ID:0:6}"
    LAST_DIGIT="${SRA_ID: -1}"

    if [ ${#SRA_ID} -eq 10 ]; then
        # 10-character IDs
        SUBDIR="00${LAST_DIGIT}"
    else
        # 9-character IDs
        SUBDIR="00"
    fi

    ENA_BASE="${ENA_FTP}/${PREFIX}/${SUBDIR}/${SRA_ID}"

    echo "Trying ENA URL: ${ENA_BASE}"

    # Download R1
    if wget -c -O "${OUTPUT_DIR}/${LABEL}_1.fastq.gz" \
        "${ENA_BASE}/${SRA_ID}_1.fastq.gz" 2>&1 | tee -a "${OUTPUT_DIR}/${SRA_ID}_download.log"; then

        # Download R2
        wget -c -O "${OUTPUT_DIR}/${LABEL}_2.fastq.gz" \
            "${ENA_BASE}/${SRA_ID}_2.fastq.gz"

        echo "✓ ENA download successful!"
        return 0
    fi

    echo "Method 3 failed, trying Method 4..."

    ###########################################################################
    # Method 4: AWS S3 (if available)
    ###########################################################################

    echo "Method 4: Trying AWS S3..."

    if command -v aws &> /dev/null; then
        # SRA data is available in AWS Open Data
        # s3://sra-pub-run-odp/sra/SRR8697636/SRR8697636

        aws s3 cp --no-sign-request \
            "s3://sra-pub-run-odp/sra/${SRA_ID}/${SRA_ID}" \
            "${OUTPUT_DIR}/${SRA_ID}/${SRA_ID}.sra"

        if [ -f "${OUTPUT_DIR}/${SRA_ID}/${SRA_ID}.sra" ]; then
            # Convert to FASTQ
            fasterq-dump \
                --split-files \
                --threads ${THREADS} \
                --outdir "${OUTPUT_DIR}" \
                "${OUTPUT_DIR}/${SRA_ID}/${SRA_ID}.sra"

            mv "${OUTPUT_DIR}/${SRA_ID}_1.fastq" "${OUTPUT_DIR}/${LABEL}_1.fastq"
            mv "${OUTPUT_DIR}/${SRA_ID}_2.fastq" "${OUTPUT_DIR}/${LABEL}_2.fastq"

            pigz -p ${THREADS} "${OUTPUT_DIR}/${LABEL}_1.fastq" &
            pigz -p ${THREADS} "${OUTPUT_DIR}/${LABEL}_2.fastq" &
            wait

            echo "✓ AWS download successful!"
            return 0
        fi
    else
        echo "AWS CLI not available, skipping Method 4"
    fi

    ###########################################################################
    # All methods failed
    ###########################################################################

    echo ""
    echo "❌ ERROR: All download methods failed for ${SRA_ID}"
    echo "Check the log: ${OUTPUT_DIR}/${SRA_ID}_download.log"
    return 1
}

###############################################################################
# Download all three samples
###############################################################################

SUCCESS_COUNT=0

download_direct "${PROBAND}" "proband" && ((SUCCESS_COUNT++)) || echo "⚠ Proband download had issues"
download_direct "${MOTHER}" "mother" && ((SUCCESS_COUNT++)) || echo "⚠ Mother download had issues"
download_direct "${FATHER}" "father" && ((SUCCESS_COUNT++)) || echo "⚠ Father download had issues"

echo ""
echo "=========================================="
echo "Download Summary"
echo "=========================================="
echo "Successfully downloaded: ${SUCCESS_COUNT}/3 samples"
echo ""

if [ ${SUCCESS_COUNT} -eq 3 ]; then
    echo "✓ All trio samples downloaded successfully!"
    echo ""
    echo "FASTQ files:"
    ls -lh "${OUTPUT_DIR}"/*.fastq.gz
    echo ""
    echo "Next step: Run preprocessing"
    echo "  ./preprocess_trio_optimized.sh"
elif [ ${SUCCESS_COUNT} -gt 0 ]; then
    echo "⚠ Partial success - some samples downloaded"
    echo "Review logs and retry failed samples"
else
    echo "❌ All downloads failed!"
    echo ""
    echo "Troubleshooting steps:"
    echo "1. Check internet connection"
    echo "2. Try manual download from NCBI SRA Run Selector:"
    echo "   https://www.ncbi.nlm.nih.gov/Traces/study/?acc=${PROBAND}"
    echo "3. Contact system admin for network/firewall issues"
    echo "4. Use ENA Browser as alternative:"
    echo "   https://www.ebi.ac.uk/ena/browser/view/${PROBAND}"
fi

echo ""
echo "Completed: $(date)"
