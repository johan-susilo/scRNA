#!/bin/bash
###############################################################################
# Optimized Trio Preprocessing Script - NEW TRIO
# Uses corrected sample IDs:
#   - SRR8697683: Proband (autism)
#   - SRR8697706: Mother (unaffected)
#   - SRR8697650: Father (unaffected)
###############################################################################

date
set -euo pipefail

# Configuration
THREADS=60  # Use 60 out of 64 cores
OUTPUT_DIR="${HOME}/output/autism_new"
REF_DIR="${HOME}/johan/johan/reference/human"
REFERENCE="${REF_DIR}/hg19.fa"

# NEW Trio sample IDs
PROBAND="SRR8697683"
MOTHER="SRR8697706"
FATHER="SRR8697650"

echo "=========================================="
echo "Optimized Trio Preprocessing Pipeline"
echo "=========================================="
echo "Date: $(date)"
echo "CPUs available: $(nproc)"
echo "Using threads: ${THREADS}"
echo ""
echo "NEW Trio samples:"
echo "  Proband (autism): ${PROBAND}"
echo "  Mother:           ${MOTHER}"
echo "  Father:           ${FATHER}"
echo ""

# Create directories
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

###############################################################################
# STEP 1: Download SRA files (if not already downloaded)
###############################################################################

echo "=== STEP 1: Downloading SRA files ==="
for SRA in ${PROBAND} ${MOTHER} ${FATHER}; do
    if [ -f "${OUTPUT_DIR}/${SRA}/${SRA}.sra" ]; then
        echo "✓ ${SRA}.sra already downloaded"
    else
        echo "Downloading ${SRA}..."
        prefetch ${SRA} -O "${OUTPUT_DIR}"
    fi
done

###############################################################################
# STEP 2: Convert to FASTQ using fastq-dump (more reliable)
###############################################################################

echo ""
echo "=== STEP 2: Converting to FASTQ ==="

convert_to_fastq() {
    local SRA=$1
    local LABEL=$2
    local SRA_FILE="${OUTPUT_DIR}/${SRA}/${SRA}.sra"

    # Check if FASTQ already exists
    if [ -f "${OUTPUT_DIR}/${LABEL}_1.fastq.gz" ] && [ -f "${OUTPUT_DIR}/${LABEL}_2.fastq.gz" ]; then
        echo "✓ ${LABEL} FASTQ already exists, skipping"
        return 0
    fi

    echo "Converting ${SRA} to FASTQ (${LABEL})..."

    # Use fastq-dump (older but more reliable for aligned SRA files)
    # --split-files: separate R1 and R2
    # --skip-technical: skip technical reads
    # --readids: append read id after spot id
    # --read-filter pass: only pass reads
    # --dumpbase: formats sequence using base space
    # --clip: apply left and right clips
    fastq-dump \
        --split-files \
        --skip-technical \
        --readids \
        --read-filter pass \
        --dumpbase \
        --clip \
        --outdir "${OUTPUT_DIR}" \
        "${SRA_FILE}"

    # Rename files
    if [ -f "${OUTPUT_DIR}/${SRA}_pass_1.fastq" ]; then
        mv "${OUTPUT_DIR}/${SRA}_pass_1.fastq" "${OUTPUT_DIR}/${LABEL}_1.fastq"
        mv "${OUTPUT_DIR}/${SRA}_pass_2.fastq" "${OUTPUT_DIR}/${LABEL}_2.fastq"
    else
        mv "${OUTPUT_DIR}/${SRA}_1.fastq" "${OUTPUT_DIR}/${LABEL}_1.fastq"
        mv "${OUTPUT_DIR}/${SRA}_2.fastq" "${OUTPUT_DIR}/${LABEL}_2.fastq"
    fi

    # Compress with pigz (parallel)
    echo "Compressing ${LABEL}..."
    pigz -p ${THREADS} "${OUTPUT_DIR}/${LABEL}_1.fastq" &
    pigz -p ${THREADS} "${OUTPUT_DIR}/${LABEL}_2.fastq" &
    wait

    echo "✓ ${LABEL} complete"
    ls -lh "${OUTPUT_DIR}/${LABEL}"_*.fastq.gz
}

# Convert all three samples
convert_to_fastq "${PROBAND}" "proband"
convert_to_fastq "${MOTHER}" "mother"
convert_to_fastq "${FATHER}" "father"

echo ""
echo "=== FASTQ files created ==="
ls -lh "${OUTPUT_DIR}"/*.fastq.gz

###############################################################################
# STEP 3: Quality Control with FastQC (parallel)
###############################################################################

echo ""
echo "=== STEP 3: Quality Control with FastQC ==="

mkdir -p "${OUTPUT_DIR}/qc"

if command -v fastqc &> /dev/null; then
    fastqc \
        --threads ${THREADS} \
        --outdir "${OUTPUT_DIR}/qc" \
        "${OUTPUT_DIR}"/*.fastq.gz
    echo "✓ FastQC reports generated in ${OUTPUT_DIR}/qc/"
else
    echo "⚠ FastQC not found, skipping QC"
fi

###############################################################################
# STEP 4: Alignment with BWA-MEM (parallel per sample)
###############################################################################

echo ""
echo "=== STEP 4: Alignment with BWA-MEM ==="

# Check if reference is indexed
if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "ERROR: BWA index not found for ${REFERENCE}"
    echo "Run: bwa index ${REFERENCE}"
    exit 1
fi

align_sample() {
    local LABEL=$1
    local FASTQ_R1="${OUTPUT_DIR}/${LABEL}_1.fastq.gz"
    local FASTQ_R2="${OUTPUT_DIR}/${LABEL}_2.fastq.gz"
    local BAM_OUT="${OUTPUT_DIR}/${LABEL}_aligned.bam"

    if [ -f "${BAM_OUT}" ]; then
        echo "✓ ${LABEL} BAM already exists, skipping alignment"
        return 0
    fi

    echo "Aligning ${LABEL}..."

    bwa mem \
        -t ${THREADS} \
        -M \
        -R "@RG\tID:${LABEL}\tSM:${LABEL}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
        "${REFERENCE}" \
        "${FASTQ_R1}" \
        "${FASTQ_R2}" | \
    samtools view -@ 8 -Sb - | \
    samtools sort -@ ${THREADS} -o "${BAM_OUT}" -

    # Index BAM
    samtools index -@ ${THREADS} "${BAM_OUT}"

    echo "✓ ${LABEL} aligned"
}

# Align all samples
align_sample "proband"
align_sample "mother"
align_sample "father"

###############################################################################
# STEP 5: BAM Statistics
###############################################################################

echo ""
echo "=== STEP 5: BAM Statistics ==="

for LABEL in proband mother father; do
    echo ""
    echo "--- ${LABEL} alignment stats ---"
    samtools flagstat -@ ${THREADS} "${OUTPUT_DIR}/${LABEL}_aligned.bam"
done

###############################################################################
# STEP 6: Mark Duplicates with GATK (parallel)
###############################################################################

echo ""
echo "=== STEP 6: Mark Duplicates ==="

# Activate varcall environment
source "${HOME}/tool/anaconda3/etc/profile.d/conda.sh"
conda activate varcall

mark_duplicates() {
    local LABEL=$1
    local BAM_IN="${OUTPUT_DIR}/${LABEL}_aligned.bam"
    local BAM_OUT="${OUTPUT_DIR}/${LABEL}_dedup.bam"
    local METRICS="${OUTPUT_DIR}/${LABEL}_dup_metrics.txt"

    if [ -f "${BAM_OUT}" ]; then
        echo "✓ ${LABEL} dedup BAM already exists"
        return 0
    fi

    echo "Marking duplicates for ${LABEL}..."

    gatk MarkDuplicates \
        --INPUT "${BAM_IN}" \
        --OUTPUT "${BAM_OUT}" \
        --METRICS_FILE "${METRICS}" \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY LENIENT

    echo "✓ ${LABEL} duplicates marked"
}

# Process all samples
mark_duplicates "proband"
mark_duplicates "mother"
mark_duplicates "father"

###############################################################################
# STEP 7: Base Quality Score Recalibration (BQSR)
###############################################################################

echo ""
echo "=== STEP 7: Base Quality Recalibration ==="

DBSNP="${REF_DIR}/dbsnp_138.hg19.vcf.gz"
MILLS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"

# Check if known sites exist
if [ ! -f "${DBSNP}" ]; then
    echo "WARNING: dbSNP file not found: ${DBSNP}"
    echo "Skipping BQSR - using dedup BAM as final BAM"
    for LABEL in proband mother father; do
        cp "${OUTPUT_DIR}/${LABEL}_dedup.bam" "${OUTPUT_DIR}/${LABEL}_final.bam"
        cp "${OUTPUT_DIR}/${LABEL}_dedup.bai" "${OUTPUT_DIR}/${LABEL}_final.bai"
    done
else
    if [ ! -f "${MILLS}" ]; then
        echo "WARNING: Mills indels not found, using dbSNP only"
        MILLS=""
    fi

    recalibrate_bases() {
        local LABEL=$1
        local BAM_IN="${OUTPUT_DIR}/${LABEL}_dedup.bam"
        local RECAL_TABLE="${OUTPUT_DIR}/${LABEL}_recal_data.table"
        local BAM_OUT="${OUTPUT_DIR}/${LABEL}_final.bam"

        if [ -f "${BAM_OUT}" ]; then
            echo "✓ ${LABEL} final BAM already exists"
            return 0
        fi

        echo "Base recalibration for ${LABEL}..."

        # BaseRecalibrator
        KNOWN_SITES_ARGS="--known-sites ${DBSNP}"
        if [ -n "${MILLS}" ]; then
            KNOWN_SITES_ARGS="${KNOWN_SITES_ARGS} --known-sites ${MILLS}"
        fi

        gatk BaseRecalibrator \
            --input "${BAM_IN}" \
            --reference "${REFERENCE}" \
            ${KNOWN_SITES_ARGS} \
            --output "${RECAL_TABLE}"

        # ApplyBQSR
        gatk ApplyBQSR \
            --input "${BAM_IN}" \
            --reference "${REFERENCE}" \
            --bqsr-recal-file "${RECAL_TABLE}" \
            --output "${BAM_OUT}"

        # Index final BAM
        samtools index -@ ${THREADS} "${BAM_OUT}"

        echo "✓ ${LABEL} recalibration complete"
    }

    # Recalibrate all samples
    recalibrate_bases "proband"
    recalibrate_bases "mother"
    recalibrate_bases "father"
fi

###############################################################################
# STEP 8: Final Statistics
###############################################################################

echo ""
echo "=========================================="
echo "Preprocessing Complete!"
echo "=========================================="
echo ""

echo "=== Final BAM files ==="
ls -lh "${OUTPUT_DIR}"/*_final.bam

echo ""
echo "=== Coverage statistics ==="
for LABEL in proband mother father; do
    echo ""
    echo "--- ${LABEL} ---"
    samtools depth "${OUTPUT_DIR}/${LABEL}_final.bam" | \
        awk '{sum+=$3; cnt++} END {print "Mean coverage: " sum/cnt "x"}'
done

echo ""
echo "=== Disk space used ==="
du -sh "${OUTPUT_DIR}"

echo ""
echo "=========================================="
echo "Ready for variant calling!"
echo "=========================================="
echo ""
echo "Final BAM files for NEW trio:"
echo "  ${OUTPUT_DIR}/proband_final.bam (${PROBAND})"
echo "  ${OUTPUT_DIR}/mother_final.bam (${MOTHER})"
echo "  ${OUTPUT_DIR}/father_final.bam (${FATHER})"
echo ""
echo "Completed: $(date)"
