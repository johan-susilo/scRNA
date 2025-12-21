#!/bin/bash
###############################################################################
# Optimized Trio Preprocessing Script for FAM92 De Novo Analysis
# Uses 64 CPU cores for maximum performance
#
# Trio samples:
#   - SRR8697636: Proband (autistic child, IQ=60, ADI=42)
#   - SRR8697627: Mother (unaffected)
#   - SRR8697645: Father (unaffected)
###############################################################################

date
set -euo pipefail

# Configuration
THREADS=60  # Use 60 out of 64 cores (leave 4 for system)
OUTPUT_DIR="${HOME}/output/autism"
REF_DIR="${HOME}/johan/johan/reference/human"
REFERENCE="${REF_DIR}/hg19.fa"

# Trio sample IDs (CORRECTED)
PROBAND="SRR8697636"
MOTHER="SRR8697627"
FATHER="SRR8697645"

echo "=========================================="
echo "Optimized Trio Preprocessing Pipeline"
echo "=========================================="
echo "Date: $(date)"
echo "CPUs available: $(nproc)"
echo "Using threads: ${THREADS}"
echo ""
echo "Trio samples:"
echo "  Proband: ${PROBAND}"
echo "  Mother:  ${MOTHER}"
echo "  Father:  ${FATHER}"
echo ""

# Create directories
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

###############################################################################
# STEP 1: Download SRA files (if not already downloaded)
###############################################################################

echo "=== STEP 1: Downloading SRA files ==="
for SRA in ${PROBAND} ${MOTHER} ${FATHER}; do
    if [ ! -f "${OUTPUT_DIR}/${SRA}/${SRA}.sra" ]; then
        echo "Downloading ${SRA}..."
        prefetch ${SRA} -O "${OUTPUT_DIR}"
    else
        echo "✓ ${SRA}.sra already downloaded"
    fi
done

###############################################################################
# STEP 2: Convert to FASTQ (parallel with max threads)
###############################################################################

echo ""
echo "=== STEP 2: Converting to FASTQ (parallel) ==="

convert_to_fastq() {
    local SRA=$1
    local LABEL=$2
    local SRA_FILE="${OUTPUT_DIR}/${SRA}/${SRA}.sra"

    # Check if FASTQ already exists
    if [ -f "${OUTPUT_DIR}/${LABEL}_1.fastq.gz" ]; then
        echo "✓ ${LABEL} FASTQ already exists, skipping"
        return 0
    fi

    echo "Converting ${SRA} to FASTQ (${LABEL})..."

    # Use fasterq-dump with multi-threading
    fasterq-dump \
        --split-files \
        --threads ${THREADS} \
        --progress \
        --outdir "${OUTPUT_DIR}" \
        "${SRA_FILE}"

    # Rename immediately
    mv "${OUTPUT_DIR}/${SRA}_1.fastq" "${OUTPUT_DIR}/${LABEL}_1.fastq"
    mv "${OUTPUT_DIR}/${SRA}_2.fastq" "${OUTPUT_DIR}/${LABEL}_2.fastq"

    # Compress with pigz (parallel gzip) for speed
    echo "Compressing ${LABEL} with pigz..."
    pigz -p ${THREADS} "${OUTPUT_DIR}/${LABEL}_1.fastq" &
    pigz -p ${THREADS} "${OUTPUT_DIR}/${LABEL}_2.fastq" &
    wait

    echo "✓ ${LABEL} complete"
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
    echo "⚠ FastQC not found, skipping QC (install with: conda install -c bioconda fastqc)"
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

# Align all samples (run sequentially since each uses all cores)
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

# Activate varcall environment if not already
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
    echo "ERROR: dbSNP file not found: ${DBSNP}"
    exit 1
fi

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
echo "Next steps:"
echo "  1. Run variant calling with GATK HaplotypeCaller"
echo "  2. Use: /home/johan/pipeline/varcall/denovo_fam92_analysis.sh"
echo "     (Skip to STEP 2 since BAM files are ready)"
echo ""
echo "Final BAM files for trio:"
echo "  ${OUTPUT_DIR}/proband_final.bam"
echo "  ${OUTPUT_DIR}/mother_final.bam"
echo "  ${OUTPUT_DIR}/father_final.bam"
echo ""
echo "Completed: $(date)"
