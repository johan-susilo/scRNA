#!/bin/bash
###############################################################################
# Re-align Autism Trio to UCSC hg19 Reference (chr naming)
#
# Purpose: Re-align FASTQ files to proper UCSC hg19 reference with chr1, chr2...
#          naming to match dbSNP and Mills VCF files
#
# Input: FASTQ files in /home/johan/output/autism/
# Output: Aligned, deduplicated, and BQSR BAM files
#
# Trio: SRR8697636 (proband), SRR8697627 (mother), SRR8697645 (father)
###############################################################################

set -euo pipefail

# Configuration
THREADS=60
OUTPUT_DIR="${HOME}/output/autism"
REF_DIR="${HOME}/johan/johan/johan/reference/human"
REFERENCE="${REF_DIR}/hg19_ucsc.fa"

# Known sites for BQSR (with chr naming)
DBSNP="${REF_DIR}/dbsnp_138.hg19.vcf.gz"
MILLS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"

# Sample information
SAMPLES=("proband" "mother" "father")

echo "=========================================="
echo "Re-alignment with UCSC hg19 Reference"
echo "=========================================="
echo "Date: $(date)"
echo "Reference: ${REFERENCE}"
echo "Threads: ${THREADS}"
echo ""

# Verify reference files exist
if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: Reference not found: ${REFERENCE}"
    exit 1
fi

if [ ! -f "${REFERENCE}.fai" ]; then
    echo "ERROR: Reference index not found: ${REFERENCE}.fai"
    exit 1
fi

if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "ERROR: BWA index not complete. Please wait for indexing to finish."
    exit 1
fi

# Activate conda environment
source "${HOME}/tool/anaconda3/etc/profile.d/conda.sh"
conda activate varcall

cd "${OUTPUT_DIR}"

###############################################################################
# Process each sample
###############################################################################

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "=========================================="
    echo "Processing: ${SAMPLE}"
    echo "=========================================="

    FASTQ1="${OUTPUT_DIR}/${SAMPLE}_1.fastq.gz"
    FASTQ2="${OUTPUT_DIR}/${SAMPLE}_2.fastq.gz"

    # Check if FASTQ files exist
    if [ ! -f "${FASTQ1}" ] || [ ! -f "${FASTQ2}" ]; then
        echo "ERROR: FASTQ files not found for ${SAMPLE}"
        echo "  Expected: ${FASTQ1} and ${FASTQ2}"
        exit 1
    fi

    ###########################################################################
    # STEP 1: Alignment with BWA-MEM
    ###########################################################################

    BAM_ALIGNED="${OUTPUT_DIR}/${SAMPLE}_aligned_ucsc.bam"

    if [ -f "${BAM_ALIGNED}" ]; then
        echo "✓ Aligned BAM already exists: ${BAM_ALIGNED}"
    else
        echo "=== STEP 1: Aligning ${SAMPLE} with BWA-MEM ==="

        # Read group information
        RGID="${SAMPLE}"
        RGSM="${SAMPLE}"
        RGLB="lib1"
        RGPL="ILLUMINA"
        RGPU="unit1"

        bwa mem \
            -t ${THREADS} \
            -R "@RG\tID:${RGID}\tSM:${RGSM}\tLB:${RGLB}\tPL:${RGPL}\tPU:${RGPU}" \
            "${REFERENCE}" \
            "${FASTQ1}" \
            "${FASTQ2}" | \
        samtools sort \
            -@ ${THREADS} \
            -o "${BAM_ALIGNED}" \
            -

        samtools index "${BAM_ALIGNED}"

        echo "✓ Alignment complete: ${BAM_ALIGNED}"
    fi

    ###########################################################################
    # STEP 2: Mark Duplicates
    ###########################################################################

    BAM_DEDUP="${OUTPUT_DIR}/${SAMPLE}_dedup_ucsc.bam"
    METRICS="${OUTPUT_DIR}/${SAMPLE}_dedup_metrics.txt"

    if [ -f "${BAM_DEDUP}" ]; then
        echo "✓ Deduplicated BAM already exists: ${BAM_DEDUP}"
    else
        echo "=== STEP 2: Marking duplicates for ${SAMPLE} ==="

        gatk MarkDuplicates \
            --INPUT "${BAM_ALIGNED}" \
            --OUTPUT "${BAM_DEDUP}" \
            --METRICS_FILE "${METRICS}" \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY LENIENT

        echo "✓ Duplicates marked: ${BAM_DEDUP}"
    fi

    ###########################################################################
    # STEP 3: Base Quality Score Recalibration (BQSR)
    ###########################################################################

    RECAL_TABLE="${OUTPUT_DIR}/${SAMPLE}_recal_data.table"
    BAM_FINAL="${OUTPUT_DIR}/${SAMPLE}_final_ucsc.bam"

    if [ -f "${BAM_FINAL}" ]; then
        echo "✓ Final BAM already exists: ${BAM_FINAL}"
    else
        echo "=== STEP 3: Base quality recalibration for ${SAMPLE} ==="

        # Generate recalibration table
        echo "Generating recalibration table..."
        gatk BaseRecalibrator \
            --input "${BAM_DEDUP}" \
            --reference "${REFERENCE}" \
            --known-sites "${DBSNP}" \
            --known-sites "${MILLS}" \
            --output "${RECAL_TABLE}"

        # Apply BQSR
        echo "Applying base recalibration..."
        gatk ApplyBQSR \
            --input "${BAM_DEDUP}" \
            --reference "${REFERENCE}" \
            --bqsr-recal-file "${RECAL_TABLE}" \
            --output "${BAM_FINAL}"

        # Index final BAM
        samtools index "${BAM_FINAL}"

        echo "✓ BQSR complete: ${BAM_FINAL}"
    fi

    echo ""
    echo "✓ ${SAMPLE} processing complete!"
    echo ""
done

###############################################################################
# Verify chromosome naming
###############################################################################

echo "=========================================="
echo "Verifying Chromosome Naming"
echo "=========================================="
echo ""

for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FINAL="${OUTPUT_DIR}/${SAMPLE}_final_ucsc.bam"
    echo "Chromosomes in ${SAMPLE}:"
    samtools view -H "${BAM_FINAL}" | grep "^@SQ" | head -5
    echo ""
done

###############################################################################
# Summary
###############################################################################

echo "=========================================="
echo "Re-alignment Complete!"
echo "=========================================="
echo ""
echo "Final BAM files (ready for variant calling):"
for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FINAL="${OUTPUT_DIR}/${SAMPLE}_final_ucsc.bam"
    SIZE=$(du -h "${BAM_FINAL}" | cut -f1)
    echo "  ${SAMPLE}: ${BAM_FINAL} (${SIZE})"
done
echo ""
echo "Next step: Run variant calling"
echo "  bash /home/johan/pipeline/varcall/varcall_trio_ucsc.sh"
echo ""
echo "Completed: $(date)"
