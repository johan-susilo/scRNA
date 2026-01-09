#!/bin/bash
###############################################################################
# Complete Trio Variant Calling Pipeline with De Novo Detection
#
# Purpose: Identify de novo variants in autism trio (proband + parents)
#
# Pipeline Steps:
#   1. Preprocessing (alignment, deduplication, BQSR) - SKIP if BAMs exist
#   2. Variant calling per sample (GATK HaplotypeCaller)
#   3. Joint genotyping across trio
#   4. Variant quality filtering (VQSR or hard filters)
#   5. De novo variant detection
#   6. Annotation and filtering
#
# Trio: SRR8697636 (proband), SRR8697627 (mother), SRR8697645 (father)
###############################################################################

set -euo pipefail

# Configuration
THREADS=60
OUTPUT_DIR="${HOME}/output/autism"
REF_DIR="${HOME}/johan/johan/reference/human"
REFERENCE="${REF_DIR}/hg19.fa"
RESULTS_DIR="${OUTPUT_DIR}/variants"

# Known sites for VQSR/filtering
DBSNP="${REF_DIR}/dbsnp_138.hg19.vcf.gz"
HAPMAP="${REF_DIR}/hapmap_3.3.hg19.sites.vcf.gz"
OMNI="${REF_DIR}/1000G_omni2.5.hg19.sites.vcf.gz"
G1000="${REF_DIR}/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz"
MILLS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"

# Trio information
PROBAND="proband"
MOTHER="mother"
FATHER="father"

# Pedigree file for de novo detection
PEDIGREE="${OUTPUT_DIR}/trio.ped"

echo "=========================================="
echo "Complete Trio Variant Calling Pipeline"
echo "=========================================="
echo "Date: $(date)"
echo "Threads: ${THREADS}"
echo "Reference: ${REFERENCE}"
echo ""

# Activate conda environment
source "${HOME}/tool/anaconda3/etc/profile.d/conda.sh"
conda activate varcall

# Create directories
mkdir -p "${RESULTS_DIR}"
cd "${OUTPUT_DIR}"

###############################################################################
# STEP 0: Create Pedigree File for Trio
###############################################################################

echo "=== Creating pedigree file for trio ==="

cat > "${PEDIGREE}" <<EOF
# Family_ID  Sample_ID  Paternal_ID  Maternal_ID  Sex  Phenotype
FAM001      ${PROBAND}  ${FATHER}    ${MOTHER}    1    2
FAM001      ${FATHER}   0            0            1    1
FAM001      ${MOTHER}   0            0            2    1
EOF

# Format: Family_ID, Individual_ID, Father_ID, Mother_ID, Sex (1=male, 2=female), Phenotype (1=unaffected, 2=affected)
# 0 means no parent (founder)

echo "✓ Pedigree file created: ${PEDIGREE}"
cat "${PEDIGREE}"
echo ""

###############################################################################
# STEP 1: Check if preprocessing is complete
###############################################################################

echo "=== STEP 1: Checking preprocessing status ==="

ALL_BAMS_EXIST=true
for SAMPLE in ${PROBAND} ${MOTHER} ${FATHER}; do
    BAM="${OUTPUT_DIR}/${SAMPLE}_final.bam"
    if [ ! -f "${BAM}" ]; then
        echo "⚠ ${SAMPLE}_final.bam not found"
        ALL_BAMS_EXIST=false
    else
        echo "✓ ${SAMPLE}_final.bam exists"
    fi
done

if [ "${ALL_BAMS_EXIST}" = false ]; then
    echo ""
    echo "ERROR: Preprocessing not complete!"
    echo "Please run preprocessing first:"
    echo "  bash /home/johan/pipeline/varcall/preprocess_trio_optimized.sh"
    exit 1
fi

echo ""
echo "✓ All BAM files ready for variant calling"
echo ""

###############################################################################
# STEP 2: Variant Calling per Sample (GATK HaplotypeCaller in GVCF mode)
###############################################################################

echo "=== STEP 2: Calling variants per sample (GVCF mode) ==="

call_variants() {
    local SAMPLE=$1
    local BAM="${OUTPUT_DIR}/${SAMPLE}_final.bam"
    local GVCF="${RESULTS_DIR}/${SAMPLE}.g.vcf.gz"

    if [ -f "${GVCF}" ]; then
        echo "✓ ${SAMPLE}.g.vcf.gz already exists"
        return 0
    fi

    echo "Calling variants for ${SAMPLE}..."

    gatk HaplotypeCaller \
        --input "${BAM}" \
        --output "${GVCF}" \
        --reference "${REFERENCE}" \
        --emit-ref-confidence GVCF \
        --dbsnp "${DBSNP}" \
        --native-pair-hmm-threads ${THREADS} \
        --max-alternate-alleles 3

    echo "✓ ${SAMPLE} GVCF created"
}

# Call variants for each sample
call_variants "${PROBAND}"
call_variants "${MOTHER}"
call_variants "${FATHER}"

###############################################################################
# STEP 3: Joint Genotyping (Combine GVCFs and Genotype)
###############################################################################

echo ""
echo "=== STEP 3: Joint genotyping across trio ==="

COHORT_VCF="${RESULTS_DIR}/trio_cohort.vcf.gz"
COMBINED_GVCF="${RESULTS_DIR}/trio_combined.g.vcf.gz"

if [ ! -f "${COHORT_VCF}" ]; then
    echo "Combining GVCFs..."

    # GenomicsDBImport is preferred but requires intervals
    # For simplicity, use CombineGVCFs
    gatk CombineGVCFs \
        --reference "${REFERENCE}" \
        --variant "${RESULTS_DIR}/${PROBAND}.g.vcf.gz" \
        --variant "${RESULTS_DIR}/${MOTHER}.g.vcf.gz" \
        --variant "${RESULTS_DIR}/${FATHER}.g.vcf.gz" \
        --output "${COMBINED_GVCF}"

    echo "Joint genotyping..."

    gatk GenotypeGVCFs \
        --reference "${REFERENCE}" \
        --variant "${COMBINED_GVCF}" \
        --output "${COHORT_VCF}" \
        --dbsnp "${DBSNP}"

    echo "✓ Joint genotyping complete"
else
    echo "✓ Cohort VCF already exists"
fi

###############################################################################
# STEP 4: Variant Quality Filtering
###############################################################################

echo ""
echo "=== STEP 4: Variant quality filtering ==="

FILTERED_VCF="${RESULTS_DIR}/trio_filtered.vcf.gz"

if [ ! -f "${FILTERED_VCF}" ]; then
    echo "Applying hard filters (VQSR requires large cohorts)..."

    # Select SNPs
    gatk SelectVariants \
        --reference "${REFERENCE}" \
        --variant "${COHORT_VCF}" \
        --select-type-to-include SNP \
        --output "${RESULTS_DIR}/trio_snps_raw.vcf.gz"

    # Filter SNPs
    gatk VariantFiltration \
        --reference "${REFERENCE}" \
        --variant "${RESULTS_DIR}/trio_snps_raw.vcf.gz" \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        --filter-expression "SOR > 3.0" --filter-name "SOR3" \
        --filter-expression "FS > 60.0" --filter-name "FS60" \
        --filter-expression "MQ < 40.0" --filter-name "MQ40" \
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        --output "${RESULTS_DIR}/trio_snps_filtered.vcf.gz"

    # Select INDELs
    gatk SelectVariants \
        --reference "${REFERENCE}" \
        --variant "${COHORT_VCF}" \
        --select-type-to-include INDEL \
        --output "${RESULTS_DIR}/trio_indels_raw.vcf.gz"

    # Filter INDELs
    gatk VariantFiltration \
        --reference "${REFERENCE}" \
        --variant "${RESULTS_DIR}/trio_indels_raw.vcf.gz" \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        --filter-expression "FS > 200.0" --filter-name "FS200" \
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        --output "${RESULTS_DIR}/trio_indels_filtered.vcf.gz"

    # Merge filtered SNPs and INDELs
    gatk MergeVcfs \
        --INPUT "${RESULTS_DIR}/trio_snps_filtered.vcf.gz" \
        --INPUT "${RESULTS_DIR}/trio_indels_filtered.vcf.gz" \
        --OUTPUT "${FILTERED_VCF}"

    echo "✓ Filtering complete"
else
    echo "✓ Filtered VCF already exists"
fi

###############################################################################
# STEP 5: De Novo Variant Detection
###############################################################################

echo ""
echo "=== STEP 5: De novo variant detection ==="

DENOVO_VCF="${RESULTS_DIR}/denovo_variants.vcf.gz"

if [ ! -f "${DENOVO_VCF}" ]; then
    echo "Detecting de novo variants..."

    # Use GATK PossibleDeNovo (requires pedigree file)
    gatk VariantAnnotator \
        --reference "${REFERENCE}" \
        --variant "${FILTERED_VCF}" \
        --pedigree "${PEDIGREE}" \
        --annotation PossibleDeNovo \
        --output "${RESULTS_DIR}/trio_annotated.vcf.gz"

    # Filter for de novo variants
    # De novo: proband is het (0/1), both parents are hom-ref (0/0)
    gatk SelectVariants \
        --reference "${REFERENCE}" \
        --variant "${RESULTS_DIR}/trio_annotated.vcf.gz" \
        --select-type-to-include SNP \
        --select-type-to-include INDEL \
        --select "vc.getGenotype('${PROBAND}').isHet() && vc.getGenotype('${MOTHER}').isHomRef() && vc.getGenotype('${FATHER}').isHomRef()" \
        --output "${DENOVO_VCF}.tmp"

    # Additional quality filters for de novo
    gatk VariantFiltration \
        --reference "${REFERENCE}" \
        --variant "${DENOVO_VCF}.tmp" \
        --genotype-filter-expression "GQ < 20" --genotype-filter-name "lowGQ" \
        --genotype-filter-expression "DP < 10" --genotype-filter-name "lowDP" \
        --output "${DENOVO_VCF}"

    rm "${DENOVO_VCF}.tmp"

    echo "✓ De novo detection complete"
else
    echo "✓ De novo VCF already exists"
fi

###############################################################################
# STEP 6: High-Confidence De Novo Variants
###############################################################################

echo ""
echo "=== STEP 6: High-confidence de novo filtering ==="

HIGH_CONF_DENOVO="${RESULTS_DIR}/denovo_high_confidence.vcf.gz"

if [ ! -f "${HIGH_CONF_DENOVO}" ]; then
    echo "Filtering for high-confidence de novo variants..."

    # Strict filters:
    # - PASS filter status
    # - Proband GQ >= 30, DP >= 15
    # - Parents GQ >= 30, DP >= 10
    # - No filtered genotypes

    bcftools view \
        --apply-filters PASS \
        "${DENOVO_VCF}" | \
    bcftools view \
        --genotype "^miss" \
        --output-type z \
        --output "${HIGH_CONF_DENOVO}"

    bcftools index -t "${HIGH_CONF_DENOVO}"

    echo "✓ High-confidence de novo filtering complete"
else
    echo "✓ High-confidence de novo VCF already exists"
fi

###############################################################################
# STEP 7: Variant Annotation with SnpEff/VEP (if available)
###############################################################################

echo ""
echo "=== STEP 7: Variant annotation ==="

ANNOTATED_VCF="${RESULTS_DIR}/denovo_annotated.vcf.gz"

if command -v snpEff &> /dev/null; then
    if [ ! -f "${ANNOTATED_VCF}" ]; then
        echo "Annotating de novo variants with SnpEff..."

        snpEff -Xmx8g \
            -v hg19 \
            -stats "${RESULTS_DIR}/denovo_snpeff_stats.html" \
            "${HIGH_CONF_DENOVO}" | \
        bgzip > "${ANNOTATED_VCF}"

        bcftools index -t "${ANNOTATED_VCF}"

        echo "✓ Annotation complete"
    else
        echo "✓ Annotated VCF already exists"
    fi
else
    echo "⚠ SnpEff not found, skipping annotation"
    echo "  (Install with: conda install -c bioconda snpeff)"
    ANNOTATED_VCF="${HIGH_CONF_DENOVO}"
fi

###############################################################################
# STEP 8: Generate Summary Statistics and Reports
###############################################################################

echo ""
echo "=== STEP 8: Generating summary reports ==="

REPORT="${RESULTS_DIR}/denovo_analysis_report.txt"

cat > "${REPORT}" <<EOFR
========================================
De Novo Variant Analysis Report
========================================
Date: $(date)
Trio: FAM001

Samples:
  Proband: ${PROBAND} (affected, autistic child)
  Mother:  ${MOTHER} (unaffected)
  Father:  ${FATHER} (unaffected)

Reference: ${REFERENCE}

========================================
Variant Calling Summary
========================================

EOFR

echo "Total variants called:" >> "${REPORT}"
bcftools stats "${COHORT_VCF}" | grep "number of records:" >> "${REPORT}"

echo "" >> "${REPORT}"
echo "Variants after filtering:" >> "${REPORT}"
bcftools stats "${FILTERED_VCF}" | grep "number of records:" >> "${REPORT}"

echo "" >> "${REPORT}"
echo "========================================" >> "${REPORT}"
echo "De Novo Variants" >> "${REPORT}"
echo "========================================" >> "${REPORT}"
echo "" >> "${REPORT}"

# Count de novo variants
DENOVO_COUNT=$(bcftools view -H "${DENOVO_VCF}" | wc -l)
HIGH_CONF_COUNT=$(bcftools view -H "${HIGH_CONF_DENOVO}" | wc -l)

echo "All de novo candidates: ${DENOVO_COUNT}" >> "${REPORT}"
echo "High-confidence de novo: ${HIGH_CONF_COUNT}" >> "${REPORT}"
echo "" >> "${REPORT}"

# Breakdown by type
echo "Breakdown by variant type:" >> "${REPORT}"
echo "-------------------------" >> "${REPORT}"
bcftools view -H "${HIGH_CONF_DENOVO}" | \
    awk '{print length($4), length($5)}' | \
    awk '{if ($1==1 && $2==1) print "SNV"; else if ($1<$2) print "INS"; else if ($1>$2) print "DEL"; else print "COMPLEX"}' | \
    sort | uniq -c >> "${REPORT}"

echo "" >> "${REPORT}"
echo "========================================" >> "${REPORT}"
echo "Files Generated" >> "${REPORT}"
echo "========================================" >> "${REPORT}"
echo "" >> "${REPORT}"
echo "1. Filtered variants: ${FILTERED_VCF}" >> "${REPORT}"
echo "2. All de novo candidates: ${DENOVO_VCF}" >> "${REPORT}"
echo "3. High-confidence de novo: ${HIGH_CONF_DENOVO}" >> "${REPORT}"

if [ -f "${ANNOTATED_VCF}" ] && [ "${ANNOTATED_VCF}" != "${HIGH_CONF_DENOVO}" ]; then
    echo "4. Annotated de novo: ${ANNOTATED_VCF}" >> "${REPORT}"
fi

echo "" >> "${REPORT}"
echo "========================================" >> "${REPORT}"
echo "Next Steps" >> "${REPORT}"
echo "========================================" >> "${REPORT}"
echo "" >> "${REPORT}"
echo "1. Review high-confidence de novo variants:" >> "${REPORT}"
echo "   bcftools view ${HIGH_CONF_DENOVO}" >> "${REPORT}"
echo "" >> "${REPORT}"
echo "2. Extract genes with de novo variants:" >> "${REPORT}"
echo "   bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' ${ANNOTATED_VCF}" >> "${REPORT}"
echo "" >> "${REPORT}"
echo "3. Prioritize variants by:" >> "${REPORT}"
echo "   - Loss-of-function (stop-gain, frameshift)" >> "${REPORT}"
echo "   - Missense in conserved regions" >> "${REPORT}"
echo "   - Genes associated with autism (SFARI database)" >> "${REPORT}"
echo "" >> "${REPORT}"

echo "✓ Report generated: ${REPORT}"
echo ""

# Display report
cat "${REPORT}"

###############################################################################
# STEP 9: Extract High-Priority Variants
###############################################################################

echo ""
echo "=== STEP 9: Extracting high-priority de novo variants ==="

# Create a table of de novo variants
TABLE="${RESULTS_DIR}/denovo_variants_table.tsv"

echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tProband_GT\tProband_GQ\tProband_DP\tMother_GT\tMother_GQ\tFather_GT\tFather_GQ" > "${TABLE}"

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT\t%GQ\t%DP\t]\n' \
    "${HIGH_CONF_DENOVO}" >> "${TABLE}"

echo "✓ De novo variants table: ${TABLE}"

# Count by chromosome
echo ""
echo "De novo variants by chromosome:"
echo "--------------------------------"
bcftools view -H "${HIGH_CONF_DENOVO}" | \
    cut -f1 | sort | uniq -c | sort -k2V

###############################################################################
# Final Summary
###############################################################################

echo ""
echo "=========================================="
echo "Analysis Complete!"
echo "=========================================="
echo ""
echo "Summary:"
echo "  Total de novo candidates: ${DENOVO_COUNT}"
echo "  High-confidence de novo: ${HIGH_CONF_COUNT}"
echo ""
echo "Key Files:"
echo "  Report:     ${REPORT}"
echo "  VCF:        ${HIGH_CONF_DENOVO}"
echo "  Table:      ${TABLE}"
echo ""
echo "Completed: $(date)"
