#!/bin/bash
###############################################################################
# Trio Variant Calling with UCSC hg19 Reference (chr naming)
#
# Purpose: Call variants, identify de novo mutations, and annotate with SnpEff
#
# Pipeline Steps:
#   1. Variant calling per sample (GATK HaplotypeCaller in GVCF mode)
#   2. Joint genotyping across trio
#   3. Variant quality filtering
#   4. De novo variant detection
#   5. SnpEff functional annotation
#   6. Extraction of high-priority variants (LoF, autism genes)
#   7. Final report generation
#
# Trio: proband (affected), mother (unaffected), father (unaffected)
###############################################################################

set -euo pipefail

# Configuration
THREADS=60
OUTPUT_DIR="${HOME}/output/autism"
REF_DIR="${HOME}/johan/johan/johan/reference/human"
REFERENCE="${REF_DIR}/hg19_ucsc.fa"
RESULTS_DIR="${OUTPUT_DIR}/variants"
ANNOTATION_DIR="${OUTPUT_DIR}/annotation"

# Known sites for filtering (with chr naming)
DBSNP="${REF_DIR}/dbsnp_138.hg19.vcf.gz"
MILLS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"

# Trio information
PROBAND="proband"
MOTHER="mother"
FATHER="father"

# Pedigree file for de novo detection
PEDIGREE="${OUTPUT_DIR}/trio.ped"

# Autism gene list
AUTISM_GENES="${OUTPUT_DIR}/known_autism_genes.txt"

echo "=========================================="
echo "Trio Variant Calling Pipeline (UCSC hg19)"
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
mkdir -p "${ANNOTATION_DIR}"
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

echo "✓ Pedigree file created: ${PEDIGREE}"
cat "${PEDIGREE}"
echo ""

###############################################################################
# STEP 1: Check if preprocessing is complete
###############################################################################

echo "=== STEP 1: Checking preprocessing status ==="

ALL_BAMS_EXIST=true
for SAMPLE in ${PROBAND} ${MOTHER} ${FATHER}; do
    BAM="${OUTPUT_DIR}/${SAMPLE}_final_ucsc.bam"
    if [ ! -f "${BAM}" ]; then
        echo "⚠ ${SAMPLE}_final_ucsc.bam not found"
        ALL_BAMS_EXIST=false
    else
        echo "✓ ${SAMPLE}_final_ucsc.bam exists"
    fi
done

if [ "${ALL_BAMS_EXIST}" = false ]; then
    echo ""
    echo "ERROR: Preprocessing not complete!"
    echo "Please run re-alignment first:"
    echo "  bash /home/johan/pipeline/varcall/realign_with_ucsc_hg19.sh"
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
    local BAM="${OUTPUT_DIR}/${SAMPLE}_final_ucsc.bam"
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

COHORT_VCF="${RESULTS_DIR}/trio_joint_called.vcf.gz"
COMBINED_GVCF="${RESULTS_DIR}/trio_combined.g.vcf.gz"

if [ ! -f "${COHORT_VCF}" ]; then
    echo "Combining GVCFs..."

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
    echo "Applying hard filters..."

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

    # Use GATK PossibleDeNovo
    gatk VariantAnnotator \
        --reference "${REFERENCE}" \
        --variant "${FILTERED_VCF}" \
        --pedigree "${PEDIGREE}" \
        --annotation PossibleDeNovo \
        --output "${RESULTS_DIR}/trio_annotated.vcf.gz"

    # Filter for de novo variants
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

# Count de novo variants
DENOVO_COUNT=$(bcftools view -H "${DENOVO_VCF}" 2>/dev/null | wc -l || echo "0")
HIGH_CONF_COUNT=$(bcftools view -H "${HIGH_CONF_DENOVO}" 2>/dev/null | wc -l || echo "0")

echo "  All de novo candidates: ${DENOVO_COUNT}"
echo "  High-confidence de novo: ${HIGH_CONF_COUNT}"

###############################################################################
# STEP 7: SnpEff Functional Annotation
###############################################################################

echo ""
echo "=== STEP 7: SnpEff functional annotation ==="

ANNOTATED_VCF="${ANNOTATION_DIR}/denovo_annotated.vcf.gz"
SNPEFF_HTML="${ANNOTATION_DIR}/snpeff_summary.html"
SNPEFF_CSV="${ANNOTATION_DIR}/snpeff_summary.csv"

if [ ! -f "${ANNOTATED_VCF}" ]; then
    echo "Running SnpEff annotation (this may take 10-20 minutes)..."

    # Check if SnpEff is available
    if ! command -v snpEff &> /dev/null; then
        echo "ERROR: SnpEff not found in conda environment"
        echo "Installing SnpEff..."
        conda install -y -c bioconda snpeff=5.1
    fi

    # Run SnpEff with increased memory
    snpEff ann \
        -Xmx16g \
        -v \
        -stats "${SNPEFF_HTML}" \
        -csvStats "${SNPEFF_CSV}" \
        hg19 \
        "${HIGH_CONF_DENOVO}" \
        2>&1 | tee "${ANNOTATION_DIR}/snpeff.log" | \
        grep -E "Reading|done|Annotating|Creating|variants" || true

    # The output goes to stdout, capture it
    snpEff ann \
        -Xmx16g \
        -stats "${SNPEFF_HTML}" \
        -csvStats "${SNPEFF_CSV}" \
        hg19 \
        "${HIGH_CONF_DENOVO}" \
        > "${ANNOTATION_DIR}/denovo_annotated.vcf" 2>/dev/null

    # Compress and index
    bgzip -f "${ANNOTATION_DIR}/denovo_annotated.vcf"
    tabix -p vcf "${ANNOTATED_VCF}"

    echo "✓ SnpEff annotation complete"
    echo "  HTML report: ${SNPEFF_HTML}"
else
    echo "✓ Annotated VCF already exists"
fi

###############################################################################
# STEP 8: Extract High-Impact Variants (Loss-of-Function)
###############################################################################

echo ""
echo "=== STEP 8: Extracting loss-of-function variants ==="

LOF_FILE="${ANNOTATION_DIR}/lof_variants.txt"

if [ -f "${ANNOTATED_VCF}" ]; then
    echo "Extracting HIGH impact variants..."

    # Extract variants with HIGH impact
    bcftools view -H "${ANNOTATED_VCF}" | grep "|HIGH|" | \
    awk '{
        chr=$1; pos=$2; ref=$4; alt=$5; qual=$6
        for(i=1;i<=NF;i++){
            if($i ~ /^ANN=/){
                ann=$i
                sub(/.*;ANN=/,"",ann)
                print chr "\t" pos "\t" ref "\t" alt "\t" qual "\t" ann
            }
        }
    }' | while IFS=$'\t' read chr pos ref alt qual ann; do
        echo "$ann" | tr ',' '\n' | grep "|HIGH|" | while IFS='|' read -r allele effect impact gene rest; do
            hgvs=$(echo "$rest" | cut -d'|' -f7)
            echo -e "$gene\t$chr\t$pos\t$ref\t$alt\t$qual\t$effect\t$hgvs"
        done
    done > "${LOF_FILE}.tmp"

    if [ -s "${LOF_FILE}.tmp" ]; then
        echo -e "Gene\tChrom\tPos\tRef\tAlt\tQual\tEffect\tHGVS_Protein" > "${LOF_FILE}"
        cat "${LOF_FILE}.tmp" >> "${LOF_FILE}"
        rm "${LOF_FILE}.tmp"
        LOF_COUNT=$(tail -n +2 "${LOF_FILE}" | wc -l)
        echo "  ✓ Found ${LOF_COUNT} loss-of-function variants"
    else
        echo "  No HIGH impact variants found"
        echo -e "Gene\tChrom\tPos\tRef\tAlt\tQual\tEffect\tHGVS_Protein" > "${LOF_FILE}"
        rm -f "${LOF_FILE}.tmp"
    fi
else
    echo "  ⚠ Skipping - annotated VCF not found"
fi

###############################################################################
# STEP 9: Extract Variants in Autism Genes
###############################################################################

echo ""
echo "=== STEP 9: Identifying variants in autism genes ==="

AUTISM_VARIANTS="${ANNOTATION_DIR}/autism_gene_variants.txt"

if [ -f "${ANNOTATED_VCF}" ] && [ -f "${AUTISM_GENES}" ]; then
    echo "Extracting all gene names..."

    bcftools view -H "${ANNOTATED_VCF}" | \
    awk '{for(i=1;i<=NF;i++){if($i ~ /^ANN=/){print $i}}}' | \
    sed 's/.*;ANN=//' | tr ',' '\n' | cut -d'|' -f4 | \
    sort -u | grep -v "^$" > "${ANNOTATION_DIR}/all_genes.txt"

    TOTAL_GENES=$(wc -l < "${ANNOTATION_DIR}/all_genes.txt")
    echo "  Total genes affected: ${TOTAL_GENES}"

    # Match with autism genes
    grep -Ff "${AUTISM_GENES}" "${ANNOTATION_DIR}/all_genes.txt" > "${ANNOTATION_DIR}/autism_genes_found.txt" || true

    if [ -s "${ANNOTATION_DIR}/autism_genes_found.txt" ]; then
        AUTISM_COUNT=$(wc -l < "${ANNOTATION_DIR}/autism_genes_found.txt")
        echo "  ✓ Found variants in ${AUTISM_COUNT} autism genes"

        # Extract details for autism gene variants
        echo -e "Gene\tChrom\tPos\tRef\tAlt\tQual\tEffect\tImpact\tHGVS" > "${AUTISM_VARIANTS}"

        while read gene; do
            bcftools view -H "${ANNOTATED_VCF}" | grep -w "$gene" | \
            awk -v g="$gene" '{
                chr=$1; pos=$2; ref=$4; alt=$5; qual=$6
                for(i=1;i<=NF;i++){
                    if($i ~ /^ANN=/){
                        ann=$i
                        sub(/.*;ANN=/,"",ann)
                        print chr "\t" pos "\t" ref "\t" alt "\t" qual "\t" ann
                    }
                }
            }' | while IFS=$'\t' read chr pos ref alt qual ann; do
                echo "$ann" | tr ',' '\n' | grep "|$gene|" | head -1 | \
                awk -F'|' -v chr="$chr" -v pos="$pos" -v ref="$ref" -v alt="$alt" -v qual="$qual" \
                    '{print $4 "\t" chr "\t" pos "\t" ref "\t" alt "\t" qual "\t" $2 "\t" $3 "\t" $11}'
            done
        done < "${ANNOTATION_DIR}/autism_genes_found.txt" >> "${AUTISM_VARIANTS}"

        echo "  ✓ Autism gene details: ${AUTISM_VARIANTS}"
    else
        echo "  No variants found in known autism genes"
        echo -e "Gene\tChrom\tPos\tRef\tAlt\tQual\tEffect\tImpact\tHGVS" > "${AUTISM_VARIANTS}"
    fi
else
    echo "  ⚠ Skipping - missing annotated VCF or autism gene list"
fi

###############################################################################
# STEP 10: Identify Top Candidates
###############################################################################

echo ""
echo "=== STEP 10: Creating prioritized candidate list ==="

TOP_CANDIDATES="${ANNOTATION_DIR}/TOP_CANDIDATES.txt"

{
    echo "=========================================="
    echo "TOP CANDIDATE VARIANTS FOR AUTISM"
    echo "=========================================="
    echo ""
    echo "Generated: $(date)"
    echo "Pipeline: Trio WGS variant calling with SnpEff annotation"
    echo ""
    echo "Total de novo variants: ${HIGH_CONF_COUNT}"
    echo ""

    if [ -f "${LOF_FILE}" ] && [ $(wc -l < "${LOF_FILE}") -gt 1 ]; then
        LOF_COUNT=$(($(wc -l < "${LOF_FILE}") - 1))
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "PRIORITY 1: LOSS-OF-FUNCTION VARIANTS (${LOF_COUNT} total)"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""

        # Check if any LoF variants are in autism genes
        if [ -f "${AUTISM_GENES}" ]; then
            cut -f1 "${LOF_FILE}" | tail -n +2 | sort -u | \
            grep -Ff "${AUTISM_GENES}" > /tmp/lof_autism_genes.txt 2>/dev/null || true

            if [ -s /tmp/lof_autism_genes.txt ]; then
                echo "⭐ LOSS-OF-FUNCTION VARIANTS IN AUTISM GENES:"
                echo ""
                head -1 "${LOF_FILE}"
                while read gene; do
                    grep "^$gene" "${LOF_FILE}"
                done < /tmp/lof_autism_genes.txt
                echo ""
                rm /tmp/lof_autism_genes.txt
            fi
        fi

        echo "All loss-of-function variants (top 20):"
        echo ""
        head -21 "${LOF_FILE}" | column -t

        if [ $(wc -l < "${LOF_FILE}") -gt 21 ]; then
            echo ""
            echo "... (showing first 20 of ${LOF_COUNT})"
        fi
        echo ""
    fi

    if [ -f "${AUTISM_VARIANTS}" ] && [ $(wc -l < "${AUTISM_VARIANTS}") -gt 1 ]; then
        AUTISM_VAR_COUNT=$(($(wc -l < "${AUTISM_VARIANTS}") - 1))
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "PRIORITY 2: ALL VARIANTS IN AUTISM GENES (${AUTISM_VAR_COUNT} total)"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""

        # Show HIGH and MODERATE impact first
        echo "HIGH and MODERATE impact variants in autism genes:"
        echo ""
        head -1 "${AUTISM_VARIANTS}"
        grep -E "HIGH|MODERATE" "${AUTISM_VARIANTS}" | head -20 || echo "None found"
        echo ""
    fi

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "NEXT STEPS"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    echo "1. Review this file for top candidates"
    echo "2. Open SnpEff HTML report: ${SNPEFF_HTML}"
    echo "3. Check population frequencies (gnomAD)"
    echo "4. Validate top candidates with Sanger sequencing"
    echo "5. Clinical geneticist consultation"
    echo ""

} > "${TOP_CANDIDATES}"

echo "✓ Top candidates list: ${TOP_CANDIDATES}"
cat "${TOP_CANDIDATES}"

###############################################################################
# STEP 11: Generate Final Summary Report
###############################################################################

echo ""
echo "=== STEP 11: Generating final summary report ==="

FINAL_REPORT="${ANNOTATION_DIR}/FINAL_ANALYSIS_REPORT.txt"

cat > "${FINAL_REPORT}" <<EOFR
========================================
AUTISM TRIO ANALYSIS - FINAL REPORT
========================================
Date: $(date)
Reference: ${REFERENCE} (UCSC hg19 with chr naming)
Trio: FAM001

Samples:
  Proband: ${PROBAND} (affected, autistic child)
  Mother:  ${MOTHER} (unaffected)
  Father:  ${FATHER} (unaffected)

========================================
VARIANT CALLING SUMMARY
========================================

Total variants called: $(bcftools view -H "${COHORT_VCF}" 2>/dev/null | wc -l || echo "N/A")
Variants after filtering: $(bcftools view -H "${FILTERED_VCF}" 2>/dev/null | wc -l || echo "N/A")

De novo candidates: ${DENOVO_COUNT}
High-confidence de novo: ${HIGH_CONF_COUNT}

========================================
SNPEFF ANNOTATION SUMMARY
========================================

EOFR

if [ -f "${LOF_FILE}" ]; then
    LOF_COUNT=$(($(wc -l < "${LOF_FILE}") - 1))
    echo "Loss-of-function (HIGH impact) variants: ${LOF_COUNT}" >> "${FINAL_REPORT}"
else
    echo "Loss-of-function (HIGH impact) variants: 0" >> "${FINAL_REPORT}"
fi

if [ -f "${ANNOTATION_DIR}/all_genes.txt" ]; then
    TOTAL_GENES=$(wc -l < "${ANNOTATION_DIR}/all_genes.txt")
    echo "Total genes affected: ${TOTAL_GENES}" >> "${FINAL_REPORT}"
fi

if [ -f "${ANNOTATION_DIR}/autism_genes_found.txt" ]; then
    AUTISM_COUNT=$(wc -l < "${ANNOTATION_DIR}/autism_genes_found.txt")
    echo "Autism genes with variants: ${AUTISM_COUNT}" >> "${FINAL_REPORT}"
fi

cat >> "${FINAL_REPORT}" <<EOFR

========================================
OUTPUT FILES
========================================

Variant Calling Results:
  ${COHORT_VCF}
  ${HIGH_CONF_DENOVO}

SnpEff Annotations:
  ${ANNOTATED_VCF}
  ${SNPEFF_HTML} (⭐ OPEN IN BROWSER)

Analysis Results:
  ${LOF_FILE}
  ${AUTISM_VARIANTS}
  ${TOP_CANDIDATES} (⭐ REVIEW THIS FIRST)

Reports:
  ${FINAL_REPORT} (this file)

========================================
NEXT STEPS
========================================

1. Review TOP_CANDIDATES.txt
2. Open SnpEff HTML report in web browser
3. Check population frequencies in gnomAD
4. Validate top candidates with Sanger sequencing
5. Clinical geneticist consultation

========================================
COMPLETED: $(date)
========================================
EOFR

echo "✓ Final report: ${FINAL_REPORT}"
echo ""

###############################################################################
# Display Final Summary
###############################################################################

cat "${FINAL_REPORT}"

echo ""
echo "=========================================="
echo "PIPELINE COMPLETE!"
echo "=========================================="
echo ""
echo "⭐ KEY OUTPUT FILES:"
echo ""
echo "  1. TOP_CANDIDATES.txt"
echo "     ${TOP_CANDIDATES}"
echo ""
echo "  2. SnpEff HTML Report (open in browser)"
echo "     ${SNPEFF_HTML}"
echo ""
echo "  3. Loss-of-function variants"
echo "     ${LOF_FILE}"
echo ""
echo "  4. Autism gene variants"
echo "     ${AUTISM_VARIANTS}"
echo ""
echo "  5. Annotated VCF"
echo "     ${ANNOTATED_VCF}"
echo ""
echo "=========================================="
echo "Completed: $(date)"
echo "=========================================="
echo ""
