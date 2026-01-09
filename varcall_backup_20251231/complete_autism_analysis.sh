#!/bin/bash
###############################################################################
# Complete Autism Trio Exome Analysis Pipeline
#
# Purpose: Comprehensive variant analysis following published methodology
#
# Analysis Components:
#   1. De Novo Variants - New mutations in proband
#   2. Rare Inherited Variants - Rare variants from parents
#
# Methodology based on:
#   - Alignment: BWA-MEM → hg19
#   - Variant Calling: GATK HaplotypeCaller
#   - Quality: Phred ≥30, Coverage ≥20x (all trio members)
#   - De Novo: Proband het, both parents hom-ref
#   - Rare Inherited: MAF ≤0.01, nonsynonymous, present in parent(s)
#
# Trio: SRR8697636 (proband), SRR8697627 (mother), SRR8697645 (father)
#
# Usage:
#   bash complete_autism_analysis.sh [options]
#
# Options:
#   --skip-download      Skip FASTQ download
#   --skip-preprocess    Skip alignment/preprocessing
#   --skip-varcall       Skip variant calling
#   --denovo-only        Only run de novo analysis
#   --inherited-only     Only run inherited variant analysis
###############################################################################

set -euo pipefail

# Configuration
OUTPUT_DIR="${HOME}/output/autism"
RESULTS_DIR="${OUTPUT_DIR}/variants"

# Parse arguments
SKIP_DOWNLOAD=false
SKIP_PREPROCESS=false
SKIP_VARCALL=false
DENOVO_ONLY=false
INHERITED_ONLY=false

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
        --skip-varcall)
            SKIP_VARCALL=true
            shift
            ;;
        --denovo-only)
            DENOVO_ONLY=true
            shift
            ;;
        --inherited-only)
            INHERITED_ONLY=true
            shift
            ;;
        --help)
            head -30 "$0" | grep "^#" | sed 's/^# //'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo "=========================================="
echo "Complete Autism Trio Exome Analysis"
echo "=========================================="
echo "Date: $(date)"
echo ""
echo "This comprehensive pipeline identifies:"
echo "  1. De novo variants (new mutations in child)"
echo "  2. Rare inherited variants (from parents)"
echo ""
echo "Trio:"
echo "  Proband: SRR8697636 (affected, autistic child)"
echo "  Mother:  SRR8697627 (unaffected)"
echo "  Father:  SRR8697645 (unaffected)"
echo ""

# Create results directory
mkdir -p "${RESULTS_DIR}"

###############################################################################
# PHASE 1: Data Preparation
###############################################################################

echo "=========================================="
echo "PHASE 1: Data Preparation"
echo "=========================================="
echo ""

# Step 1: Download FASTQ files
if [ "${SKIP_DOWNLOAD}" = false ]; then
    echo "--- Step 1.1: Downloading FASTQ from ENA ---"
    bash /home/johan/pipeline/varcall/convert_sra_fixed.sh
    echo ""
    echo "✓ FASTQ download complete"
    echo ""
else
    echo "Skipping FASTQ download (--skip-download specified)"
    echo ""
fi

# Step 2: Preprocessing
if [ "${SKIP_PREPROCESS}" = false ]; then
    echo "--- Step 1.2: Preprocessing (Alignment, Deduplication, BQSR) ---"
    echo ""
    echo "This will take 4-6 hours with 60 threads..."
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
# PHASE 2: Variant Calling
###############################################################################

if [ "${SKIP_VARCALL}" = false ]; then
    echo "=========================================="
    echo "PHASE 2: Variant Calling"
    echo "=========================================="
    echo ""

    # Check if variant calling already done
    COHORT_VCF="${RESULTS_DIR}/trio_cohort.vcf.gz"
    if [ -f "${COHORT_VCF}" ]; then
        echo "✓ Variant calling already complete"
        echo ""
    else
        echo "--- Running GATK HaplotypeCaller & Joint Genotyping ---"
        echo ""
        echo "This will take 2-4 hours..."
        echo ""

        read -p "Start variant calling? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            # Run variant calling from complete_trio_pipeline.sh
            # (We'll extract just the variant calling steps)
            bash /home/johan/pipeline/varcall/complete_trio_pipeline.sh
            echo ""
            echo "✓ Variant calling complete"
            echo ""
        else
            echo "Variant calling skipped by user"
            exit 0
        fi
    fi
else
    echo "Skipping variant calling (--skip-varcall specified)"
    echo ""
fi

###############################################################################
# PHASE 3: De Novo Variant Analysis
###############################################################################

if [ "${INHERITED_ONLY}" = false ]; then
    echo "=========================================="
    echo "PHASE 3: De Novo Variant Analysis"
    echo "=========================================="
    echo ""
    echo "Identifying new mutations in the proband..."
    echo "Criteria:"
    echo "  - Proband: Heterozygous (0/1)"
    echo "  - Mother: Homozygous reference (0/0)"
    echo "  - Father: Homozygous reference (0/0)"
    echo "  - Quality: GQ≥30, DP≥20x"
    echo ""

    DENOVO_VCF="${RESULTS_DIR}/denovo_high_confidence.vcf.gz"
    if [ -f "${DENOVO_VCF}" ]; then
        echo "✓ De novo analysis already complete"
        DENOVO_COUNT=$(bcftools view -H "${DENOVO_VCF}" | wc -l)
        echo "  High-confidence de novo variants: ${DENOVO_COUNT}"
        echo ""
    else
        read -p "Run de novo analysis? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            # De novo analysis is part of complete_trio_pipeline.sh
            # If not run yet, run it now
            if [ ! -f "${DENOVO_VCF}" ]; then
                bash /home/johan/pipeline/varcall/complete_trio_pipeline.sh
            fi
            echo ""
            echo "✓ De novo analysis complete"
            echo ""
        else
            echo "De novo analysis skipped by user"
        fi
    fi
fi

###############################################################################
# PHASE 4: Rare Inherited Variant Analysis
###############################################################################

if [ "${DENOVO_ONLY}" = false ]; then
    echo "=========================================="
    echo "PHASE 4: Rare Inherited Variant Analysis"
    echo "=========================================="
    echo ""
    echo "Identifying rare variants inherited from parents..."
    echo "Criteria:"
    echo "  - Quality: Phred≥30, DP≥20x (all members)"
    echo "  - Rarity: MAF≤0.01 (1% in population)"
    echo "  - Function: Nonsynonymous variants"
    echo "  - Inheritance: Present in proband AND parent(s)"
    echo ""

    INHERITED_VCF="${RESULTS_DIR}/rare_inherited/rare_inherited_variants.vcf.gz"
    if [ -f "${INHERITED_VCF}" ]; then
        echo "✓ Rare inherited analysis already complete"
        INHERITED_COUNT=$(bcftools view -H "${INHERITED_VCF}" | wc -l)
        echo "  Rare inherited variants: ${INHERITED_COUNT}"
        echo ""
    else
        read -p "Run rare inherited variant analysis? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            bash /home/johan/pipeline/varcall/rare_inherited_variants_pipeline.sh
            echo ""
            echo "✓ Rare inherited analysis complete"
            echo ""
        else
            echo "Rare inherited analysis skipped by user"
        fi
    fi
fi

###############################################################################
# PHASE 5: Comprehensive Analysis Report
###############################################################################

echo "=========================================="
echo "PHASE 5: Generating Comprehensive Report"
echo "=========================================="
echo ""

FINAL_REPORT="${RESULTS_DIR}/complete_analysis_summary.txt"

cat > "${FINAL_REPORT}" <<EOFR
========================================
Complete Autism Trio Exome Analysis
========================================
Date: $(date)

Trio Information:
  Family ID: FAM001
  Proband:   SRR8697636 (affected, autistic child, IQ=60, ADI=42)
  Mother:    SRR8697627 (unaffected)
  Father:    SRR8697645 (unaffected)

Reference Genome: hg19 (GRCh37)

========================================
Analysis Methodology
========================================

Based on published autism exome sequencing protocols:

1. Alignment:
   - Tool: BWA-MEM
   - Reference: hg19
   - Quality: Phred ≥30

2. Variant Calling:
   - Tool: GATK HaplotypeCaller (GVCF mode)
   - Joint genotyping across trio
   - Hard filtering for quality

3. De Novo Detection:
   - Proband: Heterozygous (0/1)
   - Parents: Homozygous reference (0/0)
   - Minimum coverage: 20x in all members
   - Minimum genotype quality: GQ≥30

4. Rare Inherited Variants:
   - Maximum MAF: 0.01 (1%)
   - Databases: 1000G, ESP6500, ExAC
   - Functional impact: Nonsynonymous only
   - Inherited from ≥1 parent

========================================
Results Summary
========================================

EOFR

# Add de novo results
if [ -f "${RESULTS_DIR}/denovo_high_confidence.vcf.gz" ]; then
    DENOVO_COUNT=$(bcftools view -H "${RESULTS_DIR}/denovo_high_confidence.vcf.gz" | wc -l)
    echo "De Novo Variants:" >> "${FINAL_REPORT}"
    echo "  High-confidence: ${DENOVO_COUNT}" >> "${FINAL_REPORT}"

    # Breakdown by type
    echo "" >> "${FINAL_REPORT}"
    echo "  By variant type:" >> "${FINAL_REPORT}"
    bcftools view -H "${RESULTS_DIR}/denovo_high_confidence.vcf.gz" | \
        awk '{print length($4), length($5)}' | \
        awk '{if ($1==1 && $2==1) print "SNV"; else if ($1<$2) print "INS"; else if ($1>$2) print "DEL"; else print "COMPLEX"}' | \
        sort | uniq -c | awk '{print "    "$2": "$1}' >> "${FINAL_REPORT}"
    echo "" >> "${FINAL_REPORT}"
fi

# Add rare inherited results
if [ -f "${RESULTS_DIR}/rare_inherited/rare_inherited_variants.vcf.gz" ]; then
    INHERITED_COUNT=$(bcftools view -H "${RESULTS_DIR}/rare_inherited/rare_inherited_variants.vcf.gz" | wc -l)
    echo "Rare Inherited Variants:" >> "${FINAL_REPORT}"
    echo "  Total: ${INHERITED_COUNT}" >> "${FINAL_REPORT}"

    # Inheritance patterns
    if [ -f "${RESULTS_DIR}/rare_inherited/inherited_maternal.vcf.gz" ]; then
        MATERNAL=$(bcftools view -H "${RESULTS_DIR}/rare_inherited/inherited_maternal.vcf.gz" | wc -l)
        PATERNAL=$(bcftools view -H "${RESULTS_DIR}/rare_inherited/inherited_paternal.vcf.gz" | wc -l)
        COMPOUND=$(bcftools view -H "${RESULTS_DIR}/rare_inherited/inherited_compound.vcf.gz" | wc -l)
        RECESSIVE=$(bcftools view -H "${RESULTS_DIR}/rare_inherited/inherited_recessive.vcf.gz" | wc -l)

        echo "" >> "${FINAL_REPORT}"
        echo "  By inheritance pattern:" >> "${FINAL_REPORT}"
        echo "    Maternal: ${MATERNAL}" >> "${FINAL_REPORT}"
        echo "    Paternal: ${PATERNAL}" >> "${FINAL_REPORT}"
        echo "    Compound heterozygous: ${COMPOUND}" >> "${FINAL_REPORT}"
        echo "    Homozygous recessive: ${RECESSIVE}" >> "${FINAL_REPORT}"
    fi
    echo "" >> "${FINAL_REPORT}"
fi

cat >> "${FINAL_REPORT}" <<EOFR

========================================
Key Output Files
========================================

Variant Calling:
  ${RESULTS_DIR}/trio_cohort.vcf.gz
  ${RESULTS_DIR}/trio_filtered.vcf.gz

De Novo Analysis:
  ${RESULTS_DIR}/denovo_high_confidence.vcf.gz
  ${RESULTS_DIR}/denovo_variants_table.tsv
  ${RESULTS_DIR}/denovo_analysis_report.txt

Rare Inherited Analysis:
  ${RESULTS_DIR}/rare_inherited/rare_inherited_variants.vcf.gz
  ${RESULTS_DIR}/rare_inherited/rare_inherited_variants_table.tsv
  ${RESULTS_DIR}/rare_inherited/genes_with_rare_variants.tsv
  ${RESULTS_DIR}/rare_inherited/rare_inherited_analysis_report.txt

By Inheritance Pattern:
  ${RESULTS_DIR}/rare_inherited/inherited_maternal.vcf.gz
  ${RESULTS_DIR}/rare_inherited/inherited_paternal.vcf.gz
  ${RESULTS_DIR}/rare_inherited/inherited_compound.vcf.gz
  ${RESULTS_DIR}/rare_inherited/inherited_recessive.vcf.gz

========================================
Recommended Next Steps
========================================

1. Prioritization:
   - Review de novo variants (highest priority)
   - Check autism gene databases (SFARI Gene, AutDB)
   - Assess variant pathogenicity (ClinVar, OMIM)

2. Functional Analysis:
   - Loss-of-function variants (stop-gain, frameshift)
   - Missense in conserved regions (GERP, PhyloP)
   - Splice-site variants

3. Gene-Level Analysis:
   - Genes with multiple variants
   - Compound heterozygous pairs
   - Known autism risk genes

4. Validation:
   - Sanger sequencing for high-priority variants
   - Check parental origin
   - Family segregation analysis

5. Further Investigation:
   - Gene expression in brain tissues
   - Protein structure modeling
   - Functional studies

========================================
Quality Metrics
========================================

Filtering Criteria Applied:
  ✓ Variant quality (QUAL) ≥ 30
  ✓ Genotype quality (GQ) ≥ 30
  ✓ Read depth (DP) ≥ 20x (all trio members)
  ✓ PASS filter status
  ✓ No missing genotypes

For Rare Variants:
  ✓ MAF ≤ 0.01 in population databases
  ✓ Nonsynonymous impact only

========================================
References
========================================

Methodology based on:
- GATK Best Practices for Germline Variant Discovery
- Autism exome sequencing consortium guidelines
- De novo mutation detection protocols
- Rare variant association study methods

========================================
Analysis completed: $(date)
========================================
EOFR

echo "✓ Comprehensive report generated: ${FINAL_REPORT}"
echo ""

# Display the report
cat "${FINAL_REPORT}"

###############################################################################
# Final Summary
###############################################################################

echo ""
echo "=========================================="
echo "ANALYSIS COMPLETE!"
echo "=========================================="
echo ""
echo "All results are in: ${RESULTS_DIR}/"
echo ""
echo "Summary:"

if [ -f "${RESULTS_DIR}/denovo_high_confidence.vcf.gz" ]; then
    DENOVO_COUNT=$(bcftools view -H "${RESULTS_DIR}/denovo_high_confidence.vcf.gz" | wc -l)
    echo "  De novo variants: ${DENOVO_COUNT}"
fi

if [ -f "${RESULTS_DIR}/rare_inherited/rare_inherited_variants.vcf.gz" ]; then
    INHERITED_COUNT=$(bcftools view -H "${RESULTS_DIR}/rare_inherited/rare_inherited_variants.vcf.gz" | wc -l)
    echo "  Rare inherited variants: ${INHERITED_COUNT}"
fi

echo ""
echo "Review the comprehensive report:"
echo "  ${FINAL_REPORT}"
echo ""
echo "Completed: $(date)"
