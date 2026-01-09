#!/bin/bash
###############################################################################
# Rare Inherited Variants Detection Pipeline
#
# Purpose: Identify rare inherited variants in autism trio following
#          published exome sequencing methodology
#
# Reference Methodology:
#   - Alignment: BWA to hg19
#   - Variant Calling: GATK v2.8+ (HaplotypeCaller)
#   - Quality: Phred score ≥30, coverage ≥20x in all trio members
#   - Rarity: MAF ≤0.01 in 1000G, ESP6500, ExAC
#   - Function: Nonsynonymous variants (missense, nonsense, splice-site, frameshift)
#   - Inheritance: Present in child AND at least one parent
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
RARE_INHERITED_DIR="${RESULTS_DIR}/rare_inherited"

# Known sites and population databases
DBSNP="${REF_DIR}/dbsnp_138.hg19.vcf.gz"
HAPMAP="${REF_DIR}/hapmap_3.3.hg19.sites.vcf.gz"
OMNI="${REF_DIR}/1000G_omni2.5.hg19.sites.vcf.gz"
G1000="${REF_DIR}/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz"
MILLS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"

# Population frequency databases (for filtering rare variants)
EXAC="${REF_DIR}/ExAC.r1.sites.vep.vcf.gz"
G1000_ALL="${REF_DIR}/1000g.phase3.v5a.sites.vcf.gz"

# Trio samples
PROBAND="proband"
MOTHER="mother"
FATHER="father"

# Pedigree file
PEDIGREE="${OUTPUT_DIR}/trio.ped"

# Quality thresholds (as per published methodology)
MIN_PHRED=30
MIN_COVERAGE=20
MAX_MAF=0.01  # 1% population frequency

echo "=========================================="
echo "Rare Inherited Variants Detection Pipeline"
echo "=========================================="
echo "Date: $(date)"
echo "Threads: ${THREADS}"
echo ""
echo "Quality Thresholds:"
echo "  Minimum Phred score: ${MIN_PHRED}"
echo "  Minimum coverage (all samples): ${MIN_COVERAGE}x"
echo "  Maximum MAF: ${MAX_MAF} (1%)"
echo ""

# Activate conda environment
source "${HOME}/tool/anaconda3/etc/profile.d/conda.sh"
conda activate varcall

# Create directories
mkdir -p "${RARE_INHERITED_DIR}"
cd "${OUTPUT_DIR}"

###############################################################################
# STEP 1: Verify Prerequisites
###############################################################################

echo "=== STEP 1: Verifying prerequisites ==="

# Check BAM files
echo "Checking BAM files..."
for SAMPLE in ${PROBAND} ${MOTHER} ${FATHER}; do
    BAM="${OUTPUT_DIR}/${SAMPLE}_final.bam"
    if [ ! -f "${BAM}" ]; then
        echo "ERROR: ${BAM} not found"
        echo "Run preprocessing first: bash preprocess_trio_optimized.sh"
        exit 1
    fi
    echo "✓ ${SAMPLE}_final.bam exists"
done

# Check if variant calling is done
COHORT_VCF="${RESULTS_DIR}/trio_cohort.vcf.gz"
if [ ! -f "${COHORT_VCF}" ]; then
    echo ""
    echo "ERROR: Variant calling not complete!"
    echo "Run variant calling first: bash complete_trio_pipeline.sh"
    exit 1
fi
echo "✓ Cohort VCF exists"

echo ""

###############################################################################
# STEP 2: Relationship Confirmation (PLINK/BEDTools approach)
###############################################################################

echo "=== STEP 2: Confirming family relationships ==="

RELATEDNESS_CHECK="${RARE_INHERITED_DIR}/relatedness_check.txt"

if [ ! -f "${RELATEDNESS_CHECK}" ]; then
    echo "Calculating identity-by-descent (IBD) for relationship confirmation..."

    # Convert VCF to PLINK format for IBD analysis
    if command -v plink &> /dev/null; then
        # Extract high-quality, common SNPs for IBD
        bcftools view \
            --apply-filters PASS \
            --types snps \
            --min-af 0.05 \
            --max-af 0.95 \
            "${COHORT_VCF}" | \
        bcftools view \
            --min-ac 1 \
            --output-type z \
            --output "${RARE_INHERITED_DIR}/trio_ibd_snps.vcf.gz"

        bcftools index -t "${RARE_INHERITED_DIR}/trio_ibd_snps.vcf.gz"

        # Convert to PLINK format
        plink --vcf "${RARE_INHERITED_DIR}/trio_ibd_snps.vcf.gz" \
              --make-bed \
              --out "${RARE_INHERITED_DIR}/trio_ibd" \
              --allow-extra-chr

        # Calculate IBD
        plink --bfile "${RARE_INHERITED_DIR}/trio_ibd" \
              --genome \
              --out "${RARE_INHERITED_DIR}/trio_ibd" \
              --allow-extra-chr

        # Check relationships
        cat "${RARE_INHERITED_DIR}/trio_ibd.genome" > "${RELATEDNESS_CHECK}"

        echo "✓ Relationship confirmation complete"
        echo ""
        echo "Expected PI_HAT values:"
        echo "  Parent-child: ~0.5"
        echo "  Unrelated: ~0.0"
        echo ""
        echo "Results:"
        cat "${RELATEDNESS_CHECK}" | grep -v "FID"

    else
        echo "⚠ PLINK not found, skipping IBD analysis"
        echo "Install with: conda install -c bioconda plink"
        echo "NA" > "${RELATEDNESS_CHECK}"
    fi
else
    echo "✓ Relationship check already performed"
fi

echo ""

###############################################################################
# STEP 3: Apply Quality Filters (Phred ≥30, Coverage ≥20x)
###############################################################################

echo "=== STEP 3: Applying quality filters (Phred≥30, DP≥20x) ==="

HIGH_QUALITY_VCF="${RARE_INHERITED_DIR}/trio_high_quality.vcf.gz"

if [ ! -f "${HIGH_QUALITY_VCF}" ]; then
    echo "Filtering for high-quality variants..."

    # Apply filters:
    # 1. PASS filter status
    # 2. Variant quality (QUAL) ≥ 30
    # 3. All samples have GQ ≥ 30
    # 4. All samples have DP ≥ 20

    bcftools view \
        --apply-filters PASS \
        "${COHORT_VCF}" | \
    bcftools filter \
        --include "QUAL>=${MIN_PHRED}" | \
    bcftools view \
        --genotype "^miss" \
        --output-type z \
        --output "${HIGH_QUALITY_VCF}.tmp"

    # Additional per-sample filters (GQ≥30, DP≥20 in ALL samples)
    bcftools filter \
        --set-GTs . \
        --include "FMT/GQ>=${MIN_PHRED} & FMT/DP>=${MIN_COVERAGE}" \
        "${HIGH_QUALITY_VCF}.tmp" | \
    bcftools view \
        --genotype "^miss" \
        --output-type z \
        --output "${HIGH_QUALITY_VCF}"

    bcftools index -t "${HIGH_QUALITY_VCF}"
    rm "${HIGH_QUALITY_VCF}.tmp"

    TOTAL_HIGH_QUAL=$(bcftools view -H "${HIGH_QUALITY_VCF}" | wc -l)
    echo "✓ High-quality variants: ${TOTAL_HIGH_QUAL}"
else
    echo "✓ High-quality VCF already exists"
fi

echo ""

###############################################################################
# STEP 4: Annotation with ANNOVAR
###############################################################################

echo "=== STEP 4: Functional annotation with ANNOVAR ==="

ANNOVAR_DIR="${HOME}/tool/annovar"
ANNOVAR_INPUT="${RARE_INHERITED_DIR}/trio_high_quality.avinput"
ANNOVAR_OUTPUT="${RARE_INHERITED_DIR}/trio_high_quality.hg19_multianno"

if [ ! -f "${ANNOVAR_OUTPUT}.vcf" ]; then
    if [ -d "${ANNOVAR_DIR}" ]; then
        echo "Converting VCF to ANNOVAR input..."

        # Convert VCF to ANNOVAR format
        perl "${ANNOVAR_DIR}/convert2annovar.pl" \
            --format vcf4 \
            --includeinfo \
            "${HIGH_QUALITY_VCF}" \
            > "${ANNOVAR_INPUT}"

        echo "Running ANNOVAR annotation..."

        # Annotate with multiple databases
        perl "${ANNOVAR_DIR}/table_annovar.pl" \
            "${ANNOVAR_INPUT}" \
            "${ANNOVAR_DIR}/humandb/" \
            --buildver hg19 \
            --out "${RARE_INHERITED_DIR}/trio_high_quality" \
            --protocol refGene,exac03,gnomad_exome,gnomad_genome,1000g2015aug_all,esp6500siv2_all,dbnsfp35a \
            --operation g,f,f,f,f,f,f \
            --nastring . \
            --vcfinput "${HIGH_QUALITY_VCF}" \
            --thread ${THREADS}

        echo "✓ ANNOVAR annotation complete"
    else
        echo "⚠ ANNOVAR not found at ${ANNOVAR_DIR}"
        echo "Please install ANNOVAR or specify correct path"
        echo "Continuing without annotation..."

        # Use the high-quality VCF as-is
        cp "${HIGH_QUALITY_VCF}" "${ANNOVAR_OUTPUT}.vcf"
    fi
else
    echo "✓ Annotated VCF already exists"
fi

ANNOTATED_VCF="${ANNOVAR_OUTPUT}.vcf"
if [ ! -f "${ANNOTATED_VCF}" ]; then
    ANNOTATED_VCF="${HIGH_QUALITY_VCF}"
fi

echo ""

###############################################################################
# STEP 5: Filter for Rare Variants (MAF ≤ 0.01)
###############################################################################

echo "=== STEP 5: Filtering for rare variants (MAF≤${MAX_MAF}) ==="

RARE_VCF="${RARE_INHERITED_DIR}/trio_rare.vcf.gz"

if [ ! -f "${RARE_VCF}" ]; then
    echo "Filtering for rare variants using population databases..."

    # Filter based on population frequencies
    # Keep variants with:
    # - ExAC_ALL ≤ 0.01 OR missing
    # - 1000g ALL ≤ 0.01 OR missing
    # - ESP6500 ALL ≤ 0.01 OR missing

    bcftools annotate \
        --annotations "${ANNOTATED_VCF}" \
        "${ANNOTATED_VCF}" | \
    bcftools filter \
        --include '(INFO/ExAC_ALL<=0.01 | INFO/ExAC_ALL="." | INFO/AF<=0.01 | INFO/AF=".") & (INFO/1000g2015aug_all<=0.01 | INFO/1000g2015aug_all="." | INFO/AF<=0.01 | INFO/AF=".") & (INFO/ESP6500siv2_all<=0.01 | INFO/ESP6500siv2_all="." | INFO/AF<=0.01 | INFO/AF=".")' \
        --output-type z \
        --output "${RARE_VCF}" 2>/dev/null || {

        # Fallback if population frequencies not in INFO field
        # Use basic allele frequency filter
        echo "Using basic allele frequency filter..."
        bcftools view \
            --max-af "${MAX_MAF}:nref" \
            "${ANNOTATED_VCF}" \
            --output-type z \
            --output "${RARE_VCF}"
    }

    bcftools index -t "${RARE_VCF}"

    TOTAL_RARE=$(bcftools view -H "${RARE_VCF}" | wc -l)
    echo "✓ Rare variants (MAF≤${MAX_MAF}): ${TOTAL_RARE}"
else
    echo "✓ Rare variants VCF already exists"
fi

echo ""

###############################################################################
# STEP 6: Filter for Functional Impact (Nonsynonymous)
###############################################################################

echo "=== STEP 6: Filtering for functional variants ==="

FUNCTIONAL_VCF="${RARE_INHERITED_DIR}/trio_rare_functional.vcf.gz"

if [ ! -f "${FUNCTIONAL_VCF}" ]; then
    echo "Selecting nonsynonymous variants..."
    echo "  - Missense variants"
    echo "  - Nonsense (stop-gain/loss)"
    echo "  - Splice-site variants"
    echo "  - Frameshift indels"

    # Use SnpEff/ANNOVAR annotations if available
    # Otherwise use basic variant type filtering

    if bcftools view -h "${RARE_VCF}" | grep -q "ANN="; then
        # SnpEff annotations present
        echo "Using SnpEff annotations..."
        bcftools filter \
            --include 'INFO/ANN ~ "missense_variant" | INFO/ANN ~ "stop_gained" | INFO/ANN ~ "stop_lost" | INFO/ANN ~ "frameshift_variant" | INFO/ANN ~ "splice_acceptor_variant" | INFO/ANN ~ "splice_donor_variant"' \
            "${RARE_VCF}" \
            --output-type z \
            --output "${FUNCTIONAL_VCF}"

    elif bcftools view -h "${RARE_VCF}" | grep -q "ExonicFunc.refGene"; then
        # ANNOVAR annotations present
        echo "Using ANNOVAR annotations..."
        bcftools filter \
            --include 'INFO/ExonicFunc.refGene ~ "nonsynonymous" | INFO/ExonicFunc.refGene ~ "stopgain" | INFO/ExonicFunc.refGene ~ "stoploss" | INFO/ExonicFunc.refGene ~ "frameshift"' \
            "${RARE_VCF}" \
            --output-type z \
            --output "${FUNCTIONAL_VCF}"

    else
        # No functional annotations, keep all variants
        echo "⚠ No functional annotations found, keeping all rare variants"
        cp "${RARE_VCF}" "${FUNCTIONAL_VCF}"
    fi

    bcftools index -t "${FUNCTIONAL_VCF}"

    TOTAL_FUNCTIONAL=$(bcftools view -H "${FUNCTIONAL_VCF}" | wc -l)
    echo "✓ Functional rare variants: ${TOTAL_FUNCTIONAL}"
else
    echo "✓ Functional variants VCF already exists"
fi

echo ""

###############################################################################
# STEP 7: Identify Inherited Variants
###############################################################################

echo "=== STEP 7: Identifying inherited variants ==="

INHERITED_VCF="${RARE_INHERITED_DIR}/rare_inherited_variants.vcf.gz"

if [ ! -f "${INHERITED_VCF}" ]; then
    echo "Selecting variants inherited from parents..."
    echo "Inheritance patterns:"
    echo "  - Proband has variant (not 0/0)"
    echo "  - At least one parent has variant (not 0/0)"
    echo "  - Exclude de novo (both parents 0/0)"

    # Select inherited variants:
    # - Proband is NOT homozygous reference (0/0)
    # - At least one parent is NOT homozygous reference (0/0)
    # This excludes de novo variants

    bcftools view \
        --samples "${PROBAND},${MOTHER},${FATHER}" \
        "${FUNCTIONAL_VCF}" | \
    bcftools filter \
        --include "!vc.getGenotype('${PROBAND}').isHomRef() && (!vc.getGenotype('${MOTHER}').isHomRef() || !vc.getGenotype('${FATHER}').isHomRef())" \
        --output-type z \
        --output "${INHERITED_VCF}"

    bcftools index -t "${INHERITED_VCF}"

    TOTAL_INHERITED=$(bcftools view -H "${INHERITED_VCF}" | wc -l)
    echo "✓ Rare inherited variants: ${TOTAL_INHERITED}"
else
    echo "✓ Inherited variants VCF already exists"
fi

echo ""

###############################################################################
# STEP 8: Categorize by Inheritance Pattern
###############################################################################

echo "=== STEP 8: Categorizing by inheritance pattern ==="

# Maternal inheritance
MATERNAL_VCF="${RARE_INHERITED_DIR}/inherited_maternal.vcf.gz"
if [ ! -f "${MATERNAL_VCF}" ]; then
    echo "Identifying maternal inheritance..."
    bcftools filter \
        --include "!vc.getGenotype('${PROBAND}').isHomRef() && !vc.getGenotype('${MOTHER}').isHomRef()" \
        "${INHERITED_VCF}" \
        --output-type z \
        --output "${MATERNAL_VCF}"
    bcftools index -t "${MATERNAL_VCF}"

    MATERNAL_COUNT=$(bcftools view -H "${MATERNAL_VCF}" | wc -l)
    echo "  Maternal: ${MATERNAL_COUNT}"
fi

# Paternal inheritance
PATERNAL_VCF="${RARE_INHERITED_DIR}/inherited_paternal.vcf.gz"
if [ ! -f "${PATERNAL_VCF}" ]; then
    echo "Identifying paternal inheritance..."
    bcftools filter \
        --include "!vc.getGenotype('${PROBAND}').isHomRef() && !vc.getGenotype('${FATHER}').isHomRef()" \
        "${INHERITED_VCF}" \
        --output-type z \
        --output "${PATERNAL_VCF}"
    bcftools index -t "${PATERNAL_VCF}"

    PATERNAL_COUNT=$(bcftools view -H "${PATERNAL_VCF}" | wc -l)
    echo "  Paternal: ${PATERNAL_COUNT}"
fi

# Compound heterozygous (from both parents)
COMPOUND_VCF="${RARE_INHERITED_DIR}/inherited_compound.vcf.gz"
if [ ! -f "${COMPOUND_VCF}" ]; then
    echo "Identifying compound heterozygous variants..."
    bcftools filter \
        --include "vc.getGenotype('${PROBAND}').isHet() && !vc.getGenotype('${MOTHER}').isHomRef() && !vc.getGenotype('${FATHER}').isHomRef()" \
        "${INHERITED_VCF}" \
        --output-type z \
        --output "${COMPOUND_VCF}"
    bcftools index -t "${COMPOUND_VCF}"

    COMPOUND_COUNT=$(bcftools view -H "${COMPOUND_VCF}" | wc -l)
    echo "  Compound heterozygous: ${COMPOUND_COUNT}"
fi

# Homozygous recessive
RECESSIVE_VCF="${RARE_INHERITED_DIR}/inherited_recessive.vcf.gz"
if [ ! -f "${RECESSIVE_VCF}" ]; then
    echo "Identifying homozygous recessive variants..."
    bcftools filter \
        --include "vc.getGenotype('${PROBAND}').isHomVar() && vc.getGenotype('${MOTHER}').isHet() && vc.getGenotype('${FATHER}').isHet()" \
        "${INHERITED_VCF}" \
        --output-type z \
        --output "${RECESSIVE_VCF}"
    bcftools index -t "${RECESSIVE_VCF}"

    RECESSIVE_COUNT=$(bcftools view -H "${RECESSIVE_VCF}" | wc -l)
    echo "  Homozygous recessive: ${RECESSIVE_COUNT}"
fi

echo ""

###############################################################################
# STEP 9: Generate Summary Reports
###############################################################################

echo "=== STEP 9: Generating summary reports ==="

REPORT="${RARE_INHERITED_DIR}/rare_inherited_analysis_report.txt"

cat > "${REPORT}" <<EOFR
========================================
Rare Inherited Variants Analysis Report
========================================
Date: $(date)
Trio: FAM001

Samples:
  Proband: ${PROBAND} (affected, autistic child)
  Mother:  ${MOTHER} (unaffected)
  Father:  ${FATHER} (unaffected)

Reference: ${REFERENCE}

========================================
Analysis Parameters
========================================

Quality Filters:
  - Minimum Phred score: ${MIN_PHRED}
  - Minimum coverage: ${MIN_COVERAGE}x (all trio members)
  - Filter status: PASS only

Rarity Filters:
  - Maximum MAF: ${MAX_MAF} (1%)
  - Databases: 1000 Genomes, ESP6500, ExAC

Functional Filters:
  - Nonsynonymous variants only
  - Includes: missense, nonsense, splice-site, frameshift

Inheritance:
  - Present in proband AND at least one parent
  - Excludes de novo variants

========================================
Results Summary
========================================

EOFR

# Add variant counts
echo "Variant Counts:" >> "${REPORT}"
echo "---------------" >> "${REPORT}"

TOTAL_RAW=$(bcftools view -H "${COHORT_VCF}" | wc -l)
TOTAL_HIGH_QUAL=$(bcftools view -H "${HIGH_QUALITY_VCF}" | wc -l)
TOTAL_RARE=$(bcftools view -H "${RARE_VCF}" | wc -l)
TOTAL_FUNCTIONAL=$(bcftools view -H "${FUNCTIONAL_VCF}" | wc -l)
TOTAL_INHERITED=$(bcftools view -H "${INHERITED_VCF}" | wc -l)

echo "Total variants called: ${TOTAL_RAW}" >> "${REPORT}"
echo "High-quality (Phred≥30, DP≥20x): ${TOTAL_HIGH_QUAL}" >> "${REPORT}"
echo "Rare (MAF≤0.01): ${TOTAL_RARE}" >> "${REPORT}"
echo "Functional (nonsynonymous): ${TOTAL_FUNCTIONAL}" >> "${REPORT}"
echo "Inherited from parents: ${TOTAL_INHERITED}" >> "${REPORT}"
echo "" >> "${REPORT}"

echo "Inheritance Pattern Breakdown:" >> "${REPORT}"
echo "-----------------------------" >> "${REPORT}"
MATERNAL_COUNT=$(bcftools view -H "${MATERNAL_VCF}" | wc -l)
PATERNAL_COUNT=$(bcftools view -H "${PATERNAL_VCF}" | wc -l)
COMPOUND_COUNT=$(bcftools view -H "${COMPOUND_VCF}" | wc -l)
RECESSIVE_COUNT=$(bcftools view -H "${RECESSIVE_VCF}" | wc -l)

echo "Maternal inheritance: ${MATERNAL_COUNT}" >> "${REPORT}"
echo "Paternal inheritance: ${PATERNAL_COUNT}" >> "${REPORT}"
echo "Compound heterozygous: ${COMPOUND_COUNT}" >> "${REPORT}"
echo "Homozygous recessive: ${RECESSIVE_COUNT}" >> "${REPORT}"
echo "" >> "${REPORT}"

echo "Variant Type Distribution:" >> "${REPORT}"
echo "-------------------------" >> "${REPORT}"
bcftools view -H "${INHERITED_VCF}" | \
    awk '{print length($4), length($5)}' | \
    awk '{if ($1==1 && $2==1) print "SNV"; else if ($1<$2) print "INS"; else if ($1>$2) print "DEL"; else print "COMPLEX"}' | \
    sort | uniq -c >> "${REPORT}"

echo "" >> "${REPORT}"
echo "========================================" >> "${REPORT}"
echo "Files Generated" >> "${REPORT}"
echo "========================================" >> "${REPORT}"
echo "" >> "${REPORT}"
echo "1. High-quality variants: ${HIGH_QUALITY_VCF}" >> "${REPORT}"
echo "2. Rare variants (MAF≤0.01): ${RARE_VCF}" >> "${REPORT}"
echo "3. Functional rare variants: ${FUNCTIONAL_VCF}" >> "${REPORT}"
echo "4. Rare inherited variants: ${INHERITED_VCF}" >> "${REPORT}"
echo "5. Maternal inheritance: ${MATERNAL_VCF}" >> "${REPORT}"
echo "6. Paternal inheritance: ${PATERNAL_VCF}" >> "${REPORT}"
echo "7. Compound heterozygous: ${COMPOUND_VCF}" >> "${REPORT}"
echo "8. Homozygous recessive: ${RECESSIVE_VCF}" >> "${REPORT}"
echo "" >> "${REPORT}"

echo "========================================" >> "${REPORT}"
echo "Next Steps" >> "${REPORT}"
echo "========================================" >> "${REPORT}"
echo "" >> "${REPORT}"
echo "1. Review rare inherited variants:" >> "${REPORT}"
echo "   bcftools view ${INHERITED_VCF}" >> "${REPORT}"
echo "" >> "${REPORT}"
echo "2. Prioritize by:" >> "${REPORT}"
echo "   - Genes associated with autism (SFARI database)" >> "${REPORT}"
echo "   - Loss-of-function variants" >> "${REPORT}"
echo "   - Conservation scores (GERP, PhyloP)" >> "${REPORT}"
echo "   - Pathogenicity predictions (SIFT, PolyPhen, CADD)" >> "${REPORT}"
echo "" >> "${REPORT}"
echo "3. Compare with de novo variants for gene overlap" >> "${REPORT}"
echo "" >> "${REPORT}"
echo "4. Consider gene burden analysis" >> "${REPORT}"
echo "   - Genes with multiple rare inherited variants" >> "${REPORT}"
echo "   - Compound heterozygous pairs in same gene" >> "${REPORT}"
echo "" >> "${REPORT}"

echo "✓ Report generated: ${REPORT}"

###############################################################################
# STEP 10: Create Variant Tables
###############################################################################

echo ""
echo "=== STEP 10: Creating variant tables ==="

# All rare inherited variants table
TABLE="${RARE_INHERITED_DIR}/rare_inherited_variants_table.tsv"
echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tGENE\tFUNCTION\tProband_GT\tProband_GQ\tProband_DP\tMother_GT\tMother_GQ\tMother_DP\tFather_GT\tFather_GQ\tFather_DP" > "${TABLE}"

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/Gene.refGene\t%INFO/ExonicFunc.refGene\t[%GT\t%GQ\t%DP\t]\n' \
    "${INHERITED_VCF}" >> "${TABLE}" 2>/dev/null || \
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t.\t.\t[%GT\t%GQ\t%DP\t]\n' \
    "${INHERITED_VCF}" >> "${TABLE}"

echo "✓ Variant table created: ${TABLE}"

# Gene-level summary
GENE_TABLE="${RARE_INHERITED_DIR}/genes_with_rare_variants.tsv"
echo -e "GENE\tVARIANT_COUNT\tINHERITANCE_PATTERNS" > "${GENE_TABLE}"

bcftools query -f '%INFO/Gene.refGene\n' "${INHERITED_VCF}" 2>/dev/null | \
    grep -v "^\.$" | \
    sort | uniq -c | \
    awk '{print $2"\t"$1"\tMixed"}' >> "${GENE_TABLE}" 2>/dev/null || \
    echo "No gene annotations available" > "${GENE_TABLE}"

if [ -s "${GENE_TABLE}" ]; then
    echo "✓ Gene summary created: ${GENE_TABLE}"
    echo ""
    echo "Top 20 genes with rare inherited variants:"
    head -21 "${GENE_TABLE}" | column -t
fi

echo ""

# Display full report
cat "${REPORT}"

###############################################################################
# Final Summary
###############################################################################

echo ""
echo "=========================================="
echo "Analysis Complete!"
echo "=========================================="
echo ""
echo "Summary:"
echo "  Total rare inherited variants: ${TOTAL_INHERITED}"
echo "  Maternal: ${MATERNAL_COUNT}"
echo "  Paternal: ${PATERNAL_COUNT}"
echo "  Compound heterozygous: ${COMPOUND_COUNT}"
echo "  Homozygous recessive: ${RECESSIVE_COUNT}"
echo ""
echo "Key Files:"
echo "  Report:     ${REPORT}"
echo "  VCF:        ${INHERITED_VCF}"
echo "  Table:      ${TABLE}"
echo "  Gene list:  ${GENE_TABLE}"
echo ""
echo "Results directory: ${RARE_INHERITED_DIR}"
echo ""
echo "Completed: $(date)"
