#!/bin/bash
#######################################################################
# De Novo FAM92 Variant Analysis for Severe ASD
# Main Question: Do de novo variants in FAM92 explain the proband's
#                severe ASD phenotype (IQ=60)?
#
# Pipeline Steps:
# 1. Run GATK variant calling on trio → FAM92_trio.vcf.gz
# 2. Filter for de novo variants (Proband=variant, Parents=normal)
# 3. Annotate variants (brain expression, pathogenicity)
# 4. Assess clinical relevance to severe ASD (IQ=60)
#######################################################################

set -e
set -o pipefail

# Activate conda environment
source /home/johan/tool/anaconda3/etc/profile.d/conda.sh
conda activate varcall

# Configuration
PIPELINE_DIR="/home/johan/pipeline/varcall"
OUTPUT_DIR="/home/johan/output/autism"
REF_DIR="/home/johan/johan/johan/reference/human"
WORK_DIR="${OUTPUT_DIR}/trio_analysis"

# Reference files
REFERENCE="${REF_DIR}/hg19.fa"
DBSNP="${REF_DIR}/dbsnp_138.hg19.vcf.gz"
MILLS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"

# Sample information (MODIFY THESE WITH YOUR TRIO SRA IDs)
PROBAND_SRA="SRR8697636"  # The affected child
FATHER_SRA="SRR8697637"   # Father (CHANGE THIS)
MOTHER_SRA="SRR8697638"   # Mother (CHANGE THIS)

# Create working directory
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

echo "=========================================="
echo "De Novo FAM92 Variant Analysis Pipeline"
echo "=========================================="
echo "Start time: $(date)"
echo ""

#######################################################################
# STEP 1: Preprocessing and Alignment (if BAM files don't exist)
#######################################################################

process_sample() {
    local SRA_ID=$1
    local SAMPLE_NAME=$2
    local SRA_DIR="${OUTPUT_DIR}/${SRA_ID}"

    echo "Processing ${SAMPLE_NAME} (${SRA_ID})..."

    # Check if BAM already exists
    if [ -f "${WORK_DIR}/${SAMPLE_NAME}_final.bam" ]; then
        echo "  → Final BAM exists, skipping preprocessing"
        return 0
    fi

    # Download SRA if needed
    if [ ! -f "${SRA_DIR}/${SRA_ID}.sra" ]; then
        echo "  → Downloading ${SRA_ID}..."
        mkdir -p "${SRA_DIR}"
        cd "${SRA_DIR}"
        prefetch ${SRA_ID}
    fi

    # Convert to FASTQ
    if [ ! -f "${SRA_DIR}/${SRA_ID}_1.fastq.gz" ]; then
        echo "  → Converting to FASTQ..."
        fasterq-dump --split-files -O "${SRA_DIR}" "${SRA_DIR}/${SRA_ID}.sra"
        gzip "${SRA_DIR}/${SRA_ID}"_*.fastq
    fi

    # Alignment with BWA
    echo "  → Aligning with BWA-MEM..."
    bwa mem -t 8 -R "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA" \
        "${REFERENCE}" \
        "${SRA_DIR}/${SRA_ID}_1.fastq.gz" \
        "${SRA_DIR}/${SRA_ID}_2.fastq.gz" | \
        samtools view -Sb - > "${WORK_DIR}/${SAMPLE_NAME}_raw.bam"

    # Sort BAM
    echo "  → Sorting BAM..."
    samtools sort -@ 8 -o "${WORK_DIR}/${SAMPLE_NAME}_sorted.bam" \
        "${WORK_DIR}/${SAMPLE_NAME}_raw.bam"
    rm "${WORK_DIR}/${SAMPLE_NAME}_raw.bam"

    # Mark duplicates
    echo "  → Marking duplicates..."
    gatk MarkDuplicates \
        -I "${WORK_DIR}/${SAMPLE_NAME}_sorted.bam" \
        -O "${WORK_DIR}/${SAMPLE_NAME}_dedup.bam" \
        -M "${WORK_DIR}/${SAMPLE_NAME}_metrics.txt" \
        --CREATE_INDEX true
    rm "${WORK_DIR}/${SAMPLE_NAME}_sorted.bam"

    # Base recalibration
    echo "  → Base quality recalibration..."
    gatk BaseRecalibrator \
        -I "${WORK_DIR}/${SAMPLE_NAME}_dedup.bam" \
        -R "${REFERENCE}" \
        --known-sites "${DBSNP}" \
        --known-sites "${MILLS}" \
        -O "${WORK_DIR}/${SAMPLE_NAME}_recal.table"

    gatk ApplyBQSR \
        -I "${WORK_DIR}/${SAMPLE_NAME}_dedup.bam" \
        -R "${REFERENCE}" \
        --bqsr-recal-file "${WORK_DIR}/${SAMPLE_NAME}_recal.table" \
        -O "${WORK_DIR}/${SAMPLE_NAME}_final.bam"

    rm "${WORK_DIR}/${SAMPLE_NAME}_dedup.bam"*

    echo "  → ${SAMPLE_NAME} preprocessing complete!"
}

# Process all trio members
echo "=== STEP 1: Preprocessing Trio Samples ==="
process_sample "${PROBAND_SRA}" "proband"
process_sample "${FATHER_SRA}" "father"
process_sample "${MOTHER_SRA}" "mother"

#######################################################################
# STEP 2: Variant Calling with GATK HaplotypeCaller
#######################################################################

echo ""
echo "=== STEP 2: GATK Variant Calling ==="

# Call variants per sample
for SAMPLE in proband father mother; do
    if [ ! -f "${WORK_DIR}/${SAMPLE}.g.vcf.gz" ]; then
        echo "Calling variants for ${SAMPLE}..."
        gatk HaplotypeCaller \
            -R "${REFERENCE}" \
            -I "${WORK_DIR}/${SAMPLE}_final.bam" \
            -O "${WORK_DIR}/${SAMPLE}.g.vcf.gz" \
            -ERC GVCF \
            --dbsnp "${DBSNP}"
    fi
done

# Joint genotyping
echo "Performing joint genotyping..."
gatk GenomicsDBImport \
    -V "${WORK_DIR}/proband.g.vcf.gz" \
    -V "${WORK_DIR}/father.g.vcf.gz" \
    -V "${WORK_DIR}/mother.g.vcf.gz" \
    --genomicsdb-workspace-path "${WORK_DIR}/trio_db" \
    -L chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
    --overwrite-existing-genomicsdb-workspace true

gatk GenotypeGVCFs \
    -R "${REFERENCE}" \
    -V gendb://"${WORK_DIR}/trio_db" \
    -O "${WORK_DIR}/FAM92_trio.vcf.gz"

echo "All variants called: FAM92_trio.vcf.gz"

#######################################################################
# STEP 3: Filter for De Novo Variants
#######################################################################

echo ""
echo "=== STEP 3: De Novo Variant Filtering ==="

# Extract variants where proband has variant and both parents are reference
bcftools view -H "${WORK_DIR}/FAM92_trio.vcf.gz" | \
    awk 'BEGIN{OFS="\t"} {
        # Parse genotypes (assuming order: proband, father, mother)
        split($10, p, ":");  # proband GT
        split($11, f, ":");  # father GT
        split($12, m, ":");  # mother GT

        proband_gt = p[1];
        father_gt = f[1];
        mother_gt = m[1];

        # De novo: proband has variant (0/1 or 1/1), parents are 0/0
        if ((proband_gt == "0/1" || proband_gt == "1/1") &&
            father_gt == "0/0" && mother_gt == "0/0") {
            print $0;
        }
    }' > "${WORK_DIR}/denovo_variants.vcf"

# Add VCF header
bcftools view -h "${WORK_DIR}/FAM92_trio.vcf.gz" > "${WORK_DIR}/denovo_variants_with_header.vcf"
cat "${WORK_DIR}/denovo_variants.vcf" >> "${WORK_DIR}/denovo_variants_with_header.vcf"
bgzip -f "${WORK_DIR}/denovo_variants_with_header.vcf"
tabix -f "${WORK_DIR}/denovo_variants_with_header.vcf.gz"

DENOVO_COUNT=$(bcftools view -H "${WORK_DIR}/denovo_variants_with_header.vcf.gz" | wc -l)
echo "De novo variants identified: ${DENOVO_COUNT}"

#######################################################################
# STEP 4: Annotate with VEP (Variant Effect Predictor)
#######################################################################

echo ""
echo "=== STEP 4: Functional Annotation ==="

# Download VEP cache if not present
VEP_CACHE="${HOME}/.vep"
if [ ! -d "${VEP_CACHE}/homo_sapiens" ]; then
    echo "Downloading VEP cache (this may take a while)..."
    vep_install -a cf -s homo_sapiens -y GRCh37 -c "${VEP_CACHE}" --CONVERT
fi

# Annotate with VEP
vep \
    --input_file "${WORK_DIR}/denovo_variants_with_header.vcf.gz" \
    --output_file "${WORK_DIR}/denovo_annotated.vcf" \
    --vcf \
    --cache \
    --dir_cache "${VEP_CACHE}" \
    --assembly GRCh37 \
    --species homo_sapiens \
    --fork 4 \
    --everything \
    --gene_phenotype \
    --pubmed \
    --regulatory \
    --force_overwrite

#######################################################################
# STEP 5: Filter for FAM92 Gene Variants
#######################################################################

echo ""
echo "=== STEP 5: FAM92 Gene-Specific Analysis ==="

# Extract FAM92 family genes (FAM92A1, FAM92A2, FAM92B)
bcftools view -H "${WORK_DIR}/denovo_annotated.vcf" | \
    grep -E "FAM92A1|FAM92A2|FAM92B" > "${WORK_DIR}/FAM92_denovo.txt" || true

FAM92_COUNT=$(cat "${WORK_DIR}/FAM92_denovo.txt" | wc -l)
echo "De novo variants in FAM92 genes: ${FAM92_COUNT}"

#######################################################################
# STEP 6: Clinical Interpretation
#######################################################################

echo ""
echo "=== STEP 6: Clinical Interpretation Report ==="

cat > "${WORK_DIR}/clinical_report.txt" << 'REPORT'
========================================================================
         DE NOVO FAM92 VARIANT ANALYSIS - CLINICAL REPORT
========================================================================

CLINICAL QUESTION:
  Do de novo variants in FAM92 explain the proband's severe ASD
  phenotype (IQ=60)?

PATIENT PHENOTYPE:
  - Diagnosis: Severe Autism Spectrum Disorder (ASD)
  - IQ: 60 (Mild Intellectual Disability)
  - Inheritance pattern: Sporadic (parents unaffected)

------------------------------------------------------------------------
ANALYSIS RESULTS:
------------------------------------------------------------------------

Total de novo variants identified: DENOVO_COUNT_PLACEHOLDER
De novo variants in FAM92 genes: FAM92_COUNT_PLACEHOLDER

FAM92 FAMILY GENES ANALYZED:
  - FAM92A1: Chromosome Xq21.31
  - FAM92A2: Chromosome Xq28
  - FAM92B: Chromosome 19p13.3

FAM92 GENE FUNCTIONS:
  - Brain expression: FAM92 genes show expression in brain tissue
  - Protein function: Involved in intracellular signaling
  - ASD relevance: Limited direct evidence in literature
  - Neurodevelopmental role: Emerging evidence

------------------------------------------------------------------------
VARIANT DETAILS:
------------------------------------------------------------------------

REPORT

# Add FAM92 variant details
if [ ${FAM92_COUNT} -gt 0 ]; then
    echo "" >> "${WORK_DIR}/clinical_report.txt"
    echo "FAM92 DE NOVO VARIANTS FOUND:" >> "${WORK_DIR}/clinical_report.txt"
    echo "" >> "${WORK_DIR}/clinical_report.txt"

    cat "${WORK_DIR}/FAM92_denovo.txt" | while read line; do
        CHR=$(echo "$line" | awk '{print $1}')
        POS=$(echo "$line" | awk '{print $2}')
        REF=$(echo "$line" | awk '{print $4}')
        ALT=$(echo "$line" | awk '{print $5}')

        echo "Variant: ${CHR}:${POS} ${REF}>${ALT}" >> "${WORK_DIR}/clinical_report.txt"
        echo "$line" | grep -oP 'SYMBOL=[^;]+|Consequence=[^;]+|SIFT=[^;]+|PolyPhen=[^;]+|CADD_PHRED=[^;]+' >> "${WORK_DIR}/clinical_report.txt"
        echo "" >> "${WORK_DIR}/clinical_report.txt"
    done
else
    echo "" >> "${WORK_DIR}/clinical_report.txt"
    echo "NO de novo variants found in FAM92 genes." >> "${WORK_DIR}/clinical_report.txt"
    echo "" >> "${WORK_DIR}/clinical_report.txt"
fi

# Replace placeholders
sed -i "s/DENOVO_COUNT_PLACEHOLDER/${DENOVO_COUNT}/g" "${WORK_DIR}/clinical_report.txt"
sed -i "s/FAM92_COUNT_PLACEHOLDER/${FAM92_COUNT}/g" "${WORK_DIR}/clinical_report.txt"

cat >> "${WORK_DIR}/clinical_report.txt" << 'REPORT'
------------------------------------------------------------------------
PATHOGENICITY ASSESSMENT:
------------------------------------------------------------------------

Criteria for damaging variants:
  ✓ SIFT: deleterious (< 0.05)
  ✓ PolyPhen: probably/possibly damaging (> 0.5)
  ✓ CADD score: pathogenic (> 20)
  ✓ Consequence: missense, nonsense, frameshift, splice site
  ✓ Gene constraint: pLI > 0.9 (loss-of-function intolerant)

Brain expression evidence:
  - Check GTEx database for brain tissue expression
  - Review autism gene databases (SFARI, AutDB)

------------------------------------------------------------------------
CLINICAL INTERPRETATION:
------------------------------------------------------------------------

Genotype-Phenotype Correlation:
  - IQ=60 suggests moderate impact on cognition
  - Severe ASD indicates disrupted neurodevelopment
  - De novo variants consistent with sporadic inheritance

Alternative genetic explanations:
  1. Check other ASD risk genes (CHD8, SCN2A, SHANK3, etc.)
  2. Consider copy number variants (CNVs)
  3. Evaluate polygenic risk score

Recommendations:
  [ ] Validate variants by Sanger sequencing
  [ ] Parental mosaicism testing
  [ ] Functional studies (if FAM92 variants found)
  [ ] Clinical genetic counseling
  [ ] Whole exome sequencing for other candidate genes

------------------------------------------------------------------------
CONCLUSION:
------------------------------------------------------------------------

REPORT

# Final conclusion
if [ ${FAM92_COUNT} -gt 0 ]; then
    echo "De novo variant(s) in FAM92 genes were identified and require" >> "${WORK_DIR}/clinical_report.txt"
    echo "further investigation to determine causality for the severe ASD" >> "${WORK_DIR}/clinical_report.txt"
    echo "phenotype. Functional validation and segregation analysis are" >> "${WORK_DIR}/clinical_report.txt"
    echo "recommended." >> "${WORK_DIR}/clinical_report.txt"
else
    echo "No de novo variants in FAM92 genes were identified. This suggests" >> "${WORK_DIR}/clinical_report.txt"
    echo "that FAM92 variants do NOT explain the proband's severe ASD" >> "${WORK_DIR}/clinical_report.txt"
    echo "phenotype. Alternative genetic causes should be investigated," >> "${WORK_DIR}/clinical_report.txt"
    echo "including other ASD candidate genes, CNVs, and polygenic factors." >> "${WORK_DIR}/clinical_report.txt"
fi

echo "" >> "${WORK_DIR}/clinical_report.txt"
echo "=======================================================================" >> "${WORK_DIR}/clinical_report.txt"
echo "Analysis completed: $(date)" >> "${WORK_DIR}/clinical_report.txt"
echo "=======================================================================" >> "${WORK_DIR}/clinical_report.txt"

# Display report
cat "${WORK_DIR}/clinical_report.txt"

echo ""
echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - All variants: ${WORK_DIR}/FAM92_trio.vcf.gz"
echo "  - De novo variants: ${WORK_DIR}/denovo_variants_with_header.vcf.gz"
echo "  - Annotated variants: ${WORK_DIR}/denovo_annotated.vcf"
echo "  - FAM92 variants: ${WORK_DIR}/FAM92_denovo.txt"
echo "  - Clinical report: ${WORK_DIR}/clinical_report.txt"
echo ""
echo "End time: $(date)"
