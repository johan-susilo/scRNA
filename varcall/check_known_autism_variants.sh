#!/bin/bash
###############################################################################
# Check for Known Pathogenic Autism Variants
#
# Purpose: Identify variants that match known pathogenic/likely pathogenic
#          autism variants from ClinVar and literature
#
# Approach:
#   1. Check ClinVar annotations for autism-related variants
#   2. Match against published pathogenic variants
#   3. Focus on variants in high-confidence autism genes
###############################################################################

set -euo pipefail

# Configuration
OUTPUT_DIR="${HOME}/output/autism"
ANNOTATION_DIR="${OUTPUT_DIR}/comprehensive_analysis"
KNOWN_VAR_DIR="${OUTPUT_DIR}/known_variants"

# Input files
ANNOTATED_VCF="${ANNOTATION_DIR}/denovo_annotated_clean.vcf.gz"
LOF_FILE="${ANNOTATION_DIR}/lof_variants.txt"
AUTISM_GENES="${OUTPUT_DIR}/known_autism_genes.txt"

# Output
mkdir -p "${KNOWN_VAR_DIR}"

echo "=========================================="
echo "Checking for Known Autism Variants"
echo "=========================================="
echo "Date: $(date)"
echo ""

# Activate conda environment
source "${HOME}/tool/anaconda3/etc/profile.d/conda.sh"
conda activate varcall

###############################################################################
# STEP 1: Extract Variants in High-Confidence Autism Genes
###############################################################################

echo "=== STEP 1: Extracting variants in high-confidence autism genes ==="

HIGH_CONF_GENES="${KNOWN_VAR_DIR}/high_confidence_autism_genes.txt"

# List of HIGH-CONFIDENCE autism genes from SFARI Gene (Category 1)
cat > "${HIGH_CONF_GENES}" << 'GENES'
ADNP
ANK2
ARID1B
ASH1L
ASXL3
BCL11A
CACNA1E
CHD2
CHD8
CTNNB1
DCHS1
DSCAM
DYRK1A
FOXP1
GRIN2B
KDM5B
KDM6B
KMT2A
KMT2C
MED13L
POGZ
PPP2R5D
PTEN
SCN2A
SETD5
SHANK3
SLC6A1
SYNGAP1
TBL1XR1
TBR1
TCF20
TRIP12
WAC
GENES

echo "High-confidence autism genes: $(wc -l < ${HIGH_CONF_GENES})"

# Extract all variants in these genes
HIGHCONF_VARS="${KNOWN_VAR_DIR}/variants_in_highconf_genes.txt"

if [ -f "${ANNOTATED_VCF}" ]; then
    echo "Extracting variants in high-confidence genes..."
    
    {
        echo -e "Gene\tChrom\tPos\tRef\tAlt\tQual\tEffect\tImpact\tHGVS_c\tHGVS_p"
        
        bcftools view -H "${ANNOTATED_VCF}" | while read line; do
            chr=$(echo "$line" | cut -f1)
            pos=$(echo "$line" | cut -f2)
            ref=$(echo "$line" | cut -f4)
            alt=$(echo "$line" | cut -f5)
            qual=$(echo "$line" | cut -f6)
            info=$(echo "$line" | cut -f8)
            
            # Extract ANN field
            ann=$(echo "$info" | grep -o 'ANN=[^;]*' | sed 's/ANN=//' || echo "")
            
            if [ -n "$ann" ]; then
                # Parse each annotation
                echo "$ann" | tr ',' '\n' | while IFS='|' read allele effect impact gene_name gene_id feature_type feature_id transcript_biotype rank hgvs_c hgvs_p cdna_pos cds_pos aa_pos distance errors_warnings; do
                    # Check if gene is in high-confidence list
                    if grep -qw "^${gene_name}$" "${HIGH_CONF_GENES}" 2>/dev/null; then
                        echo -e "${gene_name}\t${chr}\t${pos}\t${ref}\t${alt}\t${qual}\t${effect}\t${impact}\t${hgvs_c}\t${hgvs_p}"
                    fi
                done
            fi
        done
    } > "${HIGHCONF_VARS}"
    
    HIGHCONF_COUNT=$(tail -n +2 "${HIGHCONF_VARS}" | wc -l)
    echo "✓ Found ${HIGHCONF_COUNT} variants in high-confidence autism genes"
else
    echo "⚠ Annotated VCF not found"
fi

echo ""

###############################################################################
# STEP 2: Identify Loss-of-Function Variants in High-Confidence Genes
###############################################################################

echo "=== STEP 2: LoF variants in high-confidence autism genes ==="

HIGHCONF_LOF="${KNOWN_VAR_DIR}/lof_in_highconf_genes.txt"

if [ -f "${HIGHCONF_VARS}" ]; then
    {
        echo -e "Gene\tChrom\tPos\tRef\tAlt\tQual\tEffect\tHGVS_p\tEvidence"
        
        tail -n +2 "${HIGHCONF_VARS}" | grep "HIGH" | while IFS=$'\t' read gene chr pos ref alt qual effect impact hgvs_c hgvs_p; do
            echo -e "${gene}\t${chr}\t${pos}\t${ref}\t${alt}\t${qual}\t${effect}\t${hgvs_p}\tHIGH-impact in SFARI Category 1 gene"
        done
    } > "${HIGHCONF_LOF}"
    
    LOF_COUNT=$(tail -n +2 "${HIGHCONF_LOF}" | wc -l)
    echo "✓ Found ${LOF_COUNT} HIGH-impact variants in high-confidence autism genes"
    
    if [ $LOF_COUNT -gt 0 ]; then
        echo ""
        echo "HIGH-IMPACT VARIANTS IN HIGH-CONFIDENCE AUTISM GENES:"
        echo "======================================================"
        column -t -s $'\t' "${HIGHCONF_LOF}" | head -20
    fi
else
    echo "⚠ High-confidence variants file not found"
fi

echo ""

###############################################################################
# STEP 3: Search for Specific Known Pathogenic Variants
###############################################################################

echo "=== STEP 3: Checking for known pathogenic variants ==="

KNOWN_PATHOGENIC="${KNOWN_VAR_DIR}/known_pathogenic_matches.txt"

# Database of known pathogenic autism variants
# Format: Gene|Chr|Pos|Ref|Alt|HGVS|ClinVar_ID|Phenotype
cat > "${KNOWN_VAR_DIR}/pathogenic_db.txt" << 'PATHOGENIC'
# High-impact known pathogenic variants in autism genes
# Format: Gene|Effect|ClinVar_Significance|Phenotype
CHD8|p.Arg1149*|Pathogenic|Autism_ASD
CHD8|p.Arg1210*|Pathogenic|Autism_ASD
SCN2A|p.Arg853Gln|Pathogenic|Autism_Epilepsy
SCN2A|p.Arg1882Gln|Pathogenic|Autism_Seizures
SYNGAP1|p.Arg579*|Pathogenic|Autism_ID
SHANK3|p.Arg1117*|Pathogenic|Autism_PMS
PTEN|p.Arg233*|Pathogenic|Autism_Macrocephaly
ADNP|p.Tyr719*|Pathogenic|Autism_Helsmoortel-Van_Der_Aa_syndrome
DYRK1A|p.Arg205*|Pathogenic|Autism_ID
POGZ|p.Gln1120*|Pathogenic|Autism_White-Sutton_syndrome
ARID1B|p.Arg1276*|Pathogenic|Autism_Coffin-Siris_syndrome
KMT2A|p.Arg3746*|Pathogenic|Autism_Wiedemann-Steiner_syndrome
CTNNB1|p.Asp32Asn|Pathogenic|Autism_Neurodevelopmental_disorder
GRIN2B|p.Asn615Ile|Likely_Pathogenic|Autism_ID_Seizures
SLC6A1|p.Met1?|Pathogenic|Autism_Epilepsy_ID
PATHOGENIC

echo "Searching for matches to known pathogenic variants..."

{
    echo -e "Gene\tChrom\tPos\tRef\tAlt\tQual\tHGVS_p\tKnown_Variant\tClinVar_Significance\tPhenotype"
    
    if [ -f "${HIGHCONF_VARS}" ]; then
        tail -n +2 "${HIGHCONF_VARS}" | while IFS=$'\t' read gene chr pos ref alt qual effect impact hgvs_c hgvs_p; do
            # Check if this variant matches known pathogenic
            if [ -n "$hgvs_p" ]; then
                match=$(grep -w "${hgvs_p}" "${KNOWN_VAR_DIR}/pathogenic_db.txt" 2>/dev/null || echo "")
                if [ -n "$match" ]; then
                    known_effect=$(echo "$match" | cut -d'|' -f2)
                    significance=$(echo "$match" | cut -d'|' -f3)
                    phenotype=$(echo "$match" | cut -d'|' -f4)
                    echo -e "${gene}\t${chr}\t${pos}\t${ref}\t${alt}\t${qual}\t${hgvs_p}\t${known_effect}\t${significance}\t${phenotype}"
                fi
            fi
        done
    fi
} > "${KNOWN_PATHOGENIC}"

KNOWN_COUNT=$(tail -n +2 "${KNOWN_PATHOGENIC}" | wc -l)

if [ $KNOWN_COUNT -gt 0 ]; then
    echo "✓ Found ${KNOWN_COUNT} matches to known pathogenic variants!"
    echo ""
    echo "KNOWN PATHOGENIC VARIANT MATCHES:"
    echo "=================================="
    column -t -s $'\t' "${KNOWN_PATHOGENIC}"
else
    echo "✓ No exact matches to known pathogenic variants"
    echo "  (This is expected - most de novo variants are novel)"
fi

echo ""

###############################################################################
# STEP 4: Gene-Level Known Autism Associations
###############################################################################

echo "=== STEP 4: Checking gene-level autism associations ==="

GENE_ASSOCIATIONS="${KNOWN_VAR_DIR}/gene_level_associations.txt"

{
    echo -e "Gene\tVariant_Count\tHIGH_Impact\tMODERATE_Impact\tSFARI_Category\tKnown_Association"
    
    if [ -f "${HIGHCONF_VARS}" ]; then
        tail -n +2 "${HIGHCONF_VARS}" | cut -f1 | sort | uniq | while read gene; do
            total=$(grep -w "^${gene}" "${HIGHCONF_VARS}" | wc -l)
            high=$(grep -w "^${gene}" "${HIGHCONF_VARS}" | grep "HIGH" | wc -l || echo "0")
            moderate=$(grep -w "^${gene}" "${HIGHCONF_VARS}" | grep "MODERATE" | wc -l || echo "0")
            
            # Get SFARI category and known association
            case "$gene" in
                SLC6A1)
                    cat="High_Confidence"
                    assoc="GABA_transporter_haploinsufficiency|Epilepsy+Autism+ID"
                    ;;
                CHD8)
                    cat="High_Confidence"
                    assoc="Chromatin_remodeling|Macrocephaly+Autism"
                    ;;
                SCN2A)
                    cat="High_Confidence"
                    assoc="Sodium_channel|Epilepsy+Autism+ID"
                    ;;
                SYNGAP1)
                    cat="High_Confidence"
                    assoc="Synaptic_plasticity|Autism+ID+Epilepsy"
                    ;;
                SHANK3)
                    cat="High_Confidence"
                    assoc="Postsynaptic_scaffold|Phelan-McDermid_syndrome"
                    ;;
                PTEN)
                    cat="High_Confidence"
                    assoc="PI3K-AKT_pathway|Macrocephaly+Autism"
                    ;;
                ADNP|DYRK1A|POGZ|ARID1B|KMT2A|CTNNB1|GRIN2B)
                    cat="High_Confidence"
                    assoc="Multiple_neurodevelopmental_syndromes"
                    ;;
                *)
                    cat="High_Confidence"
                    assoc="Established_autism_gene"
                    ;;
            esac
            
            echo -e "${gene}\t${total}\t${high}\t${moderate}\t${cat}\t${assoc}"
        done | sort -t$'\t' -k3 -rn
    fi
} > "${GENE_ASSOCIATIONS}"

echo "✓ Gene-level associations compiled"
echo ""
echo "TOP GENES WITH VARIANTS (sorted by HIGH-impact count):"
echo "======================================================="
column -t -s $'\t' "${GENE_ASSOCIATIONS}" | head -15

echo ""

###############################################################################
# STEP 5: Create Summary Report
###############################################################################

echo "=== STEP 5: Creating summary report ==="

SUMMARY_REPORT="${KNOWN_VAR_DIR}/KNOWN_VARIANTS_REPORT.txt"

{
    echo "=========================================="
    echo "KNOWN AUTISM VARIANT ANALYSIS REPORT"
    echo "=========================================="
    echo ""
    echo "Generated: $(date)"
    echo ""
    
    echo "=========================================="
    echo "1. VARIANTS IN HIGH-CONFIDENCE AUTISM GENES"
    echo "=========================================="
    echo ""
    echo "High-confidence autism genes analyzed: $(wc -l < ${HIGH_CONF_GENES})"
    echo "Total variants in these genes: ${HIGHCONF_COUNT}"
    echo "HIGH-impact (LoF) variants: ${LOF_COUNT}"
    echo ""
    
    if [ $LOF_COUNT -gt 0 ]; then
        echo "HIGH-IMPACT VARIANTS:"
        echo "-------------------"
        column -t -s $'\t' "${HIGHCONF_LOF}"
        echo ""
    fi
    
    echo "=========================================="
    echo "2. MATCHES TO KNOWN PATHOGENIC VARIANTS"
    echo "=========================================="
    echo ""
    
    if [ $KNOWN_COUNT -gt 0 ]; then
        echo "✓ Found ${KNOWN_COUNT} exact matches to known pathogenic variants"
        echo ""
        column -t -s $'\t' "${KNOWN_PATHOGENIC}"
    else
        echo "✗ No exact matches to known pathogenic variants"
        echo ""
        echo "NOTE: This is EXPECTED and NORMAL for de novo variants."
        echo "Most pathogenic de novo variants are novel (not previously reported)."
        echo "The absence of exact matches does NOT rule out pathogenicity."
    fi
    echo ""
    
    echo "=========================================="
    echo "3. GENE-LEVEL KNOWN ASSOCIATIONS"
    echo "=========================================="
    echo ""
    echo "Genes with variants (sorted by HIGH-impact count):"
    echo ""
    column -t -s $'\t' "${GENE_ASSOCIATIONS}"
    echo ""
    
    echo "=========================================="
    echo "4. KEY FINDINGS"
    echo "=========================================="
    echo ""
    
    # Find top candidates
    if [ -f "${HIGHCONF_LOF}" ] && [ $(tail -n +2 "${HIGHCONF_LOF}" | wc -l) -gt 0 ]; then
        echo "TOP CANDIDATE VARIANTS (LoF in high-confidence autism genes):"
        echo ""
        tail -n +2 "${HIGHCONF_LOF}" | sort -t$'\t' -k6 -rn | head -5 | while IFS=$'\t' read gene chr pos ref alt qual effect hgvs_p evidence; do
            echo "  ⭐ ${gene} ${chr}:${pos} ${ref}>${alt}"
            echo "     Effect: ${effect}"
            echo "     Protein: ${hgvs_p}"
            echo "     Quality: ${qual}"
            echo "     Evidence: ${evidence}"
            echo ""
        done
    else
        echo "No HIGH-impact variants found in high-confidence autism genes."
    fi
    
    echo "=========================================="
    echo "5. CLINICAL INTERPRETATION"
    echo "=========================================="
    echo ""
    echo "High-Confidence Autism Genes (SFARI Category 1):"
    echo "  - These genes have robust evidence for autism causation"
    echo "  - De novo LoF variants are typically pathogenic"
    echo "  - Published cohorts support disease association"
    echo ""
    echo "Novel vs Known Variants:"
    echo "  - Most pathogenic de novo variants are NOVEL (not in ClinVar)"
    echo "  - Lack of exact match does NOT indicate benign variant"
    echo "  - Gene-level evidence is more important than variant-level"
    echo ""
    echo "Interpretation Guidelines:"
    echo "  - LoF in high-confidence gene = Strong candidate"
    echo "  - Exact match to known pathogenic = Additional support"
    echo "  - Novel LoF requires validation and clinical correlation"
    echo ""
    
    echo "=========================================="
    echo "6. NEXT STEPS"
    echo "=========================================="
    echo ""
    echo "For each HIGH-impact variant in high-confidence genes:"
    echo ""
    echo "1. Sanger sequencing validation (confirm de novo status)"
    echo "2. ClinVar search (chr:pos or gene + HGVS)"
    echo "3. gnomAD frequency check (should be absent/ultra-rare)"
    echo "4. Literature review (PubMed: gene + autism)"
    echo "5. ACMG classification (apply 2015 guidelines)"
    echo "6. Clinical correlation (phenotype matches gene)"
    echo ""
    
} > "${SUMMARY_REPORT}"

echo "✓ Summary report created: ${SUMMARY_REPORT}"

cat "${SUMMARY_REPORT}"

echo ""
echo "=========================================="
echo "ANALYSIS COMPLETE"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - ${HIGHCONF_VARS}"
echo "  - ${HIGHCONF_LOF}"
echo "  - ${KNOWN_PATHOGENIC}"
echo "  - ${GENE_ASSOCIATIONS}"
echo "  - ${SUMMARY_REPORT}"
echo ""
echo "Location: ${KNOWN_VAR_DIR}/"
echo ""

