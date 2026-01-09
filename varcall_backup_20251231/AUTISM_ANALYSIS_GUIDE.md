# Autism Trio Exome Analysis - Complete Guide

## Overview

This pipeline performs comprehensive variant analysis on an autism trio (affected child + unaffected parents) to identify genetic variants associated with autism spectrum disorder.

## Scientific Methodology

Based on published autism exome sequencing studies, this pipeline identifies:

1. **De Novo Variants** - New mutations in the child not present in either parent
2. **Rare Inherited Variants** - Rare functional variants transmitted from parents

### Quality Standards
- **Alignment**: BWA-MEM → hg19 reference genome
- **Variant Calling**: GATK HaplotypeCaller (Best Practices)
- **Quality Threshold**: Phred score ≥30, Coverage ≥20x in all trio members
- **Rarity Threshold**: MAF ≤0.01 (1% in population databases)
- **Functional Impact**: Nonsynonymous variants (missense, nonsense, splice-site, frameshift)

## Quick Start

### Option 1: Complete Analysis (Recommended)

Run the entire pipeline from start to finish:

```bash
bash /home/johan/pipeline/varcall/complete_autism_analysis.sh
```

This will:
1. Download FASTQ files from ENA
2. Perform alignment and preprocessing
3. Call variants with GATK
4. Identify de novo variants
5. Identify rare inherited variants
6. Generate comprehensive reports

**Expected Runtime**: 8-12 hours (with 60 CPU cores)

### Option 2: Step-by-Step Analysis

If you prefer to run each step separately:

```bash
# Step 1: Download FASTQ files (6-10 minutes)
bash /home/johan/pipeline/varcall/convert_sra_fixed.sh

# Step 2: Preprocessing (4-6 hours)
bash /home/johan/pipeline/varcall/preprocess_trio_optimized.sh

# Step 3: De novo variant detection (2-3 hours)
bash /home/johan/pipeline/varcall/complete_trio_pipeline.sh

# Step 4: Rare inherited variant detection (1-2 hours)
bash /home/johan/pipeline/varcall/rare_inherited_variants_pipeline.sh
```

### Option 3: Skip Completed Steps

If you've already completed some steps:

```bash
# Skip download and preprocessing
bash /home/johan/pipeline/varcall/complete_autism_analysis.sh \
    --skip-download \
    --skip-preprocess

# Run only de novo analysis
bash /home/johan/pipeline/varcall/complete_autism_analysis.sh \
    --skip-download \
    --skip-preprocess \
    --skip-varcall \
    --denovo-only

# Run only rare inherited analysis
bash /home/johan/pipeline/varcall/complete_autism_analysis.sh \
    --skip-download \
    --skip-preprocess \
    --skip-varcall \
    --inherited-only
```

## Trio Information

**Family ID**: FAM001

| Sample | SRA ID | Role | Phenotype | Details |
|--------|--------|------|-----------|---------|
| Proband | SRR8697636 | Child | Affected | Autistic child, IQ=60, ADI=42 |
| Mother | SRR8697627 | Parent | Unaffected | - |
| Father | SRR8697645 | Parent | Unaffected | - |

## Pipeline Workflow

```
┌─────────────────────────────────────────────────────────────┐
│                  PHASE 1: Data Preparation                  │
├─────────────────────────────────────────────────────────────┤
│ 1. Download FASTQ from ENA                                  │
│    → proband_1.fastq.gz, proband_2.fastq.gz                │
│    → mother_1.fastq.gz, mother_2.fastq.gz                  │
│    → father_1.fastq.gz, father_2.fastq.gz                  │
│                                                              │
│ 2. Quality Control (FastQC)                                 │
│                                                              │
│ 3. Alignment (BWA-MEM → hg19)                              │
│    → proband_aligned.bam                                    │
│    → mother_aligned.bam                                     │
│    → father_aligned.bam                                     │
│                                                              │
│ 4. Mark Duplicates (GATK)                                   │
│    → *_dedup.bam                                            │
│                                                              │
│ 5. Base Quality Recalibration (BQSR)                        │
│    → proband_final.bam                                      │
│    → mother_final.bam                                       │
│    → father_final.bam                                       │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                  PHASE 2: Variant Calling                   │
├─────────────────────────────────────────────────────────────┤
│ 1. HaplotypeCaller (GVCF mode)                              │
│    → proband.g.vcf.gz                                       │
│    → mother.g.vcf.gz                                        │
│    → father.g.vcf.gz                                        │
│                                                              │
│ 2. Joint Genotyping                                         │
│    → trio_cohort.vcf.gz                                     │
│                                                              │
│ 3. Quality Filtering                                        │
│    → trio_filtered.vcf.gz                                   │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│              PHASE 3: De Novo Variant Analysis              │
├─────────────────────────────────────────────────────────────┤
│ 1. Create Pedigree File                                     │
│                                                              │
│ 2. Annotate with PossibleDeNovo                             │
│                                                              │
│ 3. Filter for De Novo Pattern                               │
│    Proband: 0/1 (heterozygous)                             │
│    Mother:  0/0 (homozygous reference)                     │
│    Father:  0/0 (homozygous reference)                     │
│                                                              │
│ 4. Quality Filtering (GQ≥30, DP≥20x)                       │
│    → denovo_high_confidence.vcf.gz                          │
│                                                              │
│ 5. Functional Annotation (SnpEff)                           │
│    → denovo_annotated.vcf.gz                                │
│                                                              │
│ 6. Generate Reports                                         │
│    → denovo_analysis_report.txt                             │
│    → denovo_variants_table.tsv                              │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│          PHASE 4: Rare Inherited Variant Analysis           │
├─────────────────────────────────────────────────────────────┤
│ 1. Quality Filtering (Phred≥30, DP≥20x)                    │
│    → trio_high_quality.vcf.gz                               │
│                                                              │
│ 2. Functional Annotation (ANNOVAR)                          │
│    → trio_high_quality.hg19_multianno.vcf                  │
│                                                              │
│ 3. Filter for Rare Variants (MAF≤0.01)                     │
│    Databases: 1000G, ESP6500, ExAC                         │
│    → trio_rare.vcf.gz                                       │
│                                                              │
│ 4. Filter for Functional Impact                             │
│    Nonsynonymous: missense, nonsense, splice, frameshift   │
│    → trio_rare_functional.vcf.gz                            │
│                                                              │
│ 5. Identify Inherited Pattern                               │
│    Proband has variant AND ≥1 parent has variant           │
│    → rare_inherited_variants.vcf.gz                         │
│                                                              │
│ 6. Categorize by Inheritance                                │
│    → inherited_maternal.vcf.gz                              │
│    → inherited_paternal.vcf.gz                              │
│    → inherited_compound.vcf.gz (both parents)              │
│    → inherited_recessive.vcf.gz (homozygous)               │
│                                                              │
│ 7. Generate Reports                                         │
│    → rare_inherited_analysis_report.txt                     │
│    → rare_inherited_variants_table.tsv                      │
│    → genes_with_rare_variants.tsv                           │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│              PHASE 5: Comprehensive Reporting               │
├─────────────────────────────────────────────────────────────┤
│ → complete_analysis_summary.txt                             │
│ → Combined statistics and recommendations                   │
└─────────────────────────────────────────────────────────────┘
```

## Output Files

### Main Results Directory
```
/home/johan/output/autism/variants/
```

### De Novo Variants
```
variants/
├── denovo_variants.vcf.gz              # All de novo candidates
├── denovo_high_confidence.vcf.gz       # High-confidence only (recommended)
├── denovo_annotated.vcf.gz             # With functional annotations
├── denovo_variants_table.tsv           # Tab-separated variant table
└── denovo_analysis_report.txt          # Summary report
```

### Rare Inherited Variants
```
variants/rare_inherited/
├── rare_inherited_variants.vcf.gz      # All rare inherited variants
├── inherited_maternal.vcf.gz           # From mother
├── inherited_paternal.vcf.gz           # From father
├── inherited_compound.vcf.gz           # From both parents
├── inherited_recessive.vcf.gz          # Homozygous recessive
├── rare_inherited_variants_table.tsv   # Tab-separated table
├── genes_with_rare_variants.tsv        # Gene-level summary
└── rare_inherited_analysis_report.txt  # Summary report
```

### Supporting Files
```
variants/
├── trio_cohort.vcf.gz                  # All called variants
├── trio_filtered.vcf.gz                # Quality-filtered variants
├── trio.ped                            # Pedigree file
└── complete_analysis_summary.txt       # Comprehensive report
```

## Interpreting Results

### De Novo Variants

De novo mutations are the **highest priority** for autism genetics because:
- They occur spontaneously in the child
- Not inherited from parents
- Higher pathogenic potential
- Direct link to sporadic autism cases

**Typical counts**: 50-100 de novo variants per exome

**How to prioritize**:
1. **Loss-of-function** (stop-gain, frameshift) - Highest priority
2. **Splice-site** variants - High priority
3. **Missense in conserved regions** - Medium priority
4. **Check against autism databases**:
   - SFARI Gene (https://gene.sfari.org/)
   - AutDB (http://autism.mindspec.org/autdb/)
   - ClinVar for known pathogenic variants

### Rare Inherited Variants

These contribute to autism risk through:
- Compound heterozygous effects (both parents contribute)
- Homozygous recessive inheritance
- Additive effects of multiple rare variants
- Gene burden (multiple variants in same gene)

**Typical counts**: 100-500 rare inherited variants

**How to prioritize**:
1. **Compound heterozygous** - Two different variants in same gene
2. **Homozygous recessive** - Both alleles affected
3. **Gene burden** - Multiple rare variants in one gene
4. **Known autism genes**

## Viewing Variants

### Using bcftools

```bash
# View de novo variants
bcftools view /home/johan/output/autism/variants/denovo_high_confidence.vcf.gz | less

# Extract specific chromosome
bcftools view -r chr16 denovo_high_confidence.vcf.gz

# Count variants
bcftools view -H denovo_high_confidence.vcf.gz | wc -l

# Get variant at specific position
bcftools view -r chr16:29650000-29700000 denovo_high_confidence.vcf.gz
```

### Using the TSV tables

```bash
# View in spreadsheet format
column -t denovo_variants_table.tsv | less -S

# Search for specific gene
grep "SCN2A" rare_inherited_variants_table.tsv

# Sort by quality
sort -k5 -nr denovo_variants_table.tsv | head -20
```

### Using IGV (Integrative Genomics Viewer)

1. Load BAM files: `proband_final.bam`, `mother_final.bam`, `father_final.bam`
2. Load VCF file: `denovo_high_confidence.vcf.gz`
3. Navigate to variants of interest
4. Visually inspect read support

## Quality Control Checkpoints

### After Preprocessing
```bash
# Check alignment statistics
samtools flagstat /home/johan/output/autism/proband_final.bam

# Check coverage
samtools depth proband_final.bam | awk '{sum+=$3; cnt++} END {print "Mean coverage:", sum/cnt "x"}'

# Expected: 30-60x coverage for exome
```

### After Variant Calling
```bash
# Check variant counts
bcftools stats trio_cohort.vcf.gz | grep "number of records"

# Expected: 50,000-100,000 variants in trio
```

### Relationship Verification
```bash
# Check PLINK output (if available)
cat /home/johan/output/autism/variants/rare_inherited/relatedness_check.txt

# Expected PI_HAT:
#   Parent-child: ~0.5
#   Unrelated: ~0.0
```

## Troubleshooting

### Issue: Low Coverage
**Symptom**: Mean coverage <20x
**Solution**: This is exome data; 20-30x is acceptable for high-quality variants. If <10x, review alignment quality.

### Issue: Too Many De Novo Variants
**Symptom**: >200 high-confidence de novo variants
**Solution**: Increase quality thresholds (GQ≥40, DP≥30) or check for sample mix-up

### Issue: No Variants After Filtering
**Symptom**: 0 rare inherited variants
**Solution**:
- Check population database annotations
- Verify ANNOVAR installation
- Try relaxing MAF threshold to 0.05

### Issue: ANNOVAR Not Found
**Symptom**: "ANNOVAR not found" warning
**Solution**:
```bash
# Download ANNOVAR from:
# https://annovar.openbioinformatics.org/

# Install to: /home/johan/tool/annovar/
# Download databases:
cd /home/johan/tool/annovar/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad_exome humandb/
```

## Advanced Analysis

### Gene Burden Analysis

Identify genes with multiple rare variants:

```bash
# Genes with >1 rare inherited variant
cut -f6 rare_inherited_variants_table.tsv | sort | uniq -c | sort -rn | head -20
```

### Compound Heterozygous Detection

```bash
# Variants where proband is 0/1 and each parent contributes one allele
bcftools view inherited_compound.vcf.gz
```

### Pathway Analysis

Export gene lists for pathway enrichment:

```bash
# Extract genes with de novo variants
bcftools query -f '%INFO/Gene.refGene\n' denovo_high_confidence.vcf.gz | \
    sort -u > denovo_genes.txt

# Use with DAVID, KEGG, or Reactome
```

## Citation

If using this pipeline for publication, please cite:

- **GATK**: McKenna et al. (2010) The Genome Analysis Toolkit
- **BWA**: Li & Durbin (2009) Fast and accurate short read alignment
- **bcftools**: Li et al. (2011) The Sequence Alignment/Map format
- **ANNOVAR**: Wang et al. (2010) ANNOVAR: functional annotation of genetic variants

## Support Files Required

### Reference Genome
- `hg19.fa` (with BWA index)

### Known Sites
- `dbsnp_138.hg19.vcf.gz`
- `Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz`

### Population Databases (for ANNOVAR)
- `1000g2015aug_all`
- `esp6500siv2_all`
- `exac03`
- `gnomad_exome`

All should be in: `/home/johan/johan/johan/reference/human/`

## Contact & Questions

For issues or questions about this pipeline:
- Check log files in `/home/johan/output/autism/`
- Review error messages in terminal output
- Verify all dependencies are installed

## Workflow Summary

1. **Download**: 10 minutes
2. **Preprocessing**: 4-6 hours (alignment, deduplication, BQSR)
3. **Variant Calling**: 2-4 hours (HaplotypeCaller, joint genotyping)
4. **De Novo Analysis**: 30 minutes
5. **Rare Inherited Analysis**: 1-2 hours (annotation, filtering)

**Total Time**: 8-12 hours for complete analysis

**Disk Space**: ~100 GB (including BAM files and intermediate files)

## Next Steps After Analysis

1. **Review Reports**
   - `complete_analysis_summary.txt`
   - `denovo_analysis_report.txt`
   - `rare_inherited_analysis_report.txt`

2. **Prioritize Variants**
   - Focus on de novo loss-of-function
   - Check known autism genes
   - Review compound heterozygous variants

3. **Validate**
   - Sanger sequencing for top candidates
   - Check in IGV
   - Verify inheritance pattern

4. **Functional Studies**
   - Gene expression analysis
   - Protein modeling
   - Literature review
