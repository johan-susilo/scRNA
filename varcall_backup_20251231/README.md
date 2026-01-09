# Autism Trio Variant Calling Pipeline

Complete pipeline for identifying de novo and rare inherited variants in autism exome trios.

## Quick Start

```bash
# Run complete analysis (recommended)
bash /home/johan/pipeline/varcall/complete_autism_analysis.sh
```

## What This Pipeline Does

This pipeline analyzes exome sequencing data from an autism trio (affected child + both parents) to identify:

1. **De Novo Variants** - New mutations in the child (not in parents)
2. **Rare Inherited Variants** - Rare functional variants from parents

Following published autism genetics methodology with strict quality standards.

## Scientific Approach

### De Novo Variant Detection
- **Pattern**: Proband heterozygous (0/1), both parents homozygous reference (0/0)
- **Quality**: Phred ≥30, Coverage ≥20x in all samples
- **Validation**: High-confidence filtering, functional annotation

### Rare Inherited Variant Detection
- **Rarity**: MAF ≤0.01 (1% in 1000G, ESP6500, ExAC)
- **Function**: Nonsynonymous only (missense, nonsense, splice-site, frameshift)
- **Inheritance**: Present in child AND at least one parent
- **Patterns**: Maternal, paternal, compound heterozygous, recessive

## Pipeline Components

| Script | Purpose | Runtime |
|--------|---------|---------|
| `convert_sra_fixed.sh` | Download FASTQ from ENA | 10 min |
| `preprocess_trio_optimized.sh` | Align, deduplicate, BQSR | 4-6 hours |
| `complete_trio_pipeline.sh` | Variant calling + de novo | 2-4 hours |
| `rare_inherited_variants_pipeline.sh` | Rare inherited analysis | 1-2 hours |
| `complete_autism_analysis.sh` | **Master script (all steps)** | 8-12 hours |

## Trio Information

| Sample | SRA ID | Role | Phenotype |
|--------|--------|------|-----------|
| Proband | SRR8697636 | Child | Affected (autism, IQ=60) |
| Mother | SRR8697627 | Parent | Unaffected |
| Father | SRR8697645 | Parent | Unaffected |

## Output Structure

```
/home/johan/output/autism/
├── proband_1.fastq.gz                 # Raw sequencing data
├── proband_2.fastq.gz
├── mother_*.fastq.gz
├── father_*.fastq.gz
├── *_final.bam                        # Processed alignments
│
└── variants/                           # Analysis results
    ├── denovo_high_confidence.vcf.gz  # ★ De novo variants
    ├── denovo_variants_table.tsv      # ★ Tab-separated table
    ├── denovo_analysis_report.txt     # ★ Summary report
    │
    ├── rare_inherited/
    │   ├── rare_inherited_variants.vcf.gz      # ★ All rare inherited
    │   ├── inherited_maternal.vcf.gz           # From mother
    │   ├── inherited_paternal.vcf.gz           # From father
    │   ├── inherited_compound.vcf.gz           # Both parents
    │   ├── inherited_recessive.vcf.gz          # Homozygous
    │   ├── rare_inherited_variants_table.tsv   # ★ Tab table
    │   ├── genes_with_rare_variants.tsv        # ★ Gene summary
    │   └── rare_inherited_analysis_report.txt  # ★ Report
    │
    └── complete_analysis_summary.txt   # ★★ COMPREHENSIVE REPORT
```

**★ = Key files to review**
**★★ = Start here**

## Usage Examples

### Complete Analysis
```bash
# Run everything (first time)
bash /home/johan/pipeline/varcall/complete_autism_analysis.sh
```

### Step-by-Step
```bash
# Step 1: Download data
bash /home/johan/pipeline/varcall/convert_sra_fixed.sh

# Step 2: Preprocessing (alignment, QC)
bash /home/johan/pipeline/varcall/preprocess_trio_optimized.sh

# Step 3: De novo detection
bash /home/johan/pipeline/varcall/complete_trio_pipeline.sh

# Step 4: Rare inherited detection
bash /home/johan/pipeline/varcall/rare_inherited_variants_pipeline.sh
```

### Resume Analysis
```bash
# Already have BAM files? Skip to variant calling
bash /home/johan/pipeline/varcall/complete_autism_analysis.sh \
    --skip-download \
    --skip-preprocess

# Only de novo analysis
bash /home/johan/pipeline/varcall/complete_autism_analysis.sh \
    --skip-download \
    --skip-preprocess \
    --denovo-only

# Only inherited analysis
bash /home/johan/pipeline/varcall/complete_autism_analysis.sh \
    --skip-download \
    --skip-preprocess \
    --inherited-only
```

## Viewing Results

### Quick Summary
```bash
# Read the main report
cat /home/johan/output/autism/variants/complete_analysis_summary.txt
```

### De Novo Variants
```bash
# View VCF
bcftools view /home/johan/output/autism/variants/denovo_high_confidence.vcf.gz | less

# View table
column -t /home/johan/output/autism/variants/denovo_variants_table.tsv | less -S

# Count
bcftools view -H denovo_high_confidence.vcf.gz | wc -l
```

### Rare Inherited Variants
```bash
# View all rare inherited
cd /home/johan/output/autism/variants/rare_inherited/

# Tables
column -t rare_inherited_variants_table.tsv | less -S
column -t genes_with_rare_variants.tsv | less -S

# By inheritance pattern
bcftools view inherited_maternal.vcf.gz    # From mother
bcftools view inherited_paternal.vcf.gz    # From father
bcftools view inherited_compound.vcf.gz    # Both parents
bcftools view inherited_recessive.vcf.gz   # Homozygous recessive
```

### Search for Specific Gene
```bash
# De novo in specific gene
grep "SCN2A" /home/johan/output/autism/variants/denovo_variants_table.tsv

# Rare inherited in specific gene
grep "SCN2A" /home/johan/output/autism/variants/rare_inherited/rare_inherited_variants_table.tsv
```

## Quality Metrics

Expected values for high-quality exome data:

| Metric | Expected Value |
|--------|---------------|
| Mean coverage | 30-60x |
| Alignment rate | >95% |
| Duplicate rate | 10-25% |
| Total variants called | 50,000-100,000 |
| De novo variants | 50-100 |
| Rare inherited | 100-500 |

## Requirements

### Software (installed in conda environment)
- BWA-MEM
- GATK 4+
- samtools
- bcftools
- FastQC (optional)
- ANNOVAR (optional, for annotation)
- SnpEff (optional, for annotation)
- PLINK (optional, for relationship checking)

### Reference Files
Required in `/home/johan/johan/johan/reference/human/`:
- `hg19.fa` (with BWA index)
- `dbsnp_138.hg19.vcf.gz`
- `Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz`

### Computational Resources
- **CPU**: 60 cores recommended
- **RAM**: 64 GB minimum
- **Disk**: ~100 GB for complete analysis
- **Time**: 8-12 hours total

## Troubleshooting

### Common Issues

**Error: FASTQ files not found**
```bash
# Download FASTQ files first
bash /home/johan/pipeline/varcall/convert_sra_fixed.sh
```

**Error: BAM files not found**
```bash
# Run preprocessing first
bash /home/johan/pipeline/varcall/preprocess_trio_optimized.sh
```

**Error: Reference not found**
```bash
# Check reference path in script
ls -lh /home/johan/johan/johan/reference/human/hg19.fa
```

**Too many/few de novo variants**
- Expected: 50-100 high-confidence de novo variants
- Too many (>200): Sample mix-up or contamination
- Too few (<20): Very strict filtering or low coverage

### Getting Help

1. Check log files in `/home/johan/output/autism/`
2. Review error messages in terminal
3. Verify conda environment: `conda activate varcall`
4. Check disk space: `df -h /home/johan/output/`

## Documentation

- **Complete Guide**: [`AUTISM_ANALYSIS_GUIDE.md`](AUTISM_ANALYSIS_GUIDE.md)
- **ENA Download**: [`FASTQ_DOWNLOAD_GUIDE.md`](FASTQ_DOWNLOAD_GUIDE.md)
- **SRA vs ENA**: [`SRA_TOOLKIT_VS_ENA.md`](SRA_TOOLKIT_VS_ENA.md)

## Methodology References

This pipeline follows published best practices:

1. **GATK Best Practices** for germline variant calling
2. **Autism Sequencing Consortium** guidelines
3. **Quality standards**: Phred ≥30, DP ≥20x
4. **Rarity definition**: MAF ≤0.01 from population databases
5. **Functional impact**: Nonsynonymous variants only

## Pipeline Validation

Quality control checkpoints:

- ✓ Relationship verification (PLINK IBD)
- ✓ Coverage validation (≥20x in all samples)
- ✓ Variant quality filtering (Phred ≥30)
- ✓ Population frequency filtering (MAF ≤0.01)
- ✓ Functional impact assessment
- ✓ Inheritance pattern verification

## Citation

If you use this pipeline for research, please cite:

- **BWA**: Li H. and Durbin R. (2009) Bioinformatics
- **GATK**: McKenna A. et al. (2010) Genome Research
- **SAMtools**: Li H. et al. (2009) Bioinformatics
- **bcftools**: Danecek P. et al. (2021) GigaScience

## Contact

For questions or issues:
- Review documentation in this directory
- Check log files for error messages
- Verify all dependencies are installed

---

**Last Updated**: December 2025
**Pipeline Version**: 1.0
**Reference Genome**: hg19/GRCh37
