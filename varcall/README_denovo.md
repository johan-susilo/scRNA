# De Novo FAM92 Variant Analysis Pipeline

## Main Research Question
**"Do de novo variants in FAM92 explain the proband's severe ASD phenotype (IQ=60)?"**

## Pipeline Overview

```
[Trio Samples] → [Preprocessing] → [GATK Variant Calling] → [De Novo Filter]
    ↓
[FAM92 Gene Filter] → [Annotation] → [Clinical Interpretation]
```

## Prerequisites

### 1. Required Trio Data
You **MUST** have three samples:
- **Proband**: The affected child (SRR8697636 - already downloaded)
- **Father**: Unaffected father (SRA ID needed)
- **Mother**: Unaffected mother (SRA ID needed)

### 2. Find Parent SRA IDs
```bash
# Check the BioProject for SRR8697636
# Visit: https://www.ncbi.nlm.nih.gov/sra/SRR8697636
# Look for associated samples in the same study
```

### 3. Conda Environment
The script uses the `varcall` conda environment with:
- GATK4
- BWA
- samtools
- bcftools
- VEP (Variant Effect Predictor)

## Quick Start

### Step 1: Update Configuration
Edit the script and update parent SRA IDs:
```bash
nano /home/johan/pipeline/varcall/denovo_fam92_analysis.sh

# Change these lines:
FATHER_SRA="SRR8697637"   # ← Update with real SRA ID
MOTHER_SRA="SRR8697638"   # ← Update with real SRA ID
```

### Step 2: Make Executable
```bash
chmod +x /home/johan/pipeline/varcall/denovo_fam92_analysis.sh
```

### Step 3: Run Analysis
```bash
cd /home/johan/pipeline/varcall
./denovo_fam92_analysis.sh 2>&1 | tee denovo_analysis.log
```

## What the Pipeline Does

### 1. **Preprocessing** (per sample)
   - Download SRA files (if needed)
   - Convert to FASTQ
   - Align with BWA-MEM to hg19
   - Mark duplicates with GATK
   - Base quality recalibration

### 2. **GATK Variant Calling**
   - HaplotypeCaller per sample (GVCF mode)
   - Joint genotyping of trio
   - Output: `FAM92_trio.vcf.gz`

### 3. **De Novo Filtering**
   - Filter variants where:
     - Proband = heterozygous (0/1) or homozygous (1/1)
     - Father = reference (0/0)
     - Mother = reference (0/0)
   - Output: `denovo_variants_with_header.vcf.gz`

### 4. **Annotation with VEP**
   - Consequence prediction (missense, nonsense, etc.)
   - SIFT and PolyPhen scores
   - CADD pathogenicity scores
   - Gene and transcript information
   - Brain expression data

### 5. **FAM92 Gene Filtering**
   - Extract variants in:
     - FAM92A1 (Xq21.31)
     - FAM92A2 (Xq28)
     - FAM92B (19p13.3)
   - Output: `FAM92_denovo.txt`

### 6. **Clinical Report**
   - Summary of findings
   - Pathogenicity assessment
   - Genotype-phenotype correlation
   - Recommendations for validation
   - Output: `clinical_report.txt`

## Expected Output Files

```
/home/johan/output/autism/trio_analysis/
├── FAM92_trio.vcf.gz              # All variants in trio
├── denovo_variants_with_header.vcf.gz  # All de novo variants
├── denovo_annotated.vcf           # Annotated de novo variants
├── FAM92_denovo.txt               # FAM92-specific variants
├── clinical_report.txt            # Final clinical interpretation
├── proband_final.bam              # Processed BAM files
├── father_final.bam
└── mother_final.bam
```

## Interpreting Results

### Scenario 1: FAM92 Variants Found
✓ **Damaging variant in FAM92** → Likely causal
- Check SIFT < 0.05 (deleterious)
- Check PolyPhen > 0.5 (damaging)
- Check CADD > 20 (pathogenic)
- Validate with Sanger sequencing
- Functional studies recommended

### Scenario 2: No FAM92 Variants
✗ **No variants in FAM92** → Look elsewhere
- Check other ASD genes (SFARI database)
- Consider CNV analysis
- Whole exome sequencing for novel genes

## Troubleshooting

### Issue: Only have proband sample
**Solution**: Cannot do de novo analysis without parents
- Option A: Find trio data from SRA
- Option B: Do singleton analysis with population filtering

### Issue: VEP cache download fails
**Solution**: Manual download
```bash
cd ~/.vep
wget ftp://ftp.ensembl.org/pub/release-104/variation/indexed_vep_cache/homo_sapiens_vep_104_GRCh37.tar.gz
tar xzf homo_sapiens_vep_104_GRCh37.tar.gz
```

### Issue: Out of memory
**Solution**: Reduce threads or process by chromosome
```bash
# Add to GATK commands:
-L chr17  # FAM92B
-L chrX   # FAM92A1, FAM92A2
```

## Time Estimates
- Preprocessing (per sample): 2-4 hours
- Variant calling: 1-2 hours
- Annotation: 30 minutes
- **Total**: ~6-12 hours for complete trio

## References
- GATK Best Practices: https://gatk.broadinstitute.org/
- VEP Documentation: https://www.ensembl.org/vep
- SFARI Gene Database: https://gene.sfari.org/
- FAM92 genes: Search PubMed for recent literature

## Contact
For questions about this pipeline, check the GATK forum or VEP help pages.
