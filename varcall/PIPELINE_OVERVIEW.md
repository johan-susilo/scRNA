# Autism Trio Variant Calling Pipeline

**Complete end-to-end pipeline for identifying autism-causing variants**

---

## Overview

This pipeline performs:
1. Alignment to UCSC hg19 reference
2. Variant calling with GATK HaplotypeCaller
3. De novo variant detection
4. **SnpEff functional annotation** (integrated)
5. **Automatic extraction of loss-of-function variants**
6. **Automatic identification of autism gene variants**
7. **TOP_CANDIDATES.txt generation**

**Result:** Complete analysis from BAM files to prioritized candidate list in a single run.

---

## Quick Start

### Prerequisites
- Final BAM files: `proband_final_ucsc.bam`, `mother_final_ucsc.bam`, `father_final_ucsc.bam`
- Conda environment: `varcall` with GATK, bcftools, and SnpEff
- Reference: UCSC hg19 (`hg19_ucsc.fa`)
- Autism gene list: `known_autism_genes.txt`

### Run Pipeline

```bash
cd /home/johan/pipeline/varcall
bash varcall_trio_ucsc.sh
```

**That's it!** The pipeline will run all 11 steps automatically.

---

## Pipeline Steps

### STEP 0: Pedigree File Creation
- Creates trio pedigree file for de novo detection

### STEP 1: Preprocessing Check
- Verifies all BAM files exist

### STEP 2: Variant Calling (GATK HaplotypeCaller)
- Calls variants per sample in GVCF mode
- Runtime: ~10-14 hours per sample

### STEP 3: Joint Genotyping
- Combines GVCFs across trio
- Runtime: ~4-6 hours

### STEP 4: Quality Filtering
- Hard filters for SNPs and INDELs
- Runtime: ~1-2 hours

### STEP 5: De Novo Detection
- Identifies variants present in child, absent in parents
- Runtime: ~30 minutes

### STEP 6: High-Confidence Filtering
- Filters for PASS variants with no missing genotypes
- Runtime: ~5 minutes

### STEP 7: SnpEff Annotation (NEW!)
- Functional annotation with gene names, effects, impacts
- Automatic installation if needed
- Runtime: ~15-30 minutes
- **Output:** `annotation/denovo_annotated.vcf.gz`
- **HTML Report:** `annotation/snpeff_summary.html`

### STEP 8: Loss-of-Function Extraction (NEW!)
- Automatically extracts all HIGH impact variants
- **Output:** `annotation/lof_variants.txt`

### STEP 9: Autism Gene Identification (NEW!)
- Matches variants to known autism genes
- **Output:** `annotation/autism_gene_variants.txt`

### STEP 10: Top Candidates List (NEW!)
- Prioritizes LoF variants in autism genes
- **Output:** `annotation/TOP_CANDIDATES.txt` ⭐

### STEP 11: Final Report
- Comprehensive summary report
- **Output:** `annotation/FINAL_ANALYSIS_REPORT.txt`

---

## Output Files

### Directory Structure

```
/home/johan/output/autism/
├── variants/
│   ├── trio_joint_called.vcf.gz          # All variants
│   ├── denovo_high_confidence.vcf.gz     # De novo variants
│   ├── proband.g.vcf.gz                  # Individual GVCFs
│   ├── mother.g.vcf.gz
│   └── father.g.vcf.gz
│
└── annotation/                            # NEW! SnpEff results
    ├── denovo_annotated.vcf.gz           # Annotated variants
    ├── snpeff_summary.html               # ⭐ Visual report
    ├── lof_variants.txt                  # ⭐ Loss-of-function list
    ├── autism_gene_variants.txt          # ⭐ Autism gene variants
    ├── TOP_CANDIDATES.txt                # ⭐ Prioritized list
    └── FINAL_ANALYSIS_REPORT.txt         # Summary report
```

### Key Files to Review

**1. TOP_CANDIDATES.txt** (Most Important)
- Prioritized list of candidates
- LoF variants in autism genes highlighted
- Ready for clinical review

**2. snpeff_summary.html**
- Visual charts and statistics
- Open in web browser
- Interactive exploration

**3. lof_variants.txt**
- All loss-of-function variants
- Gene, position, effect, protein change
- ~400-500 variants expected

**4. autism_gene_variants.txt**
- All variants in known autism genes
- Annotated with impact and effect

---

## Runtime Estimates

| Step | Runtime |
|------|---------|
| Variant calling (3 samples) | 30-42 hours |
| Joint genotyping | 4-6 hours |
| Quality filtering | 1-2 hours |
| De novo detection | 30 minutes |
| **SnpEff annotation** | **15-30 minutes** |
| **Extraction & analysis** | **5-10 minutes** |
| **Total** | **~36-51 hours** |

**Note:** Most time is spent on variant calling. Once VCFs exist, the pipeline skips those steps.

---

## What's New (SnpEff Integration)

### Previous Pipeline (Manual)
1. Run `varcall_trio_ucsc.sh` → get VCF
2. Manually run SnpEff
3. Manually extract HIGH impact variants
4. Manually match with autism genes
5. Manually create candidate list

**Problem:** Multiple manual steps, error-prone, time-consuming

### New Pipeline (Automated)
1. Run `varcall_trio_ucsc.sh` → **Everything done automatically!**
2. Review `TOP_CANDIDATES.txt`

**Benefit:** Single command, complete analysis, ready for validation

---

## Example Output

### TOP_CANDIDATES.txt Sample

```
==========================================
TOP CANDIDATE VARIANTS FOR AUTISM
==========================================

Generated: 2026-01-01
Total de novo variants: 75285

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PRIORITY 1: LOSS-OF-FUNCTION VARIANTS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

⭐ LOSS-OF-FUNCTION VARIANTS IN AUTISM GENES:

Gene    Chrom  Pos       Ref  Alt  Qual    Effect        HGVS_Protein
SLC6A1  chr3   11061962  A    G    773.29  start_lost    p.Met1?
SYNE2   chr14  64560092  G    A    450.23  stop_gained   p.Trp4001*

All loss-of-function variants (top 20):
[...]
```

---

## Troubleshooting

### SnpEff not found
**Solution:** The pipeline will automatically install SnpEff if needed
```bash
conda install -y -c bioconda snpeff=5.1
```

### Out of memory error
**Solution:** Increase memory allocation in the script
```bash
# Edit line 336-337 in varcall_trio_ucsc.sh
-Xmx16g  # Change to -Xmx32g if you have more RAM
```

### Annotation is slow
**Normal:** Annotating 70,000+ variants takes 15-30 minutes
**Monitor:** Check `annotation/snpeff.log` for progress

### No autism genes found
**Check:** Verify `known_autism_genes.txt` exists
```bash
ls -lh /home/johan/output/autism/known_autism_genes.txt
```

---

## Configuration

### Edit Script Settings

Open `varcall_trio_ucsc.sh` and modify:

```bash
# Line 22: Number of threads
THREADS=60

# Line 23-27: Directories
OUTPUT_DIR="${HOME}/output/autism"
REF_DIR="${HOME}/johan/johan/johan/reference/human"

# Line 42: Autism gene list
AUTISM_GENES="${OUTPUT_DIR}/known_autism_genes.txt"
```

---

## Re-running the Pipeline

### Full Re-run (from scratch)
```bash
# Remove all output
rm -rf /home/johan/output/autism/variants/*
rm -rf /home/johan/output/autism/annotation/*

# Run pipeline
bash varcall_trio_ucsc.sh
```

### Re-run Only SnpEff (keep VCFs)
```bash
# Remove only annotation results
rm -rf /home/johan/output/autism/annotation/*

# Run pipeline (will skip variant calling)
bash varcall_trio_ucsc.sh
```

### Skip Existing Steps
The pipeline automatically skips completed steps:
- If GVCFs exist → skips HaplotypeCaller
- If annotated VCF exists → skips SnpEff
- If outputs exist → displays results

---

## Next Steps After Pipeline

1. **Review TOP_CANDIDATES.txt**
   - Identify highest priority variants
   - Focus on LoF in autism genes

2. **Open SnpEff HTML Report**
   ```bash
   firefox /home/johan/output/autism/annotation/snpeff_summary.html
   ```

3. **Check Population Frequencies**
   - Go to https://gnomad.broadinstitute.org/
   - Search top candidate positions
   - Filter out common variants (AF > 1%)

4. **Sanger Sequencing Validation**
   - Design primers for top candidates
   - Validate in all 3 family members
   - Confirm de novo inheritance

5. **Clinical Consultation**
   - Clinical geneticist review
   - ACMG classification
   - Genetic counseling

---

## Summary

**What you get:**
- Complete variant calling
- Functional annotation
- Loss-of-function identification
- Autism gene matching
- Prioritized candidate list
- Visual HTML report

**What you do:**
- Run one command: `bash varcall_trio_ucsc.sh`
- Review `TOP_CANDIDATES.txt`
- Validate top candidates

**No more manual steps!**

---

**Last updated:** January 1, 2026
**Version:** 2.0 (SnpEff integrated)
