# Autism Trio Variant Calling Pipeline

**Complete end-to-end pipeline with integrated SnpEff annotation**

---

## 🎯 What This Pipeline Does

Runs the complete analysis from aligned BAM files to prioritized candidate variants:

1. ✅ GATK variant calling (HaplotypeCaller)
2. ✅ Joint genotyping across trio
3. ✅ De novo variant detection
4. ✅ **SnpEff functional annotation** (automated)
5. ✅ **Loss-of-function extraction** (automated)
6. ✅ **Autism gene matching** (automated)
7. ✅ **TOP_CANDIDATES.txt generation** (automated)

**One command → Complete results!**

---

## 📁 Files in This Directory

| File | Description |
|------|-------------|
| **varcall_trio_ucsc.sh** | ⭐ Main pipeline (with SnpEff integrated) |
| **realign_with_ucsc_hg19.sh** | Re-alignment script (run first if needed) |
| **run_complete_pipeline_ucsc.sh** | Master orchestrator (runs both) |
| **PIPELINE_OVERVIEW.md** | Detailed documentation |
| **FINAL_PIPELINE_SUMMARY.md** | Technical documentation |
| **README.md** | This file |

---

## 🚀 Quick Start

### Option 1: Variant Calling Only (if BAMs exist)

```bash
cd /home/johan/pipeline/varcall
bash varcall_trio_ucsc.sh
```

**Requirements:**
- `proband_final_ucsc.bam`
- `mother_final_ucsc.bam`
- `father_final_ucsc.bam`

**Runtime:** ~36-48 hours (first run), ~30 minutes (if VCFs exist)

### Option 2: Complete Pipeline (from FASTQ)

```bash
cd /home/johan/pipeline/varcall
bash run_complete_pipeline_ucsc.sh
```

**Requirements:**
- FASTQ files in `/home/johan/output/autism/`

**Runtime:** ~3-4 days (alignment + variant calling)

---

## 📊 Output Files

### Main Results (in `/home/johan/output/autism/`)

**Variant Calling:**
- `variants/trio_joint_called.vcf.gz` - All variants
- `variants/denovo_high_confidence.vcf.gz` - De novo variants

**SnpEff Annotation (NEW!):**
- `annotation/TOP_CANDIDATES.txt` ⭐ **REVIEW THIS FIRST**
- `annotation/snpeff_summary.html` ⭐ Open in browser
- `annotation/lof_variants.txt` - Loss-of-function list
- `annotation/autism_gene_variants.txt` - Autism gene details
- `annotation/denovo_annotated.vcf.gz` - Annotated VCF

---

## 🔄 What's Different (Version 2.0)

### Old Workflow (Manual)
1. Run variant calling script
2. Wait for results
3. **Manually run SnpEff**
4. **Manually extract HIGH impact**
5. **Manually match autism genes**
6. **Manually create candidate list**

### New Workflow (Integrated) ✅
1. Run `bash varcall_trio_ucsc.sh`
2. **Everything done automatically!**
3. Review `annotation/TOP_CANDIDATES.txt`

**Benefit:** Zero manual steps after pipeline launch

---

## 📋 Pipeline Steps (Automated)

| Step | What It Does | Output |
|------|--------------|--------|
| 1 | Pedigree file creation | `trio.ped` |
| 2 | GATK HaplotypeCaller | Individual GVCFs |
| 3 | Joint genotyping | `trio_joint_called.vcf.gz` |
| 4 | Quality filtering | `trio_filtered.vcf.gz` |
| 5 | De novo detection | `denovo_variants.vcf.gz` |
| 6 | High-confidence filtering | `denovo_high_confidence.vcf.gz` |
| **7** | **SnpEff annotation** | **denovo_annotated.vcf.gz** |
| **8** | **LoF extraction** | **lof_variants.txt** |
| **9** | **Autism gene matching** | **autism_gene_variants.txt** |
| **10** | **Candidate prioritization** | **TOP_CANDIDATES.txt** |
| 11 | Final report | `FINAL_ANALYSIS_REPORT.txt` |

**Steps 7-10 are NEW and fully automated!**

---

## ⏱️ Runtime

| Component | First Run | Re-run (if outputs exist) |
|-----------|-----------|---------------------------|
| Variant calling | 30-42 hours | Skipped |
| Joint genotyping | 4-6 hours | Skipped |
| Filtering | 1-2 hours | Skipped |
| De novo detection | 30 min | Skipped |
| **SnpEff annotation** | **15-30 min** | **Skipped** |
| **Analysis** | **5-10 min** | **Instant** |
| **Total** | **36-51 hours** | **< 1 minute** |

**Smart skipping:** If outputs exist, steps are skipped automatically.

---

## 🎯 What You Get

### Immediate Results
- List of loss-of-function variants
- Variants in autism genes
- Prioritized candidate list
- Visual HTML report

### Example TOP_CANDIDATES.txt

```
⭐ LOSS-OF-FUNCTION VARIANTS IN AUTISM GENES:

Gene    Chrom  Pos       Effect        HGVS_Protein
SLC6A1  chr3   11061962  start_lost    p.Met1?
SYNE2   chr14  64560092  stop_gained   p.Trp4001*
```

---

## 🔧 Configuration

### Default Settings (in varcall_trio_ucsc.sh)

```bash
THREADS=60                                    # Parallel processing
OUTPUT_DIR="${HOME}/output/autism"            # Results location
REFERENCE="${REF_DIR}/hg19_ucsc.fa"          # UCSC hg19 reference
AUTISM_GENES="${OUTPUT_DIR}/known_autism_genes.txt"  # Gene list
```

### To Modify

Edit `varcall_trio_ucsc.sh` lines 22-42

---

## 📝 Next Steps After Pipeline

1. **Review Results**
   ```bash
   cat /home/johan/output/autism/annotation/TOP_CANDIDATES.txt
   ```

2. **Open HTML Report**
   ```bash
   firefox /home/johan/output/autism/annotation/snpeff_summary.html
   ```

3. **Check Population Frequencies**
   - Go to https://gnomad.broadinstitute.org/
   - Search top candidates
   - Filter common variants (AF > 1%)

4. **Validate Candidates**
   - Sanger sequencing
   - All 3 family members
   - Confirm de novo status

5. **Clinical Review**
   - Clinical geneticist
   - ACMG classification
   - Genetic counseling

---

## 🔍 Troubleshooting

### Pipeline fails at SnpEff
**Solution:** SnpEff will auto-install if missing
```bash
# Or install manually:
conda activate varcall
conda install -y -c bioconda snpeff=5.1
```

### Out of memory
**Solution:** Increase memory in script (line 336)
```bash
-Xmx16g  → -Xmx32g
```

### Need to re-run SnpEff only
```bash
# Remove annotation directory
rm -rf /home/johan/output/autism/annotation/

# Re-run pipeline (will skip variant calling)
bash varcall_trio_ucsc.sh
```

---

## 📚 Documentation

- **PIPELINE_OVERVIEW.md** - Complete guide
- **FINAL_PIPELINE_SUMMARY.md** - Technical details
- **varcall_trio_ucsc.sh** - Heavily commented code

---

## ✅ Summary

**What:** Complete trio variant calling + SnpEff annotation
**Input:** BAM files (or FASTQ files)
**Output:** Prioritized candidate list + annotated VCF
**Runtime:** 36-51 hours (first run), <1 minute (re-run)
**Manual steps:** Zero (fully automated)

**Just run:** `bash varcall_trio_ucsc.sh`

---

**Version:** 2.0 (SnpEff Integrated)
**Last updated:** January 1, 2026
**Status:** Production ready
