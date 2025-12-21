# 🚀 Complete Pipeline: FAM92 De Novo Analysis

## Current Date/Time
**2025-12-18 20:30 (Thursday)**

## System Resources
- **CPUs**: 64 cores available
- **Strategy**: Use 60 cores for processing (leave 4 for system)

---

## 📋 Corrected Trio Information

| Role | SRA ID | Phenotype |
|------|--------|-----------|
| **Proband** | SRR8697636 | Autistic child (IQ=60, ADI=42) |
| **Mother** | SRR8697627 | Unaffected |
| **Father** | SRR8697645 | Unaffected |

---

## 🎯 Complete 2-Step Pipeline

### **STEP 1: Preprocessing (NEW Optimized Script)**

Run the optimized preprocessing script:

```bash
cd /home/johan/pipeline/varcall

# Activate conda environment
conda activate varcall

# Run optimized preprocessing (6-10 hours)
./preprocess_trio_optimized.sh 2>&1 | tee preprocess_trio.log
```

**What this does:**
- ✅ Downloads all 3 SRA files (if needed)
- ✅ Converts to FASTQ using 60 threads (fast!)
- ✅ Compresses with `pigz` (parallel gzip)
- ✅ Runs FastQC quality control
- ✅ BWA-MEM alignment (60 threads per sample)
- ✅ Mark duplicates with GATK
- ✅ Base quality recalibration (BQSR)
- ✅ Outputs: `proband_final.bam`, `mother_final.bam`, `father_final.bam`

**Time estimate**: 6-10 hours (much faster with 60 cores!)

---

### **STEP 2: Variant Calling & Analysis**

After preprocessing completes, run the de novo analysis:

```bash
cd /home/johan/pipeline/varcall

# Make sure BAMs are in the right place
ls ~/output/autism/*_final.bam

# Run de novo analysis (4-6 hours)
./denovo_fam92_analysis.sh 2>&1 | tee denovo_analysis.log
```

**What this does:**
- ✅ GATK HaplotypeCaller (per sample GVCF)
- ✅ Joint genotyping → `FAM92_trio.vcf.gz`
- ✅ Filter de novo variants (proband=variant, parents=normal)
- ✅ Annotate with VEP (SIFT, PolyPhen, CADD)
- ✅ Extract FAM92 variants
- ✅ Generate clinical report

**Time estimate**: 4-6 hours

---

## 🔥 Key Improvements Over Original Script

### Original Script Issues Fixed:
❌ Mixed up SRA IDs (used SRR8697637 instead of SRR8697636)
❌ Used slow single-threaded `fastq-dump`
❌ No alignment or quality control
❌ Inconsistent sample naming

### New Script Benefits:
✅ **Correct trio IDs** (SRR8697636, SRR8697627, SRR8697645)
✅ **60-thread processing** (was: single-threaded)
✅ **Parallel gzip** with `pigz` (10x faster)
✅ **Complete pipeline** (FASTQ → aligned BAM → variant calling)
✅ **Quality metrics** at each step
✅ **Resume capability** (skips completed steps)

---

## ⏱️ Total Time Estimate

| Stage | Time | Bottleneck |
|-------|------|-----------|
| Download SRA | 30-60 min | Network speed |
| FASTQ conversion | 1-2 hours | Disk I/O |
| Alignment (BWA) | 3-4 hours | CPU + I/O |
| Mark duplicates | 1 hour | Memory |
| BQSR | 1-2 hours | I/O |
| Variant calling | 3-4 hours | CPU |
| Annotation (VEP) | 30 min | I/O |
| **TOTAL** | **10-15 hours** | Run overnight! |

---

## 💾 Disk Space Requirements

Estimate per sample:
- SRA file: ~2-3 GB
- FASTQ (gzipped): ~5-8 GB
- Aligned BAM: ~20-30 GB
- Final BAM: ~15-25 GB

**Total for trio**: ~120-180 GB

---

## 📂 Output Structure

```
~/output/autism/
├── SRR8697636/
│   └── SRR8697636.sra
├── SRR8697627/
│   └── SRR8697627.sra
├── SRR8697645/
│   └── SRR8697645.sra
├── proband_1.fastq.gz
├── proband_2.fastq.gz
├── mother_1.fastq.gz
├── mother_2.fastq.gz
├── father_1.fastq.gz
├── father_2.fastq.gz
├── proband_final.bam     ← Use for variant calling
├── mother_final.bam      ← Use for variant calling
├── father_final.bam      ← Use for variant calling
├── qc/                   (FastQC reports)
└── trio_analysis/        (Variant calling output)
    ├── FAM92_trio.vcf.gz
    ├── denovo_variants_with_header.vcf.gz
    ├── FAM92_denovo.txt  ⭐
    └── clinical_report.txt  ⭐
```

---

## 🔍 Monitoring Progress

### Check preprocessing status:
```bash
# Watch the log
tail -f ~/pipeline/varcall/preprocess_trio.log

# Check which step is running
ps aux | grep -E "bwa|gatk|fasterq"

# Check BAM files
ls -lht ~/output/autism/*.bam
```

### Check variant calling status:
```bash
# Watch the log
tail -f ~/pipeline/varcall/denovo_analysis.log

# Check for output VCF
ls -lht ~/output/autism/trio_analysis/*.vcf.gz
```

---

## ✅ Final Checklist

Before running:
- [ ] Verify 64 CPUs available: `nproc`
- [ ] Check disk space: `df -h ~`
- [ ] Conda varcall environment exists: `conda env list | grep varcall`
- [ ] Reference genome indexed: `ls ~/johan/johan/reference/human/hg19.fa.bwt`
- [ ] dbSNP available: `ls ~/johan/johan/reference/human/dbsnp_138.hg19.vcf.gz`

During run:
- [ ] Monitor CPU usage: `htop`
- [ ] Monitor disk space: `watch df -h`
- [ ] Check logs for errors: `tail -f *.log`

After completion:
- [ ] Verify BAM files: `samtools view -H *_final.bam`
- [ ] Check variant count: `bcftools view -H FAM92_trio.vcf.gz | wc -l`
- [ ] Read clinical report: `cat ~/output/autism/trio_analysis/clinical_report.txt`

---

## 🆘 Troubleshooting

### Problem: "pigz: command not found"
```bash
conda activate varcall
conda install -c conda-forge pigz
```

### Problem: Out of memory during alignment
Edit script, reduce threads to 30:
```bash
THREADS=30
```

### Problem: VEP cache download slow
Use mirror or download separately:
```bash
wget ftp://ftp.ensembl.org/pub/release-104/variation/indexed_vep_cache/homo_sapiens_vep_104_GRCh37.tar.gz
tar xzf homo_sapiens_vep_104_GRCh37.tar.gz -C ~/.vep/
```

---

## 🎯 Expected Answer to Main Question

**Question**: "Do de novo variants in FAM92 explain the proband's severe ASD phenotype (IQ=60)?"

**Possible outcomes**:

### ✅ Outcome 1: FAM92 variant found
```
Found: FAM92A1 chrX:12345678 C>T
       Consequence: missense
       SIFT: 0.01 (deleterious)
       PolyPhen: 0.95 (probably damaging)
       CADD: 28 (pathogenic)

→ YES, FAM92 variant likely explains severe ASD
→ Validate with Sanger sequencing
```

### ❌ Outcome 2: No FAM92 variants
```
De novo variants in FAM92: 0

→ NO, FAM92 does NOT explain the phenotype
→ Investigate other ASD genes (CHD8, SCN2A, SHANK3, etc.)
```

---

**Ready to run!** Start with preprocessing, then variant calling.
