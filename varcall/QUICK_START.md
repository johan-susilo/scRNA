# 🧬 Quick Start: FAM92 De Novo Variant Analysis

## Main Question
**"Do de novo variants in FAM92 explain the proband's severe ASD phenotype (IQ=60)?"**

---

## ⚠️ IMPORTANT: You Need Trio Data!

Currently you have:
- ✅ **Proband** (affected child): `SRR8697636`

You still need:
- ❌ **Father** (unaffected): SRA ID unknown
- ❌ **Mother** (unaffected): SRA ID unknown

**Without parents, you CANNOT do de novo analysis!**

---

## 🔍 Step 1: Find Parent Samples

```bash
# Run the checker script
cd /home/johan/pipeline/varcall
./check_trio_data.sh
```

OR manually:
1. Go to: https://www.ncbi.nlm.nih.gov/sra/SRR8697636
2. Click on the **BioProject** link
3. Find the father and mother SRA accessions
4. They should be labeled as "trio" or "family" samples

---

## ✏️ Step 2: Update Configuration

Edit the main script:
```bash
nano /home/johan/pipeline/varcall/denovo_fam92_analysis.sh
```

Find these lines (around line 38-40) and update:
```bash
PROBAND_SRA="SRR8697636"  # Keep this
FATHER_SRA="SRR8697637"   # ← CHANGE TO REAL SRA ID
MOTHER_SRA="SRR8697638"   # ← CHANGE TO REAL SRA ID
```

---

## 🚀 Step 3: Run the Analysis

```bash
cd /home/johan/pipeline/varcall

# Activate conda environment
conda activate varcall

# Run the pipeline (takes 6-12 hours)
./denovo_fam92_analysis.sh 2>&1 | tee denovo_analysis.log
```

---

## 📊 What the Pipeline Does

```
┌─────────────────────────────────────────────────────────────┐
│  STEP 1: Preprocess Trio Samples                           │
│  - Download SRA files (if needed)                           │
│  - Convert to FASTQ                                         │
│  - Align with BWA to hg19                                   │
│  - Mark duplicates                                          │
│  - Base recalibration                                       │
│  Output: proband_final.bam, father_final.bam, mother.bam   │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 2: GATK Variant Calling                              │
│  - HaplotypeCaller per sample                               │
│  - Joint genotyping                                         │
│  Output: FAM92_trio.vcf.gz (ALL variants)                  │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 3: De Novo Filter                                    │
│  - Proband has variant (0/1 or 1/1)                        │
│  - Father is reference (0/0)                                │
│  - Mother is reference (0/0)                                │
│  Output: denovo_variants_with_header.vcf.gz                │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 4: Functional Annotation (VEP)                       │
│  - Consequence: missense, nonsense, frameshift              │
│  - SIFT score (< 0.05 = deleterious)                       │
│  - PolyPhen score (> 0.5 = damaging)                       │
│  - CADD score (> 20 = pathogenic)                          │
│  - Brain expression data                                    │
│  Output: denovo_annotated.vcf                              │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 5: FAM92 Gene Filter                                 │
│  - Extract FAM92A1 (Xq21.31)                               │
│  - Extract FAM92A2 (Xq28)                                  │
│  - Extract FAM92B (19p13.3)                                │
│  Output: FAM92_denovo.txt                                  │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 6: Clinical Interpretation Report                    │
│  - Pathogenicity assessment                                 │
│  - Genotype-phenotype correlation                           │
│  - Answer: Does FAM92 explain ASD (IQ=60)?                 │
│  Output: clinical_report.txt                               │
└─────────────────────────────────────────────────────────────┘
```

---

## 📁 Output Files

All results will be in: `/home/johan/output/autism/trio_analysis/`

| File | Description |
|------|-------------|
| `FAM92_trio.vcf.gz` | All variants in trio |
| `denovo_variants_with_header.vcf.gz` | All de novo variants |
| `denovo_annotated.vcf` | Annotated de novo variants |
| `FAM92_denovo.txt` | **FAM92-specific variants** ⭐ |
| `clinical_report.txt` | **Clinical interpretation** ⭐ |

---

## 🎯 Interpreting Results

### ✅ Scenario 1: FAM92 variant found + damaging
```
Found: FAM92A1 chrX:12345678 C>T
       Missense, SIFT=0.01, PolyPhen=0.95, CADD=28

→ CONCLUSION: FAM92 variant likely explains severe ASD phenotype
→ NEXT STEPS: Validate with Sanger sequencing, functional studies
```

### ❌ Scenario 2: No FAM92 variants
```
De novo variants in FAM92: 0

→ CONCLUSION: FAM92 does NOT explain the ASD phenotype
→ NEXT STEPS: Check other ASD genes (CHD8, SCN2A, SHANK3)
              Consider CNV analysis or whole exome sequencing
```

---

## ⏱️ Estimated Runtime

- **Preprocessing** (per sample): 2-4 hours × 3 = 6-12 hours
- **Variant calling**: 1-2 hours
- **Annotation**: 30 minutes
- **Total**: ~8-15 hours

💡 **Tip**: Run overnight!

---

## 🆘 Troubleshooting

### Problem: "I don't have parent samples"
**Solution**: You have two options:
1. Find trio data from SRA (preferred)
2. Do singleton analysis (less powerful):
   - Filter variants against gnomAD population
   - Focus on known ASD genes only
   - Cannot prove de novo status

### Problem: "Script fails at VEP annotation"
**Solution**: VEP cache download may take time
```bash
# Check VEP cache
ls ~/.vep/homo_sapiens/

# If empty, the script will auto-download (20-30 GB)
```

### Problem: "Out of memory"
**Solution**: Process X chromosome only (FAM92A1/A2 are on X)
```bash
# Edit script, add to GATK commands:
-L chrX
```

---

## 📚 Additional Resources

- **GATK Best Practices**: https://gatk.broadinstitute.org/
- **SFARI ASD Gene Database**: https://gene.sfari.org/
- **VEP Documentation**: https://www.ensembl.org/vep
- **gnomAD Browser**: https://gnomad.broadinstitute.org/

---

## 🔬 Scientific Background

### FAM92 Gene Family
- **FAM92A1**: Xq21.31, brain-expressed
- **FAM92A2**: Xq28, neural function
- **FAM92B**: 19p13.3, intracellular signaling

### Why De Novo Analysis?
- Proband: IQ=60 (severe ASD)
- Parents: Unaffected
- Sporadic case → likely de novo mutation
- De novo variants explain ~30% of severe ASD

### Pathogenicity Criteria
1. **SIFT** < 0.05 (protein function disrupted)
2. **PolyPhen** > 0.5 (structural damage)
3. **CADD** > 20 (top 1% pathogenic)
4. **Brain expression** (GTEx database)
5. **Loss-of-function** (truncating variants)

---

## 📝 Citation

If this analysis leads to publication, consider citing:
- GATK: McKenna et al. (2010) Genome Res
- VEP: McLaren et al. (2016) Genome Biol
- SFARI Gene: Abrahams et al. (2013) Mol Autism

---

**Questions?** Check the full documentation: [README_denovo.md](README_denovo.md)
