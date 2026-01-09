# RNA-seq Pipeline Scripts Guide

**Location:** `/home/johan/pipeline/rnaseq/chicken/`

---

## Available Scripts

### 1. `preprocess.sh` - Quality Control

**What it does:** Runs FastQC on all raw FASTQ files

**How to use:**
```bash
bash preprocess.sh
```

**Output:** `${OUTPUT}/fastqc/*.html` - Quality reports

**When to use:** First step - check data quality before processing

---

### 2. `trim.sh` - Adapter Trimming

**What it does:** Removes adapters and low-quality bases using Cutadapt

**How to use:**
```bash
# Edit adapter sequences if needed
nano trim.sh

# Run trimming
bash trim.sh
```

**Output:** `${OUTPUT}/trimmed/${sample}/*.trimmed.fastq.gz`

**When to use:** After QC, before mapping

---

### 3. `map.sh` or `star.sh` - Genome Mapping

**What it does:** Aligns trimmed reads to reference genome using STAR

**How to use:**
```bash
# Make sure STAR index exists
bash map.sh
```

**Output:** `${OUTPUT}/mapped/${sample}/*Aligned.sortedByCoord.out.bam`

**When to use:** After trimming

**Note:** Both `map.sh` and `star.sh` do the same thing (genome alignment)

---

### 4. `count.sh` - Gene Quantification

**What it does:** Counts reads per gene using HTSeq-count

**How to use:**
```bash
# Set strandedness (yes/no/reverse)
export STRANDED=reverse

# Run counting
bash count.sh
```

**Output:** `${OUTPUT}/${sample}/count/*.counts.tsv`

**When to use:** After mapping

**Important:** Set correct strandedness for your library prep!

---

### 5. `combine_htseq_reports.py` - Combine Count Files

**What it does:** Merges individual HTSeq count files into a single count matrix

**How to use:**
```bash
python3 combine_htseq_reports.py --data_folder ${OUTPUT} --mode count
```

**Input:** Individual count files from count.sh

**Output:**
- `combined_counts.tsv` - Count matrix (genes × samples)
- `metadata.tsv` - Sample metadata extracted from filenames

**When to use:** After counting all samples, before differential expression

---

### 6. `differential.R` - Differential Expression

**What it does:** Runs DESeq2 for statistical analysis of gene expression

**How to use:**
```bash
# Edit paths and comparison groups in differential.R first
nano differential.R

# Run differential expression analysis
Rscript differential.R
```

**Input required:**
- Count matrix from combine_htseq_reports.py
- Sample metadata

**Output:** `${OUTPUT}/DGE/`
- `result_*.csv` - DE results for each comparison
- `plots/` - MA plot, volcano plot, PCA, heatmap

**When to use:** After combining count files

---

### 7. `go.R` - Gene Ontology Enrichment Analysis

**What it does:** Performs GO enrichment analysis on differentially expressed genes

**How to use:**
```bash
# Edit input file path in go.R first
nano go.R

# Run GO enrichment
Rscript go.R
```

**Input required:**
- DE results file (CSV from differential.R)
- GO annotation file (GAF format)
- Set thresholds: `padj_threshold` and `lfc_threshold`

**Output:** `${OUTPUT}/GO/`
- `upregulated_GO_*.csv` - GO terms for upregulated genes
- `downregulated_GO_*.csv` - GO terms for downregulated genes
- `dotplot_*.pdf` - GO enrichment dot plots
- `barplot_*.pdf` - GO enrichment bar plots

**When to use:** After differential expression analysis to find biological pathways

---

### 8. `check.sh` - Status Check

**What it does:** Checks which samples have completed each step

**How to use:**
```bash
bash check.sh
```

**Output:** Status report of pipeline progress

**When to use:** To monitor pipeline progress or find failed samples



## Typical Workflow
```bash
# Step 1: Quality control
bash preprocess.sh

# Step 2: Trimming
bash trim.sh

# Step 3: Mapping
bash map.sh

# Step 4: Counting
bash count.sh

# Step 5: Combine count files into matrix
python3 combine_htseq_reports.py --data_folder /home/johan/johan/output/chicken --mode count

# Step 6: Differential expression
Rscript differential.R

# Step 7: Gene Ontology enrichment (optional)
Rscript go.R
```

---

## Required Files

### 1. Sample Metadata: `sample.csv`

**Format:**
```csv
sample_name,condition,replicate
control_1,control,1
control_2,control,2
treatment_1,treatment,1
treatment_2,treatment,2
```

**First column:** Must match FASTQ filename prefix

### 2. Reference Files

- **Genome FASTA:** `${REF}/fasta/*.fa`
- **Gene Annotation GTF:** `${REF}/gtf/*.gtf`
- **STAR Index:** `${REF}/star_index/` (pre-built)

---

## Configuration

### Edit These Variables in Each Script:

```bash
OUTPUT="/home/johan/johan/output/chicken"          # Output directory
RAW="/mnt/2_80T/Data/TS250721001"                  # Raw FASTQ location
METADATA="/home/johan/pipeline/rnaseq/chicken/sample.csv"  # Sample info
```

### Common Parameters:

**Trimming:**
```bash
ADAPTER_R1="ACTGTCTCTTATACACATCT"  # Illumina TruSeq adapter
MIN_LENGTH=20                       # Minimum read length
```

**Counting:**
```bash
STRANDED="reverse"  # For Illumina TruSeq stranded RNA-seq
# Options: yes | no | reverse
```

---

## Input File Naming

Scripts expect FASTQ files named:
- `${sample}_1.fq.gz` and `${sample}_2.fq.gz` OR
- `${sample}_R1.fastq.gz` and `${sample}_R2.fastq.gz`

Where `${sample}` matches the first column in `sample.csv`

---

## Output Structure

```
output/
├── fastqc/                  # FastQC reports
├── trimmed/                 # Trimmed reads
│   └── ${sample}/
├── mapped/                  # BAM files
│   └── ${sample}/
├── ${sample}/count/         # Gene counts
└── differential_expression/ # DESeq2 results
```

---

## Quick Reference

| Script | Input | Output | Time |
|--------|-------|--------|------|
| `preprocess.sh` | Raw FASTQ | FastQC HTML | 5-10 min/sample |
| `trim.sh` | Raw FASTQ | Trimmed FASTQ | 10-20 min/sample |
| `map.sh` | Trimmed FASTQ | Sorted BAM | 20-40 min/sample |
| `count.sh` | BAM files | Count tables | 10-30 min/sample |
| `combine_htseq_reports.py` | Count tables | Count matrix | <1 min total |
| `differential.R` | Count matrix | DE results | 5-15 min total |
| `go.R` | DE results | GO enrichment | 2-5 min total |

---

## Helper Files

### R Scripts:
- `differential.R` - DESeq2 analysis code (see Script #6 above)
- `differential2.R` - Alternative DE analysis
- `go.R` - Gene Ontology enrichment analysis
- `go copy.R` - Backup of GO script

### Python Scripts:
- `combine_htseq_reports.py` - Merge count files into matrix (see Script #5 above)

### Config Files:
- `sample.csv` - Sample metadata
- `dge_color.yaml` - Color scheme for plots



## Software Requirements

- FastQC (v0.11+)
- Cutadapt (v3.0+)
- STAR (v2.7+)
- Samtools (v1.10+)
- HTSeq (v0.13+)
- Python 3 with pandas
- R (v4.0+) with packages:
  - DESeq2
  - ggplot2
  - clusterProfiler (for GO enrichment)
  - org.Gg.eg.db (chicken annotation)
  - enrichplot
- GNU Parallel

**Install command-line tools:**
```bash
sudo apt install fastqc samtools parallel
pip install cutadapt htseq pandas
conda install -c bioconda star
```

**Install R packages:**
```R
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "ggplot2", "clusterProfiler", "org.Gg.eg.db", "enrichplot"))
```

---

**Location:** `/home/johan/pipeline/rnaseq/chicken/`
