# Usage Examples for /home/johan/output/skin_pmh

This guide provides **exact command examples** for your existing data in `/home/johan/output/skin_pmh`.

## Your Current Data Structure

```
/home/johan/output/skin_pmh/
├── processed/              # Individual processed samples (9 samples)
│   ├── NSkin1_NC_DPP4_processed.rds
│   ├── NSkin2_NC_DPP4_processed.rds
│   ├── NSkin3_NC_DPP4_processed.rds
│   ├── NSkin4_MB_processed.rds
│   ├── TN131_HTY_processed.rds
│   ├── TN244_HTY_processed.rds
│   ├── TN258_Unaffected_processed.rds
│   ├── TN259_acute_processed.rds
│   └── TN260_chronic_processed.rds
├── plots/                  # QC and UMAP plots
├── tables/                 # Cell count tables
├── logs/                   # Processing logs
├── annotations/            # Annotation results (partial)
├── TN.combined_dim30.rds  # ✅ MAIN INTEGRATED DATA (1.29 GB)
├── TN.anchorsdim30.rds    # Integration anchors
└── samples_df.rds         # Sample metadata
```

---

## 1. Cell Annotation (START HERE)

Since you already have the integrated data (`TN.combined_dim30.rds`), you can directly run cell annotation.

### Option A: Run Complete Annotation with Consensus (RECOMMENDED)

```bash
cd /home/johan/pipeline/scRNA/skin

Rscript cell_annotation.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -s all \
  -o /home/johan/output/skin_pmh/annotations \
  --consensus \
  --tissue skin
```

**What this does:**
- Runs SingleR (HPCA + BlueprintEncode databases)
- Runs CelliD (PanglaoDB signatures)
- Runs scCATCH (skin-specific markers)
- Generates classical marker plots (18 cell types)
- **Creates consensus annotation** automatically
- Time: ~30-60 minutes depending on cell count

**Output files:**
```
/home/johan/output/skin_pmh/annotations/
├── singleR/
│   ├── SingleR_hpca.csv
│   ├── SingleR_hpca_summary.tsv        # SingleR HPCA annotations
│   ├── SingleR_bpe.csv
│   └── SingleR_bpe_summary.tsv         # SingleR BPE annotations
├── celliD/
│   ├── CelliD_PanglaoDB.csv
│   └── CelliD_PanglaoDB_summary.tsv    # CelliD annotations
├── scCATCH/
│   ├── scCATCH.csv
│   └── scCATCH_summary.tsv             # scCATCH annotations
├── markers/
│   ├── Classical_markers_Epithelial.pdf
│   ├── Classical_markers_Fibroblasts.pdf
│   ├── Classical_markers_T_cells.pdf
│   └── ... (18 total marker plots)
├── consensus/
│   └── consensus_annotation.tsv        # ✅ FINAL CELL TYPE ANNOTATIONS
└── logs/
    └── annotation_*.log
```

### Option B: Run Individual Annotation Methods

```bash
# Run only SingleR
Rscript cell_annotation.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -s singleR \
  -o /home/johan/output/skin_pmh/annotations

# Run only CelliD
Rscript cell_annotation.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -s celliD \
  -o /home/johan/output/skin_pmh/annotations

# Run only scCATCH
Rscript cell_annotation.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -s scCATCH \
  -o /home/johan/output/skin_pmh/annotations

# Generate only marker plots
Rscript cell_annotation.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -s markers \
  -o /home/johan/output/skin_pmh/annotations

# Generate consensus from existing results
Rscript cell_annotation.R \
  -s consensus \
  -o /home/johan/output/skin_pmh/annotations
```

---

## 2. Pseudotime Analysis

After you have the integrated data, you can perform trajectory analysis.

### Example 1: Analyze Fibroblast Clusters (like Yan's example)

```bash
cd /home/johan/pipeline/scRNA/skin

# First, check your cluster numbers in the UMAP plots
# Look at: /home/johan/output/skin_pmh/plots/TNcombined_umap_labelT.pdf

# Example: If fibroblasts are in clusters 0, 2, 3, 6, 15
Rscript pseudotime.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -c 0,2,3,6,15 \
  -o /home/johan/output/skin_pmh/pseudotime_fibroblasts
```

### Example 2: Analyze All Clusters

```bash
Rscript pseudotime.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  --all_clusters \
  -o /home/johan/output/skin_pmh/pseudotime_all
```

### Example 3: Analyze Specific Cell Type (e.g., T cells)

```bash
# If T cells are in clusters 5, 13
Rscript pseudotime.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -c 5,13 \
  -o /home/johan/output/skin_pmh/pseudotime_tcells
```

**Output files:**
```
/home/johan/output/skin_pmh/pseudotime_fibroblasts/
├── subset_umap.pdf                      # Selected clusters UMAP
├── monocle3_umap_clusters.pdf          # Monocle3 UMAP
├── trajectory_by_cluster.pdf           # Learned trajectory
├── pseudotime.pdf                      # ✅ Main pseudotime plot
├── trajectory_combined.pdf             # Combined visualization
├── trajectory_by_sample.pdf            # Split by sample
├── pseudotime_values.csv               # ✅ Per-cell pseudotime values
├── pseudotime_by_cluster.csv          # Summary statistics
├── pseudotime_violin_by_cluster.pdf   # Distribution
├── pseudotime_gene_fits.csv           # All gene test results
├── pseudotime_significant_genes.csv   # ✅ Genes changing along trajectory
├── top_pseudotime_genes.pdf           # Top 9 genes
├── monocle3_cds.rds                   # Monocle3 object for R
└── logs/
```

---

## 3. Re-running Preprocessing (if needed)

If you want to re-process from scratch with different parameters:

### Step 1: Prepare sample CSV

Create `/home/johan/data/PMH_scRNA-seq/samples.csv`:
```csv
sample_names,ident1,ident2
NSkin1_NC_DPP4,1_Normal_Skin,NSkin1
NSkin2_NC_DPP4,1_Normal_Skin,NSkin2
NSkin3_NC_DPP4,1_Normal_Skin,NSkin3
NSkin4_MB,1_Normal_Skin,NSkin4
TN131_HTY,2_Healthy_Control,HC131
TN244_HTY,2_Healthy_Control,HC244
TN258_Unaffected,3_Unaffected,UA258
TN259_acute,4_Acute,AC259
TN260_chronic,5_Chronic,CH260
```

### Step 2: Run preprocessing

```bash
cd /home/johan/pipeline/scRNA/skin

Rscript preprocessing.R \
  -f /home/johan/pipeline/scRNA/skin/input.csv \
  -d /home/johan/data/PMH_scRNA-seq \
  -s all \
  -o /home/johan/output/skin_pmh_v2 \
  --doublet_rate 0.08 \
  --min_features 200 \
  --max_features 5000 \
  --max_mt 30
```

### Step 3: Custom parameters example

```bash
# Stricter QC filtering
Rscript preprocessing.R \
  -f /home/johan/data/PMH_scRNA-seq/samples.csv \
  -d /home/johan/data/PMH_scRNA-seq \
  -s all \
  -o /home/johan/output/skin_pmh_strict \
  --doublet_rate 0.08 \
  --min_features 500 \
  --max_features 4000 \
  --max_mt 20
```

---

## 4. Complete Workflow from Existing Data

Here's the **recommended workflow** starting from your current data:

```bash
#!/bin/bash
# Complete analysis workflow for skin_pmh data

# Set working directory
cd /home/johan/pipeline/scRNA/skin

# Step 1: Cell Annotation with Consensus (30-60 min)
echo "===== Running Cell Annotation ====="
Rscript cell_annotation.R \
  -r /home/johan/output/skin_pmh_v2/TN.combined_dim30.rds \
  -s all \
  -o /home/johan/output/skin_pmh_v2/annotations \
  --consensus \
  --tissue skin

# Step 2: Review consensus results
echo "===== Consensus Annotation Results ====="
cat /home/johan/output/skin_pmh/annotations/consensus/consensus_annotation.tsv

# Step 3: Pseudotime for all clusters (optional)
echo "===== Running Pseudotime Analysis ====="
Rscript pseudotime.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  --all_clusters \
  -o /home/johan/output/skin_pmh/pseudotime_all

echo "===== Analysis Complete ====="
echo "Check results:"
echo "  - Annotations: /home/johan/output/skin_pmh/annotations/consensus/consensus_annotation.tsv"
echo "  - Pseudotime: /home/johan/output/skin_pmh/pseudotime_all/pseudotime.pdf"
```

Save this as `run_analysis.sh` and execute:
```bash
chmod +x run_analysis.sh
./run_analysis.sh
```

---

## 5. Checking Results

### View Consensus Annotations

```bash
# View consensus cell type assignments
cat /home/johan/output/skin_pmh/annotations/consensus/consensus_annotation.tsv

# Expected output:
# Cluster    Cell_Type           Count
# X0         Fibroblasts         4
# X1         T cells             3
# X2         Keratinocytes       4
# ...
```

### View Pseudotime Values

```bash
# View pseudotime statistics
cat /home/johan/output/skin_pmh/pseudotime_all/pseudotime_by_cluster.csv

# Expected output:
# cluster,n_cells,mean_pseudotime,median_pseudotime,min_pseudotime,max_pseudotime
# C0,1234,5.67,5.23,0.12,12.45
# ...
```

### Check Plots

```bash
# Open key plots
evince /home/johan/output/skin_pmh/plots/TNcombined_umap_labelT.pdf
evince /home/johan/output/skin_pmh/annotations/consensus/consensus_annotation.tsv
evince /home/johan/output/skin_pmh/pseudotime_all/pseudotime.pdf
```

---

## 6. Common Use Cases

### Use Case 1: Quick cell type annotation only

```bash
Rscript cell_annotation.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -s all \
  --consensus \
  -o /home/johan/output/skin_pmh/annotations
```

### Use Case 2: Fibroblast trajectory analysis

```bash
# 1. Check which clusters are fibroblasts from consensus
cat /home/johan/output/skin_pmh/annotations/consensus/consensus_annotation.tsv | grep -i fibro

# 2. Run pseudotime on those clusters (example: 0,2,6)
Rscript pseudotime.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -c 0,2,6 \
  -o /home/johan/output/skin_pmh/pseudotime_fibroblasts
```

### Use Case 3: Compare healthy vs disease

```bash
# Your data has different conditions in orig.ident1:
# - Normal_Skin (NSkin samples)
# - Healthy_Control (TN131, TN244)
# - Unaffected (TN258)
# - Acute (TN259)
# - Chronic (TN260)

# Check the split UMAP
evince /home/johan/output/skin_pmh/plots/TNcombined_umap_labelT_splitorigident1.pdf
```

---

## 7. Troubleshooting

### If annotation fails:

```bash
# Check log file
tail -100 /home/johan/output/skin_pmh/annotations/logs/annotation_*.log

# Run step by step
Rscript cell_annotation.R -r /home/johan/output/skin_pmh/TN.combined_dim30.rds -s read_rds -o /home/johan/output/skin_pmh/annotations
Rscript cell_annotation.R -s singleR -o /home/johan/output/skin_pmh/annotations
Rscript cell_annotation.R -s markers -o /home/johan/output/skin_pmh/annotations
```

### If pseudotime fails:

```bash
# Check log file
tail -100 /home/johan/output/skin_pmh/pseudotime_all/logs/pseudotime_*.log

# Try with fewer clusters first
Rscript pseudotime.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -c 0,1,2 \
  -o /home/johan/output/skin_pmh/pseudotime_test
```

---

## Summary

**For your current data (`/home/johan/output/skin_pmh/`), run these commands:**

```bash
# 1. Cell Annotation (REQUIRED - provides cell types)
Rscript /home/johan/pipeline/scRNA/skin/cell_annotation.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -s all --consensus \
  -o /home/johan/output/skin_pmh/annotations

# 2. Pseudotime Analysis (OPTIONAL - trajectory inference)
Rscript /home/johan/pipeline/scRNA/skin/pseudotime.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  --all_clusters \
  -o /home/johan/output/skin_pmh/pseudotime_all
```

**Results to check:**
- Cell types: `/home/johan/output/skin_pmh/annotations/consensus/consensus_annotation.tsv`
- Marker plots: `/home/johan/output/skin_pmh/annotations/markers/`
- Trajectory: `/home/johan/output/skin_pmh/pseudotime_all/pseudotime.pdf`
