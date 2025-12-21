# Quick Start Guide

## 1. Prepare Your Data

Create a CSV file with your sample information:

```csv
sample_names,ident1,ident2
sample1,1_Healthy,S1
sample2,2_Disease,S2
```

## 2. Run Preprocessing

```bash
Rscript preprocessing.R -f samples.csv -d /path/to/10X_data -s all -o output
```

**This will:**
- Remove doublets automatically
- Perform quality control
- Integrate all samples
- Generate QC plots

**Output:** `output/TN.combined_dim30.rds`

## 3. Run Cell Annotation

```bash
Rscript cell_annotation.R -r output/TN.combined_dim30.rds -s all --consensus
```

**This will:**
- Run SingleR, CelliD, scCATCH
- Generate marker plots
- Create consensus annotation

**Output:** `annotations/consensus/consensus_annotation.tsv`

## 4. Run Pseudotime (Optional)

```bash
# For specific clusters (e.g., fibroblasts)
Rscript pseudotime.R -r output/TN.combined_dim30.rds -c 0,1,2,3 -o output/pseudotime

# For all clusters
Rscript pseudotime.R -r output/TN.combined_dim30.rds --all_clusters -o output/pseudotime
```

**Output:** `output/pseudotime/pseudotime.pdf` and statistics

## Complete One-Liner

```bash
# Run everything
Rscript preprocessing.R -f samples.csv -d /data -s all -o out && \
Rscript cell_annotation.R -r out/TN.combined_dim30.rds -s all --consensus -o out/anno && \
Rscript pseudotime.R -r out/TN.combined_dim30.rds --all_clusters -o out/pseudo
```

## Key Files to Check

1. **After Preprocessing:**
   - `output/TN.combined_dim30.rds` - Main data
   - `output/plots/TNcombined_umap_labelT.pdf` - Cluster visualization
   - `output/tables/CellNumber_bygroup.csv` - Cell counts

2. **After Annotation:**
   - `annotations/consensus/consensus_annotation.tsv` - Cell types
   - `annotations/markers/` - Marker gene plots

3. **After Pseudotime:**
   - `output/pseudotime/pseudotime.pdf` - Trajectory
   - `output/pseudotime/pseudotime_values.csv` - Per-cell values

## Common Issues

**Problem:** Out of memory
**Solution:** Use fewer cores or process samples in batches

**Problem:** DoubletFinder fails
**Solution:** Check if cells per sample > 100

**Problem:** No consensus annotation
**Solution:** Run annotation with `--consensus` flag

## Next Steps

1. Review consensus annotations
2. Manually verify with marker plots
3. Apply annotations to Seurat object
4. Perform downstream analysis (DEG, pathway, etc.)
