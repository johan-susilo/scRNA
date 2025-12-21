# scRNA-seq Pipeline Cheat Sheet

## Quick Commands for /home/johan/output/skin_pmh

### 1️⃣ Cell Annotation (START HERE)
```bash
cd /home/johan/pipeline/scRNA/skin

Rscript cell_annotation.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -s all --consensus \
  -o /home/johan/output/skin_pmh/annotations
```
**Output:** `annotations/consensus/consensus_annotation.tsv` ✅

---

### 2️⃣ Pseudotime Analysis

**All clusters:**
```bash
Rscript pseudotime.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  --all_clusters \
  -o /home/johan/output/skin_pmh/pseudotime_all
```

**Specific clusters (e.g., fibroblasts in 0,2,6):**
```bash
Rscript pseudotime.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -c 0,2,6 \
  -o /home/johan/output/skin_pmh/pseudotime_fibroblasts
```
**Output:** `pseudotime_all/pseudotime.pdf` ✅

---

## Check Results

```bash
# View consensus annotations
cat /home/johan/output/skin_pmh/annotations/consensus/consensus_annotation.tsv

# View pseudotime statistics
cat /home/johan/output/skin_pmh/pseudotime_all/pseudotime_by_cluster.csv

# Open plots
evince /home/johan/output/skin_pmh/annotations/markers/Classical_markers_Fibroblasts.pdf
evince /home/johan/output/skin_pmh/pseudotime_all/pseudotime.pdf
```

---

## Pipeline Options

### Cell Annotation Steps
- `singleR` - SingleR annotation only
- `celliD` - CelliD annotation only
- `scCATCH` - scCATCH annotation only
- `markers` - Classical marker plots only
- `consensus` - Generate consensus from all methods
- `all` - Run everything

### Preprocessing Steps (if re-running)
- `read_csv` - Parse sample CSV
- `process` - Process individual samples
- `integrate` - Integrate all samples
- `plot` - Generate plots
- `all` - Complete pipeline

---

## File Locations

**Input:**
- Integrated data: `/home/johan/output/skin_pmh/TN.combined_dim30.rds`

**Outputs:**
- Annotations: `/home/johan/output/skin_pmh/annotations/consensus/consensus_annotation.tsv`
- Markers: `/home/johan/output/skin_pmh/annotations/markers/*.pdf`
- Pseudotime: `/home/johan/output/skin_pmh/pseudotime_*/pseudotime.pdf`
- Logs: `*/logs/`

---

## Common Parameters

### Cell Annotation
```bash
--tissue skin          # Tissue type for scCATCH (skin, blood, etc.)
--consensus           # Generate consensus annotation
-s all                # Run all methods
```

### Pseudotime
```bash
-c 0,1,2,3           # Specific clusters
--all_clusters       # Use all clusters
```

### Preprocessing (if needed)
```bash
--doublet_rate 0.08  # Doublet rate
--min_features 200   # Min features per cell
--max_features 5000  # Max features per cell
--max_mt 30          # Max MT%
```

---

## Troubleshooting

**Check logs:**
```bash
tail -100 /home/johan/output/skin_pmh/annotations/logs/annotation_*.log
tail -100 /home/johan/output/skin_pmh/pseudotime_*/logs/pseudotime_*.log
```

**Test individual steps:**
```bash
Rscript cell_annotation.R -r TN.combined_dim30.rds -s read_rds
Rscript cell_annotation.R -s singleR
Rscript cell_annotation.R -s consensus
```

---

## Complete Workflow

```bash
#!/bin/bash
cd /home/johan/pipeline/scRNA/skin

# Annotate
Rscript cell_annotation.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  -s all --consensus \
  -o /home/johan/output/skin_pmh/annotations

# Pseudotime
Rscript pseudotime.R \
  -r /home/johan/output/skin_pmh/TN.combined_dim30.rds \
  --all_clusters \
  -o /home/johan/output/skin_pmh/pseudotime_all

echo "Done! Check:"
echo "- annotations/consensus/consensus_annotation.tsv"
echo "- pseudotime_all/pseudotime.pdf"
```
