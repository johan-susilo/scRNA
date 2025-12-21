# scRNA-seq Analysis Pipeline for Skin Samples

Automated single-cell RNA-seq analysis pipeline adapted from Yan's PMH scRNA-seq workflow. This pipeline provides end-to-end analysis from raw 10X data to pseudotime trajectories with automatic cell type annotation.

## Overview

This pipeline consists of three main scripts:

1. **preprocessing.R** - Quality control, doublet removal, and data integration
2. **cell_annotation.R** - Automated cell type annotation with consensus calling
3. **pseudotime.R** - Trajectory inference and pseudotime analysis

## Key Features

### Preprocessing Pipeline
- **Automatic DoubletFinder workflow** with optimal pK parameter identification
- **Quality control** with configurable filtering parameters
- **Cell cycle regression** during integration
- **Mitochondrial percentage filtering**
- **Sequential processing** to avoid memory issues
- **Comprehensive QC plots** and statistics

### Cell Annotation Pipeline
- **Multiple annotation methods**:
  - SingleR (HumanPrimaryCellAtlas + BlueprintEncode databases)
  - CelliD (PanglaoDB signatures)
  - scCATCH (tissue-specific markers)
  - Classical marker gene plots
- **Automatic consensus annotation** from all methods
- **Cell type normalization** for consistent nomenclature
- **No manual annotation required**

### Pseudotime Analysis
- **Monocle3-based trajectory inference**
- **Flexible cluster selection** or use all clusters
- **Pseudotime-dependent gene identification**
- **Comprehensive visualizations** and statistics

## Installation

### Required R Packages

```r
# Core packages
install.packages(c("optparse", "dplyr", "tidyverse", "ggplot2", "ggpubr",
                   "cowplot", "gridExtra", "RColorBrewer", "tidyr", "presto",
                   "ggrepel", "stringr", "patchwork", "scales", "parallel",
                   "future", "future.apply"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Seurat", "DoubletFinder", "SingleR", "celldex",
                       "CelliD", "scCATCH", "clusterProfiler", "monocle3"))
```

### Monocle3 Installation

```r
# If monocle3 is not available on Bioconductor
devtools::install_github('cole-trapnell-lab/monocle3')
```

## Usage

### 1. Preprocessing

#### Input Requirements
- CSV file with columns: `sample_names`, `ident1`, `ident2`
- 10X data directories (one per sample)

Example CSV (`samples.csv`):
```csv
sample_names,ident1,ident2
TN131_HTY,1_Healthy_skin,HC131
TN244_HTY,1_Healthy_skin,HC244
TN258_Unaffected,2_Unaffected_skin,UA258
TN259_acute,3_Acute_skin,AC259
TN260_chronic,4_Chronic_skin,CH260
```

#### Running the Pipeline

```bash
# Run complete preprocessing pipeline
Rscript preprocessing.R \
  -f samples.csv \
  -d /path/to/10X_data \
  -s all \
  -o /path/to/output \
  -c 8

# Run individual steps
Rscript preprocessing.R -f samples.csv -s read_csv
Rscript preprocessing.R -f samples.csv -d /path/to/data -s process
Rscript preprocessing.R -s integrate
Rscript preprocessing.R -s plot

# Custom QC parameters
Rscript preprocessing.R \
  -f samples.csv \
  -d /path/to/data \
  -s all \
  --doublet_rate 0.08 \
  --min_features 200 \
  --max_features 5000 \
  --max_mt 30
```

#### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-f, --file` | Required | Input CSV file with sample information |
| `-d, --datadir` | `.` | Base directory containing 10X sample folders |
| `-s, --step` | Required | Pipeline step: `read_csv`, `process`, `integrate`, `plot`, `all` |
| `-o, --output` | `output` | Output directory |
| `-c, --cores` | Auto-detect | Number of cores (currently uses sequential) |
| `--doublet_rate` | 0.08 | Expected doublet formation rate |
| `--min_features` | 200 | Minimum features per cell |
| `--max_features` | 5000 | Maximum features per cell |
| `--max_mt` | 30 | Maximum mitochondrial percentage |

#### Output

```
output/
├── processed/              # Processed individual samples
├── plots/                  # QC and visualization plots
├── tables/                 # Cell counts and statistics
├── logs/                   # Pipeline logs
└── TN.combined_dim30.rds  # Integrated Seurat object
```

### 2. Cell Annotation

```bash
# Run complete annotation pipeline with consensus
Rscript cell_annotation.R \
  -r output/TN.combined_dim30.rds \
  -s all \
  -o output/annotations \
  --consensus

# Run individual annotation methods
Rscript cell_annotation.R -r TN.combined_dim30.rds -s singleR
Rscript cell_annotation.R -r TN.combined_dim30.rds -s celliD
Rscript cell_annotation.R -r TN.combined_dim30.rds -s scCATCH
Rscript cell_annotation.R -r TN.combined_dim30.rds -s markers

# Generate consensus only (after running all methods)
Rscript cell_annotation.R -s consensus

# Custom tissue type for scCATCH
Rscript cell_annotation.R \
  -r TN.combined_dim30.rds \
  -s scCATCH \
  --tissue blood
```

#### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-r, --rds` | Required | Path to integrated RDS file |
| `-s, --step` | `all` | Step: `singleR`, `markers`, `celliD`, `scCATCH`, `consensus`, `all` |
| `-o, --output` | `annotations` | Output directory |
| `--consensus` | FALSE | Generate consensus annotations |
| `--tissue` | `skin` | Tissue type for scCATCH |

#### Output

```
annotations/
├── singleR/
│   ├── SingleR_hpca_summary.tsv
│   └── SingleR_bpe_summary.tsv
├── celliD/
│   └── CelliD_PanglaoDB_summary.tsv
├── scCATCH/
│   └── scCATCH_summary.tsv
├── markers/
│   └── Classical_markers_*.pdf
├── consensus/
│   └── consensus_annotation.tsv    # Final consensus results
└── logs/
```

#### Consensus Annotation Logic

The consensus annotation aggregates results from all methods:
1. Reads annotation results from SingleR (HPCA + BPE), CelliD, and scCATCH
2. Normalizes cell type names (e.g., "Monocyte" → "Monocytes")
3. Counts votes for each cell type per cluster
4. Assigns the cell type(s) with maximum votes

Example output (`consensus_annotation.tsv`):
```
Cluster    Cell_Type           Count
X0         Fibroblasts         4
X1         T cells             3
X2         Keratinocytes       4
```

### 3. Pseudotime Analysis

```bash
# Analyze specific clusters
Rscript pseudotime.R \
  -r output/TN.combined_dim30.rds \
  -c 0,2,3,6,15 \
  -o output/pseudotime

# Analyze all clusters
Rscript pseudotime.R \
  -r output/TN.combined_dim30.rds \
  --all_clusters \
  -o output/pseudotime
```

#### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-r, --rds` | Required | Path to integrated RDS file |
| `-c, --clusters` | NULL | Comma-separated cluster IDs |
| `--all_clusters` | FALSE | Use all clusters |
| `-o, --output` | `output/pseudotime` | Output directory |
| `--method` | `monocle3` | Trajectory method |

#### Output

```
pseudotime/
├── subset_umap.pdf                      # Selected clusters
├── monocle3_umap_clusters.pdf          # Monocle3 UMAP
├── trajectory_by_cluster.pdf           # Learned trajectory
├── pseudotime.pdf                      # Pseudotime coloring
├── trajectory_combined.pdf             # Combined visualization
├── pseudotime_values.csv               # Per-cell pseudotime
├── pseudotime_by_cluster.csv          # Cluster statistics
├── pseudotime_violin_by_cluster.pdf   # Distribution
├── pseudotime_gene_fits.csv           # All gene tests
├── pseudotime_significant_genes.csv   # Significant genes
├── top_pseudotime_genes.pdf           # Top gene expression
├── monocle3_cds.rds                   # Monocle3 object
└── logs/
```

## Complete Workflow Example

```bash
# 1. Preprocessing
Rscript preprocessing.R \
  -f samples.csv \
  -d /home/johan/data/PMH_scRNA-seq \
  -s all \
  -o /home/johan/output/skin_analysis

# 2. Cell Annotation with Consensus
Rscript cell_annotation.R \
  -r /home/johan/output/skin_analysis/TN.combined_dim30.rds \
  -s all \
  -o /home/johan/output/skin_analysis/annotations \
  --consensus

# 3. Pseudotime Analysis (fibroblast clusters)
Rscript pseudotime.R \
  -r /home/johan/output/skin_analysis/TN.combined_dim30.rds \
  -c 0,2,3,6,15 \
  -o /home/johan/output/skin_analysis/pseudotime
```

## Methodology

### DoubletFinder Workflow (from Yan's approach)
1. Initial normalization and PCA
2. Optimal pK parameter identification using BCmetric
3. Homotypic doublet proportion estimation
4. Doublet classification with adjusted expectations
5. Dynamic column detection for robust filtering

### Integration Strategy
- Cell cycle scoring (S and G2M phases)
- Regression of cell cycle effects and MT%
- Canonical Correlation Analysis (CCA) integration
- 30 PCs for dimensionality reduction

### Cell Type Annotation
- **SingleR**: Reference-based annotation using curated databases
- **CelliD**: Gene signature-based prediction using PanglaoDB
- **scCATCH**: Tissue-specific marker-based identification
- **Consensus**: Majority voting across all methods with normalization

### Cell Type Normalization Rules
```
Monocyte → Monocytes
Endothelial_cells → Endothelial cells
Macrophage → Macrophages
DC → Dendritic cells
Killer Cell / Natural Killer Cell / NK cell → NK cells
```

## Key Improvements from Source Code

1. **Fully Automated**: No manual cluster annotation required
2. **Robust Error Handling**: Comprehensive try-catch blocks
3. **Dynamic Column Detection**: Handles DoubletFinder output variations
4. **Configurable Parameters**: Command-line control of all QC thresholds
5. **Consensus Annotation**: Automatic integration of multiple methods
6. **Sequential Processing**: Avoids memory deadlocks
7. **Comprehensive Logging**: Detailed logs with timestamps
8. **Cell Type Normalization**: Consistent nomenclature across methods

## Troubleshooting

### Memory Issues
- Use sequential processing (default)
- Process samples in smaller batches
- Reduce number of features with `--max_features`

### DoubletFinder Errors
- Check pK identification plots
- Adjust `--doublet_rate` based on cell count
- Verify sufficient cells per cluster

### Annotation Issues
- Check tissue type for scCATCH
- Review individual method outputs before consensus
- Verify marker gene expression in plots

## References

1. **Seurat**: Hao et al., Cell (2021)
2. **DoubletFinder**: McGinnis et al., Cell Systems (2019)
3. **SingleR**: Aran et al., Nature Immunology (2019)
4. **CelliD**: Cortal et al., Nucleic Acids Research (2021)
5. **scCATCH**: Shao et al., iScience (2020)
6. **Monocle3**: Cao et al., Nature (2019)

## Citation

If using this pipeline, please cite the original Yan et al. PMH scRNA-seq study and the relevant method papers above.

## Contact

For questions or issues, please refer to the original source code documentation or contact the pipeline maintainer.
