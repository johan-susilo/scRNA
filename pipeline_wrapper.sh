#!/bin/bash
set -e

# ==============================================================================
# UNIFIED DOWNSTREAM PIPELINE WRAPPER
# Subsets, annotates, and analyzes Fibroblasts, Macrophages, and Mast Cells
# ==============================================================================

# Input and Output Configuration
BASE_IN="/home/johan/output/skin_pmh_harmony_sctransform2/TN.combined_dim30.rds"
BASE_OUT="/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster_strict"
SCRIPT_DIR="/home/johan/pipeline/scRNA/skin"

# Target Clusters (Adjust these if clusters shift upstream)
FIBRO_CLUSTERS="1,3,8"
MACRO_CLUSTERS="6"
MAST_CLUSTERS="14" # Based on contextual hint, but parameterized for flexibility

# ==============================================================================
# 1. SUBSETTING (Universal Subsetter)
# ==============================================================================
echo "=== 1. SUBSETTING ==="

Rscript ${SCRIPT_DIR}/subset.R -i "${BASE_IN}" -o "${BASE_OUT}" -c "${FIBRO_CLUSTERS}" -p "fibroblast"
Rscript ${SCRIPT_DIR}/subset.R -i "${BASE_IN}" -o "${BASE_OUT}" -c "${MACRO_CLUSTERS}" -p "macrophage"
Rscript ${SCRIPT_DIR}/subset.R -i "${BASE_IN}" -o "${BASE_OUT}" -c "${MAST_CLUSTERS}" -p "mast_cell"

# ==============================================================================
# 2. DETAIL ANNOTATION (Only Fibroblast & Macrophage)
# ==============================================================================
echo "=== 2. DETAIL ANNOTATION ==="

# Using standard detailed annotation script for fibroblasts
Rscript ${SCRIPT_DIR}/detail_annotation.R

# Using the unified script parameterized for macrophages
Rscript ${SCRIPT_DIR}/detail_annotation_macrophage.R -i "${BASE_OUT}/macrophage/processed/macrophage_subset_processed.rds" -c "macrophage" -o "${BASE_OUT}"

# ==============================================================================
# 3. DIFFERENTIAL GENE EXPRESSION & PATHWAY ANALYSIS (All Subsets)
# ==============================================================================
echo "=== 3. DGE & GO ANALYSIS ==="

# Note: dge_unified.R automatically iterates through all cell type folders in BASE_OUT
Rscript ${SCRIPT_DIR}/dge_unified.R -d "${BASE_OUT}"
Rscript ${SCRIPT_DIR}/go_unified.R -d "${BASE_OUT}"

# ==============================================================================
# 4. PSEUDOTIME (All Subsets)
# ==============================================================================
echo "=== 4. PSEUDOTIME ==="

for celltype in "fibroblast" "macrophage" "mast_cell"; do
    # Use the detail annotated version if it exists, otherwise fallback to standard subset
    input_rds="${BASE_OUT}/${celltype}/processed/${celltype}_detailed_annotated.rds"
    if [ ! -f "$input_rds" ]; then
        input_rds="${BASE_OUT}/${celltype}/processed/${celltype}_subset_processed.rds"
    fi

    output_dir="${BASE_OUT}/${celltype}/pseudotime"
    mkdir -p "$output_dir"

    echo "Running Pseudotime for ${celltype}..."
    Rscript ${SCRIPT_DIR}/pseudotime.R -i "${input_rds}" -o "${output_dir}" --all_clusters
done

echo "=== PIPELINE COMPLETE ==="
