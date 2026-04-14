library(Seurat)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(tidyr)
library(ggrepel)
library(stringr)

# ==============================================================================
# 1. Define the Pseudobulk Function for Sub-clusters
# ==============================================================================
run_subcluster_dge <- function(subcluster_id, cell_type, base_out_dir, counts_mat, cells_per_pb, min_cells = 10, min_reps = 2) {
  
  message(paste("\n  -> Evaluating", cell_type, "Sub-cluster:", subcluster_id))
  
  # Setup directory
  out_dir <- file.path(base_out_dir, paste0("subcluster_", subcluster_id))
  
  # === Select Specific Sub-cluster Columns ===
  all_cols <- colnames(counts_mat)
  
  # Extract the cluster number from the pseudobulk column names safely
  cluster_prefixes <- str_extract(all_cols, "^[^_]+")
  cluster_numbers <- str_remove_all(cluster_prefixes, "[a-zA-Z]")
  target_cols <- all_cols[cluster_numbers == as.character(subcluster_id)]
  
  if (length(target_cols) == 0) {
    message(paste("    - SKIPPING: No pseudobulk columns found for Sub-cluster", subcluster_id))
    return(NULL)
  }

  # === SAFEGUARD 1: Filter Low-Cell Pseudobulk Samples ===
  # THE FIX: Strip the "g" from target_cols to match our dictionary keys
  pb_counts <- sapply(target_cols, function(col_name) {
    # Remove any letters at the start (e.g., "g0" becomes "0")
    clean_name <- stringr::str_remove(col_name, "^[a-zA-Z]+")
    
    if (clean_name %in% names(cells_per_pb)) {
      return(as.numeric(cells_per_pb[clean_name]))
    } else {
      return(0)
    }
  })
  
  valid_cols <- target_cols[pb_counts >= min_cells]
  dropped_cols <- setdiff(target_cols, valid_cols)
  
  if (length(dropped_cols) > 0) {
    message(paste("    - Dropped", length(dropped_cols), "sample(s) for having <", min_cells, "cells."))
  }
  
  if (length(valid_cols) == 0) {
    message("    - SKIPPING: No valid samples left after cell count filtering.")
    return(NULL)
  }
  
  counts_subset <- counts_mat[, valid_cols, drop = FALSE]
  
  # === Create Metadata ===
  meta <- data.frame(pseudobulk_id = valid_cols) %>%
    mutate(
      cluster = str_extract(pseudobulk_id, "^[^_]+"),
      condition = str_extract(pseudobulk_id, "[^_]+$"),
      sample_id = str_replace(pseudobulk_id, paste0("^", cluster, "_"), ""),
      sample_id = str_replace(sample_id, paste0("_", condition, "$"), "")
    )
  rownames(meta) <- meta$pseudobulk_id
  meta$condition <- factor(meta$condition, levels = c("Healthy", "PMH"))
  
  # === SAFEGUARD 2: Strict Biological Replicate Check ===
  pmh_count <- sum(meta$condition == "PMH")
  healthy_count <- sum(meta$condition == "Healthy")
  
  if (pmh_count < min_reps | healthy_count < min_reps) {
    message(paste("    - SKIPPING: Insufficient replicates. PMH:", pmh_count, "| Healthy:", healthy_count, "(Requires", min_reps, "each)"))
    return(NULL)
  }
  
  message(paste("    - Proceeding with DGE. Replicates -> PMH:", pmh_count, "| Healthy:", healthy_count))
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
  
  # === DESeq2 Pipeline ===
  dds <- DESeqDataSetFromMatrix(
    countData = counts_subset,
    colData = meta,
    design = ~ condition
  )
  
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "PMH", "Healthy"))
  
  # === Formatting and Saving Results ===
  res_df <- as.data.frame(res) %>% 
    filter(!is.na(padj)) %>% 
    mutate(gene = rownames(.)) %>%
    mutate(
      significance = case_when(
        padj < 0.05 & log2FoldChange > 1 ~ "Upregulated in PMH",
        padj < 0.05 & log2FoldChange < -1 ~ "Downregulated in PMH",
        TRUE ~ "Not Significant"
      )
    ) %>%
    select(gene, everything()) %>%
    arrange(padj)
  
  write.csv(res_df, file.path(out_dir, paste0("deseq2_", cell_type, "_subcluster_", subcluster_id, ".csv")), row.names = FALSE)
  
  # === Visualization ===
  top_genes <- res_df %>%
    filter(significance != "Not Significant") %>% 
    group_by(significance) %>% 
    arrange(padj) %>% 
    slice_head(n = 15) %>% 
    ungroup() 
  
  # Volcano plot
  volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = top_genes, 
                    aes(label = gene), 
                    color = "black", 
                    box.padding = 0.5, 
                    max.overlaps = Inf) +
    scale_color_manual(values = c("Upregulated in PMH" = "red", 
                                  "Downregulated in PMH" = "blue", 
                                  "Not Significant" = "grey80")) +
    theme_minimal() +
    labs(title = paste(tools::toTitleCase(cell_type), "Sub-cluster", subcluster_id, "(PMH vs Healthy)"))
  
  ggsave(filename = file.path(out_dir, paste0("volcano_plot_subcluster_", subcluster_id, ".png")),
         plot = volcano_plot, width = 6, height = 5)
  
  return(res_df)
}

# ==============================================================================
# 2. Automated Execution Loop (The Auto-Discovery Block)
# ==============================================================================

base_subset_dir <- "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster"

# Auto-discover all folders inside the subset directory
cell_type_folders <- list.dirs(base_subset_dir, recursive = FALSE, full.names = FALSE)

message(paste("Found", length(cell_type_folders), "cell type folders to process:", paste(cell_type_folders, collapse = ", ")))

for (cell_type in cell_type_folders) {
  
  message("\n==================================================================")
  message(paste("=== Starting Statistically-Clean DGE for:", toupper(cell_type), "==="))
  message("==================================================================")
  
  rds_path <- file.path(base_subset_dir, cell_type, "processed", paste0(cell_type, "_subset_processed.rds"))
  
  if (!file.exists(rds_path)) {
    message(paste("WARNING: Could not find RDS file at", rds_path, "- Skipping", cell_type))
    next
  }
  
  sub_obj <- readRDS(rds_path)
  sub_obj$Condition <- factor(sub_obj$Condition, levels = c("Healthy", "PMH"))
  
  # --- THE FIX: Create the exact same string that AggregateExpression uses, NO make.names() ---
  sub_obj$pb_id <- paste(sub_obj$seurat_clusters, sub_obj$orig.ident2, sub_obj$Condition, sep = "_")
  cells_per_pb <- table(sub_obj$pb_id)
  
  # AggregateExpression on SUB-CLUSTERS
  message(paste("Aggregating expression for", cell_type, "sub-clusters..."))
  pb_list <- AggregateExpression(
    sub_obj, 
    assays = "RNA", 
    slot = "counts",
    group.by = c("seurat_clusters", "orig.ident2", "Condition"),
    return.seurat = FALSE 
  )
  counts_matrix <- pb_list$RNA
  
  cell_type_dge_dir <- file.path(base_subset_dir, cell_type, "dge_pseudobulk")
  available_subclusters <- levels(sub_obj$seurat_clusters)
  
  for (sub_id in available_subclusters) {
    run_subcluster_dge(
      subcluster_id = sub_id, 
      cell_type = cell_type,
      base_out_dir = cell_type_dge_dir, 
      counts_mat = counts_matrix,
      cells_per_pb = cells_per_pb, # Pass the exactly matched cell counts
      min_cells = 5,              # Minimum cells a donor needs to be included
      min_reps = 2                 # Ensures at least 2 Healthy and 2 PMH donors remain
    )
  }
}

message("\n==================================================================")
message("=== All Statistically-Cleaned DGE pipelines executed! ===")
message("==================================================================")