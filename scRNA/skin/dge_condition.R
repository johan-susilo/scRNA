library(Seurat)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(tidyr)
library(ggrepel)
library(stringr)

# ==============================================================================
# Automated "All Combined" Pseudobulk DGE for ALL Cell Types
# ==============================================================================

base_dir <- "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster"

# Auto-discover all cell type folders inside the base directory
cell_types <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

message(paste("Found", length(cell_types), "cell types to process:", paste(cell_types, collapse = ", ")))

# Loop through every discovered cell type
for (cell_type in cell_types) {
  
  message(paste("\n=================================================================="))
  message(paste("=== Starting 'All Combined' DGE for:", toupper(cell_type), "==="))
  message(paste("=================================================================="))
  
  # 1. Setup paths and check if the processed RDS actually exists
  rds_path <- file.path(base_dir, cell_type, "processed", paste0(cell_type, "_subset_processed.rds"))
  
  if (!file.exists(rds_path)) {
    message(paste("  -> WARNING: Could not find RDS file for", cell_type, "at", rds_path, "- Skipping."))
    next
  }
  
  out_dir <- file.path(base_dir, cell_type, "dge_pseudobulk", "all_combined")
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
  
  # Load data
  message(paste("Loading", cell_type, "subset object..."))
  sub_obj <- readRDS(rds_path)
  sub_obj$Condition <- factor(sub_obj$Condition, levels = c("Healthy", "PMH"))
  
  # 2. Aggregate Expression
  # Notice we removed "seurat_clusters" from group.by
  message("Aggregating ALL cells per patient into a single bulk sample...")
  pb_list <- AggregateExpression(
    sub_obj, 
    assays = "RNA", 
    slot = "counts",
    group.by = c("orig.ident2", "Condition"), 
    return.seurat = FALSE 
  )
  
  counts_matrix <- pb_list$RNA
  
  # 3. Create Metadata
  # Since there are no cluster numbers, the columns are just "PatientID_Condition"
  all_cols <- colnames(counts_matrix)
  
  meta <- data.frame(pseudobulk_id = all_cols) %>%
    mutate(
      # Look for PMH or Healthy at the very end of the column name
      condition = str_extract(pseudobulk_id, "(PMH|Healthy)$"),
      sample_id = str_remove(pseudobulk_id, "_?(PMH|Healthy)$")
    )
  rownames(meta) <- meta$pseudobulk_id
  meta$condition <- factor(meta$condition, levels = c("Healthy", "PMH"))
  
  # Print replicates to terminal
  pmh_count <- sum(meta$condition == "PMH")
  healthy_count <- sum(meta$condition == "Healthy")
  message(paste("Total Biological Replicates -> PMH:", pmh_count, "| Healthy:", healthy_count))
  
  # Safety check: Ensure we have enough replicates for DESeq2 to run math
  if (pmh_count < 2 | healthy_count < 2) {
    message(paste("  -> SKIPPING DGE for", cell_type, "due to insufficient replicates (Requires >= 2 each)."))
    next
  }
  
  # 4. Run DESeq2
  message("Running DESeq2...")
  dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = meta,
    design = ~ condition
  )
  
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "PMH", "Healthy"))
  
  # 5. Format and Save Results
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
  
  csv_out <- file.path(out_dir, paste0("deseq2_ALL_COMBINED_", cell_type, ".csv"))
  write.csv(res_df, csv_out, row.names = FALSE)
  
  # 6. Visualization
  top_genes <- res_df %>%
    filter(significance != "Not Significant") %>% 
    group_by(significance) %>% 
    arrange(padj) %>% 
    slice_head(n = 20) %>% 
    ungroup() 
  
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
    labs(title = paste("All", tools::toTitleCase(cell_type), "Combined (PMH vs Healthy)"))
  
  plot_out <- file.path(out_dir, paste0("volcano_plot_ALL_COMBINED_", cell_type, ".png"))
  ggsave(filename = plot_out, plot = volcano_plot, width = 7, height = 6)
  
  message(paste("=== Success! Files saved to:", out_dir, "==="))
}

message("\n==================================================================")
message("=== All 'Combined' DGE pipelines executed successfully! ===")
message("==================================================================")

