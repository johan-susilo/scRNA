library(Seurat)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(tidyr)
library(ggrepel)
library(stringr)

# --- 1. Load Data and Set Conditions ---
pmh_obj <- readRDS("/home/johan/output/skin_pmh/TN.combined_dim30.rds")

pmh_obj$condition <- ifelse(
  grepl("HTY|UA", pmh_obj$orig.ident2, ignore.case = TRUE),
  "Healthy",
  "PMH"
)

# Verify condition assignment
print(table(pmh_obj$orig.ident2, pmh_obj$condition))

# === STEP 1: Global AggregateExpression ===
# We do this ONCE for the whole object to save time. 
# return.seurat = FALSE gives us the raw matrices.
print("Aggregating expression... This may take a moment.")
pb_list <- AggregateExpression(
  pmh_obj, 
  assays = "RNA", 
  slot = "counts",
  group.by = c("seurat_clusters", "orig.ident", "condition"),
  return.seurat = FALSE 
)

# Extract RNA count matrix
counts_matrix <- pb_list$RNA

# --- 2. Define the Reusable Pseudobulk Function ---

run_pseudobulk_dge <- function(cell_type_name, cluster_ids, base_out_dir, counts_mat) {
  
  print(paste("=========================================="))
  print(paste("Running Pseudobulk DGE for:", cell_type_name))
  print(paste("Clusters included:", paste(cluster_ids, collapse = ", ")))
  
  # Setup directory
  out_dir <- file.path(base_out_dir, cell_type_name)
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
  
  # === STEP 2: Select Specific Clusters (ROBUST FIX) ===
  all_cols <- colnames(counts_mat)
  
  # Extract the first chunk before the underscore (e.g., extracts "g6" or "6")
  cluster_prefixes <- str_extract(all_cols, "^[^_]+")
  
  # Strip any letters so "g6" becomes just "6"
  cluster_numbers <- str_remove_all(cluster_prefixes, "[a-zA-Z]")
  
  # Find columns where the cleaned number matches our target IDs
  target_cols <- all_cols[cluster_numbers %in% as.character(cluster_ids)]
  
  if (length(target_cols) == 0) {
    # If it fails again, print the actual column names so we can see what Seurat did
    print("Available columns look like this:")
    print(head(all_cols))
    stop(paste("No pseudobulk samples found for clusters:", paste(cluster_ids, collapse = ", ")))
  }
  
  # Subset the matrix
  counts_subset <- counts_mat[, target_cols, drop = FALSE]
  
  # === STEP 3: Create Metadata ===
  meta <- data.frame(pseudobulk_id = target_cols) %>%
    mutate(
      cluster = str_extract(pseudobulk_id, "^[^_]+"),
      condition = str_extract(pseudobulk_id, "[^_]+$"),
      sample_id = str_replace(pseudobulk_id, paste0("^", cluster, "_"), ""),
      sample_id = str_replace(sample_id, paste0("_", condition, "$"), "")
    )
  
  rownames(meta) <- meta$pseudobulk_id
  
  # Verify alignment
  counts_subset <- counts_subset[, rownames(meta), drop = FALSE]
  stopifnot(all(colnames(counts_subset) == rownames(meta))) 
  
  # === STEP 4: DESeq2 Pipeline ===
  meta$condition <- factor(meta$condition, levels = c("Healthy", "PMH"))
  
  # Account for multiple clusters in the design formula if necessary
  if (length(unique(meta$cluster)) > 1) {
    design_formula <- ~ cluster + condition
  } else {
    design_formula <- ~ condition
  }
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_subset,
    colData = meta,
    design = design_formula
  )
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "PMH", "Healthy"))
  
  # === STEP 5: Formatting and Saving Results ===
  res_df <- as.data.frame(res) %>% 
    filter(!is.na(padj)) %>% 
    mutate(gene = rownames(.)) %>%
    mutate(
      significance = case_when(
        padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
        padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    ) %>%
    select(gene, everything()) %>%
    arrange(padj)
  
  write.csv(res_df, file.path(out_dir, paste0("deseq2_", cell_type_name, ".csv")), row.names = FALSE)
  
  # === STEP 6: Visualization ===
  top_genes <- res_df %>%
    filter(significance != "Not Significant") %>% 
    group_by(significance) %>% 
    arrange(padj) %>% 
    slice_head(n = 20) %>% 
    ungroup() 
  
  # MA plot
  pdf(file.path(out_dir, paste0("MA_plot_", cell_type_name, ".pdf")), width = 6, height = 5)
  DESeq2::plotMA(res, main = paste("MA Plot:", cell_type_name, "(PMH vs Healthy)"), ylim = c(-8, 8))
  dev.off()
  
  # Volcano plot
  volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = top_genes, 
                    aes(label = gene), 
                    color = "black", 
                    box.padding = 0.5, 
                    max.overlaps = Inf) +
    scale_color_manual(values = c("Upregulated" = "red", 
                                  "Downregulated" = "blue", 
                                  "Not Significant" = "grey80")) +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", cell_type_name, "(PMH vs Healthy)"))
  
  ggsave(filename = file.path(out_dir, paste0("volcano_plot_", cell_type_name, ".png")),
         plot = volcano_plot,
         width = 6,
         height = 5)
  
  print(paste("Finished", cell_type_name, "- Files saved to:", out_dir))
  return(res_df)
}

# --- 3. Execute for Each Cell Type ---

base_dir <- "/home/johan/output/skin_pmh/dge"

# 1. Macrophage (Cluster 6 only)
macrophage_res <- run_pseudobulk_dge(
  cell_type_name = "macrophage", 
  cluster_ids = c("6"), 
  base_out_dir = base_dir, 
  counts_mat = counts_matrix
)

# 2. Fibroblasts (Clusters 1, 5, 7, 10 combined)
fibroblast_res <- run_pseudobulk_dge(
  cell_type_name = "fibroblast", 
  cluster_ids = c("1", "5", "7", "10"), 
  base_out_dir = base_dir, 
  counts_mat = counts_matrix
)

# 3. Keratinocytes (Clusters 0, 2, 8, 11 combined)
keratinocyte_res <- run_pseudobulk_dge(
  cell_type_name = "keratinocyte", 
  cluster_ids = c("0", "2", "8", "11"), 
  base_out_dir = base_dir, 
  counts_mat = counts_matrix
)