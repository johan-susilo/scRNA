library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel) # Required for the volcano plot labels

run_fibroblast_dge <- function(seurat_obj, ident.1, ident.2, output_prefix, out_dir, min.pct = 0.1, logfc.threshold = 0.25) {
  
  print(paste("Running DGE for:", output_prefix))
  
  # 1. Run FindMarkers
  markers <- FindMarkers(
    seurat_obj,
    ident.1 = ident.1,
    ident.2 = ident.2,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold
  )
  
  # 2. Format the dataframe
  markers$gene <- rownames(markers)
  markers <- markers[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
  
  # 3. Add significance column
  markers <- markers %>%
    mutate(
      significance = case_when(
        p_val_adj < 0.05 & avg_log2FC > 1 ~ "Upregulated",
        p_val_adj < 0.05 & avg_log2FC < -1 ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    ) 
  
  # 4. Save the CSV
  csv_path <- file.path(out_dir, paste0(output_prefix, ".csv"))
  write.csv(markers, csv_path, row.names = FALSE)
  
  # 5. Handle p-values of 0 for plotting
  markers_plot <- markers %>%
    mutate(
      p_val_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj),
      negLog10P = -log10(p_val_adj)
    )
  
  # 6. Get top genes for the volcano plot labels
  top_genes <- markers_plot %>%
    filter(significance != "Not Significant") %>%
    group_by(significance) %>%
    arrange(p_val_adj) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  # 7. Create and save the Volcano plot
  volcano_plot <- ggplot(markers_plot, aes(x = avg_log2FC, y = negLog10P, color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = top_genes, aes(label = gene), size = 3) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", output_prefix),
         x = "Log2 Fold Change",
         y = "-log10 Adjusted p-value") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") 
  
  pdf_path <- file.path(out_dir, paste0(output_prefix, "_volcano.pdf"))
  pdf(pdf_path, width = 7, height = 5)
  print(volcano_plot)
  dev.off()
  
  print(paste("Finished! Files saved to:", out_dir))
  
  # Return the dataframe so you can keep it in your R environment if needed
  return(markers)
}

pmh_obj <- readRDS("/home/johan/output/skin_pmh/TN.combined_dim30.rds")

pmh_obj$condition <- ifelse(
  grepl("HTY|UA", pmh_obj$orig.ident2, ignore.case = TRUE),
  "Healthy",
  "PMH"
)

# Setup directories and object
out_dir <- "/home/johan/output/skin_pmh/dge/fibroblast"
out_dir_krt <- "/home/johan/output/skin_pmh/dge/keratinocyte"
if (!dir.exists(out_dir_krt)) { dir.create(out_dir_krt, recursive = TRUE) }

# Assuming pmh_obj is already loaded and identities are set to "seurat_clusters"
Idents(pmh_obj) <- "seurat_clusters"

# 1. Run for all combined conditions
fibro_all <- run_fibroblast_dge(
  seurat_obj = pmh_obj, 
  ident.1 = "1", 
  ident.2 = c("5","7","10"), 
  output_prefix = "fibroblast_markers_cluster1_vs_5_7_10",
  out_dir = out_dir
)

# 2. Run for only the PMH condition
pmh_subset <- subset(pmh_obj, subset = condition == "PMH")
fibro_pmh <- run_fibroblast_dge(
  seurat_obj = pmh_subset, 
  ident.1 = "1", 
  ident.2 = c("5","7","10"), 
  output_prefix = "fibroblast_PMH_cluster1_vs_5_7_10",
  out_dir = out_dir
)

# 3. Run for only the Healthy condition (Bonus!)
healthy_subset <- subset(pmh_obj, subset = condition == "Healthy")
fibro_healthy <- run_fibroblast_dge(
  seurat_obj = healthy_subset, 
  ident.1 = "1", 
  ident.2 = c("5","7","10"), 
  output_prefix = "fibroblast_Healthy_cluster1_vs_5_7_10",
  out_dir = out_dir
)

# 1. Run for all combined conditions
krt_all <- run_fibroblast_dge(
  seurat_obj = pmh_obj, 
  ident.1 = "11", 
  ident.2 = c("0","2","8"), 
  output_prefix = "keratinocyte_markers_cluster11_vs_0_2_8",
  out_dir = out_dir_krt
)

# 2. Run for only the PMH condition
pmh_subset <- subset(pmh_obj, subset = condition == "PMH")
krt_pmh <- run_fibroblast_dge(
  seurat_obj = pmh_subset, 
  ident.1 = "11", 
  ident.2 = c("0","2","8"), 
  output_prefix = "keratinocyte_PMH_cluster11_vs_0_2_8",
  out_dir = out_dir_krt
)

# 3. Run for only the Healthy condition (Bonus!)
healthy_subset <- subset(pmh_obj, subset = condition == "Healthy")
krt_healthy <- run_fibroblast_dge(
  seurat_obj = healthy_subset, 
  ident.1 = "11", 
  ident.2 = c("0","2","8"), 
  output_prefix = "keratinocyte_Healthy_cluster11_vs_0_2_8",
  out_dir = out_dir_krt
)