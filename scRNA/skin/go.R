library(gprofiler2)
library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. Define the universal function
run_go <- function(dge_data_location, out_dir, log_change, p_val_col, plot_name) {
  
  # Ensure output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Read data
  dge_data <- read.csv(dge_data_location)
  
  # Filter for UP and DOWN genes (Log2FC threshold AND Adjusted P-value < 0.05)
  genes_up <- dge_data$gene[
    !is.na(dge_data[[log_change]]) & dge_data[[log_change]] > 1 & 
    !is.na(dge_data[[p_val_col]]) & dge_data[[p_val_col]] < 0.05
  ]
  
  genes_down <- dge_data$gene[
    !is.na(dge_data[[log_change]]) & dge_data[[log_change]] < -1 & 
    !is.na(dge_data[[p_val_col]]) & dge_data[[p_val_col]] < 0.05
  ]
  
  # Create a named list to query both sets simultaneously
  query_list <- list()
  if(length(genes_up) > 0) query_list[["Upregulated"]] <- genes_up
  if(length(genes_down) > 0) query_list[["Downregulated"]] <- genes_down
  
  if(length(query_list) == 0) {
    message(paste("No significant genes passed the thresholds for", plot_name))
    return(NULL)
  }

  # Run gProfiler
  gostres <- gost(
    query = query_list,
    organism = "hsapiens",        
    sources = c("GO:BP")          
  )

  if(is.null(gostres)) {
    message(paste("No GO terms found for", plot_name))
    return(NULL)
  }

  res_full <- gostres$result
  
  # Flatten list columns (like 'parents' or 'intersections') so write.csv doesn't crash
  for (col in names(res_full)) {
    if (is.list(res_full[[col]])) {
      res_full[[col]] <- vapply(res_full[[col]], paste, collapse = ";", FUN.VALUE = character(1))
    }
  }
  
  csv_path <- file.path(out_dir, paste0(plot_name, "_GO_results.csv"))
  write.csv(res_full, file = csv_path, row.names = FALSE)
  message("Saved full GO results CSV to: ", csv_path)


  # Extract top 10 terms PER DIRECTION
  res <- gostres$result
  
  top_terms <- res %>%
    group_by(query) %>%
    arrange(p_value) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(
      log_p = -log10(p_value),
      plot_val = ifelse(query == "Downregulated", -log_p, log_p),
      term_name = factor(term_name, levels = unique(term_name[order(plot_val)]))
    )

  # Create a diverging bar plot
  go_plot <- ggplot(top_terms, aes(x = term_name, y = plot_val, fill = query)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    coord_flip() +  
    scale_fill_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "steelblue")) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) + 
    labs(
      x = "Biological Process", 
      y = "Directional -log10(p-value) \n(Left: Down, Right: Up)", 
      title = paste("GO:BP Enrichment -", plot_name),
      fill = "Direction"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")

  # Save the plot
  save_path <- file.path(out_dir, paste0(plot_name, "_GO_BP_plot.png"))
  ggsave(filename = save_path, plot = go_plot, width = 10, height = 8, bg = "white")
  
  message("Saved plot to: ", save_path)
  return(go_plot) 
}

# --- 2. Execute Functions ---

# Macrophage (DESeq2 usually uses 'padj' for adjusted p-value)
macrophage_plot <- run_go(
  dge_data_location = "/home/johan/output/skin_pmh/dge/macrophage/deseq2_results.csv",
  out_dir = "/home/johan/output/skin_pmh/go/macrophage",
  log_change = "log2FoldChange",
  p_val_col = "padj", 
  plot_name = "Macrophage"
)

# Fibroblast (Seurat usually uses 'p_val_adj')
fibroblast_plot <- run_go(
  dge_data_location = "/home/johan/output/skin_pmh/dge/fibroblast/fibroblast_markers_cluster1_vs_5_7_10.csv",
  out_dir = "/home/johan/output/skin_pmh/go/fibroblast",
  log_change = "avg_log2FC",
  p_val_col = "p_val_adj", 
  plot_name = "Fibroblast-Combine"
)

# Keratinocyte (Seurat usually uses 'p_val_adj')
keratinocyte_plot <- run_go(
  dge_data_location = "/home/johan/output/skin_pmh/dge/keratinocyte/keratinocyte_markers_cluster1_vs_5_7_10.csv",
  out_dir = "/home/johan/output/skin_pmh/go/keratinocyte",
  log_change = "avg_log2FC",
  p_val_col = "p_val_adj",
  plot_name = "Keratinocyte-Combine"
)

# Fibroblast (Seurat usually uses 'p_val_adj')
fibroblast_plot <- run_go(
  dge_data_location = "/home/johan/output/skin_pmh/dge/fibroblast/fibroblast_PMH_cluster1_vs_5_7_10.csv",
  out_dir = "/home/johan/output/skin_pmh/go/fibroblast",
  log_change = "avg_log2FC",
  p_val_col = "p_val_adj", 
  plot_name = "Fibroblast-PMH"
)

# Keratinocyte (Seurat usually uses 'p_val_adj')
keratinocyte_plot <- run_go(
  dge_data_location = "/home/johan/output/skin_pmh/dge/keratinocyte/keratinocyte_PMH_cluster1_vs_5_7_10.csv",
  out_dir = "/home/johan/output/skin_pmh/go/keratinocyte",
  log_change = "avg_log2FC",
  p_val_col = "p_val_adj",
  plot_name = "Keratinocyte-PMH"
)