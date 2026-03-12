suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(harmony)
  library(stringr)
})

# ==============================================================================
# 0. HELPERS: Resolution Tuning & Advanced Plotting
# ==============================================================================

tune_resolution <- function(obj, target = 6) {
  message(paste("Searching for resolution to achieve ~", target, "clusters..."))
  # Seurat v5 resolution search
  best_res <- 0.4 
  for (res in seq(0.1, 1.2, by = 0.1)) {
    obj <- FindClusters(obj, resolution = res, verbose = FALSE)
    n <- length(unique(Idents(obj)))
    if (n >= target) {
      best_res <- res
      break
    }
  }
  message(paste("Optimal resolution found:", best_res))
  return(obj)
}

generate_subset_plots <- function(sub_obj, prefix, out_dirs, option_name) {
  
  safe_save_plot <- function(plot_obj, base_filename, w = 10, h = 6) {
    pdf_path <- file.path(out_dirs$plots, paste0(base_filename, ".pdf"))
    tryCatch({ pdf(pdf_path, width = w, height = h); print(plot_obj); dev.off() }, error = function(e) { if (length(dev.list()) > 0) dev.off() })
    png_path <- file.path(out_dirs$plots, paste0(base_filename, ".png"))
    tryCatch({ png(png_path, width = w, height = h, units = "in", res = 300); print(plot_obj); dev.off() }, error = function(e) { if (length(dev.list()) > 0) dev.off() })
  }

  # 1. Condition Split UMAP
  p_cond <- DimPlot(sub_obj, split.by = "Condition", label = TRUE, label.size = 4, pt.size = 0.8) +
    ggtitle(paste(prefix, "-", option_name, "- Condition Split"))
  safe_save_plot(p_cond, paste0(prefix, "_umap_condition"), w = 12, h = 6)
  
  # 2. INDIVIDUAL SAMPLE UMAPs (Per Request)
  p_samples <- DimPlot(sub_obj, split.by = "orig.ident2", label = FALSE, pt.size = 0.6, ncol = 4) +
    ggtitle(paste(prefix, "-", option_name, "- Per Sample UMAPs"))
  safe_save_plot(p_samples, paste0(prefix, "_umap_per_sample_split"), w = 16, h = 8)
  
  # 3. PROPORTIONS BY SAMPLE BARPLOT (Per Request)
  prop_samp <- as.data.frame(prop.table(table(Idents(sub_obj), sub_obj$orig.ident2), margin = 2))
  colnames(prop_samp) <- c("Cluster", "Sample", "Proportion")
  
  p_bar_samp <- ggplot(prop_samp, aes(x = Sample, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = percent) +
    labs(title = paste(prefix, "-", option_name, "- Proportions per Patient"))
  safe_save_plot(p_bar_samp, paste0(prefix, "_proportions_per_sample"), w = 10, h = 6)

  # 4. Proportions by Condition (Simplified)
  prop_cond <- as.data.frame(prop.table(table(Idents(sub_obj), sub_obj$Condition), margin = 2))
  colnames(prop_cond) <- c("Cluster", "Condition", "Proportion")
  p_bar_cond <- ggplot(prop_cond, aes(x = Condition, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity") + theme_minimal() + scale_y_continuous(labels = percent)
  safe_save_plot(p_bar_cond, paste0(prefix, "_proportions_condition"), w = 6, h = 6)
}

# ==============================================================================
# OPTION 2: Harmony (V5 Robust)
# ==============================================================================
process_option2 <- function(seurat_obj, subset_clusters, prefix, out_base_dir) {
  message(paste("\n>>> OPTION 2: Harmony for", prefix))
  out_dir <- file.path(out_base_dir, paste0(prefix, "_Opt2_Harmony"))
  out_dirs <- list(processed = file.path(out_dir, "processed"), plots = file.path(out_dir, "plots"))
  lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

  # V5 CLEANING: Strip old scaled data/layers to prevent Assay validation errors
  sub_obj <- subset(seurat_obj, idents = subset_clusters)
  sub_obj <- DietSeurat(sub_obj, layers = "counts") 
  
  DefaultAssay(sub_obj) <- "RNA"
  sub_obj <- NormalizeData(sub_obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  sub_obj <- RunHarmony(sub_obj, group.by.vars = "orig.ident2", plot_convergence = FALSE)
  
  sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = 1:20)
  sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = 1:20)
  sub_obj <- tune_resolution(sub_obj, target = 6)
  
  sub_obj$Condition <- ifelse(grepl("HTY|UA", sub_obj$orig.ident2, ignore.case = TRUE), "Healthy", "PMH")
  generate_subset_plots(sub_obj, prefix, out_dirs, "Harmony")
  saveRDS(sub_obj, file.path(out_dirs$processed, "obj.rds"))
  return(sub_obj)
}

# ==============================================================================
# OPTION 3: Seurat Native (V5 IntegrateLayers CCA)
# ==============================================================================
process_option3 <- function(seurat_obj, subset_clusters, prefix, out_base_dir) {
  message(paste("\n>>> OPTION 3: Seurat V5 CCA for", prefix))
  out_dir <- file.path(out_base_dir, paste0(prefix, "_Opt3_SeuratNative"))
  out_dirs <- list(processed = file.path(out_dir, "processed"), plots = file.path(out_dir, "plots"))
  lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

  # V5 CLEANING
  sub_obj <- subset(seurat_obj, idents = subset_clusters)
  sub_obj <- DietSeurat(sub_obj, layers = "counts")

  DefaultAssay(sub_obj) <- "RNA"
  # V5 Layered Workflow
  sub_obj <- NormalizeData(sub_obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

  # Integrate using the V5 Layered API
  sub_obj <- IntegrateLayers(
    object = sub_obj, method = CCAIntegration, 
    orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE
  )

  sub_obj <- RunUMAP(sub_obj, reduction = "integrated.cca", dims = 1:20)
  sub_obj <- FindNeighbors(sub_obj, reduction = "integrated.cca", dims = 1:20)
  sub_obj <- tune_resolution(sub_obj, target = 6)
  
  sub_obj$Condition <- ifelse(grepl("HTY|UA", sub_obj$orig.ident2, ignore.case = TRUE), "Healthy", "PMH")
  generate_subset_plots(sub_obj, prefix, out_dirs, "SeuratNative_CCA")
  saveRDS(sub_obj, file.path(out_dirs$processed, "obj.rds"))
  return(sub_obj)
}

# ==============================================================================
# EXECUTION
# ==============================================================================
base_out_dir <- "/home/johan/output/skin_pmh/subset_comparisons"
pmh_obj <- readRDS("/home/johan/output/skin_pmh/TN.combined_dim30.rds")
Idents(pmh_obj) <- "seurat_clusters"

# Process Individual Cell Types
fib_clusters <- c("1", "5", "7", "10")
krt_clusters <- c("0", "2", "8", "11")
mac_clusters <- c("6")

# Execute pipelines
message("Processing Fibroblasts...")
fib_harmony <- process_option2(pmh_obj, fib_clusters, "Fibro", base_out_dir)
fib_native  <- process_option3(pmh_obj, fib_clusters, "Fibro", base_out_dir)

message("Processing Keratinocytes...")
krt_harmony <- process_option2(pmh_obj, krt_clusters, "KRT", base_out_dir)
krt_native  <- process_option3(pmh_obj, krt_clusters, "KRT", base_out_dir)

message("Processing Macrophages...")
mac_harmony <- process_option2(pmh_obj, mac_clusters, "Macro", base_out_dir)
mac_native  <- process_option3(pmh_obj, mac_clusters, "Macro", base_out_dir)

# Execute Combined Pipeline
combo_clusters <- c(fib_clusters, krt_clusters, mac_clusters)
message("Processing Combined Subset...")
combo_harmony <- process_option2(pmh_obj, combo_clusters, "Combo", base_out_dir)
combo_native  <- process_option3(pmh_obj, combo_clusters, "Combo", base_out_dir)