suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(SingleR)
  library(SingleCellExperiment)
  library(dplyr)
})

set.seed(42)

# ==============================================================================
# 1. SETUP, FOLDER ORGANIZATION & LOAD DATA
# ==============================================================================
base_out_dir <- "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/fibroblast/haniffa_analysis/"

# Create highly organized subdirectories
dirs <- list(
  umaps_global  = file.path(base_out_dir, "1_UMAPs", "Global"),
  umaps_sample  = file.path(base_out_dir, "1_UMAPs", "Per_Sample"),
  val_vln       = file.path(base_out_dir, "2_Validation_Plots", "VlnPlots"),
  val_dot       = file.path(base_out_dir, "2_Validation_Plots", "DotPlots"),
  master_sum    = file.path(base_out_dir, "3_Summary_Plots"),
  mucin         = file.path(base_out_dir, "4_Mucin_ECM"),
  proportions   = file.path(base_out_dir, "5_Proportions")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

message("Loading mathematically cleaned fibroblast object...")
FB.subgroup <- readRDS("/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/fibroblast/processed/fibroblast_subset_processed.rds")

# Ensure Default Assay is RNA
DefaultAssay(FB.subgroup) <- "RNA"

save_dual_format <- function(plot_obj, file_path_no_ext, w = 10, h = 8) {
  pdf(paste0(file_path_no_ext, ".pdf"), width = w, height = h)
  print(plot_obj)
  dev.off()
  
  png(paste0(file_path_no_ext, ".png"), width = w, height = h, units = "in", res = 300)
  print(plot_obj)
  dev.off()
}

# ==============================================================================
# 2. THE DICTIONARIES & COLORS
# ==============================================================================
biological_dictionary <- list(
  # ==========================================
  # F1: SUPERFICIAL
  # ==========================================
  "F1_Superficial_Healthy" = c("APCDD1", "COL18A1", "COL23A1", "COL13A1", "COMP", "NKD2", "RSPO1", "AXIN2"),
  "F1_Superficial_Disease" = c("CRABP1", "CYP26B1", "TNFRSF21", "CXCL1", "WNT5A", "COL18A1", "COL23A1", "COL13A1", "NKD2", "AXIN2", "RSPO1"),

  # ==========================================
  # F2: UNIVERSAL / RETICULAR
  # ==========================================
  "F2_Reticular_Healthy" = c("CD34", "PI16", "MFAP5", "DPP4", "PCOLCE2", "LGR5", "SLPI", "CD70", "CTHRC1"),
  "F2_Reticular_Disease" = c("PI16", "DPP4", "PCOLCE2", "MFAP5", "CD70", "LGR5"),

  # ==========================================
  # F2/F3: PERIVASCULAR
  # ==========================================
  "F2_F3_Perivascular_Healthy" = c("CXCL12", "APOE", "EFEMP1", "APOC1", "C7", "PLA2G2A", "PPARG", "MYOC", "GDF10"),
  "F2_F3_Perivascular_Disease" = c("CXCL12", "APOE", "C7", "PLA2G2A", "EFEMP1", "GDF10", "MYOC"),

  # ==========================================
  # F3: FRC-LIKE (CCL19+)
  # ==========================================
  "F3_FRC_like_Healthy" = c("CCL19", "CD74", "CH25H", "TNFSF13B", "IL33", "IRF8", "IL15", "VCAM1", "HLA-DRB1", "HLA-DRA"),
  "F3_FRC_like_Disease" = c("CCL19", "CD74", "CH25H", "TNFSF13B", "IL33", "HLA-DRA", "IRF8", "COX4I2", "RBP5", "ADAMDEC1", "CXCL9", "CXCL10", "APOE", "CXCL12"),

  # ==========================================
  # F4: HAIR FOLLICLE ASSOCIATED
  # ==========================================
  "F4_HF_DPEP_Healthy" = c("DPEP1", "MYL4", "MEF2C", "COL11A1"),
  "F4_HF_DPEP_Disease" = c("MEF2C", "MYL4", "COL11A1", "POSTN", "DPEP1"),
  
  "F4_HF_TNN_Healthy" = c("TNN", "COCH", "TNMD", "MKX", "NRG3", "SLITRK6"),
  "F4_HF_TNN_Disease" = c("COCH", "CRABP1", "COL24A1", "RSPO4", "SLITRK6", "NRG3", "MKX", "TNMD"),
  
  "F4_HF_DP_Healthy" = c("CORIN", "HHIP", "BMP7", "WNT5A", "LEF1"),
  "F4_HF_DP_Disease" = c("CRABP1", "COL24A1", "RSPO4", "RSPO3", "BMP7", "WNT5A", "LEF1", "SOX18", "HHIP"),

  # ==========================================
  # F5: SCHWANN-LIKE
  # ==========================================
  "F5_Schwann_NGFR_Healthy" = c("NGFR", "ITGA6", "SCN7A", "CDH19", "CLDN1", "SFRP4", "TENM2"),
  "F5_Schwann_NGFR_Disease" = c("NGFR", "TM4SF1", "SFRP4", "ANGPTL7", "ITGA6", "CDH19", "CLDN1", "EBF2", "OLFM2", "SCN7A"),
  
  "F5_Schwann_RAMP1_Healthy" = c("RAMP1", "RELN", "PLEKHA6", "IGFBP2", "FGFBP2", "SCN7A"),
  "F5_Schwann_RAMP1_Disease" = c("RAMP1", "IGFBP2", "RELN", "COL26A1", "PLEKHA6", "FMO2", "FGFBP2"),

  # ==========================================
  # F6, F7, F8: MYOFIBROBLASTS (Disease Only)
  # ==========================================
  "F6_Inflammatory_Myo_Disease" = c("IL11", "IL24", "CXCL5", "CXCL6", "CXCL8", "MMP9", "WNT2", "COL10A1", "MMP1", "MMP3"),
  "F7_Myofibroblast" = c("ACTA2", "TAGLN", "CTHRC1", "RUNX2", "KIF26B", "SULF1", "ADAM12", "COL8A1", "LRRC15", "CCN4", "ASPN", "POSTN", "TNC", "COL3A1", "WNT2", "COL10A1"),
  "F8_Fascia_like_Myo_Disease" = c("ACAN", "ITGA10", "CDH2", "DPP4", "CCN3"),
  
  # ==========================================
  # FASCIA
  # ==========================================
  "F_Fascia_Disease" = c("ITGA10", "CCN3", "DPP4", "CDH13", "PRG4", "CRTAC1", "PCOLCE2", "LGR5")
)

validation_dictionary <- biological_dictionary
validation_dictionary[["QC_Stress_Markers"]] <- c("MT2A", "MT1M", "MT1X", "HSP90AA1", "JUNB", "GADD45B", "IER3")

f1_f8_colors <- c(
  # F1: Superficial (Blues)
  "F1_Superficial_Healthy" = "#A6CEE3", # Light Blue
  "F1_Superficial_Disease" = "#1F77B4", # Dark Blue
  
  # F2: Universal / Reticular (Greens)
  "F2_Reticular_Healthy" = "#B2DF8A", # Light Green
  "F2_Reticular_Disease" = "#33A02C", # Dark Green
  
  # F2/F3: Perivascular (Oranges)
  "F2_F3_Perivascular_Healthy" = "#FDBF6F", # Light Orange
  "F2_F3_Perivascular_Disease" = "#FF7F00", # Dark Orange
  
  # F3: FRC-like (Purples)
  "F3_FRC_like_Healthy" = "#CAB2D6", # Light Purple
  "F3_FRC_like_Disease" = "#6A3D9A", # Dark Purple
  
  # F4: Hair Follicle (Yellows/Browns/Reds)
  "F4_HF_DPEP_Healthy" = "#FFFF99", # Pale Yellow
  "F4_HF_DPEP_Disease" = "#B15928", # Brown
  "F4_HF_TNN_Healthy"  = "#FFEDA0", # Sand
  "F4_HF_TNN_Disease"  = "#D95F0E", # Rust
  "F4_HF_DP_Healthy"   = "#E31A1C", # Red 
  "F4_HF_DP_Disease"   = "#bd0026", # Dark Red
  
  # F5: Schwann-like (Pinks/Magentas)
  "F5_Schwann_NGFR_Healthy"  = "#FBB4AE", # Light Pink
  "F5_Schwann_NGFR_Disease"  = "#E7298A", # Magenta
  "F5_Schwann_RAMP1_Healthy" = "#F1EEF6", # Very pale purple/pink
  "F5_Schwann_RAMP1_Disease" = "#CE1256", # Deep Cherry
  
  # F6, F7, F8: Disease-Specific Myofibroblasts
  "F6_Inflammatory_Myo_Disease" = "#006D2C", # Dark Cyan/Green
  "F7_Myofibroblast"    = "#c73333", # Almost Black
  "F8_Fascia_like_Myo_Disease"  = "#810F7C", # Dark Eggplant
  
  # Fascia
  "F_Fascia_Disease" = "#8B0000" # Maroon
)
# ==============================================================================
# 3. THE F1-F8 SINGLE-R AUTO-ANNOTATION ENGINE (Single-Cell & Sorted)
# ==============================================================================
message("Running Single-Cell Level SingleR Auto-Annotation...")

build_custom_singler_reference <- function(marker_dictionary, gene_universe) {
  ref_mat <- sapply(names(marker_dictionary), function(label) {
    genes <- intersect(marker_dictionary[[label]], gene_universe)
    vec <- setNames(rep(0, length(gene_universe)), gene_universe)
    vec[genes] <- 1
    vec
  })
  ref_sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(ref_mat)))
  colData(ref_sce)$label <- colnames(ref_mat)
  return(ref_sce)
}

gene_universe <- rownames(FB.subgroup)
ref_sce <- build_custom_singler_reference(biological_dictionary, gene_universe)
test_sce <- as.SingleCellExperiment(FB.subgroup, assay = "RNA")

# Run SingleR on EVERY INDIVIDUAL CELL (Ignores Seurat Resolution!)
singler_res <- SingleR(
  test = test_sce,
  ref = ref_sce,
  labels = ref_sce$label,
  # clusters = FB.subgroup$seurat_clusters,  <-- REMOVED to force single-cell level annotation
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  prune = FALSE 
)

# Map labels directly back to the individual cells
FB.subgroup$SingleR_label <- singler_res$labels

# --- MATHEMATICAL SORTING ---
# Define the strict, logical biological order from the Atlas
f1_f8_order <- c(
  "F1_Superficial_Healthy", "F1_Superficial_Disease",
  "F2_Reticular_Healthy", "F2_Reticular_Disease",
  "F2_F3_Perivascular_Healthy", "F2_F3_Perivascular_Disease",
  "F3_FRC_like_Healthy", "F3_FRC_like_Disease",
  "F4_HF_DPEP_Healthy", "F4_HF_DPEP_Disease",
  "F4_HF_TNN_Healthy", "F4_HF_TNN_Disease",
  "F4_HF_DP_Healthy", "F4_HF_DP_Disease",
  "F5_Schwann_NGFR_Healthy", "F5_Schwann_NGFR_Disease",
  "F5_Schwann_RAMP1_Healthy", "F5_Schwann_RAMP1_Disease",
  "F6_Inflammatory_Myo_Disease",
  "F7_Myofibroblast",
  "F8_Fascia_like_Myo_Disease",
  "F_Fascia_Disease"
)

# Convert the text labels into an ordered factor. 
active_levels <- intersect(f1_f8_order, unique(FB.subgroup$SingleR_label))

# SAFETY CHECK: Print any unmapped cells so we know if there is a typo!
unmapped <- setdiff(unique(FB.subgroup$SingleR_label), active_levels)
if(length(unmapped) > 0) {
  message("\nWARNING: The following clusters were not mapped in f1_f8_order:")
  print(unmapped)
}

FB.subgroup$SingleR_label <- factor(FB.subgroup$SingleR_label, levels = active_levels)
Idents(FB.subgroup) <- "SingleR_label"

message("\n--- Biological Clusters Assigned & Sorted ---")
print(table(Idents(FB.subgroup)))
message("------------------------------------\n")

# --- 3A. GLOBAL UMAP ---
p_annotated_umap <- DimPlot(FB.subgroup, label = TRUE, repel = TRUE, cols = f1_f8_colors) + 
  ggtitle("Single-Cell Annotated F1-F8 Fibroblasts")
save_dual_format(p_annotated_umap, file.path(dirs$umaps_global, "UMAP_Annotated_F1_F8"), w = 9, h = 7 )

# --- ⚡ BOLT OPTIMIZATION: SPLIT ONCE ---
# Splitting the object outside the loop to avoid O(n) Seurat subset calls.
# Called immediately before the loop to ensure metadata is completely up-to-date.
FB.split_by_sample <- SplitObject(FB.subgroup, split.by = "orig.ident2")

# --- 3B. PER-SAMPLE UMAPS ---
message("Generating sample-specific UMAPs...")
for (sample_id in names(FB.split_by_sample)) {
  sample_obj <- FB.split_by_sample[[sample_id]]
  
  if (ncol(sample_obj) > 0) {
    p_sub_umap <- DimPlot(sample_obj, label = TRUE, repel = TRUE, cols = f1_f8_colors) +
      ggtitle(paste("F1-F8 Annotation - Sample:", sample_id))
    save_dual_format(p_sub_umap, file.path(dirs$umaps_sample, paste0("UMAP_Annotated_", sample_id)), w = 8, h = 6)
  }
}

# ==============================================================================
# 4. HANIFFA F1-F8 MARKER VISUALIZATION (Automated Loop)
# ==============================================================================
message("Generating Validation Plots (VlnPlots & DotPlots)...")

for (category_name in names(validation_dictionary)) {
  genes_to_plot <- validation_dictionary[[category_name]]
  valid_genes <- intersect(genes_to_plot, rownames(FB.subgroup))
  
  if (length(valid_genes) == 0) next
  
  # --- SEURAT V5 FIX: Use layer="data" instead of slot="data" ---
  exp_data <- GetAssayData(FB.subgroup, assay = "RNA", layer = "data")
  valid_genes <- valid_genes[rowSums(exp_data[valid_genes, , drop = FALSE]) > 0]
  
  if (length(valid_genes) == 0) {
    message(paste("   -> Skipping", category_name, "(All genes have 0 expression)"))
    next
  }
  
  # DotPlots -> Saved to DotPlot Folder
  tryCatch({
    p_dot <- DotPlot(FB.subgroup, features = valid_genes) + 
      RotatedAxis() + 
      ggtitle(paste("Fibroblast -", category_name))
    ggsave(file.path(dirs$val_dot, paste0("DotPlot_", category_name, ".pdf")), plot = p_dot, width = max(6, length(valid_genes)*0.5 + 2), height = 6)
  }, error = function(e) message(paste("   -> Error generating DotPlot for", category_name)))
  
  # VlnPlots -> Saved to VlnPlot Folder
  tryCatch({
    p_vln <- VlnPlot(FB.subgroup, features = valid_genes, stack = TRUE, flip = TRUE, cols = f1_f8_colors) +
      theme(legend.position = "none") +
      ggtitle(paste("Fibroblast -", category_name))
    ggsave(file.path(dirs$val_vln, paste0("VlnPlot_", category_name, ".pdf")), plot = p_vln, width = 8, height = max(6, length(valid_genes)*1.5))
  }, error = function(e) message(paste("   -> Error generating VlnPlot for", category_name)))
}

# ==============================================================================
# 5. PUBLICATION MASTER DOTPLOT
# ==============================================================================
message("Generating Publication-Ready Master DotPlot...")

paper_signature_genes <- c(
  "APCDD1", "COL18A1", "CRABP1", "CYP26B1", "WNT5A",
  "PI16", "CD34", "MFAP5", "KLF5", "DPP4", "PPARG", "CXCL12",
  "CCL19", "CD74", "CXCL9", "DPEP1", "TNN", "CORIN", "RAMP1", "NGFR",
  "IL11", "CXCL8", "MMP1", "ACTA2", "COL11A1", "COMP", "LRRC15", "ACAN", "ITGA10", "THBS4"
)
valid_paper_genes <- intersect(paper_signature_genes, rownames(FB.subgroup))

master_dotplot <- DotPlot(FB.subgroup, features = valid_paper_genes, dot.scale = 6) +
  theme_minimal() +
  RotatedAxis() +
  scale_color_gradientn(colors = c("lightgrey", "blue", "darkred")) +
  labs(
    title = "Fibroblast Subpopulation Signatures in PMH",
    x = "Key Marker Genes",
    y = "Identified Subclusters"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

save_dual_format(master_dotplot, file.path(dirs$master_sum, "Publication_Master_DotPlot"), w = 14, h = 8)

# ==============================================================================
# 6. MUCIN & ECM GENE ANALYSIS
# ==============================================================================
message("Running Mucin & Tissue Remodeling Gene Analysis...")

mucin_ecm_genes <- c("MUC1", "HAS1", "HAS2", "MMP1", "MUC12", "HAS2", "HAS1", "HAS3", "VCAN", "FN1", "CEMIP", "HYAL1", "HYAL2", "CTGF", "TGFBI", "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "SPARC", "POSTN", "ACTA2", "TAGLN", "LOX", "LOXL2")
available_mucin <- intersect(mucin_ecm_genes, rownames(FB.subgroup))

if(length(available_mucin) > 0) {
  
  num_samples <- length(unique(FB.subgroup$orig.ident2))
  dynamic_width <- max(10, num_samples * 4) 
  
  for (gene in available_mucin) {
    p_mucin_split <- FeaturePlot(
      FB.subgroup, 
      features = gene, 
      split.by = "orig.ident2", 
      pt.size = 0.5, 
      order = FALSE,              
      label = FALSE              
    ) + 
      patchwork::plot_annotation(title = paste("Expression of", gene, "Across Samples")) &
      theme(legend.position = "right")
    
    save_dual_format(p_mucin_split, file.path(dirs$mucin, paste0("Mucin_UMAP_Split_", gene)), w = dynamic_width, h = 5)
  }
  
  mucin_dotplot <- DotPlot(FB.subgroup, features = available_mucin, dot.scale = 8) +
    theme_minimal() +
    RotatedAxis() +
    scale_color_gradientn(colors = c("lightgrey", "blue", "darkred")) +
    labs(
      title = "Mucin & ECM Production by Fibroblast Subtype",
      x = "Target Genes",
      y = "Fibroblast Subcluster"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.text.x = element_text(face = "italic", color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 12)
    )
  
  save_dual_format(mucin_dotplot, file.path(dirs$mucin, "Mucin_DotPlot_Summary"), w = 10, h = 7)

  # --- ⚡ BOLT OPTIMIZATION: SPLIT ONCE ---
  FB.split_by_sample <- SplitObject(FB.subgroup, split.by = "orig.ident2")
  for (sample_id in names(FB.split_by_sample)) {
    sample_subset <- FB.split_by_sample[[sample_id]]
    if (ncol(sample_subset) > 0) {
      p_sample_dotplot <- DotPlot(sample_subset, features = available_mucin, dot.scale = 8) +
        theme_minimal() +
        RotatedAxis() +
        scale_color_gradientn(colors = c("lightgrey", "blue", "darkred")) +
        labs(
          title = paste("Mucin & ECM Production - Sample:", sample_id),
          x = "Target Genes",
          y = "Fibroblast Subcluster"
        ) +
        theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.text.x = element_text(face = "italic", color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12)
        )
      save_dual_format(p_sample_dotplot, file.path(dirs$mucin, paste0("Mucin_DotPlot_Sample_", sample_id)), w = 10, h = 7)
    }
  }
}

# ==============================================================================
# 7. PROPORTION BARPLOTS
# ==============================================================================
create_proportion_barplot <- function(seurat_obj, group_col = "Condition", label_col = "SingleR_label", plot_title = "Cell Proportions") {
  prop_df <- seurat_obj@meta.data %>%
    group_by(.data[[group_col]], .data[[label_col]]) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    group_by(.data[[group_col]]) %>%
    mutate(Proportion = Count / sum(Count))
  
  p <- ggplot(prop_df, aes(x = .data[[group_col]], y = Proportion, fill = .data[[label_col]])) +
    geom_col(position = "fill", color = "black", linewidth = 0.2) +
    theme_classic() +
    scale_fill_manual(values = f1_f8_colors) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = plot_title, x = group_col, y = "Percentage of Cells", fill = "Subtype") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  return(p)
}

create_proportion_barplot_filtered <- function(seurat_obj, group_col = "Condition", label_col = "SingleR_label", exclude_clusters = NULL, plot_title = "Cell Proportions") {
  meta_df <- seurat_obj@meta.data
  if (!is.null(exclude_clusters)) {
    meta_df <- meta_df[!meta_df[[label_col]] %in% exclude_clusters, ]
  }
  prop_df <- meta_df %>%
    group_by(.data[[group_col]], .data[[label_col]]) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    group_by(.data[[group_col]]) %>%
    mutate(Proportion = Count / sum(Count))
  
  p <- ggplot(prop_df, aes(x = .data[[group_col]], y = Proportion, fill = .data[[label_col]])) +
    geom_col(position = "fill", color = "black", linewidth = 0.2) +
    theme_classic() +
    scale_fill_manual(values = f1_f8_colors) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = plot_title, x = group_col, y = "Relative Percentage (Filtered)", fill = "Subtype") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  return(p)
}

bar_cond <- create_proportion_barplot(FB.subgroup, group_col = "Condition", plot_title = "Fibroblast Subtypes by Condition")
save_dual_format(bar_cond, file.path(dirs$proportions, "Proportion_Barplot_Condition"), w = 8, h = 6)

bar_samp <- create_proportion_barplot(FB.subgroup, group_col = "orig.ident2", plot_title = "Fibroblast Subtypes by Sample")
save_dual_format(bar_samp, file.path(dirs$proportions, "Proportion_Barplot_Sample"), w = 10, h = 6)

# FIX: Filter out the newly updated Disease state names
clusters_to_remove <- c("F6_Inflammatory_Myo_Disease", 
                        "F7_Myofibroblast", 
                        "F8_Fascia_like_Myo_Disease")

bar_samp_filtered <- create_proportion_barplot_filtered(
  FB.subgroup, 
  group_col = "orig.ident2", 
  exclude_clusters = clusters_to_remove,
  plot_title = "Sample Variances (Disease Signatures Excluded)"
)
save_dual_format(bar_samp_filtered, file.path(dirs$proportions, "Proportion_Barplot_Sample_Filtered"), w = 10, h = 6)

message("Saving final annotated RDS...")
saveRDS(FB.subgroup, file.path("/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/fibroblast_0.6/processed", "fibroblast_annotated_final.rds"))

# ==============================================================================
# 8. MACRO-LINEAGE & ACTIVATION STATE SUMMARY
# ==============================================================================
message("Grouping clusters into Macro Lineages for summary plots...")

# 1. Create a "Macro Lineage" column (Combines the micro-clusters into their parent families)
FB.subgroup$Macro_Lineage <- dplyr::case_when(
  grepl("^F1", FB.subgroup$SingleR_label) ~ "F1_Superficial",
  grepl("^F2_Reticular", FB.subgroup$SingleR_label) ~ "F2_Reticular",
  grepl("^F2_F3", FB.subgroup$SingleR_label) ~ "F2/F3_Perivascular",
  grepl("^F3", FB.subgroup$SingleR_label) ~ "F3_FRC_like",
  grepl("^F4", FB.subgroup$SingleR_label) ~ "F4_HairFollicle",
  grepl("^F5", FB.subgroup$SingleR_label) ~ "F5_Schwann",
  grepl("^F6", FB.subgroup$SingleR_label) ~ "F6_Inflammatory_Myo",
  grepl("^F7", FB.subgroup$SingleR_label) ~ "F7_Myofibroblast",
  grepl("^F8", FB.subgroup$SingleR_label) ~ "F8_Fascia_like_Myo",
  grepl("^F_Fascia", FB.subgroup$SingleR_label) ~ "F_Fascia",
  TRUE ~ "Unknown"
)

# Lock in the logical order of the lineages
macro_order <- c("F1_Superficial", "F2_Reticular", "F2/F3_Perivascular", "F3_FRC_like", 
                 "F4_HairFollicle", "F5_Schwann", "F6_Inflammatory_Myo", 
                 "F7_Myofibroblast", "F8_Fascia_like_Myo", "F_Fascia")
FB.subgroup$Macro_Lineage <- factor(FB.subgroup$Macro_Lineage, levels = macro_order)

# 2. Create an "Activation State" column
FB.subgroup$Activation_State <- ifelse(grepl("Healthy", FB.subgroup$SingleR_label), 
                                       "Resting (Healthy)", "Activated (Disease)")
FB.subgroup$Activation_State <- factor(FB.subgroup$Activation_State, levels = c("Resting (Healthy)", "Activated (Disease)"))

# Define a clean, distinct 10-color palette for the Macro Lineages
macro_colors <- c(
  "F1_Superficial"      = "#1F77B4", # Blue
  "F2_Reticular"        = "#2CA02C", # Green
  "F2/F3_Perivascular"  = "#FF7F00", # Orange
  "F3_FRC_like"         = "#9467BD", # Purple
  "F4_HairFollicle"     = "#E31A1C", # Red
  "F5_Schwann"          = "#E7298A", # Pink
  "F6_Inflammatory_Myo" = "#17BECF", # Cyan
  "F7_Myofibroblast"    = "#7F7F7F", # Grey
  "F8_Fascia_like_Myo"  = "#8C564B", # Brown
  "F_Fascia"            = "#8B0000"  # Maroon
)

# --- PLOT A: Overall Macro-Lineage Composition by Condition ---
# This plot answers: "Are PMH samples mostly Myofibroblasts, while Healthy samples are mostly F1/F2?"
prop_macro <- FB.subgroup@meta.data %>%
  group_by(Condition, Macro_Lineage) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Condition) %>%
  mutate(Proportion = Count / sum(Count))

p_macro_bar <- ggplot(prop_macro, aes(x = Condition, y = Proportion, fill = Macro_Lineage)) +
  geom_col(position = "fill", color = "black", linewidth = 0.3) +
  theme_classic() +
  scale_fill_manual(values = macro_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Major Fibroblast Lineages by Condition", x = "Disease State", y = "Percentage of Cells", fill = "Macro Lineage") +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

save_dual_format(p_macro_bar, file.path(dirs$proportions, "Summary_MacroLineage_Barplot"), w = 8, h = 6)

# --- PLOT B: Activation Ratio (Healthy vs Disease State) ---
# This plot looks ONLY at the cells that exist in both states (F1-F5) to see what % are activated.
plastic_cells <- FB.subgroup@meta.data %>%
  filter(!Macro_Lineage %in% c("F6_Inflammatory_Myo", "F7_Myofibroblast", "F8_Fascia_like_Myo", "F_Fascia"))

prop_activation <- plastic_cells %>%
  group_by(Condition, Activation_State) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Condition) %>%
  mutate(Proportion = Count / sum(Count))

p_activation <- ggplot(prop_activation, aes(x = Condition, y = Proportion, fill = Activation_State)) +
  geom_col(position = "fill", color = "black", linewidth = 0.3, width = 0.6) +
  theme_classic() +
  scale_fill_manual(values = c("Resting (Healthy)" = "#A6CEE3", "Activated (Disease)" = "#E31A1C")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Fibroblast Activation Ratio (F1-F5 Only)", 
       subtitle = "Shows the shift from Resting to Activated states across conditions",
       x = "Patient Condition", y = "Percentage of F1-F5 Cells", fill = "Cellular State") +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, face = "italic"))

save_dual_format(p_activation, file.path(dirs$proportions, "Summary_ActivationState_Barplot"), w = 7, h = 6)

# Calculate the proportions using the grouped 'Macro_Lineage' instead of the 22 micro-labels
prop_samp_macro <- FB.subgroup@meta.data %>%
  group_by(orig.ident2, Macro_Lineage) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(orig.ident2) %>%
  mutate(Proportion = Count / sum(Count))

# Generate the cleaned-up plot
p_samp_macro <- ggplot(prop_samp_macro, aes(x = orig.ident2, y = Proportion, fill = Macro_Lineage)) +
  geom_col(position = "fill", color = "black", linewidth = 0.3) +
  theme_classic() +
  scale_fill_manual(values = macro_colors) + # Uses the clean 10-color palette we defined in Section 8
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Major Fibroblast Lineages by Sample", 
    x = "Sample ID", 
    y = "Percentage of Cells", 
    fill = "Macro Lineage"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

# Save the beautiful version
save_dual_format(p_samp_macro, file.path(dirs$proportions, "Clean_MacroLineage_BySample"), w = 10, h = 6)

# ==============================================================================
# 9. MACRO-LINEAGE UMAP VISUALIZATION
# ==============================================================================
message("Generating clean Macro Lineage UMAPs...")

# --- PLOT C: Global UMAP by Macro Lineage ---
# This overrides the active identity and colors strictly by the 10 major families
p_macro_umap <- DimPlot(FB.subgroup, group.by = "Macro_Lineage", label = TRUE, repel = TRUE, cols = macro_colors) +
  ggtitle("Global UMAP: Major Fibroblast Lineages") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

save_dual_format(p_macro_umap, file.path(dirs$umaps_global, "UMAP_MacroLineage_Global"), w = 9, h = 7)

# --- PLOT D: Per-Sample UMAPs by Macro Lineage ---
message("Generating sample-specific Macro Lineage UMAPs...")
# --- ⚡ BOLT OPTIMIZATION: SPLIT ONCE ---
# Re-splitting here because Macro_Lineage metadata was added since the last split
FB.split_by_sample <- SplitObject(FB.subgroup, split.by = "orig.ident2")
for (sample_id in names(FB.split_by_sample)) {
  sample_obj <- FB.split_by_sample[[sample_id]]
  
  if (ncol(sample_obj) > 0) {
    p_sub_macro_umap <- DimPlot(sample_obj, group.by = "Macro_Lineage", label = TRUE, repel = TRUE, cols = macro_colors) +
      ggtitle(paste("Major Lineages - Sample:", sample_id)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
      
    save_dual_format(p_sub_macro_umap, file.path(dirs$umaps_sample, paste0("UMAP_MacroLineage_", sample_id)), w = 8, h = 6)
  }
}

message("=== All Summary Plots and UMAPs Complete! ===")

# ==============================================================================
# 10. HIGHLIGHT UMAPS: ISOLATING THE MYOFIBROBLASTS
# ==============================================================================
message("Generating Highlight UMAPs for F7 and F8...")

# 1. Remove F_Fascia completely to clean the plot
FB.highlight <- subset(FB.subgroup, Macro_Lineage != "F_Fascia")

# 2. Build the first grouping: F1-F6 Combined, F7 Separate, F8 Separate
FB.highlight$Highlight_Group_1 <- dplyr::case_when(
  FB.highlight$Macro_Lineage %in% c("F1_Superficial", "F2_Reticular", "F2/F3_Perivascular", "F3_FRC_like", "F4_HairFollicle", "F5_Schwann", "F6_Inflammatory_Myo") ~ "F1-F6 (Healthy & Inflammatory)",
  FB.highlight$Macro_Lineage == "F7_Myofibroblast" ~ "F7_Myofibroblastfibroblast",
  FB.highlight$Macro_Lineage == "F8_Fascia_like_Myo" ~ "F8_Fascia-like_Myofibroblast",
  TRUE ~ "Other"
)
FB.highlight$Highlight_Group_1 <- factor(FB.highlight$Highlight_Group_1, 
                                         levels = c("F1-F6 (Healthy & Inflammatory)", "F7_Myofibroblastfibroblast", "F8_Fascia-like_Myofibroblast"))

# 3. Build the second grouping: F1-F6 Combined vs F7+F8 Combined
FB.highlight$Highlight_Group_2 <- dplyr::case_when(
  FB.highlight$Highlight_Group_1 == "F1-F6 (Healthy & Inflammatory)" ~ "F1-F6 (Healthy & Inflammatory)",
  FB.highlight$Macro_Lineage %in% c("F7_Myofibroblast", "F8_Fascia_like_Myo") ~ "F7+F8 (All Terminal Myofibroblasts)",
  TRUE ~ "Other"
)
FB.highlight$Highlight_Group_2 <- factor(FB.highlight$Highlight_Group_2, 
                                         levels = c("F1-F6 (Healthy & Inflammatory)", "F7+F8 (All Terminal Myofibroblasts)"))

# Define colors (Grey for background, Bold for targets)
color_highlight_1 <- c(
  "F1-F6 (Healthy & Inflammatory)" = "#D9D9D9",  # Light Grey
  "F7_Myofibroblastfibroblast"     = "#252525",  # Almost Black
  "F8_Fascia-like_Myofibroblast"   = "#810F7C"   # Dark Eggplant
)

color_highlight_2 <- c(
  "F1-F6 (Healthy & Inflammatory)"      = "#D9D9D9",  # Light Grey
  "F7+F8 (All Terminal Myofibroblasts)" = "#E31A1C"   # Bright Red
)

# --- PLOT E: F1-F6 Combined, F7 and F8 Separate ---
p_highlight_1 <- DimPlot(FB.highlight, group.by = "Highlight_Group_1", label = FALSE, pt.size = 0.6, cols = color_highlight_1) +
  ggtitle("Myofibroblast Isolation (F7 & F8 Split)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

save_dual_format(p_highlight_1, file.path(dirs$umaps_global, "UMAP_Highlight_F7_F8_Split"), w = 9, h = 7)

# --- PLOT F: Binary Plot (F1-F6 vs F7+F8) ---
p_highlight_2 <- DimPlot(FB.highlight, group.by = "Highlight_Group_2", label = FALSE, pt.size = 0.6, cols = color_highlight_2) +
  ggtitle("Global Myofibroblast Emergence (F7+F8 Combined)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

save_dual_format(p_highlight_2, file.path(dirs$umaps_global, "UMAP_Highlight_F7_F8_Combined"), w = 9, h = 7)

message("=== Highlight UMAPs Complete! ===")

out_dir <- "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/fibroblast/processed/fibroblast_annotated_full.rds"

dir.create(dirname(out_dir), recursive = TRUE, showWarnings = FALSE)

message("Saving final annotated RDS...")
saveRDS(FB.subgroup, out_dir)

message("=== Annotation Complete! ===")