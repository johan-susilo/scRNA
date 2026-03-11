library(Seurat)

library(dplyr)

library(ggplot2)

library(clusterProfiler)

library(org.Hs.eg.db)  # for human genes; use org.Mm.eg.db for mouse

library(gprofiler2)



pmh_obj <- readRDS("/home/johan/output/skin_pmh/TN.combined_dim30.rds")

#check how many cluster there are in pmh_obj

levels(Idents(pmh_obj))

Idents(pmh_obj) <- "seurat_clusters"



#store original cluster info

pmh_obj$global_cluster <- Idents(pmh_obj)



#check dge of fibroblast cluster 1 vs all other fibroblast clusters (5,7,10)

fibro_markers <- FindMarkers(

  pmh_obj,

  ident.1 = "1",                 # fibroblast 1

  ident.2 = c("5","7","10"),  # all other fibroblast clusters

  min.pct = 0.1,  # only genes expressed in ≥10% cells

  test.use = "wilcox",

  logfc.threshold = 0.25

)

write.csv(fibro_markers, "/home/johan/output/skin_pmh/fibroblast/fibroblast_markers_cluster1_vs_5_7_10.csv")







#go analysis

dge_data <- read.csv("/home/johan/output/skin_pmh/fibroblast/fibroblast_markers_cluster1_vs_5_7_10.csv")

genes <- dge_data$X[dge_data$avg_log2FC > 1]



gostres <- gost(

  query = genes,

  organism = "hsapiens",        # human genes

  sources = c("GO:BP")          # restrict to Biological Process

)



head(gostres$result[, c("term_name","p_value","intersection_size")])

# Take top 10 terms

top_terms <- gostres$result[1:10, ]



ggplot(top_terms, aes(x = reorder(term_name, -p_value), y = -log10(p_value))) +

  geom_bar(stat = "identity", fill = "steelblue") +

  coord_flip() +  # horizontal bars

  labs(x = "Biological Process", y = "-log10(p-value)", title = "GO:BP Enrichment") +

  theme_minimal(base_size = 12)





#only take fibroblasts clusters: 1,5,7,10

fibro_obj <- subset(pmh_obj, idents = c("1","5","7","10"))

table(Idents(fibro_obj))



saveRDS(fibro_obj, "/home/johan/output/skin_pmh/fibroblasts_only.rds")



#switch back to the actual raw gene expression values so you can detect real biological differences.

DefaultAssay(fibro_obj) <- "RNA"

fibro_obj <- NormalizeData(fibro_obj)

fibro_obj <- FindVariableFeatures(fibro_obj, selection.method = "vst", nfeatures = 2000)

fibro_obj <- ScaleData(fibro_obj)

fibro_obj <- RunPCA(fibro_obj)

fibro_obj <- RunUMAP(fibro_obj, dims = 1:20)

fibro_obj <- FindNeighbors(fibro_obj, dims = 1:20)

fibro_obj <- FindClusters(fibro_obj, resolution = 0.4)





fibro_obj$Condition <- ifelse(grepl("HTY", fibro_obj$orig.ident2, ignore.case = TRUE),

                              "Healthy",

                              "PMH")



# Verify it exists now

table(fibro_obj$Condition)

table(fibro_obj$orig.ident, fibro_obj$Condition)

colnames(fibro_obj@meta.data)





# Visual Proof: Healthy vs. Disease

comparison_fibroblast <- DimPlot(fibro_obj, split.by = "Condition", label = TRUE, label.size = 5, pt.size = 0.5) +

  ggtitle("Figure 1: Identification of PMH-Specific Fibroblast Population")



# Update output paths to the fibroblast subfolder

fibro_out_dir <- "/home/johan/output/skin_pmh/fibroblast"

dir.create(fibro_out_dir, recursive = TRUE, showWarnings = FALSE)



ggsave(file.path(fibro_out_dir, "fibroblast_condition_umap.pdf"), plot = comparison_fibroblast, width = 10, height = 6)



# Visual Proof: Split by orig.ident1

comparison_fibro <- DimPlot(

  fibro_obj,

  split.by = "orig.ident1",   # changed from "Condition" to "orig.ident1"

  label = TRUE,

  label.size = 5,

  pt.size = 0.8,

  repel = TRUE,

  ncol = 2

) +

   ggtitle("Identification of PMH-Specific fibro Population")



# Save the plot

ggsave(

  file.path(fibro_out_dir, "fibro_origident1_umap.pdf"),  # updated file name

  plot = comparison_fibro,

  width = 10,

  height = 15

)





# Visual Proof: Split by orig.ident1

comparison_fibro <- DimPlot(

  fibro_obj,

  split.by = "orig.ident2",   # changed from "Condition" to "orig.ident1"

  label = TRUE,

  label.size = 5,

  pt.size = 0.8,

  repel = TRUE,

  ncol = 2

) +

  ggtitle("Identification of PMH-Specific fibro Population")



# Save the plot

ggsave(

  file.path(fibro_out_dir, "fibro_origident2_umap.pdf"),  # updated file name

  plot = comparison_fibro,

  width = 10,

  height = 15

)









# Calculate the percentage of cells in each cluster per condition

  prop_table <- prop.table(table(Idents(fibro_obj), fibro_obj$Condition), margin = 2) # calculate proportions per condition

  prop_df <- as.data.frame(prop_table)

  colnames(prop_df) <- c("Cluster", "Condition", "Proportion")



# Create the Stacked Bar Plot

stacked_comparison <- ggplot(prop_df, aes(x = Condition, y = Proportion, fill = Cluster)) +

  geom_bar(stat = "identity") +

  theme_minimal() +

  labs(title = "Figure 2: Drastic Expansion of Pathogenic Fibroblasts",

       y = "Percentage of Total Fibroblasts") +

  scale_y_continuous(labels = scales::percent)



ggsave(file.path(fibro_out_dir, "fibroblast_stacked_barplot.pdf"), plot = stacked_comparison, width = 8, height = 6)



#check original to new cluster



cluster_compare <- table(fibro_obj$global_cluster, Idents(fibro_obj))

write.csv(cluster_compare, file.path(fibro_out_dir, "fibro_cluster_comparison.csv"), row.names = TRUE)





# --- ADDITIONS FOR DIFFERENTIAL GENE EXPRESSION AND GO ANALYSIS ---



# 1. Run Differential Gene Expression (DGE)

# To compare one cluster against the other fibroblasts, you first need to identify

# which new cluster ID (after re-clustering) corresponds to the population you are interested in.

# For example, if you want to find markers for new Cluster '0':

cluster0_markers <- FindMarkers(

  object = fibro_obj,

  ident.1 = "1",

  min.pct = 0.25,

  logfc.threshold = 0.25

)



# Save the DGE results

write.csv(cluster0_markers, file.path(fibro_out_dir, "fibroblast_cluster0_vs_others.csv"))



# Optional: Run DGE comparing PMH vs Healthy within a specific cluster (e.g., Cluster 0)

# fibro_obj$cluster_condition <- paste(Idents(fibro_obj), fibro_obj$Condition, sep = "_")

# Idents(fibro_obj) <- "cluster_condition"

# pmh_vs_healthy_cluster0 <- FindMarkers(fibro_obj, ident.1 = "0_PMH", ident.2 = "0_Healthy")



# 2. Gene Ontology (GO) Analysis Preparation

# Ensure clusterProfiler is installed: BiocManager::install("clusterProfiler")

library(clusterProfiler)

library(org.Hs.eg.db)



# Filter for significantly upregulated genes in Cluster 0

sig_genes <- rownames(cluster0_markers[cluster0_markers$p_val_adj < 0.05 & cluster0_markers$avg_log2FC > 0.5, ])



# 3. Run GO Enrichment Analysis

go_results <- enrichGO(

  gene          = sig_genes,

  OrgDb         = org.Hs.eg.db,

  keyType       = 'SYMBOL',

  ont           = "BP", # Biological Process

  pAdjustMethod = "BH",

  pvalueCutoff  = 0.05,

  qvalueCutoff  = 0.05

)



# 4. Visualize GO Results

# Dotplot of the top enriched pathways

if (!is.null(go_results)) {

  go_dotplot <- dotplot(go_results, showCategory = 15) +

    ggtitle("GO Enrichment: Fibroblast Cluster 0")

 

  ggsave(file.path(fibro_out_dir, "fibroblast_cluster0_GO_dotplot.pdf"), plot = go_dotplot, width = 8, height = 6)

 

  # Save GO tabular results

  write.csv(as.data.frame(go_results), file.path(fibro_out_dir, "fibroblast_cluster0_GO_results.csv"))

} else {

  message("No significant GO terms found for the specified thresholds.")

}













#only take macrophage clusters: 6

macro_obj <- subset(pmh_obj, idents = c("6"))

table(Idents(macro_obj ))



saveRDS(macro_obj , "/home/johan/output/skin_pmh/macrophage_only.rds")



#switch back to the actual raw gene expression values so you can detect real biological differences.

DefaultAssay(macro_obj) <- "RNA"

macro_obj <- NormalizeData(macro_obj)

macro_obj <- FindVariableFeatures(macro_obj, selection.method = "vst", nfeatures = 2000)

macro_obj <- ScaleData(macro_obj)

macro_obj <- RunPCA(macro_obj)

macro_obj <- RunUMAP(macro_obj, dims = 1:20)

macro_obj <- FindNeighbors(macro_obj, dims = 1:20)

macro_obj <- FindClusters(macro_obj, resolution = 0.4)





macro_obj$Condition <- ifelse(grepl("HTY", macro_obj$orig.ident2, ignore.case = TRUE),

                              "Healthy",

                              "PMH")



# Verify it exists now

table(macro_obj$Condition)

table(macro_obj$orig.ident, macro_obj$Condition)

colnames(macro_obj@meta.data)





# Visual Proof: Healthy vs. Disease

comparison_macro <- DimPlot(macro_obj, split.by = "Condition", label = TRUE, label.size = 5, pt.size = 0.5) +

  ggtitle("Figure 1: Identification of PMH-Specific macro Population")



# Update output paths to the macro subfolder

macro_out_dir <- "/home/johan/output/skin_pmh/macro"

dir.create(macro_out_dir, recursive = TRUE, showWarnings = FALSE)



ggsave(file.path(macro_out_dir, "macro_condition_umap.pdf"), plot = comparison_macro, width = 10, height = 6)





# Visual Proof: Split by orig.ident1

comparison_macro <- DimPlot(

  macro_obj,

  split.by = "orig.ident1",   # changed from "Condition" to "orig.ident1"

  label = TRUE,

  label.size = 5,

  pt.size = 0.8,

  repel = TRUE,

  ncol = 2

) +

  ggtitle("Identification of PMH-Specific macro Population")



# Save the plot

ggsave(

  file.path(macro_out_dir, "macro_origident1_umap.pdf"),  # updated file name

  plot = comparison_macro,

  width = 10,

  height = 15

)



# Visual Proof: Split by orig.ident2

comparison_macro <- DimPlot(

  macro_obj,

  split.by = "orig.ident2",   # changed from "Condition" to "orig.ident2"

  label = TRUE,

  label.size = 5,

  pt.size = 0.8,

  repel = TRUE,

  ncol = 2

) +

  ggtitle("Identification of PMH-Specific macro Population")



# Save the plot

ggsave(

  file.path(macro_out_dir, "macro_origident2_umap.pdf"),  # updated file name

  plot = comparison_macro,

  width = 10,

  height = 15

)







# Calculate the percentage of cells in each cluster per condition

  prop_table <- prop.table(table(Idents(macro_obj), macro_obj$Condition), margin = 2) # calculate proportions per condition

  prop_df <- as.data.frame(prop_table)

  colnames(prop_df) <- c("Cluster", "Condition", "Proportion")



# Create the Stacked Bar Plot

stacked_comparison <- ggplot(prop_df, aes(x = Condition, y = Proportion, fill = Cluster)) +

  geom_bar(stat = "identity") +

  theme_minimal() +

  labs(title = "Figure 2: Drastic Expansion of Pathogenic macros",

       y = "Percentage of Total macros") +

  scale_y_continuous(labels = scales::percent)



ggsave(file.path(macro_out_dir, "macro_stacked_barplot.pdf"), plot = stacked_comparison, width = 8, height = 6)





cluster_compare <- table(macro_obj$global_cluster, Idents(macro_obj))

write.csv(cluster_compare, file.path(macro_out_dir, "macro_cluster_comparison.csv"), row.names = TRUE)













#only take macrophage + fibroblast clusters: 1,5,6,7,10

macro_fibro_obj <- subset(pmh_obj, idents = c("1","5","6","7","10"))

table(Idents(macro_fibro_obj))



saveRDS(macro_fibro_obj , "/home/johan/output/skin_pmh/macrophage_fibroblast.rds")



#switch back to the actual raw gene expression values so you can detect real biological differences.

DefaultAssay(macro_fibro_obj) <- "RNA"

macro_fibro_obj <- NormalizeData(macro_fibro_obj)

macro_fibro_obj <- FindVariableFeatures(macro_fibro_obj, selection.method = "vst", nfeatures = 2000)

macro_fibro_obj <- ScaleData(macro_fibro_obj)

macro_fibro_obj <- RunPCA(macro_fibro_obj)

macro_fibro_obj <- RunUMAP(macro_fibro_obj, dims = 1:20)

macro_fibro_obj <- FindNeighbors(macro_fibro_obj, dims = 1:20)

macro_fibro_obj <- FindClusters(macro_fibro_obj, resolution = 0.4)





macro_fibro_obj$Condition <- ifelse(grepl("HTY", macro_fibro_obj$orig.ident2, ignore.case = TRUE),

                              "Healthy",

                              "PMH")



# Verify it exists now

table(macro_fibro_obj$Condition)

table(macro_fibro_obj$orig.ident, macro_fibro_obj$Condition)

colnames(macro_fibro_obj@meta.data)







# Visual Proof: Healthy vs. Disease

comparison_macro_fibro <- DimPlot(macro_fibro_obj, split.by = "Condition", label = TRUE, label.size = 5, pt.size = 0.5) +

  ggtitle("Figure 1: Identification of PMH-Specific macro_fibro Population")



# Update output paths to the macro_fibro subfolder

macro_fibro_out_dir <- "/home/johan/output/skin_pmh/macro_fibro"

dir.create(macro_fibro_out_dir, recursive = TRUE, showWarnings = FALSE)



ggsave(file.path(macro_fibro_out_dir, "macro_fibro_condition_umap.pdf"), plot = comparison_macro_fibro, width = 10, height = 6)







# Visual Proof: Split by orig.ident1

comparison_macro_fibro <- DimPlot(

  macro_fibro_obj,

  split.by = "orig.ident1",   # changed from "Condition" to "orig.ident1"

  label = TRUE,

  label.size = 5,

  pt.size = 0.8,

  repel = TRUE,

  ncol = 2

) +

  ggtitle("Identification of PMH-Specific macrophage + fibroblast Population")



# Save the plot

ggsave(

  file.path(macro_fibro_out_dir, "macro_fibro_origident1_umap.pdf"),  # updated file name

  plot = comparison_macro_fibro,

  width = 10,

  height = 15

)









# Visual Proof: Split by orig.ident1

comparison_macro_fibro <- DimPlot(

  macro_fibro_obj,

  split.by = "orig.ident2",   # changed from "Condition" to "orig.ident1"

  label = TRUE,

  label.size = 5,

  pt.size = 0.8,

  repel = TRUE,

  ncol = 2

) +

  ggtitle("Identification of PMH-Specific macrophage + fibroblast Population")



# Save the plot

ggsave(

  file.path(macro_fibro_out_dir, "macro_fibro_origident2_umap.pdf"),  # updated file name

  plot = comparison_macro_fibro,

  width = 10,

  height = 15

)









# Calculate the percentage of cells in each cluster per condition

  prop_table <- prop.table(table(Idents(macro_fibro_obj), macro_fibro_obj$Condition), margin = 2) # calculate proportions per condition

  prop_df <- as.data.frame(prop_table)

  colnames(prop_df) <- c("Cluster", "Condition", "Proportion")



# Create the Stacked Bar Plot

stacked_comparison <- ggplot(prop_df, aes(x = Condition, y = Proportion, fill = Cluster)) +

  geom_bar(stat = "identity") +

  theme_minimal() +

  labs(title = "Figure 2: Drastic Expansion of Pathogenic macro_fibros",

       y = "Percentage of Total macro_fibros") +

  scale_y_continuous(labels = scales::percent)



ggsave(file.path(macro_fibro_out_dir, "macro_fibro_stacked_barplot.pdf"), plot = stacked_comparison, width = 8, height = 6)



cluster_compare <- table(macro_fibro_obj$global_cluster, Idents(macro_fibro_obj))

write.csv(cluster_compare, file.path(macro_fibro_out_dir, "macro_fibro_cluster_comparison.csv"), row.names = TRUE)









# -------------------- Annotation for fibroblast subset --------------------

# Summary: run SingleR (HPCA & BlueprintEncode) + scCATCH on fibro_obj,

# add per-cell SingleR labels and per-cluster majority annotations,

# create DotPlot of fibroblast markers, annotated UMAPs, save CSVs and annotated RDS.



suppressPackageStartupMessages({

  library(SingleR)

  library(celldex)

  library(scCATCH)

  library(ExperimentHub)

})



safe_save_pdf <- function(plot_obj, filename, w = 8, h = 6) {

  filepath <- file.path(fibro_out_dir, filename)

  tryCatch({

    pdf(filepath, width = w, height = h)

    on.exit(dev.off(), add = TRUE)

    print(plot_obj)

  }, error = function(e) {

    message("Failed to save plot: ", filepath, " : ", e$message)

    if (length(dev.list()) > 0) dev.off()

  })

}



message("Starting fibroblast annotation...")



# Ensure RNA assay and joined layers available

DefaultAssay(fibro_obj) <- "RNA"

joined_fib <- tryCatch({

  JoinLayers(fibro_obj)

}, error = function(e) {

  message("JoinLayers failed, proceeding with fibro_obj directly: ", e$message)

  fibro_obj

})



# Extract expression matrix (normalized data layer if available)

expr <- tryCatch({

  as.matrix(GetAssayData(joined_fib, assay = "RNA", layer = "data"))

}, error = function(e) {

  message("Failed to get 'data' layer, falling back to default GetAssayData(): ", e$message)

  as.matrix(GetAssayData(joined_fib))

})



# Ensure cell barcodes are consistent and named everywhere

cells <- colnames(expr)



singleR_results <- list()



# Helper to safely extract labels from SingleR result and align to cell barcodes

get_labels_from_singleR <- function(pred, cells) {

  labs <- NULL



  # try the common slots/fields

  if (!is.null(pred$pruned.labels)) {

    labs <- pred$pruned.labels

  } else if (!is.null(pred$labels)) {

    labs <- pred$labels

  } else if (isS4(pred) && !is.null(pred@listData[["pruned.labels"]])) {

    labs <- pred@listData[["pruned.labels"]]

  } else if (isS4(pred) && !is.null(pred@listData[["labels"]])) {

    labs <- pred@listData[["labels"]]

  }



  # Force character vector

  labs <- as.character(labs)



  # Prepare output vector named by cells (default NA)

  out <- rep(NA_character_, length(cells))

  names(out) <- cells



  if (length(labs) == 0) {

    warning("SingleR returned no labels; returning NA vector.")

    return(out)

  }



  # If labels are named and overlap with cells, align by names

  if (!is.null(names(labs))) {

    common <- intersect(names(labs), cells)

    if (length(common) > 0) {

      out[common] <- labs[common]

      return(out)

    }

  }



  # If unlabeled but same length, assume order matches cells

  if (length(labs) == length(cells)) {

    names(labs) <- cells

    out[] <- labs

    return(out)

  }



  # Fallback: partially assign first N labels and warn

  n_assigned <- min(length(labs), length(cells))

  out[seq_len(n_assigned)] <- labs[seq_len(n_assigned)]

  warning("SingleR labels length (", length(labs), ") does not match number of cells (", length(cells), "). Assigned first ",

          n_assigned, " labels to cells; remaining set to NA.")

  return(out)

}



# Run SingleR with defensive handling

try({

  message("Loading HPCA reference...")

  hpca_ref <- tryCatch(HumanPrimaryCellAtlasData(), error = function(e) {

    message("HumanPrimaryCellAtlasData() failed: ", e$message, " — retrying with ExperimentHub fallback")

    ExperimentHub::ExperimentHub()

    celldex::HumanPrimaryCellAtlasData()

  })



  message("Running SingleR (HPCA)...")

  pred_hpca <- SingleR(test = expr, ref = hpca_ref, assay.type.test = 1, labels = hpca_ref$label.main)

  labels_hpca <- get_labels_from_singleR(pred_hpca, cells)

  fibro_obj$SingleR_HPCA_label <- labels_hpca



  message("Loading BlueprintEncode reference...")

  bpe_ref <- tryCatch(BlueprintEncodeData(), error = function(e) {

    message("BlueprintEncodeData() failed: ", e$message, " — retrying with ExperimentHub fallback")

    ExperimentHub::ExperimentHub()

    celldex::BlueprintEncodeData()

  })



  message("Running SingleR (BlueprintEncode)...")

  pred_bpe <- SingleR(test = expr, ref = bpe_ref, assay.type.test = 1, labels = bpe_ref$label.main)

  labels_bpe <- get_labels_from_singleR(pred_bpe, cells)

  fibro_obj$SingleR_BPE_label <- labels_bpe



  singleR_results$hpca <- pred_hpca

  singleR_results$bpe <- pred_bpe



  # Derive per-cluster majority label mapping (unchanged logic)

  meta_df <- fibro_obj@meta.data

  meta_df$cluster <- as.character(Idents(fibro_obj))



  get_majority <- function(label_col) {

    # Protect against missing or all-NA label columns

    if (!(label_col %in% colnames(meta_df))) {

      return(setNames(rep(NA_character_, length(unique(meta_df$cluster))), unique(meta_df$cluster)))

    }

    agg <- as.data.frame(table(meta_df$cluster, meta_df[[label_col]]), stringsAsFactors = FALSE)

    colnames(agg) <- c("cluster","label","n")

    agg <- agg[order(agg$cluster, -agg$n), ]

    top <- aggregate(n ~ cluster, data = agg, FUN = max)

    merged <- merge(top, agg, by = c("cluster","n"))

    mapping <- setNames(as.character(merged$label), merged$cluster)

    return(mapping)

  }



  hpca_map <- get_majority("SingleR_HPCA_label")

  bpe_map  <- get_majority("SingleR_BPE_label")



  # Assign cluster-level annotations (will produce NA for clusters without mapping)

  fibro_obj$SingleR_HPCA_cluster <- hpca_map[as.character(Idents(fibro_obj))]

  fibro_obj$SingleR_BPE_cluster  <- bpe_map[as.character(Idents(fibro_obj))]



  # Save per-cell labels and per-cluster summaries to fibro_out_dir

  write.csv(data.frame(Cell = cells, Label = labels_hpca),

            file.path(fibro_out_dir, "singleR_hpca_per_cell.csv"), row.names = FALSE)

  write.csv(data.frame(Cell = cells, Label = labels_bpe),

            file.path(fibro_out_dir, "singleR_bpe_per_cell.csv"), row.names = FALSE)



  write.csv(table(fibro_obj$SingleR_HPCA_label, Idents(fibro_obj)),

            file.path(fibro_out_dir, "singleR_hpca_by_cluster.csv"))

  write.csv(table(fibro_obj$SingleR_BPE_label, Idents(fibro_obj)),

            file.path(fibro_out_dir, "singleR_bpe_by_cluster.csv"))



  # UMAPs grouped by SingleR cluster majority

  if ("SingleR_HPCA_cluster" %in% colnames(fibro_obj@meta.data)) {

    p_hpca <- DimPlot(fibro_obj, group.by = "SingleR_HPCA_cluster", label = TRUE, repel = TRUE) +

      ggtitle("SingleR HPCA - cluster majority")

    safe_save_pdf(p_hpca, "fibro_umap_singleR_hpca_cluster.pdf", w = 8, h = 6)

  }

  if ("SingleR_BPE_cluster" %in% colnames(fibro_obj@meta.data)) {

    p_bpe <- DimPlot(fibro_obj, group.by = "SingleR_BPE_cluster", label = TRUE, repel = TRUE) +

      ggtitle("SingleR BPE - cluster majority")

    safe_save_pdf(p_bpe, "fibro_umap_singleR_bpe_cluster.pdf", w = 8, h = 6)

  }

}, silent = FALSE)



# scCATCH annotation (marker-based)

try({

  message("Running scCATCH on fibroblast subset...")



  # Ensure Idents are aligned and named by cell barcodes

  labels_vec <- as.character(Idents(fibro_obj))

  names(labels_vec) <- Cells(fibro_obj)



  # Intersect and reorder by common cells to guarantee overlap

  common_cells <- intersect(cells, names(labels_vec))

  if (length(common_cells) == 0) {

    stop("No cell overlap between new meta data and Seurat object for scCATCH. Check that 'expr' columns and Idents(fibro_obj) share barcodes.")

  }

  # Reorder both to the same order

  labels_vec <- labels_vec[common_cells]

  sc_data <- expr[, common_cells, drop = FALSE]



  # rev_gene may fail without geneinfo; fall back if needed

  sc_data_rev <- tryCatch({

    rev_gene(data = sc_data, data_type = "data", species = "Human", geneinfo = geneinfo)

  }, error = function(e) {

    message("rev_gene failed or geneinfo missing: ", e$message, " — proceeding without rev_gene")

    sc_data

  })



  # Create scCATCH object: pass a vector ordered to match sc_data columns

  sc_obj <- createscCATCH(data = sc_data_rev, cluster = as.character(labels_vec))



  tissue_list <- c("Skin", "Dermis", "Hair follicle", "Subcutaneous adipose tissue")

  sc_obj <- findmarkergene(object = sc_obj, species = "Human", marker = cellmatch, tissue = tissue_list, use_method = "1")

  sc_obj <- findcelltype(object = sc_obj)



  write.csv(sc_obj@celltype, file = file.path(fibro_out_dir, "fibro_scCATCH_celltype.csv"), row.names = FALSE)

  # Summary TSV

  sc_df_try <- try(as.data.frame(sc_obj@celltype, stringsAsFactors = FALSE), silent = TRUE)

  if (!inherits(sc_df_try, "try-error") && ncol(sc_df_try) >= 2) {

    write.table(sc_df_try, file = file.path(fibro_out_dir, "fibro_scCATCH_summary.tsv"),

                sep = "\t", quote = FALSE, row.names = FALSE)

  }

}, silent = FALSE)



# DotPlot of canonical fibroblast markers

try({

  fibro_markers <- c("PDGFRA","DCN","LUM","POSTN","COL1A1","COL3A1","COL6A3","ACTA2")

  common_genes <- intersect(fibro_markers, rownames(fibro_obj))

  if(length(common_genes) > 0) {

    p_dot <- DotPlot(fibro_obj, features = common_genes) + RotatedAxis() + ggtitle("Fibroblast marker expression")

    safe_save_pdf(p_dot, "fibro_markers_dotplot.pdf", w = 10, h = 6)

  } else {

    message("None of the canonical fibroblast markers found in object; skipping DotPlot.")

  }

}, silent = TRUE)



# UMAPs colored by SingleR cluster annotation and scCATCH (if available)

try({

  if ("SingleR_HPCA_cluster" %in% colnames(fibro_obj@meta.data)) {

    p1 <- DimPlot(fibro_obj, group.by = "SingleR_HPCA_cluster", label = TRUE, repel = TRUE) + ggtitle("SingleR HPCA - cluster majority")

    safe_save_pdf(p1, file.path(fibro_out_dir, "fibro_umap_singleR_hpca_cluster.pdf"), w = 8, h = 6)

  }

  if ("SingleR_BPE_cluster" %in% colnames(fibro_obj@meta.data)) {

    p2 <- DimPlot(fibro_obj, group.by = "SingleR_BPE_cluster", label = TRUE, repel = TRUE) + ggtitle("SingleR BPE - cluster majority")

    safe_save_pdf(p2, file.path(fibro_out_dir, "fibro_umap_singleR_bpe_cluster.pdf"), w = 8, h = 6)

  }

}, silent = TRUE)



# Save annotated fibroblast object

annotated_rds <- file.path(fibro_out_dir, "fibroblasts_annotated.rds")

saveRDS(fibro_obj, annotated_rds)

message("Annotated fibroblast object saved to: ", annotated_rds)

