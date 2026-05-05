suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(ggplot2)
})

# Rscript /home/johan/pipeline/scRNA/skin/cellchat_analysis.R --macro "/home/johan/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/macrophage/processed/macrophage_detailed_annotated.rds"   --fibro "/home/johan/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/fibroblast/processed/fibroblast_annotated_full.rds"   --outdir "/home/johan/johan/output/skin_pmh_harmony_sctransform2/cellchat/fibro_macro"

# ==============================================================================
# 1. SETUP COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("--fibro"), type="character", help="Path to Fibroblast RDS"),
  make_option(c("--macro"), type="character", help="Path to Macrophage RDS"),
  make_option(c("--outdir"), type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. LOAD OBJECTS
# ==============================================================================
message("Loading objects...")
fibro <- readRDS(opt$fibro)
macro <- readRDS(opt$macro)

if (!"SingleR_label" %in% colnames(fibro@meta.data)) {
  stop("ERROR: Fibroblast object is missing 'SingleR_label'. Use the file ending in '_annotated_final.rds' or '_detailed_annotated.rds'.")
}

# ==============================================================================
# 3. SEURAT V5-NATIVE DATA EXTRACTION & MERGE
# ==============================================================================
message("Extracting raw data to build a clean combined object...")

safe_extract_counts <- function(seu) {
  assay_use <- ifelse("RNA" %in% Assays(seu), "RNA", DefaultAssay(seu))
  if (inherits(seu[[assay_use]], "Assay5")) {
    try({ seu <- JoinLayers(seu) }, silent = TRUE)
  }
  mat <- tryCatch({
    GetAssayData(seu, assay = assay_use, layer = "counts")
  }, error = function(e) {
    GetAssayData(seu, assay = assay_use, slot = "counts")
  })
  return(mat)
}

counts_fibro <- safe_extract_counts(fibro)
counts_macro <- safe_extract_counts(macro)

colnames(counts_fibro) <- paste0("FB_", colnames(counts_fibro))
colnames(counts_macro) <- paste0("MAC_", colnames(counts_macro))

common_genes <- intersect(rownames(counts_fibro), rownames(counts_macro))
counts_fibro <- counts_fibro[common_genes, ]
counts_macro <- counts_macro[common_genes, ]
combined_counts <- cbind(counts_fibro, counts_macro)

meta_fibro <- fibro@meta.data
meta_macro <- macro@meta.data
rownames(meta_fibro) <- paste0("FB_", rownames(meta_fibro))
rownames(meta_macro) <- paste0("MAC_", rownames(meta_macro))

cols_to_keep <- c("orig.ident1", "orig.ident2", "Condition", "SingleR_label")
meta_fibro <- meta_fibro[, cols_to_keep, drop = FALSE]
meta_macro <- meta_macro[, cols_to_keep, drop = FALSE]
combined_meta <- rbind(meta_fibro, meta_macro)

message("Building new combined Seurat object...")
combined <- CreateSeuratObject(counts = combined_counts, meta.data = combined_meta)
combined <- NormalizeData(combined, verbose = FALSE)

Idents(combined) <- "SingleR_label"

rm(fibro, macro, counts_fibro, counts_macro)
gc()

# ==============================================================================
# 4. CELLCHAT WORKFLOW (BYPASSING INTERNAL BUGS)
# ==============================================================================
message("Initializing CellChat (bypassing Seurat v5 internal conflicts)...")

# Manually extract the normalized data matrix using v5 syntax
data.input <- tryCatch({
  GetAssayData(combined, assay = "RNA", layer = "data")
}, error = function(e) {
  GetAssayData(combined, assay = "RNA", slot = "data")
})
meta.data <- combined@meta.data

# Create CellChat using the raw matrix and metadata directly!
cellchat <- createCellChat(object = data.input, meta = meta.data, group.by = "SingleR_label")

cellchat@DB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

message("Computing probabilities (this takes time)...")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# ==============================================================================
# 5. VISUALIZATIONS
# ==============================================================================
message("Saving plots to: ", opt$outdir)

fibro_groups <- as.character(unique(combined_meta$SingleR_label[grepl("^F", combined_meta$SingleR_label)]))
macro_groups <- as.character(unique(combined_meta$SingleR_label[grepl("^M", combined_meta$SingleR_label)]))

pdf(file.path(opt$outdir, "Fibro_to_Macro_Crosstalk.pdf"), width = 12, height = 8)
print(netVisual_bubble(cellchat, sources.use = fibro_groups, targets.use = macro_groups) + 
      ggtitle("Signals: Fibroblasts -> Macrophages"))
dev.off()

pdf(file.path(opt$outdir, "Macro_to_Fibro_Crosstalk.pdf"), width = 12, height = 8)
print(netVisual_bubble(cellchat, sources.use = macro_groups, targets.use = fibro_groups) + 
      ggtitle("Signals: Macrophages -> Fibroblasts"))
dev.off()

pdf(file.path(opt$outdir, "Interaction_Network_Circle.pdf"), width = 10, height = 10)
netVisual_circle(cellchat@net$count, weight.scale = T, label.edge = F, title.name = "Number of interactions")
dev.off()

saveRDS(cellchat, file.path(opt$outdir, "fibro_macro_cellchat.rds"))
message("Analysis Complete!")