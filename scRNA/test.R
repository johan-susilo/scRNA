#setwd("./pipeline/scRNA")

#Load required library
library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggsci)
library(Matrix)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(clusterProfiler)
library(gplots)
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)
library(tidyr)

output_base <- "/mnt/80T/johan/output/full_liver_2021"
TN_combined <- readRDS(file.path(output_base, "TN.combined_dim30.rds"))

View(TN_combined@meta.data)

# Print the summary of the Seurat object
print(TN_combined)

# Specifically check the layers available within the "RNA" assay
print(Layers(TN_combined, assay = "RNA"))

# 1. Check if "EPCAM" exists as a gene name in the RNA assay
# This should print TRUE
print("EPCAM" %in% rownames(GetAssay(TN_combined, assay = "RNA")))

# 2. Pull the raw count data for EPCAM
epcam_raw_counts <- FetchData(TN_combined, vars = "EPCAM", assay = "RNA", layer = "counts")

# 3. Get a summary of the raw counts
# This will show the minimum, maximum, and average counts
print(summary(epcam_raw_counts$EPCAM))

# 4. Count how many cells have raw counts > 0
# This is the most critical check. It will tell us the exact number of positive cells.
print(sum(epcam_raw_counts$EPCAM > 0))





# Initialize variables
tsv_file <- "./input.tsv"


# Read TSV file
if (!is.null(tsv_file)) {
  samples_df <- read.delim(tsv_file, header = FALSE, stringsAsFactors = FALSE, col.names = c("sample_names", "ident1", "ident2"))
} else {
  stop("TSV file not provided. Use -f to specify.")
}

print(sample_names)
# Process each row in the TSV
sample_list <- list()

output_dir <- "./output/processed/"
rds_files <- list.files(path = output_dir, pattern = "_processed\\.rds$", full.names = TRUE)


# Loop through each .rds file
for (rds_file in rds_files) {
  # Extract the sample name from the filename
  # Assuming filenames are like "SAMPLE_ID_seurat_object.rds"
  sample_name <- gsub("_processed\\.rds$", "", basename(rds_file))

  # Create a variable name by converting the sample name to lowercase
  variable_name <- tolower(sample_name)
  print(variable_name)	
  # Read the .rds file (Seurat object)
  seur_obj <- readRDS(rds_file)

  # Assign the Seurat object to a variable with the created name in the global environment
  assign(variable_name, seur_obj, envir = .GlobalEnv)

  # Optional: Print a message to confirm
  cat("Loaded Seurat object from:", rds_file, "and assigned to variable:", variable_name, "\n")
}

View(col18_crc@meta.data)

sample_names <- samples_df$sample_names
n_samples <- nrow(samples_df)
print(sample_names)

  # Loop over each sample name
 for (i in seq_len(nrow(samples_df))) {
      sample_name <- samples_df$sample_names[i]
      sample_ident1 <- samples_df$ident1[i]
      sample_ident2 <- samples_df$ident2[i]
      
      print(sample_name)
      #save sample name as variable
      sample_id <- sample_name
      #Define the path for reading the data
      data_dir <- paste0("./",sample_name)
      print(paste("Processing sample:", sample_name))  
      
      dir.create("./output/processed/")
      output_rds <- paste0("./output/processed/", sample_id, "_processed.rds")

      ##Quality control
      # Extract orig.ident1 and orig.ident2 from the Seurat object's metadata
      sample_ident1
      sample_ident2
      sample_ident <- paste0(sample_ident1,"_",sample_ident2)
      sample_list[[sample_ident]] <- seur_obj
 
 }






x = sample_list[[1]]
y = sample_list[-1]


##scRNA-seq sample integration
#Set merge file
# Load the Seurat objects dynamically using `get()`
print(sample_names)
#seurat_objects <- lapply(sample_list, get)
print(sapply(sample_list, class))
#print (class(sample_list))
# Merge all Seurat objects
TN.C <- do.call(merge, args = list(x = sample_list[[1]], y = sample_list[-1], add.cell.ids = sample_names, project = "groups"))
#View the first few rows of the metadata
print(head(TN.C[[]]))
# Check the number of cells from each original identity
print(table(TN.C$orig.ident))

#Set splitobject file
TN.list<- SplitObject(TN.C, split.by = "orig.ident")
TN.list


# cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# NormalizeData,FindVariableFeatures, ScaleData, CellCycleScoring
TN.list <- lapply(X = TN.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = rownames(x), verbose = TRUE)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
})

#select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = TN.list)
#Perform integration
#reference.list <- TN.list[c("HC131","HC244","UA258","AC259","CH260")]
TN.anchors <- FindIntegrationAnchors(object.list = TN.list, dims = 1:30)
saveRDS(TN.anchors, file = "./output/TN.anchorsdim30.rds")

#this command creates an 'integrated' data assay
TN.combined <- IntegrateData(anchorset = TN.anchors, dims = 1:30)
#Perform an integrated analysis
#specify that we will perform downstream analysis on the corrected data note that the
#original unmodified data still resides in the 'RNA' assay
DefaultAssay(TN.combined) <- "integrated"
#Run the standard workflow for visualization and clustering
#TN.combined <- ScaleData(TN.combined, verbose = TRUE)
TN.combined <- ScaleData(TN.combined, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), features = rownames(TN.combined), verbose = TRUE)
TN.combined <- RunPCA(TN.combined, npcs = 50, verbose = TRUE)
VizDimLoadings(TN.combined, dims = 1:2)
pdf ("./output/TNcombined_elbow.pdf", width = 15, height = 15)
print(ElbowPlot(TN.combined, ndims= 50))
dev.off()
# Cluster the cells (1:30)
TN.combined <- FindNeighbors(TN.combined, reduction = "pca", dims = 1:30)
TN.combined <- FindClusters(TN.combined, resolution = 0.4)
print(table(TN.combined@active.ident))
print(table(Idents(TN.combined), TN.combined$orig.ident1))


#CellNumber
CellNumber <- table(Idents(TN.combined), TN.combined$orig.ident1)
write.csv(CellNumber , file = "./output/CellNumber.csv")
#record cluster amount for later heatmap plotting
cluster_count <- nrow(CellNumber)


#TN.combined <- RunUMAP(TN.combined, reduction="pca", dims = 1:30)
#DimPlot(TN.combined,reduction = "umap",group.by = "orig.ident")
##UMAP (1:30)
TN.combined <- RunUMAP(TN.combined, reduction = "pca", dims = 1:30)
#output these plots
pdf ("./output/TNcombined_umap_labelF.pdf", width = 15, height = 15)
print(DimPlot(TN.combined, reduction = "umap",label = FALSE, pt.size = 0.8))
dev.off()
pdf ("./output/TNcombined_umap_labelT.pdf", width = 15, height = 15)
print(DimPlot(TN.combined, reduction = "umap",label = TRUE, pt.size = 0.8))
dev.off()
pdf ("./output/TNcombined_umap_groupbyorigident.pdf", width = 15, height = 15)
print(DimPlot(TN.combined, group.by = "orig.ident", pt.size = 0.8))
dev.off()
pdf ("./output/TNcombined_umap_labelT_splitorigident1.pdf", width = 15, height = 15)
print(DimPlot(TN.combined, reduction = "umap",label = TRUE, split.by = "orig.ident1", pt.size = 0.8, ncol = 2))
dev.off()

Heatmapall <- subset(TN.combined, idents = 0:(cluster_count -1))
Heatmapall.markers <- FindAllMarkers(Heatmapall, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
head(Heatmapall.markers)
write.csv(Heatmapall.markers , file = "./output/Findallmarkers.csv", col.names = TRUE)
top10 <- Heatmapall.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf ("./output/heatmap_top10.pdf", width = 25, height = 25)
print(DoHeatmap(Heatmapall, features = top10$gene))
dev.off()

##Save the merged Seurat object
saveRDS(TN.combined, file = "./output/TN.combined_dim30.rds")

#Cellproportion
Cellproportion <- table(Idents(TN.combined), TN.combined$orig.ident1)
Cellproportion <- round(sweep(Cellproportion,MARGIN=2, STATS=colSums(Cellproportion), FUN = "/")*100,2)
write.csv(Cellproportion, file = "./output/Cellproportion.csv", row.names = T)
Cellproportion <- as.data.frame(Cellproportion)

nb.cols <- nrow(CellNumber)
mycolors <- colorRampPalette(brewer.pal(nb.cols, "Paired"))(nb.cols)

proportion_plot <- ggplot(Cellproportion, aes(x = Var2, y = Freq, fill = Var1)) + theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.6) + xlab("Sample") + ylab("Proportion") +
  theme(legend.title = element_blank())+ scale_fill_manual(values = mycolors)

pdf ("./output/proportion_plot.pdf", width = 15, height = 15)
print(proportion_plot)
dev.off()


