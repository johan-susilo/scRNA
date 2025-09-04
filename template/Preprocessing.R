#!/usr/bin/env Rscript
# usage : Rscript ./Preprocessing.R -s data1_folder,data2_folder,... -ident1 data1_nature,data2_nature,... -identi2 data1_ID,data2_ID,...
# ex    : Rscript ./Preprocessing.R -s sample1,sample2,sample3,... -ident1 1_Healthy_skin,2_Acute_skin,3_Chronic_skin -ident2 HC131,HC244,HC258...
#above commandline is just example, modified it due to your usage!


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

#create dir for output
dir.create("output")

args <- commandArgs(trailingOnly = TRUE)

# Check if any arguments are provided
if (length(args) < 1) {
  stop("Please provide sample names and identifiers as arguments.")
}

# Initialize variables and parse arguments
sample_names <- NULL
ident1 <- NULL
ident2 <- NULL


if (length(args) > 0) {
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "-s") {
      # Check if the next argument exists and is not another flag
      if (i + 1 <= length(args) && substr(args[i + 1], 1, 1) != "-") {
        sample_names <- strsplit(args[i + 1], ",")[[1]]
        options(sample_names = sample_names)
        i <- i + 1
      } else {
        stop("Please provide a comma-separated list of sample names after -s")
      }
    } else if (args[i] == "-ident1") {
      # Check if the next argument exists and is not another flag
      if (i + 1 <= length(args) && substr(args[i + 1], 1, 1) != "-") {
        ident1 <- strsplit(args[i + 1], ",")[[1]]
        options(ident1 = ident1)
        i <- i + 1
      } else {
        stop("Please provide a comma-separated list for -ident1")
      }
    } else if (args[i] == "-ident2") {
      # Check if the next argument exists and is not another flag
      if (i + 1 <= length(args) && substr(args[i + 1], 1, 1) != "-") {
        ident2 <- strsplit(args[i + 1], ",")[[1]]
        options(ident2 = ident2)
        i <- i + 1
      } else {
        stop("Please provide a comma-separated list for -ident2")
      }
    }
    i <- i + 1
  }
}

# Use getOption to retrieve the options set via command-line
sample_names <- getOption("sample_names")
ident1 <- getOption("ident1")
ident2 <- getOption("ident2")
sample_list <- list()

# Check if sample names are set, then process each
if (!is.null(sample_names)) {
  # Loop over each sample name
  for (i in seq_along(sample_names)) {
      sample_name <- sample_names[i]
      sample_ident1 <- ident1[i]
      sample_ident2 <- ident2[i]

      
      #save sample name as variable
      sample_id <- sample_name
      #Define the path for reading the data
      data_dir <- paste0("./",sample_name)
      print(data_dir)
      print(paste("Processing sample:", sample_name))  
      

      #Read the 10X data
      sample_data <- Read10X(data.dir = data_dir)
      
      
      #Create a Seurat object
      seur_obj <- CreateSeuratObject(counts = sample_data, project = sample_name, min.cells = 3, min.features = 10)
      print(seur_obj)

      ##Remove doublet cells by DoubletFinder package
      #Run general flow of scRNA-seq by Seurat package
      seur_obj <- NormalizeData(object = seur_obj)
      seur_obj <- FindVariableFeatures(object = seur_obj)
      seur_obj <- ScaleData(object = seur_obj)
      seur_obj <- RunPCA(object = seur_obj)
      #Plot elbow plot
      pdf (paste0("./output/",sample_id,"_elbow1.pdf"), width = 15, height = 15)
      print(ElbowPlot(seur_obj))
      dev.off()

      #print a notification of done plotting!
      system('echo "-----------------------------------ElbowPlot Plotting completed!-----------------------------------"')

      seur_obj <- FindNeighbors(object = seur_obj, dims = 1:50)
      seur_obj <- FindClusters(object = seur_obj)
      seur_obj <- RunUMAP(object = seur_obj, dims = 1:30)
      #Plot UMAP plot
      pdf (paste0("./output/",sample_id,"_umap1.pdf"), width = 15, height = 15)
      print(DimPlot(seur_obj, reduction = "umap", label = T))
      dev.off()

      #print a notification of done plotting!
      system('echo "-----------------------------------UMAP Plot Plotting completed!-----------------------------------"')

      #pK Identification (no ground-truth)
      sweep.res.list <- paramSweep(seur_obj, PCs = 1:20, sct = FALSE)
      sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
      bcmvn <- find.pK(sweep.stats)
      pdf (paste0("./output/",sample_id,"_pkplot.pdf"), width = 15, height = 15)
      print(ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line())
      dev.off()

      #print a notification of done plotting!
      system('echo "-----------------------------------pk Plot Plotting completed!-----------------------------------"')

      pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
      pK <- as.numeric(as.character(pK[[1]]))
      print(head(pK))

      #Homotypic Doublet Proportion Estimate
      annotations <- seur_obj@meta.data$seurat_clusters
      homotypic.prop <- modelHomotypic(annotations)
      nExp_poi <- round(0.08*nrow(seur_obj@meta.data)) ## 8% doublet formation rate from 10X Genomics form if 10,000 cells are obtained.
      nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
      print(head(nExp_poi.adj))

      #Identify doublet cells
      seur_obj <- doubletFinder(seur_obj, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
      print(head(seur_obj))
      DF_classification <- colnames(seur_obj@meta.data)[7]

      #plot unfiltered doublet plot for examination!
      pdf (paste0("./output/",sample_id,"_unfiltered_doublet.pdf"), width = 15, height = 15)
      print(DimPlot(seur_obj, reduction = 'umap', group.by = DF_classification))
      dev.off()

      #print a notification of done plotting!
      system('echo "-----------------------------------unfiltered doublet Plot Plotting completed!-----------------------------------"')

      #check doublet and singlet info before filtering!
      print(table(seur_obj@meta.data[7]))

      #Remove doublet cells
      #Get row indices for Singlet
      singlet_indices <- which(seur_obj@meta.data[[DF_classification]] == "Singlet")
      #Filtered the Seurat object retained only Singlet sample
      seur_obj <- seur_obj[,singlet_indices]

      #plot filtered doublet plot for examination!
      pdf (paste0("./output/",sample_id,"_filtered_doublet.pdf"), width = 15, height = 15)
      print(DimPlot(seur_obj, reduction = "umap", group.by= DF_classification))
      dev.off()
      
      #print a notification of done plotting!
      system('echo "-----------------------------------filtered doublet Plot Plotting completed!-----------------------------------"')

      #check sample info after filtering doublet!
      print(seur_obj)
      print(table(seur_obj@meta.data[7]))

      ##Quality control
      #The $ operator can add columns to object metadata.
      seur_obj$orig.ident1 <- sample_ident1
      seur_obj$orig.ident2 <- sample_ident2
      sample_ident <- paste0(sample_ident1,"_",sample_ident2)
      print(head(seur_obj))
      print(sample_ident)
      
      #The "[[" operator can add columns to object metadata. This is a great place to stash QC stats
      seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^MT-")
      #Show QC metrics for the first 5 cells
      head(seur_obj@meta.data, 5)
      #Visualize QC metrics as a violin plot
      pdf (paste0("./output/",sample_id,"_violin.pdf"), width = 15, height = 15)
      print(VlnPlot(seur_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
      dev.off()

      #print a notification of done plotting!
      system('echo "-----------------------------------Violin Plot Plotting completed!-----------------------------------"')

      #Filter out for Feature and mitochondria pct.
      seur_obj <- subset(seur_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)
      print(seur_obj)

      sample_list[[sample_ident]] <- seur_obj
      }     
      
      ##scRNA-seq sample integration
      #Set merge file
      # Load the Seurat objects dynamically using `get()`
      #seurat_objects <- lapply(sample_list, get) 
      #print(sample_list)
      #print(seurat_objects)
      # Merge all Seurat objects
      #TN.C <- merge(seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = sample_list, project = "groups")
      #View the first few rows of the metadata
      #print(TN.C)
      # Check the number of cells from each original identity
      #print(table(TN.C$orig.ident))
      } else {
  print("No sample names set. Please use -s to specify sample names.")
}

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


