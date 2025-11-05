#!/usr/bin/env Rscript

.libPaths("/home/johan/R/x86_64-pc-linux-gnu-library/4.2")

library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)
library(tibble)
library(vsn)
library(pheatmap)
library(DESeq2)

counts <- read.table("/home/johan/johan/output/chicken/combined.HTseq_report", header = TRUE, sep = "\t", row.names = 1)
metadata <- read.csv("/home/johan/johan/output/chicken/metadata.csv")


# Attempt to find a matching column for sample IDs
sample_col <- purrr::detect(c("sample", "sample_id", "sampleID", "sampleName", "run", "library", "Library", "X"), 
                            ~ . %in% colnames(metadata) && any(as.character(metadata[[.]]) %in% colnames(counts)))

# Assign rownames based on found column or use existing rownames
rownames(metadata) <- as.character(metadata[[sample_col]])


# Subset metadata by tissue and then subset/reorder counts to match
ovary <- metadata %>% filter(tissue == "ovary")
uterus <- metadata %>% filter(tissue == "uterus")

# Ensure we have ovary samples and that they exist in the count matrix
ovary_samples <- intersect(colnames(counts), rownames(ovary))
uterus_samples <- intersect(colnames(counts), rownames(uterus))

# Subset and reorder counts to match ovary_samples
counts_ovary <- counts[, ovary_samples, drop = FALSE]
ovary <- ovary[ovary_samples, , drop = FALSE]  # reorder metadata rows to match counts columns

counts_uterus <- counts[, uterus_samples, drop = FALSE]
uterus <- uterus[uterus_samples, , drop = FALSE]  # reorder metadata rows

ovary$egg_production <- factor(ovary$egg_production)
ovary$egg_production <- relevel(ovary$egg_production, ref = "low")

uterus$egg_production <- factor(uterus$egg_production)
uterus$egg_production <- relevel(uterus$egg_production, ref = "low")

head(counts_uterus)

dds <- DESeqDataSetFromMatrix(
  countData = counts_ovary, 
  colData = ovary, 
  design = ~egg_production
)

dds <- DESeq(dds)
ddsc <- counts(dds, normalized = TRUE)
res_egg <- results(dds, contrast = c("egg_production", "high", "low"))

resOrdered <- res_egg[order(res_egg$padj), ]


write.csv(ddsc, file = "/home/johan/johan/output/chicken/DGE/DEseq2_normalized.csv", row.names = TRUE)
write.csv(resOrdered, file = "/home/johan/johan/output/chicken/DGE/result_egg.csv", row.names = TRUE)

vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("tissue","egg_production"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


pcaPlot <- ggplot(pcaData, aes(PC1, PC2, color = tissue, shape = egg_production)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), max.overlaps = 50, size = 3) +
  labs(x = paste0("PC1: ", percentVar[1], "%"), y = paste0("PC2: ", percentVar[2], "%"),
       title = "PCA of samples (VST)") +
  theme_bw()

ggsave("/home/johan/johan/output/chicken/DGE/PCA.png", width = 10, height = 9, plot = pcaPlot)


## volcano plot

# 1. Convert the results object to a data frame for plotting
res_df <- as.data.frame(res_egg) %>%
  rownames_to_column("gene") # Keep the gene names as a column

# 2. Add a column to classify genes as upregulated, downregulated, or not significant
#    You can adjust the thresholds for padj and log2FoldChange.
lfc_threshold <- 1.0
padj_threshold <- 0.05

res_df <- res_df %>%
  mutate(regulation = case_when(
    padj < padj_threshold & log2FoldChange > lfc_threshold  ~ "Upregulated",
    padj < padj_threshold & log2FoldChange < -lfc_threshold ~ "Downregulated",
    TRUE                                                    ~ "Not Significant"
  ))

# 3. Create a smaller data frame for labeling the top genes
top_genes <- res_df %>%
  filter(regulation != "Not Significant") %>%
  group_by(regulation) %>%
  slice_min(order_by = padj, n = 7) # Select top 20 genes by adjusted p-value

# 4. Create the volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  # Add threshold lines
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed") +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed") +
  # Set custom colors
  scale_color_manual(values = c(
    "Upregulated" = "#e41a1c",
    "Downregulated" = "#377eb8",
    "Not Significant" = "grey80"
  )) +
  # Add labels for the top genes
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3.5,
    box.padding = 0.5,
    max.overlaps = Inf, # Allow all labels to be plotted
    color = "black"
  ) +
  labs(
    title = "Volcano Plot: High vs. Low Egg Production",
    x = "log2 Fold Change",
    y = "-log10(Adjusted P-value)",
    color = "Regulation"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# 5. Save the plot to a file
ggsave("/home/johan/johan/output/chicken/DGE/Volcano_Plot_Egg.png", width = 10, height = 9, plot = volcano_plot)
