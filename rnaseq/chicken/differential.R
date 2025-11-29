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

# ============================================================================
# CHICKEN EGG PRODUCTION MARKER GENES - SEPARATED BY TISSUE
# Based on peer-reviewed research papers (2020-2025)
# ============================================================================

# OVARY-SPECIFIC MARKER GENES
# ============================================================================
ovary_marker_genes <- c(
  # Steroidogenesis Pathway (Ovary-specific)
  "ENSGALG00010005802",  # STAR - Cholesterol transport
  "ENSGALG00010009244",  # CYP11A1 - Cholesterol → pregnenolone
  "ENSGALG00010005887",  # HSD3B1 - Pregnenolone → progesterone
  "ENSGALG00010004177",  # FSHR - Follicle selection/maintenance
  "ENSGALG00010012742",  # CYP17A1 - Androgen synthesis
  "ENSGALG00010028748",  # CYP19A1 - Estrogen synthesis (aromatase)
  "ENSGALG00010003002",  # HSD11B2 - Steroid hormone metabolism
  
  # Follicle Development (Ovary-specific)
  "ENSGALG00010037139",  # BMP15 - Primordial follicle development
  "ENSGALG00010015635",  # COL4A2 - Follicular ECM structure
  "ENSGALG00010008060",  # RAC1 - Granulosa cell proliferation
  "ENSGALG00010009716",  # SLC5A5 - Granulosa cell differentiation
  "ENSGALG00010046237",  # OSTN - Cell growth/steroidogenesis
  "ENSGALG00010009597",  # GADD45B - Cell growth/steroidogenesis
  "ENSGALG00010046085",  # NFXL1 - Cell growth/steroidogenesis
  "ENSGALG00010010468",  # ADAMTS17 - Oocyte maturation
  
  # Ovarian Development
  "ENSGALG00010006319",  # GRB14 - Ovarian development
  "ENSGALG00010045800",  # GALNT1 - Ovarian function
  "ENSGALG00010000864",  # VLDLR - Lipid metabolism/egg weight
  
  # Ovary Egg Production Markers
  "ENSGALG00010015093",  # CALM1 - Egg production traits
  "ENSGALG00010008468",  # DRD1 - Egg production traits
  "ENSGALG00010008933",  # MMP13 - Egg production traits
  
  # Hub Marker Gene for Egg Production (Ovary)
  "ENSGALG00010003819",  # CYP21A1 - Steroid biosynthesis hub gene
  "ENSGALG00010003186",  # PHIP - Granulosa cell proliferation (functional variant chr3:79510218A>T)
  "ENSGALG00010016722"   # Additional ovary marker
)

# UTERUS (SHELL GLAND)-SPECIFIC MARKER GENES
# ============================================================================
uterus_marker_genes <- c(
  # Eggshell Matrix Proteins (Uterus-specific)
  "ENSGALG00010010927",  # MEPE (OC-116) - Major eggshell matrix protein
  "ENSGALG00010006662",  # BPIFB3 (OVX-36) - Cuticle protein, antimicrobial
  "ENSGALG00010009594",  # RARRES1 (OVX-32) - Cuticle protein
  "ENSGALG00010020972",  # WAP (OVX-25) - Cuticle protein, Kunitz-like protease inhibitor
  "ENSGALG00010023498",  # SPP1 (Osteopontin) - Mineralization, increases at end of shell formation
  
  # Ion Transport & Shell Formation (Uterus)
  "ENSGALG00010005317",  # CA2 - Carbonic anhydrase 2, bicarbonate production for mineralization
  "ENSGALG00010015003",  # ATP1A1 - Na+/K+ ATPase, active ion transport
  "ENSGALG00010005520",  # ATP1B1 - Na+/K+ ATPase beta subunit
  "ENSGALG00010205535",  # SPP1 - Secreted phosphoprotein (duplicate entry, eggshell organic matrix)
  "ENSGALG00010414880",  # SCNN1G - Sodium channel, ion transport
  
  # Circadian Clock Genes (Uterus - timing of cuticle deposition)
  "ENSGALG00010005521",  # PER2 - Period circadian regulator 2
  "ENSGALG00010008436",  # CRY2 - Cryptochrome circadian regulator 2
  "ENSGALG00010012638",  # CRY1 - Cryptochrome circadian regulator 1
  "ENSGALG00010013793",  # CLOCK - Clock circadian regulator
  "ENSGALG00010005378",  # ARNTL (BMAL1) - Aryl hydrocarbon receptor nuclear translocator-like
  
  # Immediate Early Genes (Uterus - oviposition control)
  "ENSGALG00010011Jun",  # JUN - Transcription factor, oviposition & cuticle deposition
  "ENSGALG00010011Fos",  # FOS - Transcription factor, oviposition & cuticle deposition
  
  # Energy Metabolism & Oxidative Phosphorylation (Uterus)
  "ENSGALG00010018373",  # COX1 - Cytochrome c oxidase subunit I
  "ENSGALG00010018367",  # COX3 - Cytochrome c oxidase subunit III
  "ENSGALG00010018370",  # COX2 - Cytochrome c oxidase subunit II
  "ENSGALG00010018360",  # CYTB - Cytochrome b
  "ENSGALG00010018368",  # ATP6 - ATP synthase membrane subunit 6
  "ENSGALG00010018361",  # ND5 - NADH ubiquinone oxidoreductase core subunit 5
  "ENSGALG00010018364",  # ND4 - NADH-ubiquinone oxidoreductase chain 4
  "ENSGALG00010018382",  # ND1 - NADH-ubiquinone oxidoreductase chain 1
  "ENSGALG00010018378",  # ND2 - NADH-ubiquinone oxidoreductase chain 2
  
  # Glycosaminoglycan Binding (Uterus - cuticle glycosylation)
  "ENSGALG00010HBEGF",   # HBEGF - Heparin-binding EGF-like growth factor
  "ENSGALG00010CYR61",   # CYR61 - Cysteine-rich angiogenic inducer 61
  "ENSGALG00010THBS1",   # THBS1 - Thrombospondin 1
  "ENSGALG00010CEMIP",   # CEMIP - Cell migration inducing hyaluronidase 1
  "ENSGALG00010REG4",    # REG4 - Regenerating family member 4
  
  # Apoptosis & Tissue Maintenance (Uterus)
  "ENSGALG00010TNFSF10", # TNFSF10 (TRAIL) - TNF superfamily member 10
  
  # Other Structural/Functional (Uterus)
  "ENSGALG00010015914",  # CALB (Calbindin 1) - Calcium binding
  "ENSGALG00010003578",  # FN1 - Fibronectin 1
  "ENSGALG00010009621",  # ACTB - Actin beta
  "ENSGALG00010014442",  # GAPDH - Glyceraldehyde-3-phosphate dehydrogenase
  "ENSGALG00010GKN2"     # GKN2 (Ovocalyxin 21) - Eggshell specific protein with Brichos domain
)

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
  countData = counts_uterus, 
  colData = uterus, 
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

# 3. Create a hybrid data frame combining custom marker genes and top genes
# First, identify custom marker genes from the results
custom_markers_df <- res_df %>%
  filter(gene %in% uterus_marker_genes)

# Then, select top genes by p-value, excluding any already in custom markers
top_auto_genes <- res_df %>%
  filter(regulation != "Not Significant") %>%
  filter(!gene %in% uterus_marker_genes) %>%  # Exclude custom markers to avoid duplication
  group_by(regulation) %>%
  slice_min(order_by = padj, n = 7)

# Combine both datasets
top_genes <- bind_rows(custom_markers_df, top_auto_genes) %>%
  distinct(gene, .keep_all = TRUE) %>%  # Ensure no duplicates
  mutate(is_custom_marker = gene %in% uterus_marker_genes)  # Mark custom markers for visual distinction

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
  # Add labels for the top genes with visual distinction for custom markers
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3.5,
    box.padding = 0.5,
    max.overlaps = Inf, # Allow all labels to be plotted
    color = ifelse(top_genes$is_custom_marker, "red", "black"),  # Red for custom markers, black for auto-selected
    fontface = ifelse(top_genes$is_custom_marker, "bold", "plain")  # Bold for custom markers
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
