#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  ok <- c(
    requireNamespace("optparse", quietly = TRUE),
    requireNamespace("DESeq2",   quietly = TRUE),
    requireNamespace("ggplot2",  quietly = TRUE),
    requireNamespace("ggrepel",  quietly = TRUE),
    requireNamespace("readr",    quietly = TRUE),
    requireNamespace("dplyr",    quietly = TRUE),
    requireNamespace("tibble",   quietly = TRUE),
    requireNamespace("vsn",      quietly = TRUE)
  )
  if (!all(ok)) {
    stop("Install required packages first:\n",
         "  BiocManager::install(c('DESeq2','vsn'))\n",
         "  install.packages(c('optparse','ggplot2','ggrepel','readr','dplyr','tibble'))")
  }
})

library(optparse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)
library(tibble)
library(vsn)

# ---------- helpers ----------
make_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

first_format <- function(x, default = "png") {
  if (is.null(x) || !nzchar(x)) return(default)
  v <- tolower(trimws(unlist(strsplit(x, ",", fixed = TRUE))))
  v <- v[nzchar(v)]
  if (length(v) == 0) default else v[[1]]
}

save_plot_multi <- function(plot, base_path, width = 8, height = 6, formats = c("png","pdf"),
                            transparent_png = FALSE, dpi = 300) {
  formats <- tolower(trimws(unlist(strsplit(formats, ","))))
  formats <- formats[nzchar(formats)]
  if (length(formats) == 0) formats <- "png"

  make_transparent <- function(p) {
    p + theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      legend.background      = element_rect(fill = "transparent", colour = NA),
      legend.box.background  = element_rect(fill = "transparent", colour = NA)
    )
  }

  for (fmt in formats) {
    out <- paste0(base_path, ".", fmt)
    if (fmt == "png") {
      bg <- if (isTRUE(transparent_png)) "transparent" else "white"
      p2 <- if (identical(bg, "transparent")) make_transparent(plot) else plot
      ggsave(out, plot = p2, width = width, height = height, dpi = dpi, bg = bg)
    } else if (fmt == "pdf") {
      p2 <- make_transparent(plot)
      ggsave(out, plot = p2, width = width, height = height, device = grDevices::cairo_pdf, bg = "transparent")
    } else if (fmt == "svg" && requireNamespace("svglite", quietly = TRUE)) {
      p2 <- make_transparent(plot)
      svglite::svglite(out, width = width, height = height, bg = "transparent"); print(p2); dev.off()
    } else {
      # fallback to png
      bg <- if (isTRUE(transparent_png)) "transparent" else "white"
      p2 <- if (identical(bg, "transparent")) make_transparent(plot) else plot
      ggsave(paste0(base_path, ".png"), plot = p2, width = width, height = height, dpi = dpi, bg = bg)
    }
  }
}

setup_fonts <- function(font_family, font_size, font_dir) {
  if (!is.null(font_dir) && !is.na(font_dir) && nzchar(font_dir) && dir.exists(font_dir)) {
    if (requireNamespace("sysfonts", quietly = TRUE) && requireNamespace("showtext", quietly = TRUE)) {
      reg <- list.files(font_dir, pattern = "(?i)\\.(ttf|otf)$", full.names = TRUE)
      if (length(reg) > 0) {
        try({
          sysfonts::font_add(family = font_family, regular = reg[1])
          showtext::showtext_auto()
          message(sprintf("Registered font '%s' from %s", font_family, font_dir))
        }, silent = TRUE)
      }
    }
  }
  theme_set(theme_minimal(base_family = font_family, base_size = font_size))
}

# ---------- CLI ----------
option_list <- list(
  make_option(c("--counts"), type="character", help="Path to combined count matrix (tsv) with GeneID first"),
  make_option(c("--metadata"), type="character", help="Path to metadata csv"),
  make_option(c("--outdir"), type="character", default="./results", help="Output directory"),

  make_option(c("--sampleColumn"), type="character", default="sample", help="Metadata column with sample IDs"),
  make_option(c("--assayTypeColumn"), type="character", default="assay_type", help="Column for assay (e.g., RNA/RPF)"),
  make_option(c("--assayTypeToSelect"), type="character", default="RNA", help="Assay type to analyze"),

  make_option(c("--contrastFactor"), type="character", help="Factor name to contrast (e.g., tissue or egg_production)"),
  make_option(c("--contrastNumerator"), type="character", help="Numerator level"),
  make_option(c("--contrastDenominator"), type="character", help="Denominator level"),
  make_option(c("--batchColumn"), type="character", default=NULL, help="[Optional] Block/adjustment column (e.g., subject or egg_production or tissue)"),

  make_option(c("--orgDb"), type="character", default=NULL, help="[Optional] Organism annotation package"),
  make_option(c("--keyType"), type="character", default=NULL, help="[Optional] Key type for orgDb"),

  make_option(c("--cores"), type="integer", default=1, help="Cores for DESeq2 parallel"),
  make_option(c("--padj"), type="double", default=0.05, help="padj threshold"),
  make_option(c("--lfc"), type="double", default=1.0, help="|log2FC| threshold"),
  make_option(c("--topNLabels"), type="integer", default=10, help="Top N to label on volcano"),
  make_option(c("--topNHeatmap"), type="integer", default=50, help="[reserved]"),
  make_option(c("--filterMethod"), type="character", default="intersection", help="[reserved]"),

  make_option(c("--configFile"), type="character", default="config.yaml", help="[optional config]"),
  make_option(c("--transparentPNG"), action="store_true", default=FALSE, help="Transparent PNGs"),
  make_option(c("--outputFormats"), type="character", default="png,pdf", help="Formats to save plots"),

  make_option(c("--fontFamily"), type="character", default="Source Sans 3", help="Plot font family"),
  make_option(c("--fontSize"), type="integer", default=18, help="Plot base font size"),
  make_option(c("--fontDir"), type="character", default=NA, help="Directory of font files (optional)")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ---------- setup ----------
make_dir(opt$outdir)
setup_fonts(opt$fontFamily, opt$fontSize, opt$fontDir)

message("==== Parameters ====")
for (nm in names(opt)) message(sprintf("  %-20s : %s", nm, paste(opt[[nm]], collapse=",")))

# ---------- load data ----------
if (!file.exists(opt$counts))   stop("Counts file not found: ", opt$counts)
if (!file.exists(opt$metadata)) stop("Metadata file not found: ", opt$metadata)

counts_df <- readr::read_tsv(opt$counts, progress = FALSE, show_col_types = FALSE)
if (!"GeneID" %in% colnames(counts_df)) {
  colnames(counts_df)[1] <- "GeneID"
  warning("Counts file did not have 'GeneID' header; assuming the first column is GeneID.")
}

meta <- readr::read_csv(opt$metadata, show_col_types = FALSE)

# assay filter
if (!is.null(opt$assayTypeColumn) && !is.null(opt$assayTypeToSelect) &&
    opt$assayTypeColumn %in% colnames(meta)) {
  meta <- meta[ meta[[opt$assayTypeColumn]] == opt$assayTypeToSelect, , drop=FALSE]
  if (nrow(meta) == 0) stop("No rows left in metadata after assay type filtering.")
}

# align by sample column
if (!(opt$sampleColumn %in% colnames(meta))) {
  stop("sampleColumn '", opt$sampleColumn, "' not found in metadata.")
}

samples_in_meta <- meta[[opt$sampleColumn]]
col_keep <- c("GeneID", intersect(colnames(counts_df), samples_in_meta))
counts_df <- counts_df[, col_keep, drop = FALSE]

if (ncol(counts_df) < 2) stop("Counts matrix has <2 columns after subsetting; check inputs.")
if (length(setdiff(samples_in_meta, colnames(counts_df))) > 0) {
  missing <- setdiff(samples_in_meta, colnames(counts_df))
  stop("These metadata samples have no columns in counts: ", paste(missing, collapse=", "))
}

# order metadata rows to match counts column order
meta <- meta[ match(colnames(counts_df)[colnames(counts_df)!="GeneID"], meta[[opt$sampleColumn]]), , drop=FALSE]
rownames(meta) <- meta[[opt$sampleColumn]]

# ---------- DESeq2 dataset ----------
count_mat <- counts_df %>% tibble::column_to_rownames("GeneID") %>% as.matrix()
mode(count_mat) <- "integer"

lib_sizes <- colSums(count_mat)
message("Library sizes: ", paste(paste0(names(lib_sizes), "=", lib_sizes), collapse="; "))
if (any(lib_sizes == 0)) {
  zero <- names(lib_sizes)[lib_sizes == 0]
  stop("Zero-count libraries detected: ", paste(zero, collapse=", "))
}

# build design: ~ batch + contrast OR ~ contrast
if (is.null(opt$contrastFactor) || is.null(opt$contrastNumerator) || is.null(opt$contrastDenominator)) {
  stop("Provide --contrastFactor, --contrastNumerator, --contrastDenominator.")
}
if (!is.null(opt$batchColumn) && nzchar(opt$batchColumn) && opt$batchColumn %in% colnames(meta)) {
  design_formula <- as.formula(paste("~", opt$batchColumn, "+", opt$contrastFactor))
  meta[[opt$batchColumn]]     <- as.factor(meta[[opt$batchColumn]])
  meta[[opt$contrastFactor]]  <- as.factor(meta[[opt$contrastFactor]])
  message("Design: ", deparse(design_formula))
} else {
  design_formula <- as.formula(paste("~", opt$contrastFactor))
  meta[[opt$contrastFactor]]  <- as.factor(meta[[opt$contrastFactor]])
  message("Design: ", deparse(design_formula))
}

dds <- DESeqDataSetFromMatrix(countData = round(count_mat), colData = meta, design = design_formula)

# prefilter: keep genes with at least 10 total counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
message("Kept ", sum(keep), " genes with rowSums >= 10.")

# run DESeq
if (opt$cores > 1 && requireNamespace("BiocParallel", quietly = TRUE)) {
  if (.Platform$OS.type == "windows") {
    BPPARAM <- BiocParallel::SnowParam(workers = opt$cores)
  } else {
    BPPARAM <- BiocParallel::MulticoreParam(workers = opt$cores)
  }
  dds <- DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM)
} else {
  dds <- DESeq(dds)
}

# ---------- results & shrink ----------
contrast_vec <- c(opt$contrastFactor, opt$contrastNumerator, opt$contrastDenominator)
res <- results(dds, contrast = contrast_vec)

# LFC shrink (apeglm if available, else normal)
res_shrunk <- NULL
if (requireNamespace("apeglm", quietly = TRUE)) {
  res_shrunk <- tryCatch(
    { lfcShrink(dds, contrast = contrast_vec, type = "apeglm") },
    error = function(e) {
      message("apeglm shrinkage failed, using 'normal'. Reason: ", e$message)
      lfcShrink(dds, contrast = contrast_vec, type = "normal")
    }
  )
} else {
  res_shrunk <- lfcShrink(dds, contrast = contrast_vec, type = "normal")
}

res_df <- as.data.frame(res_shrunk) %>% rownames_to_column("GeneID")
res_df$padj[is.na(res_df$padj)]     <- 1
res_df$pvalue[is.na(res_df$pvalue)] <- 1

# write tables
readr::write_tsv(res_df %>% arrange(padj), file.path(opt$outdir, "DE_results_all.tsv"))
sig_df <- res_df %>%
  filter(!is.na(padj), padj < opt$padj, !is.na(log2FoldChange), abs(log2FoldChange) >= opt$lfc) %>%
  arrange(padj)
readr::write_tsv(sig_df, file.path(opt$outdir, "DE_results_significant.tsv"))

# normalized counts
norm_counts <- counts(dds, normalized = TRUE) %>% as.data.frame() %>% rownames_to_column("GeneID")
readr::write_tsv(norm_counts, file.path(opt$outdir, "normalized_counts.tsv"))

# ---------- QC: VST + PCA ----------
vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)
pca    <- prcomp(t(mat), center = TRUE)
var_ex <- (pca$sdev^2)/sum(pca$sdev^2)
pc_df  <- data.frame(sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2], meta[rownames(pca$x), , drop=FALSE])

color_col <- if ("tissue" %in% colnames(pc_df)) "tissue" else if ("egg_production" %in% colnames(pc_df)) "egg_production" else NULL
shape_col <- if (!is.null(color_col) && color_col == "tissue" && "egg_production" %in% colnames(pc_df)) "egg_production" else NULL

aes_base <- aes(PC1, PC2)
if (!is.null(color_col)) aes_base <- modifyList(aes_base, aes(color = .data[[color_col]]))
if (!is.null(shape_col)) aes_base <- modifyList(aes_base, aes(shape = .data[[shape_col]]))

p_pca <- ggplot(pc_df, aes_base) +
  geom_point(size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3, max.overlaps = 100) +
  labs(title = "PCA (VST)", x = sprintf("PC1 (%.1f%%)", 100*var_ex[1]), y = sprintf("PC2 (%.1f%%)", 100*var_ex[2]),
       color = color_col, shape = shape_col) +
  theme_minimal(base_size = opt$fontSize)

save_plot_multi(p_pca, file.path(opt$outdir, "PCA_samples"),
                width = 8, height = 6, formats = opt$outputFormats, transparent_png = opt$transparentPNG)

# ---------- Volcano ----------
thr_p <- opt$padj
thr_l <- opt$lfc
volc_df <- res_df %>%
  filter(is.finite(padj), padj > 0, is.finite(log2FoldChange)) %>%
  mutate(sig = ifelse(padj < thr_p & abs(log2FoldChange) >= thr_l, "Significant", "NS"))

volc <- ggplot(volc_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c("Significant" = "#d62728", "NS" = "grey70")) +
  geom_vline(xintercept = c(-thr_l, thr_l), linetype = 2, linewidth = 0.4) +
  geom_hline(yintercept = -log10(thr_p), linetype = 2, linewidth = 0.4) +
  labs(title = sprintf("Volcano: %s (%s vs %s)", opt$contrastFactor, opt$contrastNumerator, opt$contrastDenominator),
       x = "log2 fold change", y = "-log10(padj)", color = NULL) +
  theme_minimal(base_size = opt$fontSize) +
  theme(legend.position = "bottom")

lab_df <- volc_df %>% arrange(padj) %>% head(opt$topNLabels)
if (nrow(lab_df) > 0) {
  volc <- volc + ggrepel::geom_text_repel(data = lab_df, aes(label = GeneID), size = 3, max.overlaps = 100)
}
save_plot_multi(volc, file.path(opt$outdir, "Volcano"),
                width = 8, height = 6, formats = opt$outputFormats, transparent_png = opt$transparentPNG)

# ---------- MA plot ----------
ma <- plotMA(res, ylim = c(-5, 5), alpha = thr_p, main = "MA plot")
fmt1  <- first_format(opt$outputFormats, "png")
bgcol <- if (identical(fmt1, "png") && isTRUE(opt$transparentPNG)) "transparent" else "white"
ma_path <- file.path(opt$outdir, paste0("MA.", fmt1))
ggsave(ma_path, plot = ma, width = 7, height = 5, bg = bgcol)

# ---------- session info ----------
try({
  writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "R_sessionInfo.txt"))
}, silent = TRUE)

message("\nDone. Results written to: ", opt$outdir)
