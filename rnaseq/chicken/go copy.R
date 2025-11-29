#!/usr/bin/env Rscript

.libPaths("/home/johan/R/x86_64-pc-linux-gnu-library/4.2")


suppressPackageStartupMessages({
  library(ggplot2)
  library(clusterProfiler)
  library(org.Gg.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  
})

# Inputs
input_file <- "/home/johan/johan/output/chicken/DGE/ovary/result_egg.csv"  # Changed to DE results
gaf_file   <- "/home/johan/johan/reference/chicken/go/goa_chicken.gaf"

# Define significance thresholds
padj_threshold <- 0.05
lfc_threshold <- 0.7  # Standard threshold (lowered from 1.0 to get more genes for enrichment)

# Create output directory with threshold info to avoid overwriting previous results
base_output_dir <- "/home/johan/johan/output/chicken/GO/ovary"
output_dir <- "/home/johan/johan/output/chicken/GO/ovary/padj_0.05_lfc_0.7"
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1) Load differential expression results (expect Ensembl gene IDs in rownames)
de_results <- read.csv(input_file, row.names = 1, check.names = FALSE)
if (nrow(de_results) == 0) stop("DE results table has no rows -> cannot perform GO enrichment.")

# Filter for significantly differentially expressed genes
sig_genes <- de_results[!is.na(de_results$padj) &
                        de_results$padj < padj_threshold &
                        abs(de_results$log2FoldChange) > lfc_threshold, ]

cat("Total genes in DE results:", nrow(de_results), "\n")
cat("Significant DE genes (padj <", padj_threshold, ", |log2FC| >", lfc_threshold, "):", nrow(sig_genes), "\n")

# Separate into upregulated and downregulated
upregulated <- sig_genes[sig_genes$log2FoldChange > lfc_threshold, ]
downregulated <- sig_genes[sig_genes$log2FoldChange < -lfc_threshold, ]

cat("Upregulated genes:", nrow(upregulated), "\n")
cat("Downregulated genes:", nrow(downregulated), "\n")

# Extract gene lists
genes_ensem_all <- rownames(sig_genes)
genes_ensem_up <- rownames(upregulated)
genes_ensem_down <- rownames(downregulated)

if (length(genes_ensem_all) == 0) stop("No significant DE genes -> cannot perform GO enrichment.")

# 2) Prepare ID mappings for ALL significant genes, upregulated, and downregulated
# Function to map gene IDs
map_gene_ids <- function(ensembl_genes, label) {
  # Ensembl -> ENTREZ
  map_entrez <- AnnotationDbi::mapIds(
    org.Gg.eg.db,
    keys = ensembl_genes,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  map_entrez <- map_entrez[!is.na(map_entrez)]
  genes_entrez <- unique(unname(map_entrez))

  # Ensembl -> UniProt
  suppressWarnings({
    sel <- AnnotationDbi::select(org.Gg.eg.db, keys = ensembl_genes,
                                 keytype = "ENSEMBL", columns = c("UNIPROT"))
  })
  sel <- sel[!is.na(sel$UNIPROT) & nzchar(sel$UNIPROT), , drop = FALSE]
  genes_uniprot <- unique(sel$UNIPROT)

  cat(label, "- Mapped to ENTREZ:", length(genes_entrez), "| to UniProt:", length(genes_uniprot),
      "| from Ensembl:", length(ensembl_genes), "\n")

  return(list(entrez = genes_entrez, uniprot = genes_uniprot))
}

# Map all three gene sets
mapped_all <- map_gene_ids(genes_ensem_all, "All significant genes")
mapped_up <- map_gene_ids(genes_ensem_up, "Upregulated genes")
mapped_down <- map_gene_ids(genes_ensem_down, "Downregulated genes")

# For backward compatibility with original code
genes_entrez <- mapped_all$entrez
genes_uniprot <- mapped_all$uniprot

# Helper to safely save plots
save_plot <- function(plot_obj, filepath, width = 9, height = 7) {
  tryCatch({
    ggsave(filename = filepath, plot = plot_obj, width = width, height = height, limitsize = FALSE)
    cat("Saved:", filepath, "\n")
  }, error = function(e) {
    cat("Failed to save", basename(filepath), ":", conditionMessage(e), "\n")
  })
}

# Build plots from an enrichResult
generate_plots <- function(eresult, label, subdir = "") {
  if (is.null(eresult) || nrow(eresult) == 0) {
    cat("No enriched GO terms found for:", label, "\n")
    return(invisible(NULL))
  }

  cat("Top terms for", label, ":\n"); print(utils::head(eresult))

  # Create subdirectory if specified
  plot_dir <- output_dir
  if (nzchar(subdir)) {
    plot_dir <- file.path(output_dir, subdir)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Compute similarity matrix where possible
  es <- tryCatch({ pairwise_termsim(eresult) }, error = function(e) eresult)

  # Build plots one by one to be robust against missing deps
  # Dotplot
  p <- tryCatch({ dotplot(eresult, showCategory = min(20, nrow(eresult))) }, error = function(e) NULL)
  if (!is.null(p)) save_plot(p, file.path(plot_dir, paste0("GO_", label, "_dotplot.pdf")))

  # Barplot
  p <- tryCatch({ barplot(eresult, showCategory = min(20, nrow(eresult))) }, error = function(e) NULL)
  if (!is.null(p)) save_plot(p, file.path(plot_dir, paste0("GO_", label, "_barplot.pdf")))

  # Cnetplot
  p <- tryCatch({ cnetplot(es, showCategory = min(10, nrow(eresult)), circular = FALSE, colorEdge = TRUE) }, error = function(e) NULL)
  if (!is.null(p)) save_plot(p, file.path(plot_dir, paste0("GO_", label, "_cnetplot.pdf")))

  # Emapplot
  p <- tryCatch({ emapplot(es, layout = "kk") }, error = function(e) NULL)
  if (!is.null(p)) save_plot(p, file.path(plot_dir, paste0("GO_", label, "_emapplot.pdf")))

  # Treeplot
  p <- tryCatch({ treeplot(eresult) }, error = function(e) NULL)
  if (!is.null(p)) save_plot(p, file.path(plot_dir, paste0("GO_", label, "_treeplot.pdf")))

  # Save enrichment results as CSV
  result_file <- file.path(plot_dir, paste0("GO_", label, "_results.csv"))
  write.csv(as.data.frame(eresult), result_file, row.names = FALSE)
  cat("Saved enrichment results to:", result_file, "\n")
}

# Strategy A: Use GAF if available (build custom TERM2GENE and use enricher)
use_gaf <- file.exists(gaf_file)
if (use_gaf && length(genes_uniprot) > 0) {
  cat("\nUsing GAF:", gaf_file, "\n")
  gaf <- tryCatch(
    read.delim(gaf_file, header = FALSE, comment.char = "!", quote = "", stringsAsFactors = FALSE, fill = TRUE),
    error = function(e) NULL
  )
  if (!is.null(gaf) && ncol(gaf) >= 9) {
    # Columns (GAF 2.x): V2 = DB Object ID (UniProt), V5 = GO ID, V9 = Aspect (P/F/C)
    gaf <- gaf[!is.na(gaf$V2) & nzchar(gaf$V2) & !is.na(gaf$V5) & nzchar(gaf$V5), ]
    term2gene <- unique(gaf[, c(5, 2)])
    colnames(term2gene) <- c("GO", "GENE")

    # Optional: GO term names if GO.db installed
    term2name <- NULL
    if (requireNamespace("GO.db", quietly = TRUE)) {
      suppressPackageStartupMessages(library(GO.db))
      term2name <- AnnotationDbi::select(GO.db, keys = unique(term2gene$GO), columns = "TERM", keytype = "GOID")
      colnames(term2name) <- c("GO", "NAME")
      term2name <- term2name[!is.na(term2name$NAME), , drop = FALSE]
    }

    # Per-ontology mapping using Aspect (P/BP, F/MF, C/CC)
    gaf$V9[is.na(gaf$V9)] <- ""
    onto_map <- list(
      BP = unique(gaf[gaf$V9 == "P", c(5, 2)]),
      MF = unique(gaf[gaf$V9 == "F", c(5, 2)]),
      CC = unique(gaf[gaf$V9 == "C", c(5, 2)])
    )

    # Run enrichment for All, Upregulated, and Downregulated gene sets
    gene_sets <- list(
      all = list(genes = mapped_all$uniprot, subdir = "all_significant", label_suffix = "all"),
      up = list(genes = mapped_up$uniprot, subdir = "upregulated", label_suffix = "up"),
      down = list(genes = mapped_down$uniprot, subdir = "downregulated", label_suffix = "down")
    )

    for (set_name in names(gene_sets)) {
      set_info <- gene_sets[[set_name]]
      if (length(set_info$genes) == 0) {
        cat("Skipping", set_name, "- no genes mapped to UniProt\n")
        next
      }

      cat("\n=== Processing", toupper(set_name), "genes ===\n")
      for (nm in names(onto_map)) {
        df <- onto_map[[nm]]
        colnames(df) <- c("GO", "GENE")
        if (nrow(df) == 0) {
          cat("No mapping rows in GAF for", nm, "\n"); next
        }
        cat("Running GAF-based enrichment for", nm, "(", set_name, "genes:", length(set_info$genes), ")\n")
        eres <- tryCatch({
          enricher(gene = set_info$genes, TERM2GENE = df, TERM2NAME = term2name, pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2)
        }, error = function(e) NULL)
        generate_plots(eres, paste0(nm, "_GAF_", set_info$label_suffix), set_info$subdir)
      }
    }
  } else {
    cat("Failed to read/parse GAF. Falling back to OrgDb enrichGO.\n")
    use_gaf <- FALSE
  }
} else if (use_gaf && length(genes_uniprot) == 0) {
  cat("GAF found but no UniProt IDs mapped from your Ensembl list. Falling back to OrgDb enrichGO.\n")
  use_gaf <- FALSE
}

# Strategy B: Fallback to standard enrichGO with OrgDb
if (!use_gaf) {
  if (length(genes_entrez) == 0) stop("No genes mapped to ENTREZ -> cannot perform GO enrichment.")

  # Run enrichment for All, Upregulated, and Downregulated gene sets
  gene_sets_entrez <- list(
    all = list(genes = mapped_all$entrez, subdir = "all_significant", label_suffix = "all"),
    up = list(genes = mapped_up$entrez, subdir = "upregulated", label_suffix = "up"),
    down = list(genes = mapped_down$entrez, subdir = "downregulated", label_suffix = "down")
  )

  for (set_name in names(gene_sets_entrez)) {
    set_info <- gene_sets_entrez[[set_name]]
    if (length(set_info$genes) == 0) {
      cat("Skipping", set_name, "- no genes mapped to ENTREZ\n")
      next
    }

    cat("\n=== Processing", toupper(set_name), "genes with OrgDb ===\n")
    for (ont in c("BP", "CC", "MF")) {
      cat("Running enrichGO for", ont, "(", set_name, "ENTREZ IDs:", length(set_info$genes), ")...\n")
      ego <- tryCatch({
        enrichGO(
          gene = set_info$genes,
          OrgDb = org.Gg.eg.db,
          keyType = "ENTREZID",
          ont = ont,
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        )
      }, error = function(e) NULL)
      generate_plots(ego, paste0(ont, "_", set_info$label_suffix), set_info$subdir)
    }
  }
}

cat("\nGO enrichment workflow complete. Plots saved to:\n", output_dir, "\n")
