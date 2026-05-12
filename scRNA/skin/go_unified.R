suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db) # Human database for gene translation
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(optparse)
})

# ==============================================================================
# 0. HELPER FUNCTION: SAVE PDF AND HIGH-RES PNG
# ==============================================================================
save_dual_format <- function(plot_obj, file_path_no_ext, w = 12, h = 8) {
  pdf(paste0(file_path_no_ext, ".pdf"), width = w, height = h)
  print(plot_obj)
  dev.off()
  
  png(paste0(file_path_no_ext, ".png"), width = w, height = h, units = "in", res = 300)
  print(plot_obj)
  dev.off()
}

# ==============================================================================
# 1. Define the Upgraded Comparative Pathway Analysis
# ==============================================================================
run_pathway_analysis <- function(csv_path, out_base_dir) {
  
  # Extract names for titling and saving
  file_name <- basename(csv_path)
  clean_name <- str_remove(file_name, "^deseq2_") %>% str_remove("\\.csv$")
  cell_type <- str_extract(clean_name, "^[^_]+")
  
  message(paste("\n-> Running Comparative Pathway Analysis for:", clean_name))
  
  out_dir <- file.path(out_base_dir, cell_type, "pathways")
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
  
  # Load DESeq2 DGE data
  dge_data <- read.csv(csv_path)
  
  # --- STEP 1: Filter and Split (Up vs Down) ---
  # We only want significant genes, separated by their direction of change
  sig_genes <- dge_data %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5)
  
  up_genes <- sig_genes %>% filter(log2FoldChange > 0) %>% pull(gene)
  down_genes <- sig_genes %>% filter(log2FoldChange < 0) %>% pull(gene)
  
  if (length(up_genes) < 10 && length(down_genes) < 10) {
    message("   [SKIP] Not enough significant genes for meaningful pathway analysis.")
    return(NULL)
  }
  
  # --- STEP 2: Gene ID Translation ---
  up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  
  # Create a named list for compareCluster
  gene_list <- list()
  if(length(up_entrez) > 0) gene_list[["Upregulated_in_PMH"]] <- up_entrez
  if(length(down_entrez) > 0) gene_list[["Downregulated_in_PMH"]] <- down_entrez

  # --- STEP 3: RUN GO (Biological Process) COMPARISON ---
  message("   ... Running GO (BP)")
  comp_go <- compareCluster(geneCluster = gene_list, fun = "enrichGO", 
                            OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)
  
  if (!is.null(comp_go) && nrow(comp_go@compareClusterResult) > 0) {
    p_go <- dotplot(comp_go, showCategory = 10, title = paste("GO Biological Processes:", clean_name)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12),
            plot.title = element_text(face="bold", hjust = 0.5))
            
    save_dual_format(p_go, file.path(out_dir, paste0(clean_name, "_GO_Comparative")), w = 8, h = 12)
  }

  # --- STEP 4: RUN KEGG COMPARISON (With Noise Filtering) ---
  message("   ... Running KEGG")
  comp_kegg <- compareCluster(geneCluster = gene_list, fun = "enrichKEGG", 
                              organism = "hsa", pvalueCutoff = 0.05)
  
  if (!is.null(comp_kegg) && nrow(comp_kegg@compareClusterResult) > 0) {
    
    # FILTERING NOISE: Remove confusing infection/viral terms requested by professors
    noisy_terms <- c("infection", "Malaria", "disease", "virus", "viral", "Salmonella", "COVID-19", "Tuberculosis")
    clean_kegg_results <- comp_kegg@compareClusterResult %>%
      filter(!grepl(paste(noisy_terms, collapse="|"), Description, ignore.case = TRUE))
    
    comp_kegg@compareClusterResult <- clean_kegg_results
    
    if(nrow(comp_kegg@compareClusterResult) > 0) {
      p_kegg <- dotplot(comp_kegg, showCategory = 10, title = paste("Cleaned KEGG Pathways:", clean_name)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=12),
              plot.title = element_text(face="bold", hjust = 0.5))
              
      save_dual_format(p_kegg, file.path(out_dir, paste0(clean_name, "_KEGG_Comparative")), w = 8, h = 12)
    }
  }
}

# ==============================================================================
# 2. RUN THE LOOP
# ==============================================================================
option_list <- list(
  make_option(c("-d", "--dir"), type="character", help="Base directory containing cell type subfolders")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$dir)) {
  stop("Missing required argument: --dir. Use --help for options.")
}

base_subset_dir <- opt$dir

message("\n==================================================================")
message("=== Starting Comparative Pathway Analysis (Up vs Down) ===")
message("==================================================================")

cell_type_folders <- list.dirs(base_subset_dir, recursive = FALSE, full.names = FALSE)
cell_type_folders <- cell_type_folders[cell_type_folders != "logs"]

for (cell_type in cell_type_folders) {
  
  dge_dir <- file.path(base_subset_dir, cell_type, "dge_pseudobulk")
  if (!dir.exists(dge_dir)) { next }
  
  dge_files <- list.files(dge_dir, pattern = "^deseq2_.*\\.csv$", full.names = TRUE, recursive = TRUE)
  
  for (csv_file in dge_files) {
    tryCatch({
      run_pathway_analysis(csv_file, base_subset_dir)
    }, error = function(e) {
      message(paste("   [ERROR] Failed to process", basename(csv_file), ":", conditionMessage(e)))
    })
  }
}

message("\n==================================================================")
message("=== All Comparative Pathway Analyses Completed Successfully! ===")