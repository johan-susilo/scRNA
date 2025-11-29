args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  # supports "--name value" and "--name=value"
  val <- NULL
  idx <- which(args == name)
  if (length(idx) > 0 && length(args) >= idx + 1) {
    val <- args[idx + 1]
  } else {
    # check --name=value form
    pat <- paste0("^", gsub("^--", "", name), "=")
    eq <- grep(pat, args)
    if (length(eq) > 0) {
      val <- sub(pat, "", args[eq[1]])
    }
  }
  if (is.null(val) || val == "") return(default)
  return(val)
}

output_dir_default <- "data/data_separated"
output_base_default <- "data"
ref_features_default <- file.path(output_base_default, "features.tsv")

input_dir <- get_arg("--input-dir", output_dir_default)
output_base <- get_arg("--output-dir", output_base_default)
ref_features_path <- get_arg("--ref-features", ref_features_default)

# ensure directories exist
if (!dir.exists(input_dir)) {
  stop("Input directory not found: ", input_dir)
}

# List all .rds files in the input directory
rds_files <- list.files(path = input_dir, pattern = "_seurat_object\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) {
  cat("No Seurat .rds files found in", input_dir, "\n")
}

# Load necessary libraries
library(Seurat)
library(Matrix)

for (rds_file in rds_files) {
  sample_name <- gsub("_seurat_object\\.rds$", "", basename(rds_file))
  sample_dir <- file.path(output_base, sample_name)
  if (!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)

  seurat_object <- readRDS(rds_file)
  raw_counts_matrix <- GetAssayData(seurat_object, slot = "counts", assay = "RNA")
  if (!is(raw_counts_matrix, "sparseMatrix")) {
    raw_counts_matrix <- as(raw_counts_matrix, "sparseMatrix")
  }
  barcodes <- colnames(raw_counts_matrix)
  writeMM(raw_counts_matrix, file.path(sample_dir, "matrix.mtx"))

  # If reference features file is missing, create one from gene names
  if (file.exists(ref_features_path)) {
    file.copy(from = ref_features_path, to = file.path(sample_dir, "features.tsv"), overwrite = TRUE)
  } else {
    warning("Reference features file not found: ", ref_features_path, ". Generating features.tsv from Seurat object.")
    genes <- rownames(raw_counts_matrix)
    write.table(data.frame(genes), file = file.path(sample_dir, "features.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  write.table(barcodes, file.path(sample_dir, "barcodes.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("Loaded Seurat object from:", rds_file, "and saved 10X files in:", sample_dir, "\n")
}

cat("All done.\n")
