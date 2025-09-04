output_dir <- "data/data_separated"
# List all .rds files in the output directory
rds_files <- list.files(path = output_dir, pattern = "_seurat_object\\.rds$", full.names = TRUE)

# Load necessary libraries
library(Seurat)
library(Matrix)

# Path to the reference features.tsv file
ref_features_path <- "data/features.tsv"

# Loop through each .rds file
for (rds_file in rds_files) {
  # Extract the sample name from the filename
  # Assuming filenames are like "SAMPLE_ID_seurat_object.rds"
  sample_name <- gsub("_seurat_object\\.rds$", "", basename(rds_file))

  # Create a directory for the sample if it doesn't exist
  sample_dir <- file.path("data", sample_name)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir)
  }

  # Read the .rds file (Seurat object)
  seurat_object <- readRDS(rds_file)

  # Extract raw counts matrix
  raw_counts_matrix <- GetAssayData(seurat_object, slot = "counts", assay = "RNA")

  # Convert to sparse matrix if necessary
  if (!is(raw_counts_matrix, "sparseMatrix")) {
    raw_counts_matrix <- as(raw_counts_matrix, "sparseMatrix")
  }

  # Extract barcodes (cell names)
  barcodes <- colnames(raw_counts_matrix)

  # Save raw counts matrix to matrix.mtx
  writeMM(raw_counts_matrix, file.path(sample_dir, "matrix.mtx"))

  # Copy the reference features.tsv into the sample directory
  file.copy(from = ref_features_path, to = file.path(sample_dir, "features.tsv"), overwrite = TRUE)

  # Save barcodes to barcodes.tsv
  write.table(barcodes, file.path(sample_dir, "barcodes.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Optional: Print a message to confirm
  cat("Loaded Seurat object from:", rds_file, "and saved as 10X files in directory:", sample_dir, "\n")
}

cat("All Seurat objects converted to 10X files and saved in individual directories.\n")
cat("Reference features.tsv copied into each sample directory.\n")
