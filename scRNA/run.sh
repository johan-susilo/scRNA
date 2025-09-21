#!/usr/bin/env bash
set -euo pipefail

PIPELINE_DIR="/home/johan/pipeline/scRNA"
LOG_DIR="$PIPELINE_DIR/log"
mkdir -p "$LOG_DIR"
cd "$PIPELINE_DIR" || exit 1

Rscript_bin="$(command -v Rscript || true)"
if [ -z "$Rscript_bin" ]; then
  echo "Rscript not found in PATH. Install R and ensure Rscript is available." >&2
  exit 1
fi

# Scripts / inputs (adjust paths if you moved files)
SEURAT_SCRIPT="seurat_to_matrix.R"
PREPROCESS_SCRIPT="preprocessing.R"
CELL_ANN_SCRIPT="cell_annotation.R"
CNV_SCRIPT="cnv.R"
# changed to CSV input
INPUT_CSV="input.csv"

# Expected RDS output from preprocessing/integration
OUTPUT_BASE="/mnt/80T/johan/output/full_liver_2021_test"
INTEGRATED_RDS="$OUTPUT_BASE/TN.combined_dim30.rds"
ANNOT_OUTPUT="$OUTPUT_BASE/annotations"

SEURAT_INPUT_DIR="/mnt/80T/johan/data/liver_R/data/data_separated"
SEURAT_OUTPUT_BASE="$OUTPUT_BASE/seurat"      # seurat_to_matrix will write per-sample 10X folders here
SEURAT_REF_FEATURES="$OUTPUT_BASE/features.tsv"


log_and_run() {
  local name="$1"; shift
  local logfile="$LOG_DIR/${name}.log"
  echo "---- [$(date '+%F %T')] START $name ----" | tee -a "$logfile"
  # run command, append both stdout/stderr to logfile and also print to console
  ("$@" 2>&1) | tee -a "$logfile"
  echo "---- [$(date '+%F %T')] END   $name ----" | tee -a "$logfile"
}

run_seurat_to_matrix() {
  if [ ! -f "$SEURAT_SCRIPT" ]; then
    echo "Missing $SEURAT_SCRIPT in $PIPELINE_DIR" >&2; return 1
  fi
  mkdir -p "$SEURAT_OUTPUT_BASE"
  log_and_run "seurat_to_matrix" "$Rscript_bin" "$SEURAT_SCRIPT" \
    --input-dir "$SEURAT_INPUT_DIR" --output-dir "$SEURAT_OUTPUT_BASE" --ref-features "$SEURAT_REF_FEATURES"
}

run_preprocessing() {
  if [ ! -f "$PREPROCESS_SCRIPT" ]; then
    echo "Missing $PREPROCESS_SCRIPT in $PIPELINE_DIR" >&2; return 1
  fi
  if [ ! -f "$INPUT_CSV" ]; then
    echo "Missing $INPUT_CSV in $PIPELINE_DIR" >&2; return 1
  fi
  # preprocessing.R now accepts -o for output base
  mkdir -p "$OUTPUT_BASE"
  log_and_run "preprocessing" "$Rscript_bin" "$PREPROCESS_SCRIPT" -f "$INPUT_CSV" -s all -o "$OUTPUT_BASE"
}

run_cell_annotation() {
  if [ ! -f "$CELL_ANN_SCRIPT" ]; then
    echo "Missing $CELL_ANN_SCRIPT in $PIPELINE_DIR" >&2; return 1
  fi
  if [ ! -f "$INTEGRATED_RDS" ]; then
    echo "Integrated RDS not found at $INTEGRATED_RDS. Run preprocessing/integrate first." >&2; return 1
  fi
  mkdir -p "$ANNOT_OUTPUT"
  log_and_run "cell_annotation" "$Rscript_bin" "$CELL_ANN_SCRIPT" --rds "$INTEGRATED_RDS" --output "$ANNOT_OUTPUT"
}

run_cnv() {
  if [ ! -f "$CNV_SCRIPT" ]; then
    echo "Missing $CNV_SCRIPT in $PIPELINE_DIR" >&2; return 1
  fi
  if [ ! -f "$INTEGRATED_RDS" ]; then
    echo "Integrated RDS not found at $INTEGRATED_RDS. Run preprocessing/integrate first." >&2; return 1
  fi
  mkdir -p "$ANNOT_OUTPUT"
  log_and_run "cnv" "$Rscript_bin" "$CNV_SCRIPT" --rds "$INTEGRATED_RDS" --output "$ANNOT_OUTPUT"
}

usage() {
  cat <<EOF
Usage: $0 <step>
Steps:
  all                 run seurat_to_matrix -> preprocessing -> cell_annotation -> cnv (sequential)
  seurat_to_matrix
  preprocess          run preprocessing (reads input.csv)
  annotate            run cell_annotation (requires integrated RDS)
  cnv                 run CNV.R (requires integrated RDS)
  help
EOF
  exit 1
}


if [ $# -lt 1 ]; then
  usage
fi

case "$1" in
  all)
    run_seurat_to_matrix
    run_preprocessing
    # wait briefly for filesystem to settle (if using network mounts)
    sleep 2
    run_cell_annotation
    run_cnv
    ;;
  seurat_to_matrix) run_seurat_to_matrix ;;
  preprocess) run_preprocessing ;;
  annotate) run_cell_annotation ;;
  cnv) run_cnv ;;
  help|--help|-h) usage ;;
  *)
    echo "Unknown command: $1" >&2
    usage
    ;;
esac

echo "Done. Logs are in $LOG_DIR"