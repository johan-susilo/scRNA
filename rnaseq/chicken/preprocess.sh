#!/bin/bash
################################################################################
# Script: preprocess.sh
# Purpose: Quality control of raw FASTQ files using FastQC
# Step: 1 of RNA-seq pipeline
################################################################################

# Exit on error, undefined variable, or pipe failure
set -euo pipefail

# Print current date/time for logging
date

# Print script header for log visibility
echo "#################"
echo "# PREPROCESS STARTING #"
echo "#################"

################################################################################
# CONFIGURATION - Edit these paths for your system
################################################################################

OUTPUT="/mnt/80T/johan/output/chicken"          # Main output directory
REF="/mnt/80T/johan/reference/chicken"          # Reference genome directory (not used in this script)
LOG="/mnt/80T/johan/output/chicken/log"         # Log files directory
METADATA="/home/johan/pipeline/rnaseq/chicken/sample.csv"  # Sample metadata file (CSV format)

RAW="/mnt/2_80T/Data/TS250721001"               # Directory containing raw FASTQ files
FASTQC_THREADS=1                                # Number of threads per FastQC job (1-4 recommended)
MAX_CONCURRENT=8                                # Maximum number of FastQC jobs to run simultaneously

################################################################################
# CREATE OUTPUT DIRECTORIES
################################################################################

mkdir -p "${OUTPUT}/fastqc"  # Create directory for FastQC HTML reports
mkdir -p "${LOG}"            # Create directory for log files

# Print configuration for verification
echo "Using RAW dir: ${RAW}"
echo "FastQC threads per job: ${FASTQC_THREADS}"
echo "Max concurrent FastQC jobs: ${MAX_CONCURRENT}"

################################################################################
# READ SAMPLE LIST FROM METADATA
################################################################################

# Read sample names from metadata CSV file
# - awk -F,: Split lines by comma (CSV format)
# - $1 !~ /^\/\//: Skip lines starting with // (comments)
# - $1 !~ /^#/: Skip lines starting with # (comments)
# - NF>0: Skip empty lines
# - {print $1}: Print first column (sample names)
# - mapfile -t files: Store results in array called 'files'
mapfile -t files < <(awk -F, '$1 !~ /^\/\// && $1 !~ /^#/ && NF>0 {print $1}' "${METADATA}")

# Check if any samples were found
if [ ${#files[@]} -eq 0 ]; then
    echo "No files found in metadata (${METADATA}). Exiting."
    exit 1  # Exit with error code 1
fi

################################################################################
# BUILD FULL PATHS TO INPUT FILES
################################################################################

# Initialize empty array for input file paths
inputs=()

# For each sample name, create full path to FASTQ file
for f in "${files[@]}"; do
    inputs+=("${RAW}/${f}")  # Append full path to inputs array
done

################################################################################
# CHECK IF FASTQC IS INSTALLED
################################################################################

# Check if fastqc command exists in system PATH
if ! command -v fastqc >/dev/null 2>&1; then
    echo "fastqc not found in PATH. Please install FastQC."
    exit 1  # Exit if FastQC is not installed
fi

################################################################################
# RUN FASTQC (PARALLEL OR SEQUENTIAL)
################################################################################

# Check if GNU parallel is available (recommended for faster processing)
if command -v parallel >/dev/null 2>&1; then
    echo "Running FastQC with GNU parallel..."

    # Run FastQC using GNU parallel for efficient parallelization
    # - parallel: GNU parallel command
    # - -j "${MAX_CONCURRENT}": Run up to MAX_CONCURRENT jobs simultaneously
    # - --lb: Use load balancing for better CPU utilization
    # - fastqc: The FastQC command
    # - -t "${FASTQC_THREADS}": Threads per FastQC job
    # - -o "${OUTPUT}/fastqc": Output directory for results
    # - {}: Placeholder for input file (replaced by parallel)
    # - ::: "${inputs[@]}": Input file list to process
    parallel -j "${MAX_CONCURRENT}" --lb fastqc -t "${FASTQC_THREADS}" -o "${OUTPUT}/fastqc" {} ::: "${inputs[@]}"

else
    # Fallback method if GNU parallel is not available
    echo "GNU parallel not found. Falling back to background-job semaphore loop..."

    # Process files using background jobs with manual semaphore
    for fq in "${inputs[@]}"; do

        # Wait until a job slot is available
        # - jobs -rp: List process IDs of running background jobs
        # - wc -l: Count number of running jobs
        # - while loop: Wait if we've reached MAX_CONCURRENT jobs
        while [ "$(jobs -rp | wc -l)" -ge "${MAX_CONCURRENT}" ]; do
            sleep 1  # Wait 1 second before checking again
        done

        # Start FastQC in background (&) once a slot is available
        echo "Starting FastQC on ${fq}"
        fastqc -t "${FASTQC_THREADS}" -o "${OUTPUT}/fastqc" "${fq}" &
    done

    # Wait for all background jobs to complete
    wait
fi

################################################################################
# COMPLETION MESSAGE
################################################################################

echo "FastQC done. Outputs in ${OUTPUT}/fastqc"
date  # Print completion time

################################################################################
# OUTPUT FILES
################################################################################
# This script generates:
# - ${OUTPUT}/fastqc/<sample>_fastqc.html  - HTML quality report
# - ${OUTPUT}/fastqc/<sample>_fastqc.zip   - Detailed QC data
#
# Review HTML files to check:
# - Per base sequence quality (should be >28)
# - Adapter content (should be low)
# - GC content distribution
# - Sequence duplication levels (expected for RNA-seq)
################################################################################
