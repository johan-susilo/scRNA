#!/bin/bash

pipeline set -euo

date
echo "#################"
echo "# PREPROCESS STARTING #"
echo "#################"

OUTPUT="/mnt/80T/johan/output/chicken"
REF="/mnt/80T/johan/reference/chicken"
LOG="/mnt/80T/johan/output/chicken/log"
METADATA="/home/johan/pipeline/rnaseq/chicken/sample.csv"

RAW="/mnt/2_80T/Data/TS250721001"    # directory containing the .fastq.gz files; adjust if needed
FASTQC_THREADS=1                      # threads per FastQC run (fastqc -t); set to 1-4 usually
MAX_CONCURRENT=8                      # how many FastQC jobs to run at the same time

mkdir -p "${OUTPUT}/fastqc"
mkdir -p "${LOG}"

echo "Using RAW dir: ${RAW}"
echo "FastQC threads per job: ${FASTQC_THREADS}"
echo "Max concurrent FastQC jobs: ${MAX_CONCURRENT}"

# Build list of files from METADATA (first column), skip lines starting with '/' or '#'
mapfile -t files < <(awk -F, '$1 !~ /^\/\// && $1 !~ /^#/ && NF>0 {print $1}' "${METADATA}")

if [ ${#files[@]} -eq 0 ]; then
    echo "No files found in metadata (${METADATA}). Exiting."
    exit 1
fi

# Full paths to input files (assume they live in $RAW)
inputs=()
for f in "${files[@]}"; do
    inputs+=("${RAW}/${f}")
done

# Check fastqc exists
if ! command -v fastqc >/dev/null 2>&1; then
    echo "fastqc not found in PATH. Please install FastQC."
    exit 1
fi

# If GNU parallel is available, use it (recommended)
if command -v parallel >/dev/null 2>&1; then
    echo "Running FastQC with GNU parallel..."
    parallel -j "${MAX_CONCURRENT}" --lb fastqc -t "${FASTQC_THREADS}" -o "${OUTPUT}/fastqc" {} ::: "${inputs[@]}"
else
    echo "GNU parallel not found. Falling back to background-job semaphore loop..."
    # start background jobs, but never exceed MAX_CONCURRENT
    for fq in "${inputs[@]}"; do
        # wait until we have a free slot
        while [ "$(jobs -rp | wc -l)" -ge "${MAX_CONCURRENT}" ]; do
            sleep 1
        done
        echo "Starting FastQC on ${fq}"
        fastqc -t "${FASTQC_THREADS}" -o "${OUTPUT}/fastqc" "${fq}" &
    done
    # wait for remaining background jobs
    wait
fi

echo "FastQC done. Outputs in ${OUTPUT}/fastqc"
date
