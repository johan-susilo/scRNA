#!/bin/bash

# ensure output dir exists and is writable
OUTDIR="/mnt/80T/johan/output/chicken/RH1_1_R1_001/trim"
mkdir -p "${OUTDIR}" || { echo "Failed to create output dir ${OUTDIR}" >&2; exit 1; }

fastqc -o "${OUTDIR}" "/mnt/80T/johan/output/chicken/RH1_1_R1_001/trim/RH1_1_R1_001.trim.fastq.gz"

