#!/bin/bash
set -euo pipefail

date
echo "#################"
echo "# TRIM STARTING #"
echo "#################"

# --- paths you already use
OUTPUT="/home/johan/johan/output/chicken-test3"
RAW="/mnt/2_80T/Data/TS250721001"
LOG="${OUTPUT}/log"
METADATA="/home/johan/pipeline/rnaseq/chicken/sample.csv"

# Replace single ADAPTER with pair-specific adapters:
ADAPTER_R1="CTGTCTCTTATACACATCT"
ADAPTER_R2="CTGTCTCTTATACACATCT"
MIN_LENGTH=20

less "${METADATA}"
exit 0