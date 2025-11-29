#!/bin/bash
set -euo pipefail

date
echo "#################"
echo "# HTSEQ-COUNT   #"
echo "#################"

# paths consistent with map.sh
OUTPUT="/home/johan/johan/output/chicken"
LOG="${OUTPUT}/log"
METADATA="/home/johan/pipeline/rnaseq/chicken/sample.csv"

# annotation (same reference root as star.sh/map.sh)
Ref_Ara="/home/johan/johan/reference/chicken"
GTF="${Ref_Ara}/gtf/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf"

# configuration
STRANDED="${STRANDED:-reverse}"     # one of: yes | no | reverse
MAX_CONCURRENT="${MAX_CONCURRENT:-4}"

# Detect optional flags supported by your htseq-count version
HTSEQ_OPTIONAL_FLAGS=""
if htseq-count --help 2>&1 | grep -q -- '--secondary-alignments'; then
  HTSEQ_OPTIONAL_FLAGS+=" --secondary-alignments ignore"
fi
if htseq-count --help 2>&1 | grep -q -- '--supplementary-alignments'; then
  HTSEQ_OPTIONAL_FLAGS+=" --supplementary-alignments ignore"
fi
if htseq-count --help 2>&1 | grep -q -- '--nonunique'; then
  HTSEQ_OPTIONAL_FLAGS+=" --nonunique none"
fi

# CPU and per-job threads (used for samtools sort -n)
TOTAL_CPUS="${SLURM_CPUS_ON_NODE:-$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN || echo 8)}"
if (( MAX_CONCURRENT < 1 )); then MAX_CONCURRENT=1; fi
SORT_THREADS=$(( TOTAL_CPUS / MAX_CONCURRENT ))
(( SORT_THREADS < 1 )) && SORT_THREADS=1

mkdir -p "${LOG}"

# tool checks
command -v samtools >/dev/null 2>&1 || { echo "samtools not in PATH"; exit 1; }
command -v htseq-count >/dev/null 2>&1 || { echo "htseq-count not in PATH"; exit 1; }
[ -s "${GTF}" ] || { echo "GTF not found: ${GTF}"; exit 1; }

echo "Using strandedness=${STRANDED}, MAX_CONCURRENT=${MAX_CONCURRENT}, SORT_THREADS=${SORT_THREADS}"

# ===========================================
# CLEANUP OLD TEMP FILES AT START
# ===========================================
echo "Cleaning up old samtools temp files..."
find "${OUTPUT}" -type f -name "*.name.sorted.bam.tmp.*.bam" -delete 2>/dev/null || true
find "${OUTPUT}" -type f -name "*.tmp.*.bam" -delete 2>/dev/null || true
echo "Cleanup complete."

# collect unique samples
declare -A seen
samples=()
while IFS=, read -r col_file col_base col_sample _rest; do
  [[ -z "${col_file}" ]] && continue
  [[ "${col_file}" == "file" ]] && continue
  [[ "${col_file}" =~ ^# ]] && continue
  s="$(echo "${col_sample}" | xargs)"
  [[ -n "${s}" ]] || continue
  if [[ -z "${seen[$s]:-}" ]]; then
    samples+=("$s")
    seen["$s"]=1
  fi
done < "${METADATA}"

# build jobs (sample<TAB>bam<TAB>countdir)
jobs=()
for s in "${samples[@]}"; do
  bam="${OUTPUT}/${s}/map/${s}.Aligned.sortedByCoord.out.bam"
  countdir="${OUTPUT}/${s}/count"
  if [[ -s "${bam}" ]]; then
    jobs+=("${s}\t${bam}\t${countdir}")
  else
    echo "WARN: missing BAM for ${s}: ${bam}" | tee -a "${LOG}/htseq.warn.log"
  fi
done

(( ${#jobs[@]} > 0 )) || { echo "No BAMs to quantify. Exit."; exit 0; }

run_one() {
  local s="$1" bam="$2" countdir="$3"
  mkdir -p "${countdir}"

  # Create unique temp directory for this sample to avoid conflicts
  local temp_dir="${countdir}/samtools_tmp_${s}_$$"
  mkdir -p "${temp_dir}"
  
  # Clean up any existing temp files for this sample
  rm -f "${countdir}/${s}.name.sorted.bam.tmp."*.bam 2>/dev/null || true
  rm -f "${countdir}"/*.tmp.*.bam 2>/dev/null || true

  # name-sort BAM for HTSeq paired-end counting
  name_bam="${countdir}/${s}.name.sorted.bam"
  temp_prefix="${temp_dir}/${s}.tmp"
  
  echo "[$(date '+%H:%M:%S')] Name-sorting BAM for ${s} with ${SORT_THREADS} threads..."
  
  # Use -T to specify unique temp prefix for this sample
  samtools sort -n -@ "${SORT_THREADS}" \
    -T "${temp_prefix}" \
    -o "${name_bam}" \
    "${bam}" 2>&1 | tee -a "${LOG}/${s}.samtools.log"

  # htseq-count (paired-end fragments inferred from name-sorted input)
  out_counts="${countdir}/${s}.htseq.${STRANDED}.counts.tsv"
  out_log="${countdir}/${s}.htseq.${STRANDED}.log"

  echo "[$(date '+%H:%M:%S')] HTSeq-count for ${s} -> ${out_counts}"
  htseq-count \
    -f bam \
    -r name \
    -s "${STRANDED}" \
    -t exon \
    -i gene_id \
    -m union \
    ${HTSEQ_OPTIONAL_FLAGS} \
    "${name_bam}" \
    "${GTF}" > "${out_counts}" 2> "${out_log}"

  # Clean up temp directory
  rm -rf "${temp_dir}" 2>/dev/null || true

  # optional: remove intermediate name-sorted BAM to save space
  if [[ "${KEEP_NAME_BAM:-0}" != "1" ]]; then
    rm -f "${name_bam}" || true
  fi
  
  echo "[$(date '+%H:%M:%S')] Completed: ${s}"
}

export -f run_one
export OUTPUT LOG GTF STRANDED SORT_THREADS KEEP_NAME_BAM HTSEQ_OPTIONAL_FLAGS

# parallelize across samples
if command -v parallel >/dev/null 2>&1; then
  echo "Running with GNU parallel (${MAX_CONCURRENT} concurrent jobs)..."
  printf "%b\n" "${jobs[@]}" | parallel -j "${MAX_CONCURRENT}" --colsep '\t' \
    'run_one {1} {2} {3}'
else
  echo "GNU parallel not found; running with up to ${MAX_CONCURRENT} background jobs."
  running=0
  for line in "${jobs[@]}"; do
    while (( running >= MAX_CONCURRENT )); do wait -n || true; running=$((running-1)); done
    IFS=$'\t' read -r s bam countdir <<< "${line}"
    ( run_one "${s}" "${bam}" "${countdir}" ) &
    running=$((running+1))
  done
  wait
fi

# ===========================================
# FINAL CLEANUP
# ===========================================
echo ""
echo "Performing final cleanup of temp files..."
find "${OUTPUT}" -type f -name "*.tmp.*.bam" -delete 2>/dev/null || true
find "${OUTPUT}" -type d -name "samtools_tmp_*" -exec rm -rf {} + 2>/dev/null || true

# ===========================================
# SUMMARY
# ===========================================
echo ""
echo "================================"
echo "HTSEQ-COUNT COMPLETED"
echo "================================"

success_count=0
for s in "${samples[@]}"; do
  counts_file="${OUTPUT}/${s}/count/${s}.htseq.${STRANDED}.counts.tsv"
  if [[ -s "${counts_file}" ]]; then
    ((success_count++))
  fi
done

echo "Successful: ${success_count}/${#samples[@]} samples"
echo "Output: ${OUTPUT}/<sample>/count/*.htseq.${STRANDED}.counts.tsv"
echo "Logs: ${LOG}/"
echo ""

date