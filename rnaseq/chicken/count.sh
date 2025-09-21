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
STRANDED="${STRANDED:-no}"     # one of: yes | no | reverse
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

  # name-sort BAM for HTSeq paired-end counting
  name_bam="${countdir}/${s}.name.sorted.bam"
  echo "Name-sorting BAM for ${s} with ${SORT_THREADS} threads..."
  samtools sort -n -@ "${SORT_THREADS}" -o "${name_bam}" "${bam}"

  # htseq-count (paired-end fragments inferred from name-sorted input; no '-p' here)
  out_counts="${countdir}/${s}.htseq.${STRANDED}.counts.tsv"
  out_log="${countdir}/${s}.htseq.${STRANDED}.log"

  echo "HTSeq-count for ${s} -> ${out_counts}"
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

  # optional: remove intermediate to save space (keep if KEEP_NAME_BAM=1)
  if [[ "${KEEP_NAME_BAM:-0}" != "1" ]]; then
    rm -f "${name_bam}" || true
  fi
}

export -f run_one
export OUTPUT LOG GTF STRANDED SORT_THREADS KEEP_NAME_BAM HTSEQ_OPTIONAL_FLAGS

# parallelize across samples
if command -v parallel >/dev/null 2>&1; then
  printf "%b\n" "${jobs[@]}" | parallel -j "${MAX_CONCURRENT}" --colsep '\t' --lb \
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

echo "HTSeq-count done. Outputs in ${OUTPUT}/<sample>/count/*.htseq.${STRANDED}.counts.tsv"
date
