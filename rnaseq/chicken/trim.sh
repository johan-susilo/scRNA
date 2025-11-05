#!/bin/bash
set -euo pipefail

date
echo "#################"
echo "# TRIM STARTING #"
echo "#################"

# --- paths you already use
OUTPUT="/home/johan/johan/output/chicken"
RAW="/mnt/2_80T/Data/TS250721001"
LOG="${OUTPUT}/log"
METADATA="/home/johan/pipeline/rnaseq/chicken/sample.csv"

AVAILABLE_CORES=$(nproc)
CUTADAPT_THREADS=$((AVAILABLE_CORES / 2 > 0 ? AVAILABLE_CORES / 2 : 1))
FASTQC_THREADS=1
MAX_CONCURRENT=$AVAILABLE_CORES

# Replace single ADAPTER with pair-specific adapters:
ADAPTER_R1="ACTGTCTCTTATACACATCT"
ADAPTER_R2="ACTGTCTCTTATACACATCT"
MIN_LENGTH=20

# ensure base dirs exist and are writable
mkdir -p "${OUTPUT}" "${LOG}"
umask 002

echo "RAW: ${RAW}"
echo "OUTPUT: ${OUTPUT}"

# --- build per-sample R1/R2 from METADATA (skip header / comments)
declare -A r1_map r2_map s_seen
while IFS=, read -r col_file col_base col_sample _spp _subj _animal _tissue _tcode _rep _egg col_readpair; do
  [[ -z "${col_file}" ]] && continue
  [[ "${col_file}" == "file" ]] && continue
  [[ "${col_file}" =~ ^# ]] && continue
  file="$(echo "${col_file}" | xargs)"
  sample="$(echo "${col_sample}" | xargs)"
  rp="$(echo "${col_readpair:-}" | xargs)"
  [[ -n "${file}" && -n "${sample}" ]] || continue
  s_seen["$sample"]=1
  if [[ "$rp" == "1" ]]; then
    r1_map["$sample"]="$file"
  elif [[ "$rp" == "2" ]]; then
    r2_map["$sample"]="$file"
  else
    echo "WARN: unknown read_pair '${rp}' for ${file} (sample ${sample}); will treat as SE if used" >&2
  fi
done < "${METADATA}"

(( ${#s_seen[@]} > 0 )) || { echo "No samples found in ${METADATA}"; exit 1; }

# --- build job lists
pairs=()    # fields: in1<<<in2<<<out1<<<out2<<<outdir<<<a1<<<a2
singles=()  # fields: in<<<out<<<outdir<<<adapter
for sample in "${!s_seen[@]}"; do
  in1="${r1_map[$sample]:-}"
  in2="${r2_map[$sample]:-}"
  outdir="${OUTPUT}/${sample}/trim"
  raw_link_dir="${OUTPUT}/${sample}/raw"
  mkdir -p "${raw_link_dir}"

  if [[ -n "$in1" ]]; then
    inp1="${RAW}/${in1}"
    # Create symlink for R1 raw data in sample-specific raw folder
    ln -sf "$inp1" "${raw_link_dir}/${in1##*/}"
  fi

  if [[ -n "$in2" ]]; then
    inp2="${RAW}/${in2}"
    # Create symlink for R2 raw data in sample-specific raw folder
    ln -sf "$inp2" "${raw_link_dir}/${in2##*/}"
  fi

  if [[ -n "$inp1" && -n "$inp2" && -e "$inp1" && -e "$inp2" ]]; then
    out1="${outdir}/${sample}_R1.trim.fastq.gz"
    out2="${outdir}/${sample}_R2.trim.fastq.gz"
    pairs+=("${inp1}<<<${inp2}<<<${out1}<<<${out2}<<<${outdir}<<<${ADAPTER_R1}<<<${ADAPTER_R2}")
  else
    [[ -e "$inp1" ]] && singles+=("${inp1}<<<${outdir}/${sample}_R1.trim.fastq.gz<<<${outdir}<<<${ADAPTER_R1}")
    [[ -e "$inp2" ]] && singles+=("${inp2}<<<${outdir}/${sample}_R2.trim.fastq.gz<<<${outdir}<<<${ADAPTER_R2}")
    [[ -e "$inp1" || -e "$inp2" ]] || echo "WARN: missing both mates for ${sample} (${inp1}, ${inp2})"
  fi
done


if (( ${#pairs[@]} + ${#singles[@]} == 0 )); then
  echo "No valid inputs. Exit."
  exit 1
fi

# --- tool checks
command -v cutadapt >/dev/null || { echo "cutadapt not in PATH"; exit 1; }
command -v fastqc   >/dev/null || { echo "fastqc not in PATH";   exit 1; }

# --- run paired-end jobs
if (( ${#pairs[@]} > 0 )); then
  if command -v parallel >/dev/null 2>&1; then
    printf "%s\n" "${pairs[@]}" | parallel -j "${MAX_CONCURRENT}" --colsep '<<<' --lb '
      in1={1}; in2={2}; out1={3}; out2={4}; outdir={5}; a1={6}; a2={7};
      mkdir -p "$outdir" && chmod u+rwx,go+rx "$outdir" || exit 1
      cutadapt -j '"${CUTADAPT_THREADS}"' -a "$a1" -A "$a2" -m '"${MIN_LENGTH}"' \
        -o "$out1" -p "$out2" "$in1" "$in2" 2>>"'"${LOG}"'/cutadapt.log"
      echo "Finished cutadapt (PE) for $out1 and $out2, starting FastQC..." >>"'"${LOG}"'/fastqc.log"
      fastqc -t '"${FASTQC_THREADS}"' -o "$outdir" "$out1" "$out2" >>"'"${LOG}"'/fastqc.log" 2>&1 \
        || echo "FastQC failed for $out1/$out2" >>"'"${LOG}"'/fastqc.log"
    '
  else
    running=0
    for job in "${pairs[@]}"; do
      while (( running >= MAX_CONCURRENT )); do wait -n || true; running=$((running-1)); done
      (
        IFS='<<<' read -r in1 in2 out1 out2 outdir a1 a2 <<< "${job}"
        mkdir -p "$outdir" && chmod u+rwx,go+rx "$outdir"
        cutadapt -j "${CUTADAPT_THREADS}" -a "${a1}" -A "${a2}" -m "${MIN_LENGTH}" \
          -o "${out1}" -p "${out2}" --trim-n "${in1}" "${in2}" 2>>"${LOG}/cutadapt.log"
        echo "Finished cutadapt (PE) for ${out1}/${out2}, starting FastQC..." >>"${LOG}/fastqc.log"
        fastqc -t "${FASTQC_THREADS}" -o "${outdir}" "${out1}" "${out2}" >>"${LOG}/fastqc.log" 2>&1 \
          || echo "FastQC failed for ${out1}/${out2}" >> "${LOG}/fastqc.log"
      ) &
      running=$((running+1))
    done
    wait
  fi
fi

# --- run single-end jobs
if (( ${#singles[@]} > 0 )); then
  if command -v parallel >/dev/null 2>&1; then
    printf "%s\n" "${singles[@]}" | parallel -j "${MAX_CONCURRENT}" --colsep '<<<' --lb '
      inp={1}; out={2}; outdir={3}; adapter={4};
      mkdir -p "$outdir" && chmod u+rwx,go+rx "$outdir" || exit 1
      cutadapt -j '"${CUTADAPT_THREADS}"' -a "$adapter" -m '"${MIN_LENGTH}"' -o "$out" "$inp" 2>>"'"${LOG}"'/cutadapt.log"
      echo "Finished cutadapt (SE) for $out, starting FastQC..." >>"'"${LOG}"'/fastqc.log"
      fastqc -t '"${FASTQC_THREADS}"' -o "$outdir" "$out" >>"'"${LOG}"'/fastqc.log" 2>&1 \
        || echo "FastQC failed for $out" >>"'"${LOG}"'/fastqc.log"
    '
  else
    running=0
    for job in "${singles[@]}"; do
      while (( running >= MAX_CONCURRENT )); do wait -n || true; running=$((running-1)); done
      (
        IFS='<<<' read -r inp out outdir adapter <<< "${job}"
        mkdir -p "$outdir" && chmod u+rwx,go+rx "$outdir"
        cutadapt -j "${CUTADAPT_THREADS}" -a "${adapter}" -m "${MIN_LENGTH}" -o "${out}" "${inp}" 2>>"${LOG}/cutadapt.log"
        echo "Finished cutadapt (SE) for ${out}, starting FastQC..." >>"${LOG}/fastqc.log"
        fastqc -t "${FASTQC_THREADS}" -o "${outdir}" "${out}" >>"${LOG}/fastqc.log" 2>&1 \
          || echo "FastQC failed for ${out}" >> "${LOG}/fastqc.log"
      ) &
      running=$((running+1))
    done
    wait
  fi
fi

echo "Done. Trimmed reads + FastQC reports under ${OUTPUT}/<sample>/trim"
date
