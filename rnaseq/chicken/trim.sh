#!/bin/bash
# Trim adapters from raw FASTQ files using Cutadapt + FastQC

set -euo pipefail

date
echo "#################"
echo "# TRIM STARTING #"
echo "#################"

# Configuration - Edit these paths for your system
OUTPUT="/home/johan/johan/output/chicken"     # Output directory
RAW="/mnt/2_80T/Data/TS250721001"             # Raw FASTQ directory
LOG="${OUTPUT}/log"                            # Log files directory
METADATA="/home/johan/pipeline/rnaseq/chicken/sample.csv"  # Sample metadata

# Thread configuration
AVAILABLE_CORES=$(nproc)                       # Detect available CPU cores
CUTADAPT_THREADS=$((AVAILABLE_CORES / 2 > 0 ? AVAILABLE_CORES / 2 : 1))  # Half cores for cutadapt
FASTQC_THREADS=1                               # FastQC threads per job
MAX_CONCURRENT=$AVAILABLE_CORES                # Max parallel jobs

# Trimming parameters
ADAPTER_R1="ACTGTCTCTTATACACATCT"              # Illumina TruSeq adapter R1
ADAPTER_R2="ACTGTCTCTTATACACATCT"              # Illumina TruSeq adapter R2
MIN_LENGTH=20                                  # Minimum read length after trimming

# Create directories
mkdir -p "${OUTPUT}" "${LOG}"                  # Ensure output directories exist
umask 002                                      # Set permissions

echo "RAW: ${RAW}"
echo "OUTPUT: ${OUTPUT}"

# Parse metadata CSV to build R1/R2 file pairs per sample
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

# Build job lists for paired-end and single-end files
pairs=()    # Paired-end jobs: in1<<<in2<<<out1<<<out2<<<outdir<<<a1<<<a2
singles=()  # Single-end jobs: in<<<out<<<outdir<<<adapter
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

# Check required tools
command -v cutadapt >/dev/null || { echo "cutadapt not in PATH"; exit 1; }
command -v fastqc   >/dev/null || { echo "fastqc not in PATH";   exit 1; }

# Run paired-end trimming (R1 + R2)
if (( ${#pairs[@]} > 0 )); then
  if command -v parallel >/dev/null 2>&1; then
    printf "%s\n" "${pairs[@]}" | parallel -j "${MAX_CONCURRENT}" --colsep '<<<' --lb '
      in1={1}; in2={2}; out1={3}; out2={4}; outdir={5}; a1={6}; a2={7};
      mkdir -p "$outdir" && chmod u+rwx,go+rx "$outdir" || exit 1
      cutadapt -j '"${CUTADAPT_THREADS}"' -a "$a1" -A "$a2" -m '"${MIN_LENGTH}"' \
        -o "$out1" -p "$out2" "$in1" "$in2" 2>>"'"${LOG}"'/cutadapt.log"  # Trim both R1 and R2
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
          -o "${out1}" -p "${out2}" --trim-n "${in1}" "${in2}" 2>>"${LOG}/cutadapt.log"  # Trim both R1 and R2
        echo "Finished cutadapt (PE) for ${out1}/${out2}, starting FastQC..." >>"${LOG}/fastqc.log"
        fastqc -t "${FASTQC_THREADS}" -o "${outdir}" "${out1}" "${out2}" >>"${LOG}/fastqc.log" 2>&1 \
          || echo "FastQC failed for ${out1}/${out2}" >> "${LOG}/fastqc.log"
      ) &
      running=$((running+1))
    done
    wait
  fi
fi

# Run single-end trimming (orphaned R1 or R2)
if (( ${#singles[@]} > 0 )); then
  if command -v parallel >/dev/null 2>&1; then
    printf "%s\n" "${singles[@]}" | parallel -j "${MAX_CONCURRENT}" --colsep '<<<' --lb '
      inp={1}; out={2}; outdir={3}; adapter={4};
      mkdir -p "$outdir" && chmod u+rwx,go+rx "$outdir" || exit 1
      cutadapt -j '"${CUTADAPT_THREADS}"' -a "$adapter" -m '"${MIN_LENGTH}"' -o "$out" "$inp" 2>>"'"${LOG}"'/cutadapt.log"  # Trim single file
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
        cutadapt -j "${CUTADAPT_THREADS}" -a "${adapter}" -m "${MIN_LENGTH}" -o "${out}" "${inp}" 2>>"${LOG}/cutadapt.log"  # Trim single file
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
