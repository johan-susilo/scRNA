#!/bin/bash
set -euo pipefail

date
echo "#################"
echo "# MAP STARTING  #"
echo "#################"

OUTPUT="/home/johan/johan/output/chicken"
RAW="/mnt/2_80T/Data/TS250721001"
LOG="${OUTPUT}/log"
METADATA="/home/johan/pipeline/rnaseq/chicken/sample.csv"
STAR_INDEX="/home/johan/johan/reference/chicken/star_index/STAR_index_chicken_150"

THREADS="40"
mkdir -p "${LOG}"

# (A) build STAR index if missing  (sjdbOverhang choice is OK; STAR authors note 100–149 is generally fine)
Ref_Ara="/home/johan/johan/reference/chicken"
FASTA="${Ref_Ara}/fasta/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"
GTF="${Ref_Ara}/gtf/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf"

if [ ! -s "${STAR_INDEX}/SA" ]; then
  mkdir -p "${STAR_INDEX}"
  echo "Building STAR index at ${STAR_INDEX}"
  STAR --runThreadN "${THREADS}" \
       --runMode genomeGenerate \
       --genomeDir "${STAR_INDEX}" \
       --genomeFastaFiles "${FASTA}" \
       --sjdbGTFfile "${GTF}" \
       --sjdbOverhang 150 \
       --genomeSAindexNbases 13
else
  echo "STAR index found at ${STAR_INDEX}"
fi
# (Notes on sjdbOverhang: ideally read_len-1, but using a larger generic value is typically fine. )  :contentReference[oaicite:1]{index=1}

# (B) collect unique samples
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

# (C) map each sample
for s in "${samples[@]}"; do
  r1="${OUTPUT}/${s}/trim/${s}_R1.trim.fastq.gz"
  r2="${OUTPUT}/${s}/trim/${s}_R2.trim.fastq.gz"
  mapdir="${OUTPUT}/${s}/map"
  mkdir -p "${mapdir}"

  # Ensure write access & clean stale STAR temp dir
  chmod u+rwx,go+rx "${mapdir}" || true
  rm -rf "${mapdir}/${s}._STARtmp" || true   # STAR wants it removed if it exists :contentReference[oaicite:2]{index=2}

  if [[ -s "${r1}" && -s "${r2}" ]]; then
    echo "STAR mapping (PE) for ${s}"
    STAR --runThreadN "${THREADS}" \
         --genomeDir "${STAR_INDEX}" \
         --readFilesIn "${r1}" "${r2}" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${mapdir}/${s}." \
         --outTmpDir "${mapdir}/${s}._STARtmp" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattrRGline ID:${s} SM:${s} PL:ILLUMINA \
         --quantMode GeneCounts
    # outFileNamePrefix controls where outputs go; prefix can include a path. :contentReference[oaicite:3]{index=3}
  else
    echo "WARN: Missing trimmed mates for ${s}; expected ${r1} and ${r2}. Skipping." | tee -a "${LOG}/star.warn.log"
  fi
done

echo "STAR mapping done. BAMs under ${OUTPUT}/<sample>/map"
date