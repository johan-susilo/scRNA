#!/bin/bash
set -euo pipefail

date
echo "#################"
echo "# MAP STARTING  #"
echo "#################"

# ============================================
# CONFIGURATION (Edit these paths only)
# ============================================
OUTPUT="/home/johan/johan/output/chicken"
RAW="/mnt/2_80T/Data/TS250721001"
LOG="${OUTPUT}/log"
METADATA="/home/johan/pipeline/rnaseq/chicken/sample.csv"
STAR_INDEX="/home/johan/johan/reference/chicken/star_index/STAR_index_chicken_150"
Ref_Ara="/home/johan/johan/reference/chicken"
FASTA="${Ref_Ara}/fasta/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"
GTF="${Ref_Ara}/gtf/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf"

# ============================================
# AUTO-DETECT SYSTEM RESOURCES
# ============================================
TOTAL_THREADS=$(nproc)
TOTAL_RAM_GB=$(free -g | awk '/^Mem:/ {print $2}')

echo "Detected system resources:"
echo "  - CPU cores: ${TOTAL_THREADS}"
echo "  - Total RAM: ${TOTAL_RAM_GB} GB"

# Calculate optimal parallel jobs based on available RAM
if (( TOTAL_RAM_GB >= 180 )); then
    PARALLEL_JOBS=3
elif (( TOTAL_RAM_GB >= 90 )); then
    PARALLEL_JOBS=2
else
    PARALLEL_JOBS=1
fi

THREADS_PER_JOB=$((TOTAL_THREADS / PARALLEL_JOBS))
SORT_THREADS=$((THREADS_PER_JOB / 4))
[[ ${SORT_THREADS} -lt 2 ]] && SORT_THREADS=2
BAM_SORT_RAM=$(awk "BEGIN {printf \"%.0f\", (${TOTAL_RAM_GB} * 0.6 / ${PARALLEL_JOBS}) * 1000000000}")

echo "Optimized settings:"
echo "  - Parallel jobs: ${PARALLEL_JOBS}"
echo "  - Threads per job: ${THREADS_PER_JOB}"
echo "  - BAM sort threads: ${SORT_THREADS}"
echo "  - RAM per job: $((BAM_SORT_RAM / 1000000000)) GB"
echo ""

mkdir -p "${LOG}"

# ============================================
# CLEANUP OLD TEMP DIRECTORIES
# ============================================
echo "Cleaning up any old STAR temp directories..."
find "${OUTPUT}" -type d -name "*_STARtmp" -exec rm -rf {} + 2>/dev/null || true
echo "Cleanup complete."
echo ""

# ============================================
# BUILD STAR INDEX IF MISSING
# ============================================
if [ ! -s "${STAR_INDEX}/SA" ]; then
  mkdir -p "${STAR_INDEX}"
  echo "Building STAR index at ${STAR_INDEX}"
  STAR --runThreadN "${TOTAL_THREADS}" \
       --runMode genomeGenerate \
       --genomeDir "${STAR_INDEX}" \
       --genomeFastaFiles "${FASTA}" \
       --sjdbGTFfile "${GTF}" \
       --sjdbOverhang 150 \
       --genomeSAindexNbases 13
else
  echo "STAR index found at ${STAR_INDEX}"
fi

# ============================================
# COLLECT UNIQUE SAMPLES
# ============================================
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

echo "Found ${#samples[@]} samples to process"
echo ""

# ============================================
# LOAD GENOME INTO SHARED MEMORY
# ============================================
if [[ ${PARALLEL_JOBS} -gt 1 ]]; then
  echo "Loading genome into shared memory for parallel processing..."
  STAR --genomeDir "${STAR_INDEX}" \
       --genomeLoad LoadAndExit \
       2>&1 | tee "${LOG}/genome_load.log" || echo "Note: Shared memory loading failed, will proceed normally"
  USE_SHARED_MEMORY=true
else
  USE_SHARED_MEMORY=false
fi

# ============================================
# MAPPING FUNCTION
# ============================================
map_sample() {
  local s="$1"
  local r1="${OUTPUT}/${s}/trim/${s}_R1.trim.fastq.gz"
  local r2="${OUTPUT}/${s}/trim/${s}_R2.trim.fastq.gz"
  local mapdir="${OUTPUT}/${s}/map"
  local tmpdir="${mapdir}/${s}._STARtmp"
  
  mkdir -p "${mapdir}"
  
  # Robust temp directory cleanup - try multiple methods
  if [[ -d "${tmpdir}" ]]; then
    echo "[$(date '+%H:%M:%S')] Cleaning existing temp dir for ${s}..."
    # Method 1: Force recursive removal with explicit permissions
    chmod -R u+rwx "${tmpdir}" 2>/dev/null || true
    rm -rf "${tmpdir}" 2>/dev/null || true
    
    # Method 2: If still exists, use find to delete files first
    if [[ -d "${tmpdir}" ]]; then
      find "${tmpdir}" -type f -delete 2>/dev/null || true
      find "${tmpdir}" -type d -delete 2>/dev/null || true
    fi
    
    # Method 3: Last resort - move it out of the way
    if [[ -d "${tmpdir}" ]]; then
      mv "${tmpdir}" "${tmpdir}.old.$$" 2>/dev/null || true
      rm -rf "${tmpdir}.old.$$" 2>/dev/null &
    fi
  fi
  
  if [[ -s "${r1}" && -s "${r2}" ]]; then
    echo "[$(date '+%H:%M:%S')] Starting: ${s}"
    
    if [[ "${USE_SHARED_MEMORY}" == "true" ]]; then
      genome_load_option="LoadAndKeep"
    else
      genome_load_option="NoSharedMemory"
    fi
    
    STAR --runThreadN "${THREADS_PER_JOB}" \
         --genomeDir "${STAR_INDEX}" \
         --genomeLoad "${genome_load_option}" \
         --readFilesIn "${r1}" "${r2}" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${mapdir}/${s}." \
         --outTmpDir "${tmpdir}" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattrRGline ID:${s} SM:${s} PL:ILLUMINA \
         --quantMode GeneCounts \
         --limitBAMsortRAM "${BAM_SORT_RAM}" \
         --outBAMsortingThreadN "${SORT_THREADS}" \
         --outBAMcompression 6 \
         > "${LOG}/${s}.star.log" 2>&1
    
    # Cleanup temp directory after completion
    if [[ -d "${tmpdir}" ]]; then
      chmod -R u+rwx "${tmpdir}" 2>/dev/null || true
      rm -rf "${tmpdir}" 2>/dev/null || true
    fi
    
    echo "[$(date '+%H:%M:%S')] Completed: ${s}"
  else
    echo "WARNING: Missing files for ${s}" | tee -a "${LOG}/warnings.log"
  fi
}

export -f map_sample
export OUTPUT STAR_INDEX THREADS_PER_JOB SORT_THREADS BAM_SORT_RAM LOG USE_SHARED_MEMORY

# ============================================
# RUN MAPPING
# ============================================
if [[ ${PARALLEL_JOBS} -eq 1 ]]; then
  echo "Running sequentially (1 sample at a time)..."
  for s in "${samples[@]}"; do
    map_sample "$s"
  done
else
  echo "Running in parallel (${PARALLEL_JOBS} samples at a time)..."
  echo ""
  
  if command -v parallel &> /dev/null; then
    printf '%s\n' "${samples[@]}" | parallel -j "${PARALLEL_JOBS}" --joblog "${LOG}/parallel.log" map_sample {}
    echo ""
    echo "Progress log saved to: ${LOG}/parallel.log"
  else
    echo "Note: Install 'parallel' for better job management (sudo apt-get install parallel)"
    job_count=0
    for s in "${samples[@]}"; do
      map_sample "$s" &
      ((job_count++))
      if (( job_count >= PARALLEL_JOBS )); then
        wait -n
        ((job_count--))
      fi
    done
    wait
  fi
fi

# ============================================
# FINAL CLEANUP
# ============================================
echo ""
echo "Performing final cleanup of temp directories..."
find "${OUTPUT}" -type d -name "*_STARtmp*" -exec rm -rf {} + 2>/dev/null || true

if [[ "${USE_SHARED_MEMORY}" == "true" ]]; then
  echo "Removing genome from shared memory..."
  STAR --genomeDir "${STAR_INDEX}" --genomeLoad Remove 2>&1 | tee -a "${LOG}/genome_load.log" || true
fi

# ============================================
# SUMMARY
# ============================================
echo ""
echo "================================"
echo "STAR MAPPING COMPLETED"
echo "================================"

success_count=0
for s in "${samples[@]}"; do
  if [[ -s "${OUTPUT}/${s}/map/${s}.Aligned.sortedByCoord.out.bam" ]]; then
    ((success_count++))
  fi
done

echo "Successful: ${success_count}/${#samples[@]} samples"
echo "Output: ${OUTPUT}"
echo "Logs: ${LOG}/"
echo ""

date