#!/bin/bash
# Build STAR genome index for chicken (ONE-TIME setup, run before pipeline)

set -euo pipefail

date

# Reference directory
Ref_Ara='/home/johan/johan/reference/chicken'

# Create index output directory
mkdir -p "${Ref_Ara}/star_index/STAR_index_chicken_150"

# Build STAR index
STAR --runThreadN "${SLURM_CPUS_PER_TASK:-8}" \  # CPU threads (8 default)
     --runMode genomeGenerate \                   # Build index mode
     --genomeDir "${Ref_Ara}/star_index/STAR_index_chicken_150" \  # Output directory
     --genomeFastaFiles "${Ref_Ara}/fasta/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa" \  # Genome FASTA
     --sjdbGTFfile "${Ref_Ara}/gtf/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf" \  # Gene annotation
     --sjdbOverhang 150 \       # Read length - 1 (for 151bp reads)
     --genomeSAindexNbases 13   # Index parameter for ~1Gb genome

echo "STAR index complete!"
date

