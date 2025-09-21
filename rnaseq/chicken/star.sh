#!/bin/bash

set -euo pipefail


date

Ref_Ara='/home/johan/johan/reference/chicken'
mkdir -p "${Ref_Ara}/star_index/STAR_index_chicken_150"

STAR --runThreadN "${SLURM_CPUS_PER_TASK:-8}" \
     --runMode genomeGenerate \
     --genomeDir "${Ref_Ara}/star_index/STAR_index_chicken_150" \
     --genomeFastaFiles "${Ref_Ara}/fasta/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa" \
     --sjdbGTFfile "${Ref_Ara}/gtf/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf" \
     --sjdbOverhang 150 \
     --genomeSAindexNbases 13