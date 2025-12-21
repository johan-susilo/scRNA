#!/bin/bash
# Quick convert for already-downloaded SRA file

set -euo pipefail

cd ~/output/autism

echo "Converting existing SRR8697636.sra to FASTQ..."

if [ -f "SRR8697636/SRR8697636.sra" ]; then
    echo "✓ Found SRR8697636.sra (2.4 GB)"

    fasterq-dump \
        --split-files \
        --threads 60 \
        --progress \
        --outdir . \
        SRR8697636/SRR8697636.sra

    # Rename
    mv SRR8697636_1.fastq proband_1.fastq
    mv SRR8697636_2.fastq proband_2.fastq

    # Compress
    echo "Compressing with pigz..."
    pigz -p 60 proband_1.fastq &
    pigz -p 60 proband_2.fastq &
    wait

    echo "✓ Proband FASTQ ready!"
    ls -lh proband_*.fastq.gz
else
    echo "❌ SRR8697636.sra not found"
fi
