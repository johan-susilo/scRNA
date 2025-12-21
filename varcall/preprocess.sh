date

set -euo pipefail

mkdir -p ~/output/autism
cd ~/output/autism


# Use prefetch (sometimes works when fasterq-dump fails)
prefetch SRR8697636 -O ./
prefetch SRR8697645 -O ./
prefetch SRR8697646 -O ./

# Convert to FASTQ
fastq-dump --split-files SRR8697636.sra
fastq-dump --split-files SRR8697645.sra
fastq-dump --split-files SRR8697646.sra


# Proband (autistic child with IQ=60, ADI=42)
fasterq-dump --split-files SRR8697637 -O ./

# Mother (unaffected)
fasterq-dump --split-files 	SRR8697627 -O ./

# Father (unaffected - your file!)
fasterq-dump --split-files 	SRR8697645 -O ./

# Rename for clarity
mv SRR8697637_*.fastq autism_*.fastq
mv SRR8697627_*.fastq mother_*.fastq
mv SRR8697645_*.fastq father_*.fastq

# Compress to save space
gzip *.fastq