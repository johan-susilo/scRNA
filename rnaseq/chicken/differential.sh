#!/bin/bash
# Program:
#   Differential expression (RNA) for chicken data with safe preprocessing.
#   - Converts legacy metadata (SampleID/Type/etc) to the schema required by differential.R
#   - Subsets counts to metadata samples (avoids DESeq2 mismatches)
#   - Compares high_egg vs low_egg while adjusting for tissue (uterus/ovary)
#
# Usage:
#   ./differential.sh RNA
#
#SBATCH --job-name=differential_chicken
#SBATCH --output=./log/DGE/%x.log
#SBATCH -n 8
#SBATCH --mem=64GB

set -euo pipefail

# ---------------- Args / Mode ----------------
if [[ -z "${1:-}" ]]; then
  echo "Error: Analysis mode not specified."
  echo "Usage: $0 RNA"
  exit 1
fi

ANALYSIS_MODE="$1"
case "$ANALYSIS_MODE" in
  RNA) OUTPUT_PREFIX="RNA_DGE_analysis" ;;
  RPF) echo "RPF mode not supported here (needs RPF metadata + counts)."; exit 1 ;;
  TE)  echo "TE mode not supported here (needs RNA+RPF + TE design)."; exit 1 ;;
  *)   echo "Error: Invalid mode '$ANALYSIS_MODE' (use RNA)."; exit 1 ;;
esac

# ---------------- Paths ----------------
read_count="/home/johan/johan/output/chicken/combined.HTseq_report"
metadata_in="/home/johan/johan/output/chicken/metadata.csv"          # legacy or modern
r_script_path="./differential.R"
config_file="./dge_color.yaml"

# Styling (optional)
FONT_FAMILY="Source Sans 3"
FONT_SIZE=24
FONT_DIR="/homes/johansusiloe/wd/References/font/static"
TRANSPARENT_PNG=true

# Output
output_dir="/home/johan/johan/output/chicken/${OUTPUT_PREFIX}"
mkdir -p "$output_dir" ./log/DGE

echo "====================================================================="
date
echo "Starting differential analysis."
echo "Mode: $ANALYSIS_MODE"
echo "Counts: $read_count"
echo "Metadata (raw): $metadata_in"
echo "Output: $output_dir"
echo

# ---------------- Sanity ----------------
[[ -f "$read_count" ]]   || { echo "Counts not found: $read_count" >&2; exit 1; }
[[ -f "$metadata_in" ]]  || { echo "Metadata not found: $metadata_in" >&2; exit 1; }

# ---------------- 1) Normalize/convert metadata to required schema ----------------
# Required columns for differential.R: sample,subject,tissue,replicate,egg_production,assay_type
META_NORM="${output_dir}/metadata.norm.csv"
python3 - "$metadata_in" "$META_NORM" <<'PY'
import re, sys, pandas as pd
src, dst = sys.argv[1], sys.argv[2]
m = pd.read_csv(src)

def modern(df):
    need = ['sample','subject','tissue','replicate','egg_production','assay_type']
    if all(c in df.columns for c in need):
        out = df[need].copy()
        out['assay_type'] = out['assay_type'].astype(str).str.upper()
        return out
    return None

def legacy(df):
    # Expect column SampleID with IDs like RH1_ovary_rep1
    if 'SampleID' not in df.columns: return None
    rows = []
    for sid in df['SampleID'].astype(str):
        mo = re.match(r'^(RH)(\d+)_(ovary|uterus)_rep(\d+)$', sid, flags=re.I)
        if mo:
            rh, num, tissue, rep = mo.group(1).upper(), int(mo.group(2)), mo.group(3).lower(), int(mo.group(4))
            egg = 'high_egg' if 1 <= num <= 5 else 'low_egg'
            rows.append(dict(sample=sid, subject=f"{rh}{num}", tissue=tissue,
                             replicate=rep, egg_production=egg, assay_type='RNA'))
        else:
            rows.append(dict(sample=sid, subject='NA', tissue='NA',
                             replicate='NA', egg_production='low_egg', assay_type='RNA'))
    return pd.DataFrame(rows, columns=['sample','subject','tissue','replicate','egg_production','assay_type'])

out = modern(m)
if out is None: out = legacy(m)
if out is None: raise SystemExit("Unrecognized metadata format.")

try: out['replicate'] = out['replicate'].astype(int)
except Exception: pass

out.to_csv(dst, index=False)
print(f"Wrote normalized metadata: {dst} (n={len(out)})")
PY

# ---------------- 2) Keep only RNA rows (future-proof) ----------------
META_RNA="${output_dir}/metadata.RNA.csv"
python3 - "$META_NORM" "$META_RNA" <<'PY'
import sys, pandas as pd
src, dst = sys.argv[1], sys.argv[2]
m = pd.read_csv(src)
m = m[m['assay_type'].astype(str).str.upper()=='RNA'].copy()
if m.empty: raise SystemExit("No RNA rows after filtering.")
m.to_csv(dst, index=False)
print(f"Wrote RNA-only metadata: {dst} (n={len(m)})")
PY

# ---------------- 3) Subset counts to metadata samples & sync metadata ----------------
COUNT_USED="${output_dir}/counts.used.tsv"
META_USED="${output_dir}/metadata.used.csv"
python3 - "$read_count" "$META_RNA" "$COUNT_USED" "$META_USED" <<'PY'
import sys, pandas as pd
cnt_fp, meta_fp, out_cnt, out_meta = sys.argv[1:]
cnt  = pd.read_csv(cnt_fp, sep="\t")
meta = pd.read_csv(meta_fp)

if "GeneID" not in cnt.columns:
    cnt = cnt.rename(columns={cnt.columns[0]: "GeneID"})

sample_set = set(meta["sample"].astype(str))
cols = ["GeneID"] + [c for c in cnt.columns if c != "GeneID" and c in sample_set]
sub = cnt.loc[:, cols].copy()
if sub.shape[1] <= 1: raise SystemExit("No matching sample columns between counts and metadata.")

meta2 = meta[meta["sample"].isin([c for c in sub.columns if c != "GeneID"])].copy()
meta2 = meta2.set_index("sample").loc[[c for c in sub.columns if c != "GeneID"]].reset_index()

sub.to_csv(out_cnt, sep="\t", index=False)
meta2.to_csv(out_meta, index=False)
print(f"Wrote counts: {out_cnt} (genes={len(sub)} cols={sub.shape[1]})")
print(f"Wrote meta  : {out_meta} (n={len(meta2)})")
PY

read_count_used="$COUNT_USED"
metadata_used="$META_USED"

# ---------------- 4) Contrast & design ----------------
# Compare high_egg vs low_egg, adjust for tissue (uterus/ovary)
CONTRAST_FACTOR="egg_production"
NUMERATOR="high_egg"
DENOMINATOR="low_egg"
BATCH_COLUMN="tissue"   # design ~ tissue + egg_production

# To flip to uterus vs ovary (adjust for egg group), switch to:
# CONTRAST_FACTOR="tissue"; NUMERATOR="uterus"; DENOMINATOR="ovary"; BATCH_COLUMN="egg_production"

# ---------------- 5) Detect keyType ----------------
first_id="$(awk -F'\t' 'NR==2{print $1}' "$read_count_used" 2>/dev/null || echo "")"
KEY_TYPE="ENSEMBL"
if [[ "$first_id" =~ ^[0-9]+$ ]]; then KEY_TYPE="ENTREZID"; fi
ORG_DB="org.Gg.eg.db"

# ---------------- 6) Run R ----------------
echo "#######################################"
echo "# Running differential.R              #"
echo "#######################################"
echo "Design: ~ ${BATCH_COLUMN:+${BATCH_COLUMN} + }${CONTRAST_FACTOR}"
echo "Contrast: ${CONTRAST_FACTOR} : ${NUMERATOR} vs ${DENOMINATOR}"
echo "OrgDb/keyType: $ORG_DB / $KEY_TYPE"
echo

CMD=( Rscript "$r_script_path"
  --counts "$read_count_used"
  --metadata "$metadata_used"
  --outdir "$output_dir"
  --sampleColumn "sample"
  --assayTypeColumn "assay_type"
  --assayTypeToSelect "RNA"
  --contrastFactor "$CONTRAST_FACTOR"
  --contrastNumerator "$NUMERATOR"
  --contrastDenominator "$DENOMINATOR"
  --orgDb "$ORG_DB"
  --keyType "$KEY_TYPE"
  --cores "${SLURM_NTASKS:-4}"
  --padj 0.05
  --lfc 1
  --topNLabels 10
  --topNHeatmap 50
  --filterMethod "intersection"
  --configFile "$config_file"
  --transparentPNG
  --outputFormats "png,pdf"
  --fontFamily "$FONT_FAMILY"
  --fontSize "$FONT_SIZE"
  --fontDir "$FONT_DIR"
)

# Add batch column if set
if [[ -n "$BATCH_COLUMN" ]]; then
  CMD+=( --batchColumn "$BATCH_COLUMN" )
fi

echo "Executing:"
printf ' %q' "${CMD[@]}"; echo
"${CMD[@]}"

echo
echo "#############################"
echo "# Analysis complete         #"
echo "#############################"
echo "Results: $output_dir"
date
echo "====================================================================="
