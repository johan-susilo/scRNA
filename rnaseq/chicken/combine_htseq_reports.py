#!/usr/bin/env python3
"""
Combine per-sample HTSeq-count outputs into a single combined matrix and
generate metadata files.

Detection modes:
  - auto  : inspect folders and pick 'count' if count/*.htseq.*.counts.tsv exist,
            otherwise pick 'star' if STAR_out/*.HTseq_report exist.
  - count : expect <data_folder>/<sample>/count/*.htseq.*.counts.tsv
  - star  : expect <data_folder>/<sample>/STAR_out/*.HTseq_report
"""

# python3 combine_htseq_reports.py --data_folder /home/johan/johan/output/chicken --mode count
import os
import sys
import re
import glob
import argparse
import pandas as pd

def parse_htseq_file(file_path, sample_name):
    # Read 2-column HTSeq-count-like file; skip lines starting with '__'
    rows = []
    try:
        with open(file_path, 'r') as fh:
            for ln, line in enumerate(fh, 1):
                line = line.rstrip("\n")
                if not line or line.startswith("__"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                gene, cnt = parts[0], parts[1]
                try:
                    rows.append((gene, int(cnt)))
                except ValueError:
                    # skip non-integer counts
                    continue
    except Exception as e:
        print(f"Error reading {file_path}: {e}", file=sys.stderr)
        return pd.DataFrame(columns=["GeneID", sample_name])
    if not rows:
        return pd.DataFrame(columns=["GeneID", sample_name])
    df = pd.DataFrame(rows, columns=["GeneID", sample_name])
    return df

def _parse_sample_to_meta(sample_label):
    """
    Expect labels like: RH303_ovary_rep1  or RH2_uterus_rep1
    Returns dict for both new (lowercase) and legacy (caps) metadata rows.
    """
    m = re.match(r'^(RH)(\d+)_(ovary|uterus)_rep(\d+)$', str(sample_label), flags=re.IGNORECASE)
    if m:
        rh_prefix = m.group(1).upper()
        rh_num    = int(m.group(2))
        tissue    = m.group(3).lower()
        rep       = int(m.group(4))
        egg = "high" if 1 <= rh_num <= 5 else "low"
        subj = f"{rh_prefix}{rh_num}"
        # New (lowercase) schema used by DESeq driver scripts
        row_new = {
            "sample": sample_label,
            "subject": subj,
            "tissue": tissue,
            "replicate": rep,
            "egg_production": egg,
            "assay_type": "RNA"
        }
        # Legacy (capitalized) schema like your pasted CSV
        row_caps = {
            "SampleID": sample_label,
            "Treatment": egg,
            "Type": "rna",
            "Replicate": rep
        }
    else:
        # Fallback for unexpected labels
        row_new = {
            "sample": sample_label,
            "subject": "NA",
            "tissue": "NA",
            "replicate": "NA",
            "egg_production": "low",
            "assay_type": "RNA"
        }
        row_caps = {
            "SampleID": sample_label,
            "Treatment": "low",
            "Type": "rna",
            "Replicate": "NA"
        }
    return row_new, row_caps

def generate_metadata_file(sample_headers, output_dir):
    """
    Emit two metadata files:
      1) metadata.csv             : sample,subject,tissue,replicate,egg_production,assay_type
      2) metadata_legacy_caps.csv : SampleID,Treatment,Type,Replicate
    """
    rows_new, rows_caps = [], []
    for s in sample_headers:
        rnew, rcaps = _parse_sample_to_meta(s)
        rows_new.append(rnew)
        rows_caps.append(rcaps)

    # Order columns explicitly for cleanliness
    df = pd.DataFrame(rows_new, columns=[
        "sample", "subject", "tissue", "replicate", "egg_production", "assay_type"
    ])
    df_caps = pd.DataFrame(rows_caps, columns=[
        "SampleID", "Treatment", "Type", "Replicate"
    ])

    out_main = os.path.join(output_dir, "metadata.csv")
    out_caps = os.path.join(output_dir, "metadata_legacy_caps.csv")
    df.to_csv(out_main, index=False)
    df_caps.to_csv(out_caps, index=False)
    print(f"Wrote metadata: {out_main}")
    print(f"Wrote legacy metadata: {out_caps}")

def combine_dataframes(dfs):
    # dfs: list of DataFrames each with columns GeneID + sample_col
    if not dfs:
        return pd.DataFrame()
    combined = dfs[0]
    for df in dfs[1:]:
        combined = pd.merge(combined, df, on="GeneID", how="outer")
    combined = combined.fillna(0)
    # ensure integers for sample columns
    for c in combined.columns:
        if c == "GeneID":
            continue
        try:
            combined[c] = combined[c].astype(int)
        except Exception:
            combined[c] = pd.to_numeric(combined[c], errors="coerce").fillna(0).astype(int)
    return combined

def discover_files(data_folder, mode):
    subdirs = [d for d in sorted(os.listdir(data_folder)) if os.path.isdir(os.path.join(data_folder, d))]
    found = []
    for sub in subdirs:
        base = os.path.join(data_folder, sub)
        if mode == "count":
            # Prefer the sample folder name as label to match metadata expectations
            files = glob.glob(os.path.join(base, "count", "*.htseq.*.counts.tsv"))
            if not files:
                files = glob.glob(os.path.join(base, "count", "*.htseq.*.tsv"))

            # Select only ONE file per sample based on priority order
            # Priority: reverse > no > yes > (any other)
            selected_file = None
            if files:
                # Sort by priority: prefer 'reverse' strandedness
                priority_order = ['reverse', 'no', 'yes']
                for priority in priority_order:
                    matching = [f for f in files if f'.htseq.{priority}.' in f]
                    if matching:
                        selected_file = sorted(matching)[0]  # Take first if multiple matches
                        break
                # If no priority match, just take the first file
                if selected_file is None:
                    selected_file = sorted(files)[0]

                lbl = sub  # <--- crucial: use folder (sample) name
                found.append((selected_file, lbl))
        elif mode == "star":
            star_dir = os.path.join(base, "STAR_out")
            expected = os.path.join(star_dir, f"{sub}.HTseq_report")
            if os.path.isfile(expected):
                found.append((expected, sub))
            else:
                alt = glob.glob(os.path.join(star_dir, "*.HTseq_report"))
                for f in sorted(alt):
                    lbl = os.path.basename(f).split(".HTseq_report")[0]
                    found.append((f, lbl))
        else:
            raise ValueError("unknown mode: " + str(mode))
    return found

def auto_detect_mode(data_folder):
    count_pat = os.path.join(data_folder, "*", "count", "*.htseq.*.counts.tsv")
    star_pat  = os.path.join(data_folder, "*", "STAR_out", "*.HTseq_report")
    if glob.glob(count_pat):
        return "count"
    if glob.glob(star_pat):
        return "star"
    return None

def main(args):
    data_folder = args.data_folder
    mode = args.mode
    if mode == "auto":
        detected = auto_detect_mode(data_folder)
        if detected is None:
            print("No recognizable HTSeq outputs found under", data_folder, file=sys.stderr)
            sys.exit(1)
        print("Auto-detected mode:", detected)
        mode = detected

    files = discover_files(data_folder, mode)
    if not files:
        print(f"No HTSeq files found (mode={mode}) under {data_folder}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(files)} files, processing...")
    dfs = []
    samples = []
    for fp, label in files:
        print(f"  -> {fp} (sample: {label})")
        df = parse_htseq_file(fp, label)
        # ensure sample column exists in df even if empty
        if df.shape[1] == 1:
            df[label] = pd.Series(dtype=int)
        dfs.append(df)
        samples.append(label)

    combined = combine_dataframes(dfs)
    if combined.empty:
        print("Combined data is empty, nothing to write.", file=sys.stderr)
        sys.exit(1)

    out_comb = os.path.join(data_folder, "combined.HTseq_report")
    combined.to_csv(out_comb, sep="\t", index=False)
    print("Wrote combined HTSeq report to:", out_comb)

    # write metadata from sample headers
    sample_cols = [c for c in combined.columns if c != "GeneID"]
    if sample_cols:
        generate_metadata_file(sample_cols, data_folder)
    else:
        print("No sample columns found to generate metadata.", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine HTSeq-count files (star/count layout).")
    parser.add_argument("--data_folder", required=True, help="Base folder with per-sample subdirectories")
    parser.add_argument("--mode", choices=["auto","star","count"], default="auto",
                        help="Layout detection mode; default 'auto'")
    args = parser.parse_args()
    if not os.path.isdir(args.data_folder):
        print("data_folder not found:", args.data_folder, file=sys.stderr)
        sys.exit(1)
    main(args)
