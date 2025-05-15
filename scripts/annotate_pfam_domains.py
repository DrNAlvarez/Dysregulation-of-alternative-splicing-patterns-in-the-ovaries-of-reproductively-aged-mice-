#!/usr/bin/env python3
"""
annotate_pfam_domains.py
------------------------

• Parses an HMMER .domtblout file                                    → DataFrame  
• Filters to best PFAM domains per transcript (E-value threshold)    → pfam_start_stops.tsv  
• Converts protein-space PFAM coords to transcript-space (with sanity
  check against cDNA translation start/stop)                         → pfam_start_stop_with_translation.tsv  
• Adds each PFAM domain as a new feature block into a splice-event
  GTF, adjusting exon labels where domains start/stop                → *_start_stop_pfam.gtf

Example
-------
python annotate_pfam_domains.py \
    --domtbl pfam/output.domtblout \
    --translation translation_start_stops.csv \
    --gtf all.flair.collapse.isoforms_event_label.gtf \
    --out-prefix pfam/isoforms_event_label \
    --evalue 1e-5 \
    --procs 32
"""
import argparse
import multiprocessing as mp
from pathlib import Path

import pandas as pd
from tqdm import tqdm, trange


# ────────────────────────────────────────────────────────────────────────────────
#  HMMER parsing and PFAM filtering
# ────────────────────────────────────────────────────────────────────────────────
def hmmer_to_dataframe(filepath: Path) -> pd.DataFrame:
    rows = []
    with filepath.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.strip().split(maxsplit=22)
            if len(fields) < 23:
                fields.append("")
            rows.append(fields)

    cols = [
        "target_name", "target_accession", "tlen",
        "query_name", "query_accession", "qlen",
        "fs_E-value", "fs_score", "fs_bias",
        "td_#", "td_of", "td_c-Evalue", "td_i-Evalue",
        "td_score", "td_bias", "hmm_from", "hmm_to",
        "ali_from", "ali_to", "env_from", "env_to",
        "acc", "description_of_target",
    ]
    df = pd.DataFrame(rows, columns=cols)

    # Numeric columns (leave accession / names as strings)
    numeric_cols = ["tlen", "qlen", "fs_E-value", "fs_score", "fs_bias",
                    "td_#", "td_of", "td_c-Evalue", "td_i-Evalue",
                    "td_score", "td_bias", "hmm_from", "hmm_to",
                    "ali_from", "ali_to", "env_from", "env_to", "acc"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def process_hmmer_results(df: pd.DataFrame, evalue_threshold: float = 1e-5) -> pd.DataFrame:
    df = df.copy()
    df = df[df["fs_E-value"] < evalue_threshold]

    # best domain hit (max td_score) per query_name
    df_best = df.loc[df.groupby("query_name")["td_score"].idxmax()].copy()

    # split query_name "transcript_gene" (requires exactly one "_")
    df_best[["transcript_id", "gene_id"]] = df_best["query_name"].str.split("_", n=1, expand=True)

    # final tidy columns
    keep = ["query_name", "transcript_id", "gene_id", "target_name", "target_accession",
            "ali_from", "ali_to", "description_of_target"]
    df_best = df_best[keep].rename(columns={
        "query_name": "id",
        "ali_from": "start",
        "ali_to": "end",
        "description_of_target": "description"
    })
    df_best[["start", "end"]] = df_best[["start", "end"]].astype("Int64")
    return df_best


# ────────────────────────────────────────────────────────────────────────────────
#  Coordinate conversion: PFAM → transcript
# ────────────────────────────────────────────────────────────────────────────────
def convert_pfam_to_transcript_space(pfam_df: pd.DataFrame, trans_df: pd.DataFrame) -> pd.DataFrame:
    pfam = pfam_df.rename(columns={"start": "pfam_start", "end": "pfam_end"})
    trans = trans_df.rename(columns={"start": "translation_start", "end": "translation_end"})

    merged = pfam.merge(trans, on="transcript_id", how="left").dropna(subset=["translation_start"])

    merged["start"] = merged["translation_start"] + (merged["pfam_start"] - 1) * 3
    merged["end"]   = merged["translation_start"] + (merged["pfam_end"] - 1) * 3

    merged["gene_id"] = merged["gene_id_x"].combine_first(merged["gene_id_y"])
    merged["id"] = merged["transcript_id"] + "_" + merged["gene_id"]

    return merged[["id", "transcript_id", "gene_id", "target_name", "target_accession",
                   "translation_start", "translation_end", "start", "end", "description"]]


# ────────────────────────────────────────────────────────────────────────────────
#  GTF helpers (mostly unchanged from your notebook)
# ────────────────────────────────────────────────────────────────────────────────
def find_individual_overlapping_ranges(main_range, ranges):
    main_start, main_end = sorted(main_range)
    overlaps = []
    for start, end in ranges:
        ov = [n for n in range(start, end + 1) if main_start <= n <= main_end]
        if ov:
            overlaps.append([min(ov), max(ov)])
    return overlaps


def get_exon_coordinates(gtf, tid):
    exons = gtf[(gtf["transcript_id"] == tid) & (gtf["feature"].str.contains("exon", case=False))]
    strand = exons.iloc[0]["strand"] if not exons.empty else "+"
    exons = exons.sort_values("start", ascending=(strand == "+"))
    return [[row["start"], row["end"]] for _, row in exons.iterrows()]


def get_genomic_codon_coordinates_with_strand(gtf, pfam_df, tid):
    codon_row = pfam_df[pfam_df["transcript_id"] == tid]
    if codon_row.empty:
        return {"error": "no PFAM coord"}
    tx_start, tx_end = codon_row.iloc[0][["start", "end"]]

    exons = gtf[(gtf["transcript_id"] == tid) & (gtf["feature"].str.contains("exon", case=False))]
    strand = exons.iloc[0]["strand"] if not exons.empty else "+"
    exons = exons.sort_values("start", ascending=(strand == "+"))

    cum = 0
    start_genomic = stop_genomic = None
    for _, exon in exons.iterrows():
        length = exon["end"] - exon["start"] + 1
        if cum <= tx_start <= cum + length:
            start_genomic = exon["start"] + (tx_start - cum) if strand == "+" else exon["end"] - (tx_start - cum)
        if cum <= tx_end <= cum + length:
            stop_genomic  = exon["start"] + (tx_end  - cum) if strand == "+" else exon["end"]  - (tx_end  - cum)
        cum += length
        if start_genomic is not None and stop_genomic is not None:
            break

    if strand == "-" and start_genomic < stop_genomic:
        start_genomic, stop_genomic = stop_genomic, start_genomic

    return {"strand": strand,
            "start_pfam_genomic": start_genomic,
            "stop_pfam_genomic": stop_genomic}


def add_pfam_segments(gtf, pfam_df, tid):
    coords = get_genomic_codon_coordinates_with_strand(gtf, pfam_df, tid)
    if coords.get("error"):
        return pd.DataFrame()

    start_g, stop_g = coords["start_pfam_genomic"], coords["stop_pfam_genomic"]
    strand = coords["strand"]
    orf_range = [start_g, stop_g] if strand == "+" else [stop_g, start_g]

    ex_ranges = get_exon_coordinates(gtf, tid)
    overlaps = find_individual_overlapping_ranges(orf_range, ex_ranges)

    exons = gtf[(gtf["transcript_id"] == tid) & (gtf["feature"].str.contains("exon", case=False))]
    acc = pfam_df[pfam_df["transcript_id"] == tid]["target_accession"].iloc[0]

    entries = []
    for ov_start, ov_end in overlaps:
        match = exons[(exons["start"] <= ov_start) & (exons["end"] >= ov_end)]
        if match.empty:
            continue
        row = match.iloc[0].copy()
        row["feature"] = f"PFAM_{acc}_{row['feature']}"
        row["start"], row["end"] = ov_start, ov_end
        entries.append(row)

    return pd.DataFrame(entries)


# ────────────────────────────────────────────────────────────────────────────────
#  Parallel wrappers
# ────────────────────────────────────────────────────────────────────────────────
def _filter_data(tid, gtf, pfam):
    return gtf[gtf["transcript_id"] == tid], pfam[pfam["transcript_id"] == tid], tid


def _process(arg_tuple):
    gtf_subset, pfam_subset, tid = arg_tuple
    original = gtf_subset.copy()
    pfam_segs = add_pfam_segments(original, pfam_subset, tid)
    return pd.concat([original, pfam_segs], ignore_index=True).sort_values("start")


# ────────────────────────────────────────────────────────────────────────────────
#  CLI
# ────────────────────────────────────────────────────────────────────────────────
def cli():
    p = argparse.ArgumentParser(description="Annotate splice-event GTF with PFAM domains.")
    p.add_argument("--domtbl", required=True, help="HMMER domtblout file.")
    p.add_argument("--translation", required=True, help="CSV from Translation_cDNA.ipynb (tab-sep).")
    p.add_argument("--gtf", required=True, help="Input splice-event GTF (tab-sep).")
    p.add_argument("--out-prefix", default="annotated", help="Prefix for output files.")
    p.add_argument("--evalue", type=float, default=1e-5, help="E-value cutoff for PFAM filter.")
    p.add_argument("--procs", type=int, default=0, help="Worker processes (0 = all cores).")
    return p.parse_args()


def main() -> None:
    args = cli()
    domtbl_path = Path(args.domtbl).resolve()
    trans_path   = Path(args.translation).resolve()
    gtf_path     = Path(args.gtf).resolve()
    prefix       = Path(args.out_prefix).resolve()

    # ── Step 1: HMMER → PFAM start/stop (protein space) ──────────────────────
    hmmer_df = hmmer_to_dataframe(domtbl_path)
    pfam_df  = process_hmmer_results(hmmer_df, evalue_threshold=args.evalue)
    pfam_tsv = prefix.with_suffix(".pfam_start_stops.tsv")
    pfam_df.to_csv(pfam_tsv, sep="\t", index=False)
    print(f"[✓] Wrote {pfam_tsv}")

    # ── Step 2: PFAM coords → transcript space ───────────────────────────────
    trans_df = pd.read_csv(trans_path, sep="\t")
    pfam_tx  = convert_pfam_to_transcript_space(pfam_df, trans_df)
    pfam_tx_tsv = prefix.with_suffix(".pfam_start_stop_with_translation.tsv")
    pfam_tx.to_csv(pfam_tx_tsv, sep="\t", index=False)
    print(f"[✓] Wrote {pfam_tx_tsv}")

    # ── Step 3: GTF augmentation in parallel ─────────────────────────────────
    gtf_df = pd.read_csv(gtf_path, sep="\t", low_memory=False)
    tids = pfam_tx["transcript_id"].unique()
    n_procs = args.procs or min(mp.cpu_count(), len(tids))

    arg_list = [_filter_data(tid, gtf_df, pfam_tx) for tid in tids]

    with mp.Pool(processes=n_procs) as pool:
        results = list(tqdm(pool.imap(_process, arg_list),
                            total=len(arg_list), desc="Annotating GTF"))

    combined = pd.concat(results, ignore_index=True)
    out_gtf = prefix.with_suffix("_start_stop_pfam.gtf")
    combined.to_csv(out_gtf, sep="\t", index=False)
    print(f"[✓] Final GTF with PFAM features → {out_gtf}")


if __name__ == "__main__":
    main()

