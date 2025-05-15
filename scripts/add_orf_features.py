#!/usr/bin/env python3
"""
add_orf_features.py
-------------------
Merge GTF exon data with start/stop-codon coordinates and build segmented ORF,
start_codon, stop_codon, and annotated exon rows for every transcript.

Example
-------
python add_orf_features.py \
       --gtf all.flair.collapse.isoforms_event_label.gtf \
       --codon translation_start_stops.csv \
       --out all.flair.collapse.isoforms_event_label_start_stop_orf.gtf \
       --procs 96
"""
import argparse
import multiprocessing as mp
from pathlib import Path

import pandas as pd
from tqdm import tqdm


# ---------------------------------------------------------------------
#  Feature-building functions (unchanged from the notebook)
# ---------------------------------------------------------------------
def add_codon_features_to_gtf_with_segmented_orf(df_gtf, df_codon_positions, transcript_id):
    codon_coordinates = get_genomic_codon_coordinates_with_strand(
        df_gtf, df_codon_positions, transcript_id
    )
    start_codon_genomic = codon_coordinates.get("start_codon_genomic")
    stop_codon_genomic = codon_coordinates.get("stop_codon_genomic")
    if start_codon_genomic is None or stop_codon_genomic is None:
        return pd.DataFrame()  # missing codon info

    strand = codon_coordinates["strand"]
    orf_range = (
        [start_codon_genomic, stop_codon_genomic]
        if strand == "+"
        else [stop_codon_genomic, start_codon_genomic]
    )

    exon_coordinates = get_exon_coordinates(df_gtf, transcript_id)
    overlapping_exon_ranges = find_individual_overlapping_ranges(orf_range, exon_coordinates)

    transcript_exons = df_gtf[
        (df_gtf["transcript_id"] == transcript_id) & (df_gtf["feature"].str.contains("exon", case=False))
    ]
    orf_entries = []
    for overlap_start, overlap_end in overlapping_exon_ranges:
        match = transcript_exons[
            (transcript_exons["start"] <= overlap_start) & (transcript_exons["end"] >= overlap_end)
        ]
        if not match.empty:
            exon_info = match.iloc[0].copy()
            exon_info["feature"] = f"orf_{exon_info['feature']}"
            exon_info["start"] = overlap_start
            exon_info["end"] = overlap_end
            orf_entries.append(exon_info)

    return pd.DataFrame(orf_entries).reset_index(drop=True)


def add_codon_features_with_updated_exons(df_gtf, df_codon_positions, transcript_id):
    orf_df = add_codon_features_to_gtf_with_segmented_orf(df_gtf, df_codon_positions, transcript_id)
    if "error" in orf_df:
        return orf_df

    codon_coordinates = get_genomic_codon_coordinates_with_strand(
        df_gtf, df_codon_positions, transcript_id
    )
    start_codon_genomic = codon_coordinates["start_codon_genomic"]
    stop_codon_genomic = codon_coordinates["stop_codon_genomic"]

    filtered_gtf = df_gtf[df_gtf["transcript_id"] == transcript_id].copy()
    last_exon_number = filtered_gtf["exon_number"].max()

    # add start_codon / stop_codon rows
    codon_rows = []
    for codon_type, codon_position in [
        ("start_codon", start_codon_genomic),
        ("stop_codon", stop_codon_genomic),
    ]:
        row = filtered_gtf.iloc[0].copy()
        row["feature"] = codon_type
        row["start"] = row["end"] = codon_position
        codon_rows.append(row)

    # annotate exon features
    for idx, exon in filtered_gtf[filtered_gtf["feature"].str.contains("exon", case=False)].iterrows():
        if exon["start"] <= start_codon_genomic <= exon["end"]:
            filtered_gtf.at[idx, "feature"] = f"{exon['feature']}_start"
        if exon["start"] <= stop_codon_genomic <= exon["end"]:
            filtered_gtf.at[idx, "feature"] = f"{exon['feature']}_stop"
        if exon["exon_number"] == last_exon_number:
            filtered_gtf.at[idx, "feature"] = f"{exon['feature']}_last"

    # annotate ORF segments
    for idx, exon in orf_df[orf_df["feature"].str.contains("exon", case=False)].iterrows():
        if exon["start"] <= start_codon_genomic <= exon["end"]:
            orf_df.at[idx, "feature"] = f"{exon['feature']}_start"
        if exon["start"] <= stop_codon_genomic <= exon["end"]:
            orf_df.at[idx, "feature"] = f"{exon['feature']}_stop"
    for idx, exon in orf_df[orf_df["feature"].str.contains("orf", case=False)].iterrows():
        if exon["exon_number"] == last_exon_number:
            orf_df.at[idx, "feature"] = f"{exon['feature']}_last"

    updated_df = pd.concat([filtered_gtf, pd.DataFrame(codon_rows), orf_df], ignore_index=True)
    return updated_df.sort_values("start").reset_index(drop=True)


def get_genomic_codon_coordinates_with_strand(df_gtf, df_codon_positions, transcript_id):
    codon_data = df_codon_positions[df_codon_positions["transcript_id"] == transcript_id]
    if codon_data.empty:
        return {"error": "No codon data found for this transcript ID."}

    transcript_start_codon = codon_data["start"].values[0]
    transcript_stop_codon = codon_data["end"].values[0]

    transcript_exons = df_gtf[
        (df_gtf["transcript_id"] == transcript_id) & (df_gtf["feature"].str.contains("exon", case=False))
    ]
    strand = transcript_exons.iloc[0]["strand"] if not transcript_exons.empty else "+"
    transcript_exons = transcript_exons.sort_values("start", ascending=(strand == "+"))

    cumulative = 0
    start_codon_genomic = stop_codon_genomic = None
    for _, exon in transcript_exons.iterrows():
        exon_len = exon["end"] - exon["start"] + 1
        # start
        if cumulative <= transcript_start_codon <= cumulative + exon_len:
            start_codon_genomic = (
                exon["start"] + (transcript_start_codon - cumulative)
                if strand == "+"
                else exon["end"] - (transcript_start_codon - cumulative)
            )
        # stop
        if cumulative <= transcript_stop_codon <= cumulative + exon_len:
            stop_codon_genomic = (
                exon["start"] + (transcript_stop_codon - cumulative)
                if strand == "+"
                else exon["end"] - (transcript_stop_codon - cumulative)
            )
        cumulative += exon_len
        if start_codon_genomic is not None and stop_codon_genomic is not None:
            break

    if strand == "-" and start_codon_genomic < stop_codon_genomic:
        start_codon_genomic, stop_codon_genomic = stop_codon_genomic, start_codon_genomic

    return {
        "strand": strand,
        "start_codon_genomic": start_codon_genomic,
        "stop_codon_genomic": stop_codon_genomic,
    }


def get_exon_coordinates(df_gtf, transcript_id):
    exons = df_gtf[
        (df_gtf["transcript_id"] == transcript_id) & (df_gtf["feature"].str.contains("exon", case=False))
    ]
    strand = exons.iloc[0]["strand"] if not exons.empty else "+"
    exons = exons.sort_values("start", ascending=(strand == "+"))
    return [[row["start"], row["end"]] for _, row in exons.iterrows()]


def find_individual_overlapping_ranges(main_range, ranges):
    main_start, main_end = sorted(main_range)
    overlaps = []
    for start, end in ranges:
        ov = [n for n in range(start, end + 1) if main_start <= n <= main_end]
        if ov:
            overlaps.append([min(ov), max(ov)])
    return overlaps


# ---------------------------------------------------------------------
#  Parallel helpers
# ---------------------------------------------------------------------
def _prepare_arg(tid, df_gtf, df_trans):
    return df_gtf[df_gtf["transcript_id"] == tid], df_trans[df_trans["transcript_id"] == tid], tid


def _process(arg_tuple):
    gtf_subset, codon_subset, tid = arg_tuple
    return add_codon_features_with_updated_exons(gtf_subset, codon_subset, tid)


# ---------------------------------------------------------------------
#  CLI entry point
# ---------------------------------------------------------------------
def main() -> None:
    p = argparse.ArgumentParser(description="Add segmented ORF + start/stop features to a GTF.")
    p.add_argument("--gtf", required=True, help="Input GTF file (tab-separated).")
    p.add_argument("--codon", required=True, help="CSV with translation start/stop coordinates.")
    p.add_argument("--out", default="output_with_orf.gtf", help="Output GTF filename.")
    p.add_argument(
        "--procs",
        type=int,
        default=0,
        help="Number of worker processes (0 = use all available cores).",
    )
    args = p.parse_args()

    df_gtf = pd.read_csv(args.gtf, sep="\t")
    df_trans = pd.read_csv(args.codon, sep="\t")

    trans_ids = df_trans["transcript_id"].unique()
    n_workers = args.procs or min(mp.cpu_count(), len(trans_ids))

    # build argument tuples (done in serial to avoid sending huge DataFrames repeatedly)
    arg_list = [_prepare_arg(tid, df_gtf, df_trans) for tid in tqdm(trans_ids, desc="Preparing")]

    with mp.Pool(processes=n_workers) as pool:
        results = list(
            tqdm(pool.imap(_process, arg_list), total=len(arg_list), desc="Processing transcripts")
        )

    combined = pd.concat(results, ignore_index=True)
    combined.to_csv(args.out, sep="\t", index=False)
    print(f"Saved {len(combined):,} rows → {args.out}")


if __name__ == "__main__":
    main()

