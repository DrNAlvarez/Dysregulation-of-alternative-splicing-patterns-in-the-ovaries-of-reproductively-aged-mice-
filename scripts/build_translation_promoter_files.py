#!/usr/bin/env python3
"""
build_translation_promoter_files.py
-----------------------------------
Create a promoter–translation merge table and a FASTA file of AA sequences.

Example
-------
python build_translation_promoter_files.py \
    --translation translation_start_stops.csv \
    --gtf all.flair.collapse.isoforms_event_label_start_stop_orf.gtf \
    --csv   translation_data_by_promoter_group.csv \
    --fasta pfam/translation_data_by_promoter_group.fa
"""
import argparse
from pathlib import Path

import pandas as pd


def make_csv(gtf_path: Path, translation_path: Path, out_csv: Path) -> pd.DataFrame:
    """Merge promoter-group info with translation table and save to CSV."""
    gtf_df = pd.read_csv(gtf_path, sep="\t")
    translation_df = pd.read_csv(translation_path, sep="\t")

    promoter_df = (
        gtf_df[["gene_id", "transcript_id", "promoter_group"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    merged = promoter_df.merge(
        translation_df, on=["gene_id", "transcript_id"], how="inner"
    )
    merged.to_csv(out_csv, sep="\t", index=False)
    return merged


def make_fasta(df: pd.DataFrame, out_fa: Path) -> None:
    """Write AA sequences to FASTA with header transcriptID_geneID."""
    out_fa.parent.mkdir(parents=True, exist_ok=True)

    with out_fa.open("w") as handle:
        for _, row in df.iterrows():
            header = f">{row['transcript_id']}_{row['gene_id']}"
            seq = row["AA sequence"]
            handle.write(f"{header}\n{seq}\n")


def cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Merge translation start/stop data with promoter groups and emit CSV + FASTA."
    )
    p.add_argument(
        "--translation",
        required=True,
        help="Tab-separated file produced by Translation_cDNA.ipynb (translation_start_stops.csv)",
    )
    p.add_argument(
        "--gtf",
        required=True,
        help="Event-label GTF file with promoter_group column.",
    )
    p.add_argument(
        "--csv",
        default="translation_data_by_promoter_group.csv",
        help="Output CSV filename (tab-separated).",
    )
    p.add_argument(
        "--fasta",
        default="translation_data_by_promoter_group.fa",
        help="Output FASTA filename.",
    )
    return p.parse_args()


def main() -> None:
    args = cli()
    csv_path = Path(args.csv)
    fasta_path = Path(args.fasta)

    merged_df = make_csv(Path(args.gtf), Path(args.translation), csv_path)
    make_fasta(merged_df, fasta_path)

    print(f"Wrote:\n  • {csv_path.resolve()}\n  • {fasta_path.resolve()}")


if __name__ == "__main__":
    main()

