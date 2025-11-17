from __future__ import annotations

import gzip
import json
import math
import os
import re
import shutil
import subprocess
import tarfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import (
    Callable,
    Dict,
    Iterable,
    Iterator,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
    TYPE_CHECKING,
)

import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
import seaborn as sns
try:
    from scipy import sparse as sp
except ImportError:  # pragma: no cover - optional dependency
    class _SparseShim:
        @staticmethod
        def issparse(matrix) -> bool:
            return False

    sp = _SparseShim()  # type: ignore[assignment]

try:  # pragma: no cover - optional dependency
    from tqdm.auto import tqdm
except ImportError:  # pragma: no cover
    tqdm = None

if TYPE_CHECKING:  # pragma: no cover
    from anndata import AnnData

DEFAULT_GTF_COLUMNS = (
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "strand",
    "gene_id",
    "transcript_id",
    "promoter_group",
)
EVENT_PREFIX = "exon_"
MINIMAP_VERSION = "2.30"
MINIMAP_ARCHIVE = f"minimap2-{MINIMAP_VERSION}_x64-linux"
MINIMAP_URL = (
    f"https://github.com/lh3/minimap2/releases/download/v{MINIMAP_VERSION}/"
    f"{MINIMAP_ARCHIVE}.tar.bz2"
)

_ANNDATA_NULL_PATCHED = False


def _progress_iter(
    iterable: Iterable,
    description: str,
    total: Optional[int] = None,
    item_formatter: Optional[Callable[[object], str]] = None,
):
    """
    Wrap an iterable with tqdm if available; otherwise emit textual progress updates.
    """
    if tqdm is not None:
        return tqdm(iterable, desc=description, total=total, leave=False)

    def _generator():
        print(f"{description} (tqdm unavailable; streaming progress below)")
        known_total = total
        if known_total is None:
            try:
                known_total = len(iterable)  # type: ignore[arg-type]
            except TypeError:
                known_total = "?"
        for idx, item in enumerate(iterable, start=1):
            label = item_formatter(item) if item_formatter else ""
            suffix = f" {label}" if label else ""
            print(f"  [{idx}/{known_total}] {suffix}".rstrip(), flush=True)
            yield item

    return _generator()


@dataclass
class SupFigure4Config:
    pfam_paths: Mapping[str, Path]
    splice_paths: Sequence[Path]
    gtf_path: Path
    fasta_path: Path
    figure_dir: Path
    window_nt: int = 300
    tail_windows: Optional[Mapping[str, Tuple[int, int]]] = None
    summary_filename: Optional[str] = None
    fasta_filename: Optional[str] = None
    plot_filename: Optional[str] = None
    make_plot: bool = True
    skip_existing_outputs: bool = True

    def __post_init__(self) -> None:
        self.pfam_paths = {str(label): Path(path) for label, path in self.pfam_paths.items()}
        self.splice_paths = [Path(path) for path in self.splice_paths]
        self.gtf_path = Path(self.gtf_path)
        self.fasta_path = Path(self.fasta_path)
        self.figure_dir = Path(self.figure_dir)
        if self.tail_windows is None:
            self.tail_windows = dict(TAIL_WINDOWS)
        else:
            self.tail_windows = {
                str(label): (int(bounds[0]), int(bounds[1]))
                for label, bounds in self.tail_windows.items()
            }


@dataclass
class SupFigure4Outputs:
    transcript_universe: pd.Index
    pfam_hits: pd.DataFrame
    events: pd.DataFrame
    summary: pd.DataFrame
    tail_uniqueness: pd.DataFrame
    final_summary: pd.DataFrame
    unique_hits: pd.Index
    expanded_transcripts: pd.Index
    summary_path: Path
    fasta_output: Optional[Path]
    plot_path: Optional[Path]


@dataclass
class MinimapAlignmentConfig:
    reference_fasta: Path
    fastq_map: Mapping[str, Iterable[Path]]
    output_dir: Path
    minimap_bin: Optional[Path] = None
    threads: int = 8
    preset: str = "sr"
    disable_secondary: bool = True
    include_cigar: bool = True
    compress_output: bool = True
    auto_install: bool = True
    install_dir: Optional[Path] = None
    skip_existing_outputs: bool = True
    output_stems: Optional[Mapping[str, str]] = None

    def __post_init__(self) -> None:
        self.reference_fasta = Path(self.reference_fasta)
        self.output_dir = Path(self.output_dir)
        if self.minimap_bin is not None:
            self.minimap_bin = Path(self.minimap_bin)
        if self.install_dir is not None:
            self.install_dir = Path(self.install_dir)
        normalized: Dict[str, List[Path]] = {}
        for label, paths in self.fastq_map.items():
            label_str = str(label)
            if isinstance(paths, (str, os.PathLike, Path)):
                iterable = [paths]  # type: ignore[list-item]
            else:
                iterable = list(paths)
            resolved = [Path(p) for p in iterable]
            if not resolved:
                raise ValueError(f"No FASTQ paths supplied for label '{label_str}'.")
            normalized[label_str] = resolved
        self.fastq_map = normalized
        if self.output_stems is None:
            self.output_stems = {}
        else:
            self.output_stems = {str(k): str(v) for k, v in self.output_stems.items()}


@dataclass
class MinimapAlignmentOutputs:
    minimap_bin: Path
    paf_paths: Dict[str, Path]
    paf_gz_paths: Dict[str, Optional[Path]]


@dataclass
class ScanpyAnalysisConfig:
    kallisto_dir: Path
    data_dir: Path
    figure_dir: Path
    expected_cells: Mapping[str, int]
    condition_map: Mapping[str, str]
    isoform_match_paths: Mapping[str, Path]
    isoform_transcripts: Sequence[str]
    isoform_ratio_pair: Optional[Tuple[str, str]] = None
    max_hvg_genes: int = 3000
    harmony_key: str = "Sample"
    leiden_resolution: float = 0.7
    isoform_max_distance: int = 2
    write_barcode_tables: bool = True
    barcode_table_name: str = "barcodes_per_celltype.csv"
    barcode_subdir_name: str = "barcodes_by_sample"
    combined_counts_filename: str = "combined_isoform_by_celltype_counts.csv"
    write_isoform_tables: bool = True
    marker_genes: Optional[Mapping[str, Sequence[str]]] = None
    marker_gene_path: Optional[Path] = None
    skip_existing_outputs: bool = True
    adata_full_path: Optional[Path] = None
    adata_hvg_path: Optional[Path] = None
    splice_result_paths: Optional[Sequence[Path]] = None
    splice_annotation_filename: str = "merged_splice_transcript_annotations.csv"
    promoter_summary_path: Optional[Path] = None
    promoter_summary_filename: str = "Sup_Figure_4_last300_events_summary_with_uniqueness.csv"
    annotation_gtf_path: Optional[Path] = None
    isoform_group_ratio_filename: str = "isoform_ratios_by_gene_promoter_celltype.csv"

    def __post_init__(self) -> None:
        self.kallisto_dir = Path(self.kallisto_dir)
        self.data_dir = Path(self.data_dir)
        self.figure_dir = Path(self.figure_dir)
        self.isoform_match_paths = {
            str(label): Path(path) for label, path in self.isoform_match_paths.items()
        }
        self.isoform_transcripts = tuple(str(tx) for tx in self.isoform_transcripts)
        if self.isoform_ratio_pair is not None:
            self.isoform_ratio_pair = (
                str(self.isoform_ratio_pair[0]),
                str(self.isoform_ratio_pair[1]),
            )
        if self.adata_full_path is not None:
            self.adata_full_path = Path(self.adata_full_path)
        if self.adata_hvg_path is not None:
            self.adata_hvg_path = Path(self.adata_hvg_path)
        if self.marker_gene_path is not None:
            self.marker_gene_path = Path(self.marker_gene_path)
        if self.splice_result_paths is not None:
            self.splice_result_paths = tuple(Path(p) for p in self.splice_result_paths)
        if self.promoter_summary_path is not None:
            self.promoter_summary_path = Path(self.promoter_summary_path)
        if self.annotation_gtf_path is not None:
            self.annotation_gtf_path = Path(self.annotation_gtf_path)


@dataclass
class ScanpyAnalysisOutputs:
    adata_full: "AnnData"
    adata_hvg: "AnnData"
    barcode_table: pd.DataFrame
    barcode_table_path: Optional[Path]
    barcode_sample_paths: Dict[str, Path]
    cluster_labels: Mapping[str, str]
    cluster_debug: Mapping[str, Mapping[str, object]]
    marker_sets: Dict[str, List[str]]
    marker_source: str
    marker_coverage: pd.DataFrame
    marker_overlaps: pd.DataFrame
    isoform_counts_by_age: Dict[str, pd.DataFrame]
    isoform_counts_paths: Dict[str, Path]
    isoform_summary_paths: Dict[str, Path]
    combined_isoform_counts: pd.DataFrame
    combined_isoform_counts_path: Optional[Path]
    isoform_ratio_table: pd.DataFrame
    isoform_proportion_matrix: pd.DataFrame
    slope_data: pd.DataFrame
    splice_annotation_table: pd.DataFrame
    splice_annotation_path: Optional[Path]
    isoform_group_ratios: pd.DataFrame
    isoform_group_ratio_path: Optional[Path]


ALPHABET: Tuple[str, ...] = ("A", "C", "G", "T")
BARCODE_SUFFIX_RE = re.compile(r"-\d+$")
DEFAULT_MARKER_GENE_SETS: Dict[str, List[str]] = {
    "Oocyte (Primordial/Early Growing)": [
        "Nobox",
        "Figla",
        "Lhx8",
        "Sohlh1",
        "Sohlh2",
        "H1foo",
        "Marf1",
    ],
    "Oocyte (Growing/Maturing)": [
        "Gdf9",
        "Bmp15",
        "Zp2",
        "Zp3",
        "Zp1",
        "Nlrp5",
        "Padi6",
        "Tle6",
        "Zar1",
    ],
    "Oocyte (Mature/MII)": ["Mos", "Ooep", "Dppa3", "WEE2", "Fbxo43"],
    "Granulosa (Primordial/Primary)": [
        "Foxl2",
        "Lgr5",
        "Frzb",
        "Zfx",
        "Tcf21",
        "Rspo1",
        "Wnt4",
        "Ihh",
    ],
    "Granulosa (Pre-antral/Early Antral)": [
        "Amh",
        "Fshr",
        "Inha",
        "Inhbb",
        "Kctd14",
        "Gas6",
        "Kitl",
        "Esr2",
        "Fst",
    ],
    "Granulosa (Antral/Estrogenic)": [
        "Cyp19a1",
        "Hsd17b1",
        "Fshr",
        "Inhbb",
        "Igfbp5",
        "Gatm",
        "Nr5a2",
        "Hsd17b7",
    ],
    "Granulosa (Mural Preovulatory)": [
        "Lhcgr",
        "Ptgs2",
        "Pgr",
        "Star",
        "Cyp11a1",
        "Edn2",
        "Runx2",
        "Adamts1",
        "Areg",
        "Ereg",
        "Btc",
        "Egr1",
        "Nr4a1",
        "Klf4",
    ],
    "Granulosa (Cumulus—Periovulatory/Expansion)": [
        "Has2",
        "Tnfaip6",
        "Ptx3",
        "Icam1",
        "Ppap2b",
        "Ptgs2",
        "Npr2",
    ],
    "Granulosa (Mural vs Cumulus—Meiotic Arrest Axis)": ["Nppc"],
    "Granulosa (Mitotic/Proliferative)": ["Top2a", "Mki67"],
    "Granulosa (Atretic/Regressing Follicle)": [
        "Ghr",
        "Pik3ip1",
        "Nupr1",
        "Gadd45a",
        "Vat1",
        "Mitf",
    ],
    "Theca (Interna, Early/Preantral)": ["Hhip", "Ptch1"],
    "Theca (Interna, Steroidogenic)": [
        "Cyp17a1",
        "Insl3",
        "Hsd3b1",
        "Lhcgr",
        "Star",
        "Cyp11a1",
    ],
    "Theca (Externa, Contractile/Smooth Muscle-like)": [
        "Acta2",
        "Tagln",
        "Mfap5",
        "Myh11",
    ],
    "Large Luteal Cells (GC-derived)": [
        "Prlr",
        "Hsd3b1",
        "Akr1c18",
        "Cyp11a1",
        "Star",
        "Pgr",
        "Lhcgr",
    ],
    "Small Luteal Cells (Theca-derived)": [
        "Vegfa",
        "Corin",
        "Lpar3",
        "Lhcgr",
        "Star",
        "Hsd3b1",
    ],
    "Corpus Luteum (Active, Steroidogenic)": [
        "Star",
        "Cyp11a1",
        "Hsd3b1",
        "Lhcgr",
        "Prlr",
    ],
    "Corpus Luteum (Regressing/Luteolysis)": ["Cdkn1a", "Fas", "Il1b", "Mmp9"],
    "Stromal Fibroblasts (Interstitial/Fibroblast-like)": [
        "Col1a1",
        "Col1a2",
        "Dcn",
        "Lum",
        "Cxcl14",
        "Kcnk2",
        "Pdgfra",
    ],
    "Stromal (Steroidogenic/Interstitial)": ["Cyp11a1", "Star", "Ptch1"],
    "Endothelial Cells": ["Pecam1", "Kdr", "Tek", "Esam", "Emcn"],
    "Pericytes": ["Rgs5", "Pdgfrb", "Kcnj8", "Acta2", "Myh11"],
    "Macrophages/Monocytes": ["Lyz2", "Cd68", "Adgre1", "Csf1r", "Mrc1"],
    "Neutrophils": ["S100a8", "S100a9", "Cxcr2", "Lcn2"],
    "T Cells": ["Cd3d", "Cd3e", "Cd3g", "Trac"],
    "NK Cells": ["Nkg7", "Klrb1c", "Prf1", "Gzmb"],
    "Dendritic Cells": ["Itgax", "Siglech", "Xcr1", "Ccr7"],
    "Schwann/Peripheral Glia": ["Mpz", "Plp1", "Sox10", "Ngfr"],
    "OSE (Non-dividing)": ["Krt8", "Krt18", "Krt19", "Msln"],
    "OSE (Mitotic during Wound Repair)": ["Mki67", "Top2a", "Birc5"],
}

PRIORITY_ORDER: Tuple[str, ...] = (
    "Oocyte (Mature/MII)",
    "Oocyte (Growing/Maturing)",
    "Oocyte (Primordial/Early Growing)",
    "Granulosa (Mural Preovulatory)",
    "Granulosa (Cumulus—Periovulatory/Expansion)",
    "Granulosa (Antral/Estrogenic)",
    "Granulosa (Pre-antral/Early Antral)",
    "Granulosa (Primordial/Primary)",
    "Granulosa (Mitotic/Proliferative)",
    "Granulosa (Atretic/Regressing Follicle)",
    "Theca (Interna, Steroidogenic)",
    "Theca (Interna, Early/Preantral)",
    "Theca (Externa, Contractile/Smooth Muscle-like)",
    "Large Luteal Cells (GC-derived)",
    "Small Luteal Cells (Theca-derived)",
    "Corpus Luteum (Active, Steroidogenic)",
    "Corpus Luteum (Regressing/Luteolysis)",
    "OSE (Non-dividing)",
    "OSE (Mitotic during Wound Repair)",
    "Endothelial Cells",
    "Pericytes",
    "Macrophages/Monocytes",
    "Neutrophils",
    "T Cells",
    "NK Cells",
    "Dendritic Cells",
    "Schwann/Peripheral Glia",
    "Stromal (Steroidogenic/Interstitial)",
    "Stromal Fibroblasts (Interstitial/Fibroblast-like)",
)

AGG_QUANTILE = 0.90
MARGIN_THRESHOLD = 0.0
MIN_DETECT_FRAC_DEFAULT = 0.05
MIN_DETECT_GENES = 2
MIN_DETECT_FRAC_BY_TYPE = {
    "Stromal Fibroblasts (Interstitial/Fibroblast-like)": 0.40,
}
def load_pfam_hits(pfam_paths: Mapping[str, Path]) -> pd.DataFrame:
    """
    Load multiple PFAM hit tables and append a source label.

    Parameters
    ----------
    pfam_paths : Mapping[str, Path]
        Dictionary-like mapping of label -> CSV path.

    Returns
    -------
    pd.DataFrame
        Concatenated PFAM hit records with a pfam_source column.
    """
    frames = []
    for label, path in pfam_paths.items():
        df = pd.read_csv(path, index_col=0)
        df = df.reset_index(drop=True)
        df["pfam_source"] = label
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    combined["transcript_id"] = combined["transcript_id"].astype(str).str.strip()
    return combined


def load_transcript_sequences(
    fasta_path: Path, transcript_ids: Optional[Iterable[str]] = None
) -> Dict[str, str]:
    """
    Load transcript sequences from a FASTA file.

    Parameters
    ----------
    fasta_path : Path
        FASTA file containing transcript cDNA sequences.
    transcript_ids : Iterable[str], optional
        If provided, restrict the output to these transcript identifiers.

    Returns
    -------
    Dict[str, str]
        Mapping of transcript_id -> uppercase nucleotide sequence.
    """
    allowed = None
    if transcript_ids is not None:
        allowed = {str(tid).strip() for tid in transcript_ids if tid}
    sequences: Dict[str, str] = {}
    for header, sequence in _iter_fasta_records(fasta_path):
        transcript_id = header.split("_", 1)[0].strip()
        if allowed is not None and transcript_id not in allowed:
            continue
        sequences[transcript_id] = sequence.upper()
    return sequences


def load_splicing_universe(splice_paths: Iterable[Path]) -> pd.Index:
    """
    Build the universe of expressed transcripts from splice result tables.

    Parameters
    ----------
    splice_paths : Iterable[Path]
        Iterable of CSV paths that contain an `ids` column where the first token
        encodes the transcript identifier.

    Returns
    -------
    pd.Index
        Sorted unique transcript identifiers observed across all tables.
    """
    transcripts: set[str] = set()
    for path in splice_paths:
        df = pd.read_csv(path, usecols=["ids"])
        ids = (
            df["ids"]
            .dropna()
            .astype(str)
            .str.strip()
            .str.split("_", n=1)
            .str[0]
            .str.strip()
        )
        transcripts.update(ids[ids.notna()])
    return pd.Index(sorted(transcripts))


def unique_transcript_ids(pfam_hits: pd.DataFrame) -> pd.Index:
    """
    Return the unique transcript identifiers present in the PFAM hits table.

    Parameters
    ----------
    pfam_hits : pd.DataFrame
        Output from load_pfam_hits or a similar table with transcript_id column.

    Returns
    -------
    pd.Index
        Sorted unique transcript identifiers.
    """
    if "transcript_id" not in pfam_hits:
        raise KeyError("pfam_hits DataFrame must include a transcript_id column.")
    unique_ids = (
        pfam_hits["transcript_id"]
        .dropna()
        .astype(str)
        .str.strip()
        .unique()
    )
    return pd.Index(sorted(unique_ids))


def load_gtf_for_transcripts(
    gtf_path: Path, transcript_ids: Iterable[str], usecols: Sequence[str] = DEFAULT_GTF_COLUMNS
) -> pd.DataFrame:
    """
    Load the GTF/TSV rows corresponding to the supplied transcript IDs.

    Parameters
    ----------
    gtf_path : Path
        Tab-delimited file produced by FLAIR with exon level annotations.
    transcript_ids : Iterable[str]
        Transcript identifiers to retain from the GTF.
    usecols : Sequence[str], optional
        Subset of columns to load from the file. Defaults to key structural columns.

    Returns
    -------
    pd.DataFrame
        Rows limited to the requested transcripts.
    """
    transcripts = {str(tid) for tid in transcript_ids if tid}
    if not transcripts:
        columns = list(usecols)
        return pd.DataFrame(columns=columns)

    gtf_df = pd.read_csv(gtf_path, sep="\t", usecols=usecols, dtype=str)
    gtf_df = gtf_df.dropna(subset=["transcript_id"]).copy()
    gtf_df["transcript_id"] = gtf_df["transcript_id"].astype(str).str.strip()
    gtf_df = gtf_df[gtf_df["transcript_id"].isin(transcripts)].copy()
    if gtf_df.empty:
        return gtf_df
    gtf_df["start"] = gtf_df["start"].astype(int)
    gtf_df["end"] = gtf_df["end"].astype(int)
    return gtf_df


# Backward-compatible alias for legacy callers.
load_gtf_subset = load_gtf_for_transcripts


def build_promoter_group_map(gtf_subset: pd.DataFrame) -> pd.DataFrame:
    """
    Construct a promoter group to transcript mapping from a GTF subset.

    Parameters
    ----------
    gtf_subset : pd.DataFrame
        GTF rows that include promoter_group and transcript_id columns.

    Returns
    -------
    pd.DataFrame
        Columns: promoter_group, transcript_id with unique combinations.
    """
    required = {"gene_id", "promoter_group", "transcript_id"}
    missing = required.difference(gtf_subset.columns)
    if missing:
        raise KeyError(f"GTF subset missing columns: {', '.join(sorted(missing))}")
    mapping = (
        gtf_subset.loc[:, ["gene_id", "promoter_group", "transcript_id"]]
        .dropna()
        .assign(
            gene_id=lambda df: df["gene_id"].astype(str).str.strip(),
            promoter_group=lambda df: df["promoter_group"].astype(str).str.strip(),
            transcript_id=lambda df: df["transcript_id"].astype(str).str.strip(),
        )
        .drop_duplicates()
    )
    mapping = mapping[mapping["gene_id"] != ""]
    mapping = mapping[mapping["promoter_group"] != ""]
    mapping = mapping[mapping["transcript_id"] != ""]
    return mapping.reset_index(drop=True)


def annotate_transcript_exons(gtf_subset: pd.DataFrame) -> pd.DataFrame:
    """
    Annotate exon features with transcript-relative coordinates.

    Parameters
    ----------
    gtf_subset : pd.DataFrame
        Subset of the GTF limited to transcripts of interest.

    Returns
    -------
    pd.DataFrame
        Exon (including alternative event) rows with additional columns:
        exon_length, transcript_offset_start, transcript_offset_end, transcript_length.
    """
    required = {"transcript_id", "feature", "start", "end", "strand"}
    missing = required.difference(gtf_subset.columns)
    if missing:
        raise KeyError(f"GTF subset missing columns: {', '.join(sorted(missing))}")

    exon_mask = gtf_subset["feature"].astype(str).str.startswith("exon")
    exons = gtf_subset.loc[exon_mask].copy()
    if exons.empty:
        columns = list(gtf_subset.columns) + [
            "exon_length",
            "transcript_offset_start",
            "transcript_offset_end",
            "transcript_length",
        ]
        return pd.DataFrame(columns=columns)

    annotated_frames = []
    for transcript_id, group in exons.groupby("transcript_id"):
        ordered = group.copy()
        strand_series = ordered["strand"].dropna()
        strand = strand_series.iloc[0] if not strand_series.empty else "+"
        ascending = strand != "-"
        ordered = ordered.sort_values(["start", "end"], ascending=ascending)
        offset = 0
        records = []
        for _, row in ordered.iterrows():
            exon_length = int(row["end"]) - int(row["start"]) + 1
            if exon_length <= 0:
                continue
            annotated_row = row.copy()
            annotated_row["exon_length"] = exon_length
            annotated_row["transcript_offset_start"] = offset
            annotated_row["transcript_offset_end"] = offset + exon_length
            offset += exon_length
            records.append(annotated_row)
        if not records:
            continue
        frame = pd.DataFrame(records)
        frame["transcript_length"] = offset
        annotated_frames.append(frame)

    if not annotated_frames:
        columns = list(gtf_subset.columns) + [
            "exon_length",
            "transcript_offset_start",
            "transcript_offset_end",
            "transcript_length",
        ]
        return pd.DataFrame(columns=columns)

    annotated = pd.concat(annotated_frames, ignore_index=True)
    return annotated


def identify_events_within_last_window(
    annotated_exons: pd.DataFrame, window_nt: int = 300
) -> pd.DataFrame:
    """
    Identify alternative exon features that overlap the last N nucleotides.

    Parameters
    ----------
    annotated_exons : pd.DataFrame
        Output from annotate_transcript_exons.
    window_nt : int, optional
        Window size measured from the 3' end (defaults to 300 nt).

    Returns
    -------
    pd.DataFrame
        Alternative exon rows with overlap metadata for the final window.
    """
    if annotated_exons.empty:
        return annotated_exons

    events = annotated_exons[
        annotated_exons["feature"].astype(str).str.startswith(EVENT_PREFIX)
    ].copy()
    if events.empty:
        return events

    window_start = (events["transcript_length"] - window_nt).clip(lower=0).to_numpy(dtype=int)
    window_end = events["transcript_length"].to_numpy(dtype=int)
    start = events["transcript_offset_start"].to_numpy(dtype=int)
    end = events["transcript_offset_end"].to_numpy(dtype=int)

    overlap = np.minimum(end, window_end) - np.maximum(start, window_start)
    overlap = np.clip(overlap, a_min=0, a_max=None)

    events["window_overlap_nt"] = overlap
    events["overlaps_last_window"] = events["window_overlap_nt"] > 0
    events["distance_to_3prime"] = (
        events["transcript_length"] - events["transcript_offset_end"]
    ).astype(int)
    return events


def summarize_events_in_window(events: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize alternative exon events that overlap the tail window.

    Parameters
    ----------
    events : pd.DataFrame
        Output from identify_events_within_last_window.

    Returns
    -------
    pd.DataFrame
        Transcript level summary with event counts and feature list.
    """
    if events.empty:
        return pd.DataFrame(columns=["transcript_id", "events_in_window", "event_features"])

    filtered = events.loc[events["overlaps_last_window"]].copy()
    if filtered.empty:
        return pd.DataFrame(columns=["transcript_id", "events_in_window", "event_features"])

    summary = (
        filtered.groupby("transcript_id")
        .agg(
            events_in_window=("feature", "count"),
            event_features=("feature", lambda vals: ";".join(sorted(set(vals)))),
            transcript_length=("transcript_length", "max"),
        )
        .reset_index()
    )
    promoter_lookup = (
        events[["transcript_id", "promoter_group"]]
        .dropna(subset=["transcript_id"])
        .drop_duplicates()
    )
    summary = summary.merge(promoter_lookup, on="transcript_id", how="left")
    gene_lookup = (
        events[["transcript_id", "gene_id"]]
        .dropna(subset=["transcript_id"])
        .drop_duplicates()
    )
    summary = summary.merge(gene_lookup, on="transcript_id", how="left")
    return summary


def subset_fasta_by_transcripts(
    fasta_path: Path, transcript_ids: Iterable[str], output_path: Optional[Path] = None
) -> Dict[str, str]:
    """
    Extract cDNA sequences for the supplied transcript identifiers from a FASTA file.

    Parameters
    ----------
    fasta_path : Path
        FASTA file containing transcript translations.
    transcript_ids : Iterable[str]
        Transcript identifiers to retain.
    output_path : Path, optional
        If provided, write the filtered FASTA to this path.

    Returns
    -------
    Dict[str, str]
        Mapping of FASTA headers to sequence strings that matched the transcript IDs.
    """
    targets = {str(tid) for tid in transcript_ids if tid}
    if not targets:
        if output_path is not None:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text("")
        return {}

    selected: Dict[str, str] = {}
    writer = None
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        writer = output_path.open("w")

    try:
        for header, sequence in _iter_fasta_records(fasta_path):
            transcript_id = header.split("_", 1)[0]
            if transcript_id in targets:
                selected[header] = sequence
                if writer is not None:
                    writer.write(f">{header}\n{sequence}\n")
    finally:
        if writer is not None:
            writer.close()

    if output_path is not None and not selected:
        output_path.write_text("")  # ensure empty file exists

    return selected


def _read_fasta_transcript_ids(fasta_path: Path) -> List[str]:
    ids: List[str] = []
    for header, _ in _iter_fasta_records(fasta_path):
        transcript_id = header.split("_", 1)[0].strip()
        if transcript_id:
            ids.append(transcript_id)
    return ids


def _iter_fasta_records(fasta_path: Path) -> Iterator[Tuple[str, str]]:
    """
    Iterate over header/sequence pairs in a FASTA file.
    """
    header: Optional[str] = None
    seq_chunks: list[str] = []
    with Path(fasta_path).open("r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def filter_table_by_transcripts(
    table: pd.DataFrame, transcripts: Iterable[str], column: str = "transcript_id"
) -> pd.DataFrame:
    """
    Filter a transcript level table to a provided transcript universe.
    """
    universe = {str(tid).strip() for tid in transcripts if tid}
    if not universe or column not in table.columns:
        return table.copy()
    return table[table[column].astype(str).isin(universe)].copy()


def expand_transcripts_by_promoter_group(
    transcripts: Iterable[str],
    promoter_map: pd.DataFrame,
    universe: Optional[Iterable[str]] = None,
) -> pd.Index:
    """
    Expand a transcript list to include all promoter group members.

    Parameters
    ----------
    transcripts : Iterable[str]
        Input transcript identifiers.
    promoter_map : pd.DataFrame
        Output from build_promoter_group_map.
    universe : Iterable[str], optional
        Optional transcript universe to filter against.

    Returns
    -------
    pd.Index
        Sorted unique transcript identifiers including promoter group members.
    """
    seeds = {str(tid).strip() for tid in transcripts if tid}
    if not seeds or promoter_map.empty:
        return pd.Index(sorted(seeds))

    promoter_map = promoter_map.copy()
    promoter_map["gene_id"] = promoter_map["gene_id"].astype(str).str.strip()
    promoter_map["promoter_group"] = promoter_map["promoter_group"].astype(str).str.strip()
    promoter_map["transcript_id"] = promoter_map["transcript_id"].astype(str).str.strip()
    promoter_map["combo"] = list(
        zip(promoter_map["gene_id"], promoter_map["promoter_group"])
    )

    transcript_to_combos = (
        promoter_map.groupby("transcript_id")["combo"]
        .agg(lambda vals: {val for val in vals if val})
        .to_dict()
    )
    combo_to_transcripts = (
        promoter_map.groupby("combo")["transcript_id"]
        .agg(lambda vals: {val for val in vals if val})
        .to_dict()
    )

    expanded: set[str] = set()
    for tid in seeds:
        expanded.add(tid)
        combos = transcript_to_combos.get(tid, set())
        for combo in combos:
            expanded.update(combo_to_transcripts.get(combo, set()))

    if universe is not None:
        allowed = {str(tid).strip() for tid in universe if tid}
        expanded.intersection_update(allowed)

    return pd.Index(sorted(expanded))


TAIL_WINDOWS: Dict[str, Tuple[int, int]] = {
    "tail_300": (-300, 0),
    "tail_first90": (-300, -210),
    "tail_last150": (-150, 0),
}


def evaluate_tail_uniqueness(
    annotated_exons: pd.DataFrame,
    window_map: Optional[Mapping[str, Tuple[int, int]]] = None,
    sequences: Optional[Mapping[str, str]] = None,
) -> pd.DataFrame:
    """
    Evaluate whether tail windows are unique within each promoter group.

    Parameters
    ----------
    annotated_exons : pd.DataFrame
        Output from annotate_transcript_exons, must include promoter_group.
    window_map : Mapping[str, Tuple[int, int]], optional
        Mapping of window label to (start_offset, end_offset) relative to the 3'
        end of the transcript. Offsets are negative numbers where 0 represents
        the transcript end. Defaults to TAIL_WINDOWS.
    sequences : Mapping[str, str], optional
        Preloaded transcript sequences keyed by transcript_id. When provided,
        uniqueness is evaluated at the sequence level; otherwise genomic
        segments inferred from the GTF coordinates are used.

    Returns
    -------
    pd.DataFrame
        Columns include transcript_id, promoter_group, window signatures, and
        boolean uniqueness flags.
    """
    default_windows = window_map or TAIL_WINDOWS
    if annotated_exons.empty:
        base_cols = ["transcript_id", "gene_id", "promoter_group", "transcript_length"]
        extra_cols = []
        for name in default_windows:
            extra_cols.extend(
                [
                    f"{name}_length",
                    f"{name}_sequence",
                    f"{name}_signature",
                    f"is_unique_{name}",
                ]
            )
        return pd.DataFrame(columns=base_cols + extra_cols)

    windows = dict(default_windows)
    per_transcript: Dict[str, Dict[str, object]] = {}

    grouped = annotated_exons.groupby("transcript_id")
    for transcript_id, rows in grouped:
        rows = rows.sort_values("transcript_offset_start")
        promoter_group = (
            rows["promoter_group"].dropna().astype(str).str.strip().iloc[0]
            if not rows["promoter_group"].dropna().empty
            else None
        )
        gene_series = rows["gene_id"].dropna().astype(str).str.strip()
        gene_id = gene_series.iloc[0] if not gene_series.empty else None
        seq = None
        if sequences is not None:
            seq = sequences.get(transcript_id)
            if seq is not None:
                seq = seq.upper().replace("\n", "").strip()
                if not seq:
                    seq = None
        if seq is not None:
            transcript_length = len(seq)
        else:
            length_series = rows["transcript_length"].dropna()
            transcript_length = int(length_series.iloc[0]) if not length_series.empty else None
        promoter_group = (
            promoter_group
        )
        signatures: Dict[str, Dict[str, object]] = {}
        for label, (offset_start, offset_end) in windows.items():
            if transcript_length is None:
                signatures[label] = {
                    "sequence": "",
                    "signature": "",
                    "length": 0,
                }
                continue
            rel_start = max(transcript_length + offset_start, 0)
            rel_end = max(transcript_length + offset_end, 0)
            if rel_end <= rel_start:
                signatures[label] = {
                    "sequence": "",
                    "signature": "",
                    "length": 0,
                }
                continue
            if seq is not None:
                snippet = seq[rel_start:rel_end]
                signature_value = snippet
                length_value = len(snippet)
                signatures[label] = {
                    "sequence": snippet,
                    "signature": signature_value,
                    "length": length_value,
                }
            else:
                segments = _compute_tail_segments(rows, rel_start, rel_end)
                length_value = rel_end - rel_start
                signatures[label] = {
                    "sequence": "",
                    "signature": _format_segments(segments),
                    "length": length_value,
                }
        per_transcript[transcript_id] = {
            "promoter_group": promoter_group,
            "gene_id": gene_id,
            "transcript_length": transcript_length,
            "signatures": signatures,
        }

    by_combo: Dict[Tuple[Optional[str], Optional[str]], Dict[str, Dict[str, list]]] = {}
    for transcript_id, data in per_transcript.items():
        promoter_group = data["promoter_group"]
        gene_id = data.get("gene_id")
        if promoter_group is None or gene_id is None:
            continue
        combo = (gene_id, promoter_group)
        promoter_bucket = by_combo.setdefault(combo, {})
        for label, sig_info in data["signatures"].items():
            promoter_bucket.setdefault(label, {}).setdefault(sig_info["signature"], []).append(
                transcript_id
            )

    records = []
    for transcript_id, data in per_transcript.items():
        promoter_group = data["promoter_group"]
        gene_id = data.get("gene_id")
        record = {
            "transcript_id": transcript_id,
            "promoter_group": promoter_group,
            "gene_id": gene_id,
            "transcript_length": data["transcript_length"],
        }
        for label, sig_info in data["signatures"].items():
            signature = sig_info.get("signature", "")
            sequence_value = sig_info.get("sequence", "")
            length_value = int(sig_info.get("length", 0) or 0)
            record[f"{label}_length"] = length_value
            record[f"{label}_sequence"] = sequence_value
            record[f"{label}_signature"] = signature
            if length_value <= 0 or promoter_group is None or gene_id is None or not signature:
                record[f"is_unique_{label}"] = pd.NA
                continue
            matches = by_combo.get((gene_id, promoter_group), {}).get(label, {}).get(
                signature, []
            )
            record[f"is_unique_{label}"] = len(matches) <= 1
        records.append(record)

    return pd.DataFrame(records)


def _compute_tail_segments(rows: pd.DataFrame, rel_start: int, rel_end: int) -> Tuple[Tuple[str, int, int, str], ...]:
    """
    Compute genomic segments covering the requested transcript window.
    """
    segments: List[Tuple[str, int, int, str]] = []
    for _, row in rows.iterrows():
        exon_start = int(row["transcript_offset_start"])
        exon_end = int(row["transcript_offset_end"])
        overlap_start = max(rel_start, exon_start)
        overlap_end = min(rel_end, exon_end)
        if overlap_end <= overlap_start:
            continue
        within_start = overlap_start - exon_start
        within_end = overlap_end - exon_start
        genomic_start, genomic_end = _map_within_exon_to_genome(
            int(row["start"]),
            int(row["end"]),
            row["strand"],
            within_start,
            within_end,
        )
        segments.append((str(row["seqname"]), int(genomic_start), int(genomic_end), str(row["strand"])))
    return tuple(segments)


def _map_within_exon_to_genome(
    exon_start: int,
    exon_end: int,
    strand: str,
    within_start: int,
    within_end: int,
) -> Tuple[int, int]:
    """
    Convert transcript-relative offsets within an exon back to genomic coordinates.
    """
    if strand == "-":
        genomic_end = exon_end - within_start
        genomic_start = exon_end - within_end + 1
    else:
        genomic_start = exon_start + within_start
        genomic_end = exon_start + within_end - 1
    if genomic_start > genomic_end:
        genomic_start, genomic_end = genomic_end, genomic_start
    return genomic_start, genomic_end


def _format_segments(segments: Tuple[Tuple[str, int, int, str], ...]) -> str:
    if not segments:
        return ""
    parts = []
    for seqname, start, end, strand in segments:
        parts.append(f"{seqname}:{start}-{end}:{strand}")
    return ";".join(parts)


def generate_sup_figure_4_outputs(config: SupFigure4Config) -> SupFigure4Outputs:
    """
    Run the Supplemental Figure 4 preprocessing pipeline and materialize outputs.
    """
    summary_filename = (
        config.summary_filename
        or f"Sup_Figure_4_last{config.window_nt}_events_summary_with_uniqueness.csv"
    )
    fasta_filename = (
        config.fasta_filename
        or f"Sup_Figure_4_last{config.window_nt}_promoter_group_transcripts.fa"
    )
    plot_filename = (
        config.plot_filename or f"Sup_Figure_4_last{config.window_nt}_event_counts.png"
    )
    summary_path = config.figure_dir / summary_filename
    fasta_output_path = config.figure_dir / fasta_filename
    plot_path = config.figure_dir / plot_filename if config.make_plot else None

    if (
        config.skip_existing_outputs
        and summary_path.exists()
        and fasta_output_path.exists()
    ):
        return _load_cached_sup_figure_outputs(summary_path, fasta_output_path, plot_path)

    config.figure_dir.mkdir(parents=True, exist_ok=True)
    missing_pfam = [label for label, path in config.pfam_paths.items() if not path.exists()]
    if missing_pfam:
        raise FileNotFoundError(f"Missing PFAM tables for: {', '.join(sorted(missing_pfam))}")

    missing_splice = [str(path) for path in config.splice_paths if not path.exists()]
    if missing_splice:
        raise FileNotFoundError(f"Missing splice universe files: {', '.join(missing_splice)}")
    if not config.gtf_path.exists():
        raise FileNotFoundError(config.gtf_path)
    if not config.fasta_path.exists():
        raise FileNotFoundError(config.fasta_path)

    transcript_universe = load_splicing_universe(config.splice_paths)
    pfam_hits = load_pfam_hits(config.pfam_paths)
    pfam_hits = filter_table_by_transcripts(pfam_hits, transcript_universe)
    pfam_transcripts = unique_transcript_ids(pfam_hits)

    gtf_universe = load_gtf_for_transcripts(config.gtf_path, transcript_universe)
    annotated_universe = annotate_transcript_exons(gtf_universe)
    promoter_map = build_promoter_group_map(gtf_universe)
    transcript_sequences = load_transcript_sequences(config.fasta_path, transcript_universe)

    events = identify_events_within_last_window(annotated_universe, window_nt=config.window_nt)
    summary_all = summarize_events_in_window(events)
    summary = filter_table_by_transcripts(summary_all, pfam_transcripts)

    tail_uniqueness_all = evaluate_tail_uniqueness(
        annotated_universe,
        window_map=config.tail_windows,
        sequences=transcript_sequences,
    )
    tail_uniqueness = filter_table_by_transcripts(tail_uniqueness_all, pfam_transcripts)
    final_summary = _merge_summary_with_uniqueness(summary, tail_uniqueness)

    final_summary.to_csv(summary_path, index=False)

    unique_hits_series = final_summary.loc[
        final_summary["is_tail_unique_all"] == True, "transcript_id"
    ].dropna()
    unique_hits = pd.Index(unique_hits_series.astype(str))

    if unique_hits.empty:
        expanded = pd.Index([])
        fasta_output_path: Optional[Path] = None
    else:
        expanded = expand_transcripts_by_promoter_group(unique_hits, promoter_map, transcript_universe)
        if expanded.empty:
            fasta_output_path = None
        else:
            subset_fasta_by_transcripts(config.fasta_path, expanded, fasta_output_path)

    if plot_path and not final_summary.empty:
        _save_event_distribution_plot(final_summary, config.window_nt, plot_path)
    elif plot_path:
        plot_path = None

    return SupFigure4Outputs(
        transcript_universe=transcript_universe,
        pfam_hits=pfam_hits,
        events=events,
        summary=summary,
        tail_uniqueness=tail_uniqueness,
        final_summary=final_summary,
        unique_hits=unique_hits,
        expanded_transcripts=expanded,
        summary_path=summary_path,
        fasta_output=fasta_output_path,
        plot_path=plot_path,
    )


def _merge_summary_with_uniqueness(
    summary: pd.DataFrame, tail_uniqueness: pd.DataFrame
) -> pd.DataFrame:
    if summary.empty:
        return summary.assign(is_tail_unique_all=False)

    uniqueness_cols = [
        col for col in tail_uniqueness.columns if col.startswith("is_unique_")
    ]
    merged = summary.merge(
        tail_uniqueness.drop(columns=["transcript_length"], errors="ignore"),
        on=["transcript_id", "gene_id", "promoter_group"],
        how="left",
    )
    merged["is_tail_unique_all"] = _resolve_tail_uniqueness_flag(merged, uniqueness_cols)
    return merged


def _resolve_tail_uniqueness_flag(
    final_summary: pd.DataFrame, uniqueness_cols: Sequence[str]
) -> pd.Series:
    if not uniqueness_cols:
        return pd.Series(False, index=final_summary.index)

    def _col_true(column: str) -> pd.Series:
        return final_summary[column].fillna(False).eq(True)

    first90 = "is_unique_tail_first90"
    last150 = "is_unique_tail_last150"
    colset = set(uniqueness_cols)
    if {first90, last150}.issubset(colset):
        return _col_true(first90) | _col_true(last150)
    if first90 in colset:
        return _col_true(first90)
    if last150 in colset:
        return _col_true(last150)

    combined = pd.Series(True, index=final_summary.index)
    for column in uniqueness_cols:
        combined &= _col_true(column)
    return combined


def _save_event_distribution_plot(
    final_summary: pd.DataFrame, window_nt: int, output_path: Path
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts = final_summary["events_in_window"].value_counts().sort_index()
    if counts.empty:
        output_path.touch()
        return
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4))
    counts.plot(kind="bar", ax=ax, title=f"Alternative events within last {window_nt} nt")
    ax.set_xlabel("Events in final window")
    ax.set_ylabel("Transcript count")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def _load_cached_sup_figure_outputs(
    summary_path: Path, fasta_path: Path, plot_path: Optional[Path]
) -> SupFigure4Outputs:
    final_summary = pd.read_csv(summary_path)
    summary = final_summary.copy()
    unique_hits = pd.Index(
        final_summary.loc[
            final_summary.get("is_tail_unique_all", False) == True, "transcript_id"
        ]
        .dropna()
        .astype(str)
        .unique()
    )
    expanded = pd.Index(sorted(set(_read_fasta_transcript_ids(fasta_path))))
    return SupFigure4Outputs(
        transcript_universe=pd.Index([]),
        pfam_hits=pd.DataFrame(),
        events=pd.DataFrame(),
        summary=summary,
        tail_uniqueness=pd.DataFrame(),
        final_summary=final_summary,
        unique_hits=unique_hits,
        expanded_transcripts=expanded,
        summary_path=summary_path,
        fasta_output=fasta_path,
        plot_path=plot_path if plot_path and plot_path.exists() else None,
    )


def _decompress_gzip_file(source: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rb") as src, dest.open("wb") as dst:
        shutil.copyfileobj(src, dst)


def _ensure_anndata_null_reader() -> None:
    """
    anndata<=0.10 does not always register a reader for datasets encoded as "null".
    Install a no-op reader so legacy H5AD files with empty scalars can load.
    """
    global _ANNDATA_NULL_PATCHED
    if _ANNDATA_NULL_PATCHED:
        return
    try:
        import h5py  # type: ignore
        from anndata._io.specs import IOSpec, registry as _anndata_registry  # type: ignore
    except ImportError:
        return

    spec = IOSpec("null", "0.1.0")
    registry = _anndata_registry._REGISTRY  # type: ignore[attr-defined]

    def _register(dataset_type: type) -> None:
        key = (dataset_type, spec, frozenset())
        if key in registry.read:
            return

        @registry.register_read(dataset_type, spec)  # type: ignore[misc]
        def _read_null(elem, *, _reader):
            return None

    _register(h5py.Dataset)
    try:
        import zarr  # type: ignore
    except ImportError:  # pragma: no cover - optional dependency
        zarr = None
    if zarr is not None:
        _register(zarr.Array)  # type: ignore[attr-defined]
    _ANNDATA_NULL_PATCHED = True


def clean_barcode(value: object) -> str:
    if not isinstance(value, str):
        return ""
    return value.strip().upper()


def _standardize_barcode(value: object) -> str:
    cleaned = clean_barcode(value)
    if not cleaned:
        return ""
    return BARCODE_SUFFIX_RE.sub("", cleaned)


def harmonize_gene_list(
    gene_symbols: Iterable[str], var_names: Sequence[str]
) -> Tuple[List[str], List[str]]:
    lookup = {g.upper(): g for g in var_names}
    mapped: List[str] = []
    missing: List[str] = []
    for symbol in gene_symbols:
        key = symbol.upper()
        if key in lookup:
            mapped.append(lookup[key])
        else:
            missing.append(symbol)
    # preserve order but drop duplicates
    seen: Set[str] = set()
    unique_mapped: List[str] = []
    for gene in mapped:
        if gene not in seen:
            unique_mapped.append(gene)
            seen.add(gene)
    return unique_mapped, missing


def deduplicate_markers(
    markers: Mapping[str, Sequence[str]],
    var_names: Sequence[str],
    priority: Sequence[str],
) -> Tuple[Dict[str, List[str]], pd.DataFrame, pd.DataFrame]:
    coverage_rows: List[Tuple[str, int, int, int]] = []
    mapped_markers: Dict[str, List[str]] = {}
    for cell_type, genes in markers.items():
        mapped, missing = harmonize_gene_list(genes, var_names)
        mapped_markers[cell_type] = mapped
        coverage_rows.append((cell_type, len(genes), len(mapped), len(missing)))
    coverage_df = pd.DataFrame(
        coverage_rows, columns=["Cell Type", "n_input", "n_present", "n_missing"]
    ).set_index("Cell Type")

    unique_markers: Dict[str, List[str]] = {key: [] for key in markers}
    used: Set[str] = set()
    for cell_type in priority:
        genes = mapped_markers.get(cell_type)
        if not genes:
            continue
        for gene in genes:
            if gene not in used:
                unique_markers[cell_type].append(gene)
                used.add(gene)

    overlap_rows = []
    gene_to_types: Dict[str, Set[str]] = {}
    for cell_type, genes in mapped_markers.items():
        for gene in genes:
            gene_to_types.setdefault(gene, set()).add(cell_type)
    for gene, cell_types in gene_to_types.items():
        if len(cell_types) <= 1:
            continue
        kept = next((ct for ct in priority if ct in cell_types), None)
        overlap_rows.append(
            {
                "gene": gene,
                "n_types": len(cell_types),
                "kept_for": kept,
                "others": sorted(ct for ct in cell_types if ct != kept),
            }
        )
    overlap_df = pd.DataFrame(overlap_rows)
    if not overlap_df.empty:
        overlap_df = overlap_df.sort_values("n_types", ascending=False)
    return unique_markers, coverage_df, overlap_df


def _normalize_marker_column(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.strip().lower())


def _load_marker_gene_sets_from_frame(frame: pd.DataFrame) -> Dict[str, List[str]]:
    if frame.empty:
        raise ValueError("Marker gene table is empty.")
    normalized = {_normalize_marker_column(col): col for col in frame.columns}
    cell_col = None
    gene_col = None
    for key in ("celltype", "celltypes", "celltypename", "celltypeage", "cluster"):
        if key in normalized:
            cell_col = normalized[key]
            break
    for key in ("gene", "genes", "genesymbol", "marker", "markers"):
        if key in normalized:
            gene_col = normalized[key]
            break
    if cell_col is None or gene_col is None:
        raise ValueError(
            "Marker gene table must contain columns for cell type and gene names."
        )
    subset = frame[[cell_col, gene_col]].dropna(subset=[gene_col])
    subset[cell_col] = subset[cell_col].astype(str).str.strip()
    subset[gene_col] = subset[gene_col].astype(str).str.strip()
    subset = subset[(subset[cell_col] != "") & (subset[gene_col] != "")]
    marker_map: Dict[str, List[str]] = {}
    for cell_type, group in subset.groupby(cell_col, sort=False):
        genes = [gene for gene in group[gene_col].tolist() if gene]
        if genes:
            marker_map[str(cell_type)] = genes
    if not marker_map:
        raise ValueError("Marker gene table did not produce any usable entries.")
    return marker_map


def _load_marker_gene_sets_from_file(marker_path: Path) -> Dict[str, List[str]]:
    if not marker_path.exists():
        raise FileNotFoundError(marker_path)
    suffix = marker_path.suffix.lower()
    if suffix in {".json", ".jsn"}:
        data = json.loads(marker_path.read_text())
        if isinstance(data, dict):
            mapping_obj = data
        else:
            raise ValueError(
                "Marker gene JSON must be a mapping of cell types to gene lists."
            )
        marker_map: Dict[str, List[str]] = {}
        for cell_type, genes in mapping_obj.items():
            cell_key = str(cell_type).strip()
            if not cell_key:
                continue
            if genes is None:
                continue
            if isinstance(genes, (list, tuple, set)):
                items = list(genes)
            else:
                items = [genes]
            cleaned = [str(g).strip() for g in items if str(g).strip()]
            if cleaned:
                marker_map[cell_key] = cleaned
        if not marker_map:
            raise ValueError(
                f"Marker gene JSON at {marker_path} did not contain usable entries."
            )
        return marker_map
    if suffix in {".tsv", ".tab"}:
        frame = pd.read_csv(marker_path, sep="\t")
        return _load_marker_gene_sets_from_frame(frame)
    if suffix in {".csv", ".txt"}:
        frame = pd.read_csv(marker_path)
        return _load_marker_gene_sets_from_frame(frame)
    raise ValueError(
        f"Unsupported marker gene file format '{marker_path.suffix}'. "
        "Use .json, .csv, or .tsv."
    )


def detect_fraction(
    adata: "AnnData",
    cell_mask: np.ndarray,
    genes: Sequence[str],
    threshold: float = 0.0,
    min_detect: int = 1,
) -> float:
    if len(genes) == 0 or cell_mask.sum() == 0:
        return 0.0
    gene_idx = adata.var_names.get_indexer(genes)
    gene_idx = gene_idx[gene_idx >= 0]
    if gene_idx.size == 0:
        return 0.0
    matrix = adata.raw[cell_mask, gene_idx].X
    if sp.issparse(matrix):
        detected = (matrix > threshold).sum(axis=1).A1
    else:
        detected = (matrix > threshold).sum(axis=1).ravel()
    need = min(min_detect, len(genes))
    return float((detected >= need).mean())


class BarcodeCorrector:
    """
    Correct raw barcodes to a whitelist with substitution-only Hamming distance.
    """

    def __init__(self, whitelist: Iterable[str], alphabet: Sequence[str] = ALPHABET):
        cleaned = [clean_barcode(b) for b in whitelist if isinstance(b, str)]
        self.whitelist: Set[str] = {b for b in cleaned if b}
        if not self.whitelist:
            raise ValueError("Whitelist is empty; cannot perform barcode correction.")
        self.alphabet = tuple(alphabet)
        self.lengths = {len(b) for b in self.whitelist}

    def correct_one(
        self, raw: str, max_distance: int = 2
    ) -> Tuple[Optional[str], Optional[int], bool]:
        candidate = clean_barcode(raw)
        if not candidate or len(candidate) not in self.lengths:
            return (None, None, False)
        if candidate in self.whitelist:
            return (candidate, 0, False)
        if max_distance >= 1:
            matches = self._neighbors(candidate, 1)
            if len(matches) == 1:
                return (matches[0], 1, False)
            if len(matches) > 1:
                return (None, 1, True)
        if max_distance >= 2:
            matches = self._neighbors(candidate, 2)
            if len(matches) == 1:
                return (matches[0], 2, False)
            if len(matches) > 1:
                return (None, 2, True)
        return (None, None, False)

    def correct_many(
        self,
        raws: Iterable[str],
        max_distance: int = 2,
        progress_callback: Optional[Callable[[int], None]] = None,
        progress_step: int = 5000,
    ) -> pd.DataFrame:
        rows = []
        seen: Set[str] = set()
        since_callback = 0
        for raw in raws:
            barcode = clean_barcode(raw)
            if not barcode or barcode in seen:
                continue
            seen.add(barcode)
            corrected, dist, ambiguous = self.correct_one(barcode, max_distance=max_distance)
            rows.append((barcode, corrected, dist, ambiguous))
            if progress_callback is not None:
                since_callback += 1
                if since_callback >= progress_step:
                    progress_callback(since_callback)
                    since_callback = 0
        if progress_callback is not None and since_callback:
            progress_callback(since_callback)
        return pd.DataFrame(rows, columns=["raw", "corrected", "distance", "ambiguous"])

    def _neighbors(self, barcode: str, distance: int) -> List[str]:
        matches: List[str] = []
        if distance == 1:
            for i, base in enumerate(barcode):
                for alt in self.alphabet:
                    if alt == base:
                        continue
                    candidate = barcode[:i] + alt + barcode[i + 1 :]
                    if candidate in self.whitelist:
                        matches.append(candidate)
        elif distance == 2:
            length = len(barcode)
            for i in range(length - 1):
                for j in range(i + 1, length):
                    for alt_i in self.alphabet:
                        if alt_i == barcode[i]:
                            continue
                        for alt_j in self.alphabet:
                            if alt_j == barcode[j]:
                                continue
                            candidate = (
                                barcode[:i]
                                + alt_i
                                + barcode[i + 1 : j]
                                + alt_j
                                + barcode[j + 1 :]
                            )
                            if candidate in self.whitelist:
                                matches.append(candidate)
        return matches


def load_raw_matches(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"cb", "umi", "transcript_id"}
    if not required.issubset(df.columns):
        raise KeyError(f"Raw match table missing columns: {required - set(df.columns)}")
    df = df.copy()
    df["cb"] = df["cb"].map(clean_barcode)
    return df


def summarize_correction(mapping_df: pd.DataFrame) -> pd.DataFrame:
    if mapping_df.empty:
        return pd.DataFrame(columns=["category", "count", "fraction"])

    def bucket(row: pd.Series) -> str:
        if row["ambiguous"]:
            return "ambiguous"
        if row["distance"] == 0:
            return "exact(0)"
        if row["distance"] == 1:
            return "1-mismatch"
        if row["distance"] == 2:
            return "2-mismatch"
        return "unmatched"

    summary = (
        mapping_df.assign(category=lambda d: d.apply(bucket, axis=1))
        .value_counts("category")
        .rename("count")
        .reset_index()
    )
    summary["fraction"] = summary["count"] / len(mapping_df)
    return summary


def compute_counts_per_isoform_celltype(
    raw_df: pd.DataFrame,
    whitelist_df: pd.DataFrame,
    max_distance: int = 2,
    progress_callback: Optional[Callable[[int], None]] = None,
    progress_step: int = 5000,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if whitelist_df.empty:
        raise ValueError("Whitelist dataframe is empty; cannot assign cell types.")
    corrector = BarcodeCorrector(whitelist_df["barcode"].tolist())
    mapping = corrector.correct_many(
        raw_df["cb"].tolist(),
        max_distance=max_distance,
        progress_callback=progress_callback,
        progress_step=progress_step,
    )
    merged = raw_df.merge(
        mapping.rename(columns={"raw": "cb", "corrected": "cb_corrected"}),
        on="cb",
        how="left",
    )
    usable = merged[
        (merged["ambiguous"] == False)
        & (merged["cb_corrected"].notna())
        & merged["cb_corrected"].ne("")
    ].copy()
    usable = usable.merge(
        whitelist_df.rename(columns={"barcode": "cb_corrected"}), on="cb_corrected", how="left"
    )
    grouped = usable.groupby(["transcript_id", "Cell Type", "Age"])
    counts = (
        grouped.agg(
            n_cells=("cb_corrected", "nunique"),
            n_molecules=("umi", "nunique"),
        )
        .reset_index()
        .sort_values(["transcript_id", "Cell Type", "Age"])
        .reset_index(drop=True)
    )
    summary = summarize_correction(mapping)
    return merged, counts, summary


def _import_scanpy_modules():
    try:
        import scanpy as sc  # type: ignore
        import scanpy.external as sce  # type: ignore
    except ImportError as exc:  # pragma: no cover - import guard
        raise ImportError(
            "run_scanpy_analysis requires scanpy and scanpy.external to be installed."
        ) from exc
    return sc, sce


def _load_kallisto_sample(sample_dir: Path, sc_module) -> "AnnData":
    mtx_path = sample_dir / "genes.mtx"
    barcodes_path = sample_dir / "genes.barcodes.txt"
    genes_path = sample_dir / "genes.genes.txt"
    missing = [str(path) for path in (mtx_path, barcodes_path, genes_path) if not path.exists()]
    if missing:
        raise FileNotFoundError(
            f"Kallisto sample {sample_dir} is missing required files: {', '.join(missing)}"
        )
    adata = sc_module.read_mtx(str(mtx_path))
    barcodes = pd.read_csv(barcodes_path, header=None)[0].astype(str)
    genes = pd.read_csv(genes_path, header=None)[0].astype(str)
    obs_matches = adata.n_obs == len(barcodes) and adata.n_vars == len(genes)
    var_matches = adata.n_obs == len(genes) and adata.n_vars == len(barcodes)
    if not (obs_matches or var_matches):
        raise ValueError(
            f"Dimension mismatch in sample {sample_dir}: "
            f"{adata.n_obs} obs vs {len(barcodes)} barcodes, "
            f"{adata.n_vars} vars vs {len(genes)} genes."
        )
    if var_matches and not obs_matches:
        adata = adata.transpose().copy()
    adata.obs_names = barcodes
    adata.var_names = genes
    adata.var_names_make_unique()
    return adata


def _compute_knee_cutoff(counts: np.ndarray, expected_cells: int) -> float:
    if expected_cells <= 0:
        raise ValueError("expected_cells must be > 0 for knee cutoff calculation.")
    if expected_cells >= counts.size:
        raise ValueError(
            f"Expected cell index {expected_cells} is >= total barcodes ({counts.size})."
        )
    sorted_counts = np.sort(counts)[::-1]
    return float(sorted_counts[expected_cells])


def _filter_cells_by_counts(adata: "AnnData", min_counts: float, sc_module) -> None:
    sc_module.pp.filter_cells(adata, min_counts=max(1, int(math.floor(min_counts))))
    sc_module.pp.filter_genes(adata, min_cells=0)


def _sort_celltype_age_index(values: Sequence[str], ages: Sequence[str]) -> List[str]:
    age_order = list(dict.fromkeys(ages))

    def key(entry: str) -> Tuple[str, int]:
        if "_" in entry:
            cell_type, age = entry.rsplit("_", 1)
            try:
                idx = age_order.index(age)
            except ValueError:
                idx = len(age_order)
            return (cell_type, idx)
        return (entry, len(age_order))

    return sorted(values, key=key)



def _filter_by_ids(df: pd.DataFrame, transcripts=None, gene_ids=None) -> pd.DataFrame:
    data = df.copy()
    if transcripts and "transcript_id" in data.columns:
        data = data[data["transcript_id"].isin(transcripts)]
    if gene_ids and "gene_id" in data.columns:
        data = data[data["gene_id"].isin(gene_ids)]
    return data


def plot_isoform_fraction_heatmap(
    ax,
    combined_counts: pd.DataFrame,
    transcripts=None,
    gene_ids=None,
    sorter: Optional[Callable[[Sequence[str]], Sequence[str]]] = None,
    value_col: str = "n_cells",
    cmap: str = "viridis",
    title: str = "Isoform fraction by cell type/age",
    ylabel: str = "Cell Type / Age",
    show_y_ticks: bool = True,
    empty_message: str = "No isoform counts available",
) -> bool:
    if combined_counts is None or combined_counts.empty:
        ax.text(0.5, 0.5, empty_message, ha="center", va="center")
        ax.axis("off")
        return False
    data = _filter_by_ids(combined_counts, transcripts, gene_ids)
    if data.empty:
        ax.text(0.5, 0.5, "Requested transcripts missing", ha="center", va="center")
        ax.axis("off")
        return False
    if value_col not in data.columns:
        value_col = "n_cells"
    data = data.copy()
    data["Age"] = data["Age"].astype(str)
    data["Cell Type"] = data["Cell Type"].astype(str)
    totals = (
        data.groupby(["Age"], observed=False)[value_col]
        .sum()
        .rename("total_value")
        .reset_index()
    )
    data = data.merge(totals, on="Age", how="left")
    data["fraction"] = data[value_col] / data["total_value"].replace(0, np.nan)
    data["Cell Type_Age"] = data["Cell Type"] + "_" + data["Age"]
    heat = data.pivot_table(
        index="Cell Type_Age",
        columns="transcript_id",
        values="fraction",
        fill_value=0.0,
        observed=False,
    )
    if heat.empty:
        ax.text(0.5, 0.5, "No isoform fractions to plot", ha="center", va="center")
        ax.axis("off")
        return False
    ordered_index = sorter(heat.index.tolist()) if sorter else heat.index.tolist()
    heat = heat.loc[ordered_index]
    sns.heatmap(
        heat,
        cmap=cmap,
        ax=ax,
        vmin=float(heat.values.min()),
        vmax=float(heat.values.max()),
        cbar_kws={"label": "Fraction of cells expressing transcript"},
    )
    ax.set_xlabel("Transcripts")
    ax.set_ylabel(ylabel if show_y_ticks else "")
    ax.set_title(title)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    if not show_y_ticks:
        ax.tick_params(axis="y", left=False, labelleft=False)
    return True


def plot_isoform_ratio_heatmap(
    ax,
    ratio_df: Optional[pd.DataFrame] = None,
    combined_counts: Optional[pd.DataFrame] = None,
    transcripts=None,
    gene_ids=None,
    sorter: Optional[Callable[[Sequence[str]], Sequence[str]]] = None,
    cmap: str = "magma",
    title: str = "Isoform ratio by cell type/age",
    show_y_ticks: bool = False,
    empty_message: str = "No isoform ratio data",
) -> bool:
    ratio_source = pd.DataFrame()
    if ratio_df is not None and not ratio_df.empty:
        ratio_source = ratio_df.copy()
    elif combined_counts is not None and not combined_counts.empty:
        ratio_source = _compute_isoform_group_ratios(combined_counts)
    if ratio_source.empty:
        ax.text(0.5, 0.5, empty_message, ha="center", va="center")
        ax.axis("off")
        return False
    ratio_source = _filter_by_ids(ratio_source, transcripts, gene_ids)
    if ratio_source.empty:
        ax.text(0.5, 0.5, "Requested transcripts missing", ha="center", va="center")
        ax.axis("off")
        return False
    ratio_source = ratio_source.copy()
    ratio_source["Age"] = ratio_source["Age"].astype(str)
    ratio_source["Cell Type"] = ratio_source["Cell Type"].astype(str)
    if "isoform_ratio" not in ratio_source.columns:
        ax.text(0.5, 0.5, "Ratio column unavailable", ha="center", va="center")
        ax.axis("off")
        return False
    ratio_source["Cell Type_Age"] = ratio_source["Cell Type"] + "_" + ratio_source["Age"]
    heat = ratio_source.pivot_table(
        index="Cell Type_Age",
        columns="transcript_id",
        values="isoform_ratio",
        fill_value=0.0,
        observed=False,
    )
    if heat.empty:
        ax.text(0.5, 0.5, "No isoform ratios to plot", ha="center", va="center")
        ax.axis("off")
        return False
    ordered_index = sorter(heat.index.tolist()) if sorter else heat.index.tolist()
    heat = heat.loc[ordered_index]
    sns.heatmap(
        heat,
        cmap=cmap,
        ax=ax,
        vmin=0,
        vmax=1,
        cbar_kws={"label": "Isoform ratio (molecule share)"},
    )
    ax.set_xlabel("Transcripts")
    if show_y_ticks:
        ax.set_ylabel("Cell Type / Age")
    else:
        ax.set_ylabel("")
        ax.tick_params(axis="y", left=False, labelleft=False)
    ax.set_title(title)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    return True


def plot_promoter_group_heatmap(
    ax,
    combined_counts: pd.DataFrame,
    promoter_meta: Optional[pd.DataFrame],
    sorter: Optional[Callable[[Sequence[str]], Sequence[str]]] = None,
    transcripts=None,
    gene_ids=None,
    cmap: str = "viridis",
    title: str = "All isoforms (scaled per promoter group)",
    empty_message: str = "Isoform metadata unavailable",
) -> bool:
    if combined_counts is None or combined_counts.empty:
        ax.text(0.5, 0.5, empty_message, ha="center", va="center")
        ax.axis("off")
        return False
    if promoter_meta is None or promoter_meta.empty:
        ax.text(0.5, 0.5, "No promoter group annotations", ha="center", va="center")
        ax.axis("off")
        return False
    promoter_meta = (
        promoter_meta[["transcript_id", "promoter_group"]]
        .dropna(subset=["transcript_id", "promoter_group"])
        .drop_duplicates()
    )
    if promoter_meta.empty:
        ax.text(0.5, 0.5, "No promoter group annotations", ha="center", va="center")
        ax.axis("off")
        return False
    data = combined_counts.merge(promoter_meta, on="transcript_id", how="inner")
    data = _filter_by_ids(data, transcripts, gene_ids)
    if data.empty:
        ax.text(0.5, 0.5, "No overlap between counts and promoter groups", ha="center", va="center")
        ax.axis("off")
        return False
    value_col = "n_molecules" if "n_molecules" in data.columns and data["n_molecules"].notna().any() else "n_cells"
    data = data.copy()
    data["Age"] = data["Age"].astype(str)
    data["Cell Type"] = data["Cell Type"].astype(str)
    totals = (
        data.groupby(["Cell Type", "Age", "promoter_group"], observed=False)[value_col]
        .sum()
        .rename("group_total")
        .reset_index()
    )
    data = data.merge(
        totals,
        on=["Cell Type", "Age", "promoter_group"],
        how="left",
    )
    data["ratio_within_group"] = data[value_col] / data["group_total"].replace(0, np.nan)
    data["Cell Type_Age"] = data["Cell Type"] + "_" + data["Age"]
    heat = data.pivot_table(
        index="Cell Type_Age",
        columns="transcript_id",
        values="ratio_within_group",
        fill_value=0.0,
        observed=False,
    )
    if heat.empty:
        ax.text(0.5, 0.5, "No isoform ratios to plot", ha="center", va="center")
        ax.axis("off")
        return False
    ordered_index = sorter(heat.index.tolist()) if sorter else heat.index.tolist()
    heat = heat.loc[ordered_index]
    transcript_groups = promoter_meta.set_index("transcript_id")["promoter_group"].to_dict()
    heat = heat[[col for col in heat.columns if col in transcript_groups]]
    if heat.empty:
        ax.text(0.5, 0.5, "No transcripts with promoter groups available", ha="center", va="center")
        ax.axis("off")
        return False
    def _pg_sort_key(col):
        value = transcript_groups[col]
        try:
            return (float(value), col)
        except (TypeError, ValueError):
            return (float("inf"), col)
    ordered_cols = sorted(heat.columns, key=_pg_sort_key)
    heat = heat[ordered_cols]
    scaled_values = heat.values.copy()
    boundaries = []
    prev_group = None
    start = 0
    for idx, col in enumerate(ordered_cols):
        group = transcript_groups[col]
        if group != prev_group:
            if prev_group is not None:
                boundaries.append((prev_group, start, idx))
            start = idx
            prev_group = group
    if prev_group is not None:
        boundaries.append((prev_group, start, len(ordered_cols)))
    for group, start, end in boundaries:
        block = scaled_values[:, start:end]
        max_val = np.nanmax(block)
        if not np.isfinite(max_val) or max_val == 0:
            scaled_values[:, start:end] = 0
        else:
            scaled_values[:, start:end] = block / max_val
    im = ax.imshow(
        scaled_values,
        aspect="auto",
        cmap=cmap,
        vmin=0,
        vmax=1,
    )
    ax.set_yticks(np.arange(len(heat.index)))
    ax.set_yticklabels(heat.index, fontsize=7)
    ax.set_xlabel("Transcripts ordered by promoter group")
    ax.set_ylabel("Cell Type / Age")
    ax.set_title(title)
    tick_positions = []
    tick_labels = []
    for group, start, end in boundaries:
        midpoint = (start + end - 1) / 2
        tick_positions.append(midpoint)
        try:
            value = float(group)
            label = f"PG {int(value)}" if value.is_integer() else f"PG {value}"
        except (TypeError, ValueError):
            label = f"PG {group}"
        tick_labels.append(f"{label}\n(n={end - start})")
        ax.axvline(start - 0.5, color="white", linewidth=0.3)
    ax.axvline(len(ordered_cols) - 0.5, color="white", linewidth=0.3)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=0, fontsize=7)
    fig = ax.get_figure()
    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cbar.set_label("Relative fraction (scaled within promoter group)")
    return True
def _build_transcript_gene_promoter_map(
    gtf_path: Optional[Path], transcript_ids: Iterable[str]
) -> pd.DataFrame:
    """
    Load gene_id and promoter_group assignments for the requested transcripts from a GTF file.
    """
    columns = ["transcript_id", "gene_id", "promoter_group"]
    transcripts = {str(tid).strip() for tid in transcript_ids if str(tid).strip()}
    if not transcripts:
        return pd.DataFrame(columns=columns)
    if gtf_path is None or not gtf_path.exists():
        print(
            f"[SupFig4][scanpy] Warning: annotation GTF path not found ({gtf_path}); "
            "gene/promoter annotations unavailable."
        )
        return pd.DataFrame(columns=columns)
    subset = load_gtf_subset(gtf_path, transcripts, usecols=DEFAULT_GTF_COLUMNS)
    if subset.empty:
        print(
            "[SupFig4][scanpy] Warning: requested transcripts missing from annotation GTF; "
            "gene/promoter annotations unavailable."
        )
        return pd.DataFrame(columns=columns)
    mapping = (
        subset.loc[:, ["transcript_id", "gene_id", "promoter_group"]]
        .dropna(subset=["transcript_id"])
        .assign(
            transcript_id=lambda df: df["transcript_id"].astype(str).str.strip(),
            gene_id=lambda df: df["gene_id"].astype(str).str.strip(),
            promoter_group=lambda df: df["promoter_group"].astype(str).str.strip(),
        )
        .drop_duplicates(subset=["transcript_id"])
    )
    mapping = mapping[mapping["transcript_id"] != ""]
    return mapping.reset_index(drop=True)


def _compile_splice_annotations(
    splice_paths: Sequence[Path],
    gtf_path: Optional[Path],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Merge multiple splice result tables into a unique transcript annotation sheet and
    return both the merged table and the transcript → gene/promoter lookup extracted from the GTF.
    """
    tables: List[pd.DataFrame] = []
    transcript_ids: Set[str] = set()
    for path in splice_paths:
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        if "ids" not in frame.columns:
            raise KeyError(f"Splice result {path} missing 'ids' column.")
        frame = frame.copy()
        frame["transcript_id"] = frame["ids"].astype(str).str.split("_").str[0]
        frame["splice_source"] = path.name
        frame = frame.dropna(subset=["transcript_id"])
        frame = frame[frame["transcript_id"].ne("")]
        if frame.empty:
            continue
        transcript_ids.update(frame["transcript_id"].astype(str))
        tables.append(frame)

    if not tables:
        return pd.DataFrame(), pd.DataFrame(columns=["transcript_id", "gene_id", "promoter_group"])

    combined = pd.concat(tables, ignore_index=True)
    combined = combined.drop_duplicates(subset=["transcript_id"])
    gene_hint = (
        combined["ids"].astype(str).str.extract(r"_(ENSMUSG\d+)$").rename(columns={0: "gene_id_hint"})
    )
    combined = pd.concat([combined, gene_hint], axis=1)

    annotation_lookup = _build_transcript_gene_promoter_map(gtf_path, transcript_ids)
    if not annotation_lookup.empty:
        combined = combined.merge(
            annotation_lookup, on="transcript_id", how="left", suffixes=("", "_annotation")
        )
    else:
        combined["gene_id"] = pd.NA
        combined["promoter_group"] = pd.NA

    combined["gene_id"] = combined["gene_id"].fillna(combined.get("gene_id_hint"))
    combined = combined.drop(columns=["gene_id_hint"], errors="ignore")
    return combined, annotation_lookup


def _filter_gene_promoter_min_transcripts(
    df: pd.DataFrame, min_transcripts: int = 2
) -> pd.DataFrame:
    """
    Keep only rows whose (gene_id, promoter_group) combination has >= min_transcripts isoforms.
    """
    required = {"transcript_id", "gene_id", "promoter_group"}
    if df.empty or not required.issubset(df.columns):
        return df

    working = df.copy()

    def _normalize(series: pd.Series) -> pd.Series:
        normalized = series.astype(str).str.strip()
        normalized = normalized.where(series.notna(), "")
        normalized = normalized.replace({"nan": "", "None": ""})
        return normalized

    working["gene_id_norm"] = _normalize(working["gene_id"])
    working["promoter_group_norm"] = _normalize(working["promoter_group"])
    mask = (working["gene_id_norm"] != "") & (working["promoter_group_norm"] != "")
    if not mask.any():
        return df.iloc[0:0].copy()

    combo_counts = (
        working.loc[mask]
        .groupby(["gene_id_norm", "promoter_group_norm"])["transcript_id"]
        .nunique()
        .reset_index(name="transcript_count")
    )
    valid_pairs = combo_counts[combo_counts["transcript_count"] >= min_transcripts]
    if valid_pairs.empty:
        return df.iloc[0:0].copy()

    filtered = working.merge(
        valid_pairs[["gene_id_norm", "promoter_group_norm"]],
        on=["gene_id_norm", "promoter_group_norm"],
        how="inner",
    )
    filtered = filtered[df.columns]
    return filtered.reset_index(drop=True)


def _compute_isoform_group_ratios(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate isoform ratios within each gene/promoter group per cell type and age.
    """
    required = {"gene_id", "promoter_group", "transcript_id", "Cell Type", "Age"}
    if counts.empty or not required.issubset(counts.columns):
        return pd.DataFrame(columns=list(required) + ["group_total", "isoform_ratio"])

    working = counts.copy()
    value_column = "n_molecules" if "n_molecules" in working.columns else "n_cells"
    if value_column not in working.columns:
        return pd.DataFrame(columns=list(required) + ["group_total", "isoform_ratio"])

    for col in ("gene_id", "promoter_group", "Cell Type", "Age"):
        working[col] = working[col].astype(str).str.strip()
    working[value_column] = working[value_column].fillna(0).astype(float)

    totals = (
        working.groupby(["gene_id", "promoter_group", "Cell Type", "Age"], observed=False)[
            value_column
        ]
        .sum()
        .rename("group_total")
        .reset_index()
    )
    ratios = working.merge(
        totals, on=["gene_id", "promoter_group", "Cell Type", "Age"], how="left"
    )
    ratios["isoform_ratio"] = ratios[value_column] / ratios["group_total"].replace(0, np.nan)
    return ratios


def _normalize_dataframe_columns(df: pd.DataFrame) -> pd.DataFrame:
    data = df.copy()
    data.columns = [c.strip() for c in data.columns]
    return data


def _ensure_transcript_id_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Guarantee a transcript_id column exists, deriving it from an 'ids' column if needed.
    """
    data = df.copy()
    if "transcript_id" in data.columns:
        return data
    if "ids" not in data.columns:
        raise ValueError("Splice result is missing both 'transcript_id' and 'ids' columns.")
    data["transcript_id"] = (
        data["ids"]
        .astype(str)
        .where(data["ids"].notna())
        .str.split("_", n=1)
        .str[0]
    )
    return data


def _fold_change_direction(df: pd.DataFrame, column: str = "fold_change_mean") -> pd.DataFrame:
    """
    Map fold_change_mean to categorical directions:
      >1 => 'up' (AMA higher)
      <1 => 'down' (AMA lower)
      otherwise => 'tie'
    """
    data = df.copy()
    data[column] = pd.to_numeric(data[column], errors="coerce")

    def _dir_map(value: float) -> str:
        if pd.isna(value):
            return "tie"
        if value > 1:
            return "up"
        if value < 1:
            return "down"
        return "tie"

    data["fc_dir"] = data[column].apply(_dir_map)
    return data


def _collapse_transcript_direction(df: pd.DataFrame, dir_col: str = "fc_dir") -> pd.DataFrame:
    """
    Collapse per-row splice directions to a single consistent direction per transcript.
    """
    keep = (
        df[["transcript_id", dir_col]]
        .dropna(subset=["transcript_id"])
        .query(f"{dir_col} in ['up', 'down']")
    )
    if keep.empty:
        return pd.DataFrame(columns=["transcript_id", "dir"])
    collapsed = (
        keep.groupby("transcript_id")[dir_col]
        .agg(lambda vals: vals.iloc[0] if vals.nunique() == 1 else None)
        .dropna()
        .reset_index()
        .rename(columns={dir_col: "dir"})
    )
    return collapsed


def _build_same_direction_transcripts(a_df: pd.DataFrame, b_df: pd.DataFrame) -> pd.DataFrame:
    """
    Return transcripts that share the same splice direction in both contrasts.
    """
    a_dirs = _collapse_transcript_direction(a_df)
    b_dirs = _collapse_transcript_direction(b_df)
    merged = a_dirs.merge(b_dirs, on="transcript_id", suffixes=("_a", "_b"))
    merged = merged[merged["dir_a"] == merged["dir_b"]].copy()
    merged = merged.rename(columns={"dir_a": "shared_dir"}).drop(columns=["dir_b"])
    if merged.empty:
        merged["fc_a"] = pd.Series(dtype=float)
        merged["fc_b"] = pd.Series(dtype=float)
        return merged
    a_fc = (
        a_df.groupby("transcript_id")["fold_change_mean"]
        .median()
        .rename("fc_a")
        .reset_index()
    )
    b_fc = (
        b_df.groupby("transcript_id")["fold_change_mean"]
        .median()
        .rename("fc_b")
        .reset_index()
    )
    merged = merged.merge(a_fc, on="transcript_id", how="left").merge(b_fc, on="transcript_id", how="left")
    return merged


def _build_isoform_pairs(iso_df: pd.DataFrame, transcript_filter: Set[str]) -> pd.DataFrame:
    """
    Create AMA vs YNG isoform ratio pairs per (transcript, cell type, promoter).
    """
    required = {"transcript_id", "promoter_group", "Cell Type", "Age", "isoform_ratio"}
    missing = required - set(iso_df.columns)
    if missing:
        raise ValueError(f"Isoform ratio table is missing required columns: {missing}")
    if not transcript_filter:
        return pd.DataFrame(
            columns=["transcript_id", "Cell Type", "promoter_group", "AMA", "YNG", "iso_dir"]
        )

    data = iso_df.copy()
    data["Age"] = data["Age"].astype(str).str.upper().str.strip()
    data = data[data["transcript_id"].isin(transcript_filter)]
    if data.empty:
        return pd.DataFrame(
            columns=["transcript_id", "Cell Type", "promoter_group", "AMA", "YNG", "iso_dir"]
        )

    wide = (
        data.pivot_table(
            index=["transcript_id", "Cell Type", "promoter_group"],
            columns="Age",
            values="isoform_ratio",
            aggfunc="first",
        )
        .reset_index()
    )

    for age in ("AMA", "YNG"):
        if age not in wide.columns:
            wide[age] = np.nan
    wide = wide.dropna(subset=["AMA", "YNG"])
    if wide.empty:
        wide["iso_dir"] = pd.Series(dtype="object")
        return wide

    def _iso_dir(row: pd.Series) -> str:
        if row["AMA"] > row["YNG"]:
            return "up"
        if row["AMA"] < row["YNG"]:
            return "down"
        return "tie"

    wide["iso_dir"] = wide.apply(_iso_dir, axis=1)
    return wide


def generate_isoform_vs_splice_summary(
    isoform_ratio_path: Path,
    splice_a_path: Path,
    splice_b_path: Path,
    output_dir: Path,
) -> Dict[str, object]:
    """
    Compare isoform ratios against two splice contrasts and write validation summaries.

    Parameters
    ----------
    isoform_ratio_path
        CSV with per-(transcript, Cell Type, Age, promoter_group) isoform ratios.
    splice_a_path
        CSV for the a_vs_c splice contrast (YNG vs AMA).
    splice_b_path
        CSV for the b_vs_d splice contrast (YNG vs AMA).
    output_dir
        Directory where summary CSV/JSON artifacts will be written.

    Returns
    -------
    dict
        Dictionary mirroring the JSON summary.
    """
    isoform_ratio_path = Path(isoform_ratio_path)
    splice_a_path = Path(splice_a_path)
    splice_b_path = Path(splice_b_path)
    output_dir = Path(output_dir)

    missing = [p for p in (isoform_ratio_path, splice_a_path, splice_b_path) if not p.exists()]
    if missing:
        missing_str = ", ".join(str(p) for p in missing)
        raise FileNotFoundError(f"Required input files not found: {missing_str}")

    output_dir.mkdir(parents=True, exist_ok=True)

    iso = _normalize_dataframe_columns(pd.read_csv(isoform_ratio_path))
    splice_a = _normalize_dataframe_columns(pd.read_csv(splice_a_path))
    splice_b = _normalize_dataframe_columns(pd.read_csv(splice_b_path))

    splice_a = _fold_change_direction(_ensure_transcript_id_column(splice_a))
    splice_b = _fold_change_direction(_ensure_transcript_id_column(splice_b))

    same_direction = _build_same_direction_transcripts(splice_a, splice_b)

    iso_unique_tx_total = int(iso["transcript_id"].dropna().nunique())
    splice_overlap = int(
        pd.merge(
            splice_a[["transcript_id"]].dropna().drop_duplicates(),
            splice_b[["transcript_id"]].dropna().drop_duplicates(),
            on="transcript_id",
            how="inner",
        )["transcript_id"].nunique()
    )
    same_direction_count = int(same_direction["transcript_id"].nunique())
    same_up = int(same_direction[same_direction["shared_dir"] == "up"]["transcript_id"].nunique())
    same_down = int(same_direction[same_direction["shared_dir"] == "down"]["transcript_id"].nunique())

    same_dir_set = set(same_direction["transcript_id"])
    iso_dir_set = set(iso["transcript_id"].dropna().unique())
    iso_presence = len(same_dir_set & iso_dir_set)

    iso_pairs = _build_isoform_pairs(iso, same_dir_set)
    if not iso_pairs.empty:
        iso_pairs = iso_pairs.merge(
            same_direction[["transcript_id", "shared_dir"]].drop_duplicates(),
            on="transcript_id",
            how="left",
        )
        iso_pairs["match_splice"] = iso_pairs["iso_dir"] == iso_pairs["shared_dir"]
    else:
        iso_pairs["shared_dir"] = pd.Series(dtype="object")
        iso_pairs["match_splice"] = pd.Series(dtype=bool)

    pairs_total = int(len(iso_pairs))
    pairs_match = int(iso_pairs["match_splice"].sum()) if "match_splice" in iso_pairs else 0

    per_tx = pd.DataFrame(columns=["transcript_id", "n_pairs", "any_match", "shared_dir"])
    if not iso_pairs.empty:
        per_tx = (
            iso_pairs.groupby("transcript_id")
            .agg(n_pairs=("iso_dir", "size"), any_match=("match_splice", "any"))
            .reset_index()
            .merge(
                same_direction[["transcript_id", "shared_dir"]].drop_duplicates(),
                on="transcript_id",
                how="left",
            )
        )

    transcripts_with_pairs = int(per_tx["transcript_id"].nunique()) if not per_tx.empty else 0
    transcripts_validated = int(per_tx["any_match"].sum()) if not per_tx.empty else 0

    up_group = per_tx[per_tx["shared_dir"] == "up"] if not per_tx.empty else pd.DataFrame()
    down_group = per_tx[per_tx["shared_dir"] == "down"] if not per_tx.empty else pd.DataFrame()
    up_with_pairs = int(up_group["transcript_id"].nunique()) if not up_group.empty else 0
    down_with_pairs = int(down_group["transcript_id"].nunique()) if not down_group.empty else 0
    up_validated = int(up_group["any_match"].sum()) if not up_group.empty else 0
    down_validated = int(down_group["any_match"].sum()) if not down_group.empty else 0

    summary = {
        "isoform_table_unique_transcripts_total": iso_unique_tx_total,
        "splice_overlap_unique_transcripts": splice_overlap,
        "same_direction_unique_transcripts": same_direction_count,
        "same_direction_up": same_up,
        "same_direction_down": same_down,
        "isoform_any_presence_within_same_direction": iso_presence,
        "transcripts_with_both_ages_in_any_celltype_promoter": transcripts_with_pairs,
        "validated_transcripts_any_celltype_promoter_match": transcripts_validated,
        "pairs_total_both_ages": pairs_total,
        "pairs_matching_splice_direction": pairs_match,
        "up_group": {
            "same_direction_up": same_up,
            "with_pairs": up_with_pairs,
            "validated": up_validated,
        },
        "down_group": {
            "same_direction_down": same_down,
            "with_pairs": down_with_pairs,
            "validated": down_validated,
        },
    }

    same_direction[["transcript_id", "shared_dir", "fc_a", "fc_b"]].drop_duplicates().to_csv(
        output_dir / "same_direction_transcripts.csv", index=False
    )
    iso_pairs.to_csv(output_dir / "isoform_pairwise_by_tx_celltype_promoter.csv", index=False)
    per_tx.to_csv(output_dir / "per_transcript_iso_match_summary.csv", index=False)

    summary_json = output_dir / "summary_stats_isoform_vs_splice.json"
    summary_csv = output_dir / "isoform_vs_splice_summary_stats.csv"
    summary_json.write_text(json.dumps(summary, indent=2))
    pd.DataFrame([summary]).to_csv(summary_csv, index=False)

    return summary


def run_scanpy_analysis(config: ScanpyAnalysisConfig) -> ScanpyAnalysisOutputs:
    """
    Execute the Scanpy + isoform counting workflow for Supplemental Figure 4.
    """

    use_cached_adata = (
        config.adata_full_path is not None
        and config.adata_hvg_path is not None
        and config.adata_full_path.exists()
        and config.adata_hvg_path.exists()
    )

    print("[SupFig4][scanpy] Importing Scanpy modules...")
    sc, sce = _import_scanpy_modules()
    _ensure_anndata_null_reader()

    adata_full: "AnnData"
    adata_hvg: "AnnData"
    age_order: List[str]

    if use_cached_adata:
        print(
            f"[SupFig4][scanpy] Reusing cached AnnData objects from "
            f"{config.adata_full_path} / {config.adata_hvg_path}"
        )
        adata_full = sc.read_h5ad(config.adata_full_path)
        adata_hvg = sc.read_h5ad(config.adata_hvg_path)
        if "Age" not in adata_full.obs:
            raise KeyError("Cached AnnData is missing 'Age' in obs.")
        age_series = adata_full.obs["Age"]
        if not isinstance(age_series.dtype, CategoricalDtype):
            ordered = list(dict.fromkeys(age_series.astype(str)))
            adata_full.obs["Age"] = pd.Categorical(age_series.astype(str), categories=ordered, ordered=True)
        age_order = list(adata_full.obs["Age"].cat.categories)
        if "leiden" not in adata_full.obs and "leiden" in adata_hvg.obs:
            adata_full.obs["leiden"] = (
                adata_hvg.obs["leiden"].reindex(adata_full.obs_names).astype("category")
            )
        if "leiden" not in adata_full.obs:
            raise KeyError("Cached AnnData is missing 'leiden' clusters.")
        if "X_umap" not in adata_full.obsm and "X_umap" in adata_hvg.obsm:
            adata_full.obsm["X_umap"] = adata_hvg.obsm["X_umap"].copy()
        print("[SupFig4][scanpy] Cached AnnData loaded successfully.")
    else:
        if not config.kallisto_dir.exists():
            raise FileNotFoundError(
                f"{config.kallisto_dir} is missing and no cached AnnData files were found."
            )
        sample_dirs = sorted(
            [path for path in config.kallisto_dir.iterdir() if path.is_dir()],
            key=lambda p: p.name,
        )
        if not sample_dirs:
            raise FileNotFoundError(
                f"No kallisto sample directories found under {config.kallisto_dir}"
            )

        adata_by_sample: Dict[str, "AnnData"] = {}
        print(f"[SupFig4][scanpy] Loading {len(sample_dirs)} kallisto samples...")
        for sample_dir in sample_dirs:
            sample_name = sample_dir.name
            if sample_name not in config.expected_cells:
                raise KeyError(
                    f"No expected cell count provided for sample '{sample_name}'. "
                    "Update ScanpyAnalysisConfig.expected_cells."
                )
            print(f"[SupFig4][scanpy]  - Reading {sample_name} matrices...")
            adata = _load_kallisto_sample(sample_dir, sc)
            counts = np.asarray(adata.X.sum(axis=1)).ravel()
            cutoff = _compute_knee_cutoff(counts, int(config.expected_cells[sample_name]))
            _filter_cells_by_counts(adata, cutoff, sc)
            print(
                f"[SupFig4][scanpy]    retained {adata.n_obs} cells after knee cutoff (min {cutoff:.1f} counts)"
            )
            adata.obs["Sample"] = sample_name
            adata_by_sample[sample_name] = adata

        print("[SupFig4][scanpy] Concatenating samples...")
        adata_all = sc.concat(
            adata_by_sample.values(),
            label="Sample",
            keys=list(adata_by_sample.keys()),
            join="outer",
            merge=None,
        )
        adata_all.var_names_make_unique()
        print(
            f"[SupFig4][scanpy]  Combined matrix shape: {adata_all.n_obs} cells × {adata_all.n_vars} genes"
        )

        condition_map = {str(k): str(v) for k, v in config.condition_map.items()}
        adata_all.obs["Age"] = adata_all.obs["Sample"].map(condition_map)
        adata_all = adata_all[adata_all.obs["Age"].notna()].copy()
        age_order = []
        for age in condition_map.values():
            if age not in age_order:
                age_order.append(age)
        observed_ages = set(adata_all.obs["Age"].unique())
        age_order = [age for age in age_order if age in observed_ages]
        adata_all.obs["Age"] = pd.Categorical(adata_all.obs["Age"], categories=age_order, ordered=True)

        print("[SupFig4][scanpy] Normalizing + log1p transforming...")
        sc.pp.normalize_total(adata_all, target_sum=1e4)
        sc.pp.log1p(adata_all)
        adata_full = adata_all.copy()
        adata_full.raw = adata_full

        print("[SupFig4][scanpy] Selecting highly variable genes...")
        sc.pp.highly_variable_genes(adata_all, batch_key="Sample", n_top_genes=config.max_hvg_genes)
        hvg_mask = adata_all.var["highly_variable"].astype(bool)
        if not hvg_mask.any():
            raise RuntimeError("No highly variable genes detected; cannot continue Scanpy pipeline.")
        adata_hvg = adata_all[:, hvg_mask].copy()
        print(f"[SupFig4][scanpy]  HVGs retained: {adata_hvg.n_vars}")

        print("[SupFig4][scanpy] Scaling, PCA, Harmony integration, neighbors, UMAP, Leiden...")
        sc.pp.scale(adata_hvg, max_value=10)
        sc.tl.pca(adata_hvg, svd_solver="arpack")
        sce.pp.harmony_integrate(adata_hvg, key=config.harmony_key)
        sc.pp.neighbors(adata_hvg, use_rep="X_pca_harmony")
        sc.tl.umap(adata_hvg)
        sc.tl.leiden(adata_hvg, resolution=config.leiden_resolution)
        print("[SupFig4][scanpy]  Dimensionality reduction complete.")

        adata_full.obs["leiden"] = (
            adata_hvg.obs["leiden"].reindex(adata_full.obs_names).astype("category")
        )
        adata_full.obsm["X_umap"] = adata_hvg.obsm["X_umap"].copy()

        if config.adata_full_path is not None:
            config.adata_full_path.parent.mkdir(parents=True, exist_ok=True)
            print(f"[SupFig4][scanpy] Saving cached adata_full to {config.adata_full_path}...")
            adata_full.write(config.adata_full_path)
        if config.adata_hvg_path is not None:
            config.adata_hvg_path.parent.mkdir(parents=True, exist_ok=True)
            print(f"[SupFig4][scanpy] Saving cached adata_hvg to {config.adata_hvg_path}...")
            adata_hvg.write(config.adata_hvg_path)

    marker_source = "DEFAULT_MARKER_GENE_SETS"
    if config.marker_genes is not None:
        marker_sets = {
            str(cell_type): [str(gene) for gene in genes]
            for cell_type, genes in config.marker_genes.items()
        }
        marker_source = "config.marker_genes"
    elif config.marker_gene_path is not None:
        marker_sets = _load_marker_gene_sets_from_file(config.marker_gene_path)
        marker_source = f"file:{config.marker_gene_path}"
    else:
        marker_sets = DEFAULT_MARKER_GENE_SETS
    markers_unique, coverage_df, overlap_df = deduplicate_markers(
        marker_sets, adata_full.var_names, PRIORITY_ORDER
    )

    print(f"[SupFig4][scanpy] Marker gene source: {marker_source}")
    print("[SupFig4][scanpy] Marker coverage after harmonization (pre-dedup):")
    coverage_sorted = coverage_df.sort_values(
        ["n_present", "n_input"], ascending=[False, False]
    )
    print(coverage_sorted.to_string(max_rows=999))
    if overlap_df.empty:
        print("[SupFig4][scanpy] No overlaps detected across marker sets.")
    else:
        print("[SupFig4][scanpy] Overlapping genes removed (kept for first-in-priority type):")
        print(overlap_df.head(50).to_string(index=False))

    print("[SupFig4][scanpy] Scoring marker gene sets...")
    score_cols: List[str] = []
    for cell_type, genes in markers_unique.items():
        if not genes:
            adata_full.obs[cell_type] = np.nan
            continue
        sc.tl.score_genes(
            adata_full,
            gene_list=genes,
            score_name=cell_type,
            use_raw=True,
        )
        score_cols.append(cell_type)

    if not score_cols:
        raise RuntimeError("Marker gene scoring produced no usable gene sets.")

    cluster_labels: Dict[str, str] = {}
    cluster_debug: Dict[str, Dict[str, object]] = {}
    if "leiden" not in adata_full.obs:
        raise KeyError("AnnData lacks 'leiden' clusters required for annotation.")
    if not isinstance(adata_full.obs["leiden"].dtype, CategoricalDtype):
        adata_full.obs["leiden"] = adata_full.obs["leiden"].astype("category")

    print("[SupFig4][scanpy] Assigning cluster annotations...")
    for cluster_id in adata_full.obs["leiden"].cat.categories:
        mask = (adata_full.obs["leiden"] == cluster_id).values
        row_scores = []
        for cell_type in score_cols:
            scores = adata_full.obs.loc[mask, cell_type].to_numpy()
            agg = float(np.nanquantile(scores, AGG_QUANTILE)) if scores.size else -np.inf
            det_frac = detect_fraction(
                adata_full,
                mask,
                markers_unique[cell_type],
                threshold=0.0,
                min_detect=MIN_DETECT_GENES,
            )
            needed = MIN_DETECT_FRAC_BY_TYPE.get(cell_type, MIN_DETECT_FRAC_DEFAULT)
            if det_frac < needed:
                agg = -np.inf
            row_scores.append((cell_type, agg, det_frac))
        row_scores.sort(key=lambda item: item[1], reverse=True)
        top = row_scores[0]
        second = row_scores[1] if len(row_scores) > 1 else ("None", -np.inf, 0.0)
        margin = top[1] - second[1]
        chosen = top[0] if np.isfinite(top[1]) and margin >= MARGIN_THRESHOLD else "Ambiguous"
        cluster_labels[cluster_id] = chosen
        cluster_debug[cluster_id] = {
            "n_cells": int(mask.sum()),
            "top": top,
            "second": second,
            "margin": float(margin),
        }

    adata_full.obs["Cell Type"] = adata_full.obs["leiden"].map(cluster_labels).astype("category")
    print("\n[SupFig4][scanpy] Cluster \u2192 cell_type mapping (with top/second & margins):")
    for cluster_id in adata_full.obs["leiden"].cat.categories:
        info = cluster_debug.get(cluster_id)
        label = cluster_labels.get(cluster_id, "Unknown")
        if info is None:
            print(f"  {cluster_id}: {label:<40s}  (no scoring info)")
            continue
        top = info["top"]
        second = info["second"]
        margin = info["margin"]
        print(
            f"  {cluster_id}: {label:<40s}  n={info['n_cells']:>5d}  "
            f"top={top[0]}({top[1]:.3f}, det={top[2]:.2f})  "
            f"second={second[0]}({second[1]:.3f}, det={second[2]:.2f})  "
            f"margin={margin:.3f}"
        )


    barcode_table_path = config.data_dir / config.barcode_table_name
    if (
        config.skip_existing_outputs
        and config.write_barcode_tables
        and barcode_table_path.exists()
    ):
        print(f"[SupFig4][scanpy] Reusing barcode table from {barcode_table_path}")
        barcode_table = pd.read_csv(barcode_table_path)
    else:
        barcode_table = (
            adata_full.obs[["Cell Type", "Age"]]
            .reset_index(names="barcode")
            .assign(barcode=lambda df: df["barcode"].map(_standardize_barcode))
            .dropna(subset=["barcode"])
        )
        barcode_table = barcode_table.drop_duplicates(subset=["barcode"])
        print(f"[SupFig4][scanpy] Barcode table prepared with {len(barcode_table)} entries.")
        if config.write_barcode_tables:
            config.data_dir.mkdir(parents=True, exist_ok=True)
            barcode_table.to_csv(barcode_table_path, index=False)

    barcode_sample_paths: Dict[str, Path] = {}
    if config.write_barcode_tables:
        sample_dir = config.data_dir / config.barcode_subdir_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        for age_label, sub_df in barcode_table.groupby("Age", observed=False):
            age_str = str(age_label)
            output_path = sample_dir / f"barcodes_{age_str}.csv"
            barcode_sample_paths[age_str] = output_path
            if config.skip_existing_outputs and output_path.exists():
                continue
            sub_df.to_csv(output_path, index=False)

    isoform_counts_by_age: Dict[str, pd.DataFrame] = {}
    isoform_counts_paths: Dict[str, Path] = {}
    isoform_summary_paths: Dict[str, Path] = {}
    combined_counts: List[pd.DataFrame] = []
    reused_conditions: List[str] = []

    isoform_items = list(config.isoform_match_paths.items())
    print("[SupFig4][scanpy] Computing isoform counts per condition...")
    total_raw_barcodes = 0
    processed_raw = 0
    last_bucket = -1
    text_total = 0
    progress_bar = (
        tqdm(
            total=0,
            desc="[SupFig4][scanpy] Correcting barcodes/UMIs",
            unit="raw barcodes",
            leave=True,
        )
        if tqdm is not None
        else None
    )

    def _make_progress_callback(label: str, step: int):
        def _callback(delta: int) -> None:
            nonlocal processed_raw, last_bucket
            processed_raw += delta
            if progress_bar is not None:
                progress_bar.update(delta)
            else:
                bucket = processed_raw // step if step > 0 else processed_raw
                if bucket != last_bucket:
                    total_display = text_total if text_total > 0 else "?"
                    print(
                        f"[SupFig4][scanpy]  -> {processed_raw}/{total_display} raw barcodes processed "
                        f"(latest: {label})"
                    )
                    last_bucket = bucket

        return _callback
    for age_label, path_tsv in isoform_items:
        if not path_tsv.exists():
            raise FileNotFoundError(path_tsv)
        label_str = str(age_label)
        counts_path = config.data_dir / f"{label_str.lower()}_isoform_by_celltype_counts.csv"
        summary_path = config.data_dir / f"{label_str.lower()}_correction_summary.csv"
        if (
            config.skip_existing_outputs
            and config.write_isoform_tables
            and counts_path.exists()
            and summary_path.exists()
        ):
            print(f"[SupFig4][scanpy]  - Reusing isoform counts for {label_str}")
            counts = pd.read_csv(counts_path)
            if "n_molecules" not in counts.columns:
                counts["n_molecules"] = np.nan
            isoform_counts_by_age[label_str] = counts
            isoform_counts_paths[label_str] = counts_path
            isoform_summary_paths[label_str] = summary_path
            combined_counts.append(counts)
            reused_conditions.append(label_str)
            continue

        whitelist = barcode_table[barcode_table["Age"] == age_label].copy()
        if whitelist.empty:
            raise ValueError(
                f"No barcodes available for age label '{age_label}'. "
                "Verify condition_map and isoform_match_paths."
            )
        print(f"[SupFig4][scanpy]  - Processing isoforms for {label_str}...")
        raw_matches = load_raw_matches(path_tsv)
        cleaned_unique = (
            pd.Series(raw_matches["cb"].map(clean_barcode), dtype="string")
            .dropna()
            .replace("", pd.NA)
            .dropna()
            .nunique()
        )
        total_raw_barcodes += cleaned_unique
        text_total = total_raw_barcodes
        if progress_bar is not None:
            progress_bar.total += cleaned_unique
            progress_bar.refresh()
        merged, counts, summary = compute_counts_per_isoform_celltype(
            raw_matches,
            whitelist,
            max_distance=config.isoform_max_distance,
            progress_callback=_make_progress_callback(
                label_str, max(cleaned_unique // 50, 500) if cleaned_unique > 0 else 500
            ),
            progress_step=max(cleaned_unique // 50, 500) if cleaned_unique > 0 else 500,
        )
        isoform_counts_by_age[label_str] = counts
        combined_counts.append(counts)
        if config.write_isoform_tables:
            counts.to_csv(counts_path, index=False)
            summary.to_csv(summary_path, index=False)
            isoform_summary_paths[label_str] = summary_path
            isoform_counts_paths[label_str] = counts_path

    if progress_bar is not None:
        progress_bar.close()
    print(
        "[SupFig4][scanpy] Isoform counting complete."
        + (
            f" Reused cached counts for: {', '.join(reused_conditions)}."
            if reused_conditions
            else ""
        )
    )

    combined_counts_path = config.data_dir / config.combined_counts_filename
    if combined_counts:
        combined_isoform_counts = pd.concat(combined_counts, ignore_index=True)
        agg_map = {"n_cells": "sum"}
        if "n_molecules" in combined_isoform_counts.columns:
            agg_map["n_molecules"] = "sum"
        combined_isoform_counts = (
            combined_isoform_counts.groupby(
                ["transcript_id", "Cell Type", "Age"], as_index=False
            )[list(agg_map.keys())]
            .sum()
        )
    else:
        combined_isoform_counts = pd.DataFrame(
            columns=["transcript_id", "Cell Type", "Age", "n_cells", "n_molecules"]
        )
        if config.skip_existing_outputs and config.write_isoform_tables and combined_counts_path.exists():
            combined_isoform_counts = pd.read_csv(combined_counts_path)
            if "n_molecules" not in combined_isoform_counts.columns:
                combined_isoform_counts["n_molecules"] = np.nan

    # Build merged splice annotations once we have isoform counts
    splice_annotation_path = config.data_dir / config.splice_annotation_filename
    splice_paths: List[Path] = []
    if config.splice_result_paths:
        splice_paths = [Path(p) for p in config.splice_result_paths]
    else:
        for candidate in ("results_a_vs_c_splice.csv", "results_b_vs_d_splice.csv"):
            candidate_path = config.data_dir / candidate
            if candidate_path.exists():
                splice_paths.append(candidate_path)

    annotation_gtf_path = config.annotation_gtf_path
    if annotation_gtf_path is None:
        candidate = config.data_dir / "all.flair.collapse.isoforms_event_label.gtf"
        if candidate.exists():
            annotation_gtf_path = candidate

    annotation_lookup = pd.DataFrame(columns=["transcript_id", "gene_id", "promoter_group"])
    splice_annotation_table = pd.DataFrame()
    if splice_paths:
        reuse_annotations = config.skip_existing_outputs and splice_annotation_path.exists()
        if reuse_annotations:
            print(
                f"[SupFig4][scanpy] Reusing splice annotation table from {splice_annotation_path}"
            )
            splice_annotation_table = pd.read_csv(splice_annotation_path)
            annotation_lookup = (
                splice_annotation_table.loc[:, ["transcript_id", "gene_id", "promoter_group"]]
                .dropna(subset=["transcript_id"])
                .assign(transcript_id=lambda df: df["transcript_id"].astype(str).str.strip())
                .drop_duplicates(subset=["transcript_id"])
            )
        else:
            print(
                f"[SupFig4][scanpy] Building splice annotation table from {len(splice_paths)} files..."
            )
            splice_annotation_table, annotation_lookup = _compile_splice_annotations(
                splice_paths, annotation_gtf_path
            )
            if not splice_annotation_table.empty and config.write_isoform_tables:
                config.data_dir.mkdir(parents=True, exist_ok=True)
                splice_annotation_table.to_csv(splice_annotation_path, index=False)
    elif splice_annotation_path.exists():
        splice_annotation_table = pd.read_csv(splice_annotation_path)
        annotation_lookup = (
            splice_annotation_table.loc[:, ["transcript_id", "gene_id", "promoter_group"]]
            .dropna(subset=["transcript_id"])
            .assign(transcript_id=lambda df: df["transcript_id"].astype(str).str.strip())
            .drop_duplicates(subset=["transcript_id"])
        )

    if annotation_lookup.empty and not combined_isoform_counts.empty:
        transcript_ids = combined_isoform_counts["transcript_id"].dropna().astype(str).unique()
        annotation_lookup = _build_transcript_gene_promoter_map(
            annotation_gtf_path, transcript_ids
        )

    if not combined_isoform_counts.empty and not annotation_lookup.empty:
        combined_isoform_counts = combined_isoform_counts.drop(
            columns=[col for col in ("gene_id", "promoter_group") if col in combined_isoform_counts.columns]
        )
        combined_isoform_counts = combined_isoform_counts.merge(
            annotation_lookup, on="transcript_id", how="left"
        )
        ordered_cols = [
            col
            for col in [
                "transcript_id",
                "gene_id",
                "promoter_group",
                "Cell Type",
                "Age",
                "n_cells",
                "n_molecules",
            ]
            if col in combined_isoform_counts.columns
        ]
        remaining_cols = [
            col for col in combined_isoform_counts.columns if col not in ordered_cols
        ]
        combined_isoform_counts = combined_isoform_counts[ordered_cols + remaining_cols]

    if not combined_isoform_counts.empty:
        before_filter = len(combined_isoform_counts)
        combined_isoform_counts = _filter_gene_promoter_min_transcripts(
            combined_isoform_counts, min_transcripts=2
        )
        after_filter = len(combined_isoform_counts)
        if after_filter != before_filter:
            print(
                "[SupFig4][scanpy] Filtered isoform counts to gene/promoter groups with "
                f">=2 isoforms ({after_filter}/{before_filter} rows retained)."
            )

    if config.write_isoform_tables and not combined_isoform_counts.empty:
        config.data_dir.mkdir(parents=True, exist_ok=True)
        combined_isoform_counts.to_csv(combined_counts_path, index=False)

    isoform_group_ratio_path = config.data_dir / config.isoform_group_ratio_filename
    isoform_group_ratios = pd.DataFrame()
    if not combined_isoform_counts.empty:
        isoform_group_ratios = _compute_isoform_group_ratios(combined_isoform_counts)
        if config.write_isoform_tables:
            config.data_dir.mkdir(parents=True, exist_ok=True)
            isoform_group_ratios.to_csv(isoform_group_ratio_path, index=False)
    elif isoform_group_ratio_path.exists():
        isoform_group_ratios = pd.read_csv(isoform_group_ratio_path)

    splice_annotation_output_path = (
        splice_annotation_path if splice_annotation_path.exists() else None
    )
    isoform_group_ratio_output_path = (
        isoform_group_ratio_path if isoform_group_ratio_path.exists() else None
    )

    ratio_pair = config.isoform_ratio_pair
    if ratio_pair is None and len(config.isoform_transcripts) >= 2:
        ratio_pair = (config.isoform_transcripts[0], config.isoform_transcripts[1])

    isoform_ratio_table = pd.DataFrame()
    slope_data = pd.DataFrame()
    value_col_ratio = (
        "n_molecules" if "n_molecules" in combined_isoform_counts.columns else "n_cells"
    )
    if ratio_pair and not combined_isoform_counts.empty:
        isoform_counts = (
            combined_isoform_counts[
                combined_isoform_counts["transcript_id"].isin(list(ratio_pair))
            ]
            .pivot_table(
                index=["Cell Type", "Age"],
                columns="transcript_id",
                values=value_col_ratio,
                fill_value=0,
                observed=False,
            )
        )
        numerator = isoform_counts.get(ratio_pair[0], 0)
        denominator = numerator + isoform_counts.get(ratio_pair[1], 0)
        denom_safe = denominator.replace(0, np.nan)
        isoform_counts["ratio"] = numerator / denom_safe
        isoform_counts = isoform_counts.reset_index()
        isoform_ratio_table = isoform_counts.rename(
            columns={"ratio": f"ratio_{ratio_pair[0]}"}
        )
        slope_data = isoform_ratio_table.pivot_table(
            index="Cell Type", columns="Age", values=f"ratio_{ratio_pair[0]}", observed=False
        )

    heatmap_matrix = pd.DataFrame(columns=list(config.isoform_transcripts))
    if not combined_isoform_counts.empty and config.isoform_transcripts:
        subset = combined_isoform_counts[
            combined_isoform_counts["transcript_id"].isin(config.isoform_transcripts)
        ].copy()
        if not subset.empty:
            totals = (
                subset.groupby(["Cell Type", "Age"], observed=False)["n_cells"]
                .sum()
                .rename("total_cells")
                .reset_index()
            )
            subset = subset.merge(totals, on=["Cell Type", "Age"], how="left")
            subset["proportion"] = subset["n_cells"] / subset["total_cells"].replace(0, np.nan)
            subset["Cell Type_Age"] = (
                subset["Cell Type"].astype(str) + "_" + subset["Age"].astype(str)
            )
            heatmap_matrix = subset.pivot_table(
                index="Cell Type_Age",
                columns="transcript_id",
                values="proportion",
                fill_value=0.0,
                observed=False,
            )
            ordered_index = _sort_celltype_age_index(
                heatmap_matrix.index.tolist(), age_order
            )
            heatmap_matrix = heatmap_matrix.loc[ordered_index]

    return ScanpyAnalysisOutputs(
        adata_full=adata_full,
        adata_hvg=adata_hvg,
        barcode_table=barcode_table,
        barcode_table_path=barcode_table_path if config.write_barcode_tables else None,
        barcode_sample_paths=barcode_sample_paths,
        cluster_labels=cluster_labels,
        cluster_debug=cluster_debug,
        marker_sets=markers_unique,
        marker_source=marker_source,
        marker_coverage=coverage_df,
        marker_overlaps=overlap_df,
        isoform_counts_by_age=isoform_counts_by_age,
        isoform_counts_paths=isoform_counts_paths,
        isoform_summary_paths=isoform_summary_paths,
        combined_isoform_counts=combined_isoform_counts,
        combined_isoform_counts_path=combined_counts_path,
        isoform_ratio_table=isoform_ratio_table,
        isoform_proportion_matrix=heatmap_matrix,
        slope_data=slope_data,
        splice_annotation_table=splice_annotation_table,
        splice_annotation_path=splice_annotation_output_path,
        isoform_group_ratios=isoform_group_ratios,
        isoform_group_ratio_path=isoform_group_ratio_output_path,
    )

def run_minimap_alignments(config: MinimapAlignmentConfig) -> MinimapAlignmentOutputs:
    """
    Align R1 FASTQ files against the promoter-group FASTA using minimap2.
    """
    config.output_dir.mkdir(parents=True, exist_ok=True)

    paf_paths: Dict[str, Path] = {}
    paf_gz_paths: Dict[str, Optional[Path]] = {}
    to_process: List[Tuple[str, Path, Path, Iterable[Path]]] = []
    missing_fastqs: Dict[str, List[Path]] = {}

    for label, fastq_paths in config.fastq_map.items():
        output_stem = config.output_stems.get(label, f"{label}_pfamhits")
        paf_path = config.output_dir / f"{output_stem}.paf"
        paf_gz = paf_path.with_suffix(paf_path.suffix + ".gz")
        outputs_exist = paf_path.exists() or paf_gz.exists()
        missing_inputs = [path for path in fastq_paths if not path.exists()]
        reuse_reason: Optional[str] = None

        if missing_inputs:
            if outputs_exist:
                reuse_reason = "missing FASTQ inputs: " + ", ".join(
                    str(path) for path in missing_inputs
                )
            else:
                missing_fastqs[label] = missing_inputs
                continue
        elif config.skip_existing_outputs and outputs_exist:
            reuse_reason = "existing output detected"

        if reuse_reason:
            if not paf_path.exists() and paf_gz.exists():
                _decompress_gzip_file(paf_gz, paf_path)
            reused_path = paf_path if paf_path.exists() else paf_gz
            if reused_path is not None:
                print(
                    f"[SupFig4][minimap] Skipping minimap for {label}: "
                    f"{reused_path} already generated ({reuse_reason})."
                )
            paf_paths[label] = paf_path
            paf_gz_paths[label] = paf_gz if paf_gz.exists() else None
            continue

        to_process.append((label, paf_path, paf_gz, fastq_paths))

    if missing_fastqs:
        missing_desc = "; ".join(
            f"{label}: {', '.join(str(path) for path in paths)}"
            for label, paths in missing_fastqs.items()
        )
        raise FileNotFoundError(
            f"Missing FASTQ inputs with no cached outputs available ({missing_desc})."
        )

    if to_process and not config.reference_fasta.exists():
        raise FileNotFoundError(config.reference_fasta)

    binary = _resolve_minimap_binary(config)

    for label, paf_path, paf_gz, fastq_paths in to_process:
        cmd = [
            str(binary),
            "-t",
            str(config.threads),
            "-x",
            config.preset,
        ]
        if config.disable_secondary:
            cmd.append("--secondary=no")
        if config.include_cigar:
            cmd.append("-c")
        cmd.append(str(config.reference_fasta))
        cmd.extend(str(path) for path in fastq_paths)

        with paf_path.open("w") as handle:
            subprocess.run(cmd, check=True, stdout=handle, stderr=subprocess.PIPE)

        paf_paths[label] = paf_path

        if config.compress_output:
            with paf_path.open("rb") as src, gzip.open(paf_gz, "wb") as dst:
                shutil.copyfileobj(src, dst)
            paf_gz_paths[label] = paf_gz
        else:
            paf_gz_paths[label] = paf_gz if paf_gz.exists() else None

    return MinimapAlignmentOutputs(
        minimap_bin=binary, paf_paths=paf_paths, paf_gz_paths=paf_gz_paths
    )


def _resolve_minimap_binary(config: MinimapAlignmentConfig) -> Path:
    if config.minimap_bin is not None:
        if not config.minimap_bin.exists():
            raise FileNotFoundError(config.minimap_bin)
        return config.minimap_bin

    found = shutil.which("minimap2")
    if found:
        return Path(found)

    if not config.auto_install:
        raise FileNotFoundError(
            "minimap2 binary was not provided and could not be located in PATH."
        )

    install_dir = config.install_dir or (config.output_dir / "tools")
    install_dir.mkdir(parents=True, exist_ok=True)
    binary = install_dir / MINIMAP_ARCHIVE / "minimap2"
    if binary.exists():
        return binary

    archive_path = install_dir / f"{MINIMAP_ARCHIVE}.tar.bz2"
    urllib.request.urlretrieve(MINIMAP_URL, archive_path)
    with tarfile.open(archive_path, mode="r:bz2") as tar:
        _safe_extract(tar, install_dir)
    binary.chmod(0o755)
    archive_path.unlink(missing_ok=True)
    return binary


def _safe_extract(tar: tarfile.TarFile, target_dir: Path) -> None:
    target_dir = target_dir.resolve()
    target_str = str(target_dir)
    for member in tar.getmembers():
        member_path = (target_dir / member.name).resolve()
        member_str = str(member_path)
        if not member_str.startswith(target_str):
            raise RuntimeError(f"Unsafe path detected in archive: {member.name}")
    tar.extractall(target_dir)
