from __future__ import annotations

import io
import json
import os
import re
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import gridspec

try:  # pragma: no cover - optional dependency
    from matplotlib_venn import venn2
except ImportError:  # pragma: no cover - optional dependency
    venn2 = None  # type: ignore[assignment]

PathLike = Union[str, os.PathLike]

@dataclass
class SupFigure3PanelStyle:
    figsize: Tuple[float, float] = (8.0, 4.5)
    width_ratios: Tuple[float, float] = (1.0, 1.0)
    wspace: float = 0.4
    panel_labels: Tuple[str, str] = ("A", "B")
    panel_label_offsets: Tuple[Tuple[float, float], Tuple[float, float]] = (
        (-0.01, 0.02),
        (-0.01, 0.02),
    )
    panel_label_fontsize: float = 13.0
    single_panel_size: Tuple[float, float] = (5.5, 5.5)
    dpi: int = 300
    font_family: str = "sans-serif"
    font_preferences: Tuple[str, ...] = ("Arial", "DejaVu Sans", "Liberation Sans")
    base_fontsize: float = 8.0
    title_size: float = 7.0
    label_size: float = 7.0
    legend_size: float = 7.0
    tick_label_size: float = 7.0
    venn_title_fontsize: Optional[float] = None
    venn_set_label_fontsizes: Tuple[Optional[float], Optional[float]] = (None, None)
    venn_colors: Tuple[str, str] = ("#4C72B0", "#DD8452")
    venn_alpha: float = 0.6
    panel_c_figsize: Tuple[float, float] = (6.0, 3.2)
    panel_c_colors: Tuple[str, str] = ("#55A868", "#C44E52")
    panel_c_hide_box: bool = False


@dataclass
class SupFigure3PanelLabels:
    panelA_title: Optional[str] = None
    panelB_title: Optional[str] = None
    panelA_left: Optional[str] = None
    panelA_right: Optional[str] = None
    panelB_left: Optional[str] = None
    panelB_right: Optional[str] = None
    panelC_title: Optional[str] = None
    panelC_detected: Optional[str] = None
    panelC_not_detected: Optional[str] = None
    panelC_xlabel: Optional[str] = None


@dataclass
class SupFigure3PanelData:
    pfam_union: Set[str]
    ms_strict: Set[str]
    ms_relaxed: Set[str]
    relaxed_label: str
    relaxed_suffix: str
    panel_titles: Tuple[str, str]
    sensitivity_preset: str
    panelC_detected: Set[str]
    panelC_not_detected: Set[str]
    panelC_variant_label: str
    panelC_suffix: str

    @property
    def union_size(self) -> int:
        return len(self.pfam_union)

    @property
    def strict_size(self) -> int:
        return len(self.ms_strict)

    @property
    def relaxed_size(self) -> int:
        return len(self.ms_relaxed)

    @property
    def strict_overlap(self) -> int:
        return len(self.pfam_union & self.ms_strict)

    @property
    def relaxed_overlap(self) -> int:
        return len(self.pfam_union & self.ms_relaxed)

    @property
    def strict_detected(self) -> Set[str]:
        return self.pfam_union & self.ms_strict

    @property
    def missed_after_strict(self) -> Set[str]:
        return self.pfam_union - self.strict_detected

    @property
    def panelC_detected_count(self) -> int:
        return len(self.panelC_detected)

    @property
    def panelC_not_detected_count(self) -> int:
        return len(self.panelC_not_detected)


# --------------------------
# Utilities and normalization
# --------------------------

def _strip_version(value: str) -> str:
    if not isinstance(value, str):
        return value
    return re.sub(r"\.\d+$", "", value.strip())


def _normalize_transcript_id(value: str) -> str:
    if not isinstance(value, str):
        return value
    return _strip_version(value.strip())


def _parse_transcript_from_leading_razor(values: str) -> Optional[str]:
    if not isinstance(values, str) or not values.strip():
        return None
    token = values.split(";", 1)[0].strip()
    if "_" in token:
        token = token.split("_", 1)[0]
    return _normalize_transcript_id(token)


def _split_proteins_to_txs(proteins: str) -> List[str]:
    if not isinstance(proteins, str) or not proteins.strip():
        return []
    out: List[str] = []
    for token in proteins.split(";"):
        token = token.strip()
        if not token:
            continue
        tx = _parse_transcript_from_leading_razor(token)
        if tx:
            out.append(tx)
    return out


def _as_path(value: PathLike) -> Path:
    if isinstance(value, Path):
        return value
    return Path(value)


def _ensure_dir(path: PathLike) -> Path:
    path = _as_path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path

def apply_sup_fig3_style(style: SupFigure3PanelStyle) -> None:
    """Apply global matplotlib rcParams similar to Sup Figure 4 styling."""
    plt.rcParams["font.family"] = style.font_family
    preferred_fonts = list(style.font_preferences)
    existing_fonts = plt.rcParams.get("font.sans-serif", [])
    plt.rcParams["font.sans-serif"] = preferred_fonts + [
        f for f in existing_fonts if f not in preferred_fonts
    ]
    plt.rcParams["font.size"] = style.base_fontsize
    plt.rcParams["axes.titlesize"] = style.title_size
    plt.rcParams["axes.labelsize"] = style.label_size
    plt.rcParams["legend.fontsize"] = style.legend_size
    plt.rcParams.update({
        "xtick.labelsize": style.tick_label_size,
        "ytick.labelsize": style.tick_label_size,
    })


def resolve_sup_fig3_panel_labels(panel_data: SupFigure3PanelData,
                                  overrides: Optional[SupFigure3PanelLabels] = None) -> Dict[str, str]:
    cfg = overrides or SupFigure3PanelLabels()
    panelA_title = cfg.panelA_title or panel_data.panel_titles[0]
    panelB_title = cfg.panelB_title or panel_data.panel_titles[1]
    panelA_left = cfg.panelA_left or "PFAM (A ∪ B)"
    panelA_right = cfg.panelA_right or "Mass spec (unique)"
    panelB_left = cfg.panelB_left or "PFAM (A ∪ B)"
    panelB_right = cfg.panelB_right or panel_data.relaxed_label
    default_c_title = f"Detection among PFAM-missed ({panel_data.panelC_variant_label})"
    panelC_title = cfg.panelC_title or default_c_title
    panelC_detected = cfg.panelC_detected or "Detected at all"
    panelC_not = cfg.panelC_not_detected or "Not detected"
    panelC_xlabel = cfg.panelC_xlabel or "PFAM (A ∪ B) transcripts missed in Panel A (count)"
    return {
        "panelA_title": panelA_title,
        "panelB_title": panelB_title,
        "panelA_left": panelA_left,
        "panelA_right": panelA_right,
        "panelB_left": panelB_left,
        "panelB_right": panelB_right,
        "panelC_title": panelC_title,
        "panelC_detected": panelC_detected,
        "panelC_not_detected": panelC_not,
        "panelC_xlabel": panelC_xlabel,
    }


def _place_panel_label(fig, ax, label: str, x_pad: float = -0.01, y_pad: float = 0.02,
                       fontsize: int = 13) -> None:
    bbox = ax.get_position()
    fig.text(
        bbox.x0 + x_pad,
        bbox.y1 + y_pad,
        label,
        fontsize=fontsize,
        fontweight="bold",
        ha="right",
        va="bottom",
    )

def _get_set_label_fontsizes(style: SupFigure3PanelStyle) -> Optional[Tuple[Optional[float], Optional[float]]]:
    fontsizes: Optional[Sequence[Optional[float]]] = getattr(style, "venn_set_label_fontsizes", None)
    if not fontsizes:
        return None
    fonts_tuple = tuple(fontsizes)
    if len(fonts_tuple) < 2:
        fonts = fonts_tuple + (None,) * (2 - len(fonts_tuple))
    else:
        fonts = fonts_tuple[:2]
    if any(size is not None for size in fonts):
        return fonts  # type: ignore[return-value]
    return None


def build_sup_fig3_gridspec(panel_data: SupFigure3PanelData,
                            style: Optional[SupFigure3PanelStyle] = None,
                            panel_labels: Optional[SupFigure3PanelLabels] = None):
    style = style or SupFigure3PanelStyle()
    label_values = resolve_sup_fig3_panel_labels(panel_data, panel_labels)
    venn_label_sizes = _get_set_label_fontsizes(style)
    fig = plt.figure(figsize=style.figsize)
    gs = gridspec.GridSpec(
        1,
        2,
        width_ratios=style.width_ratios,
        wspace=style.wspace,
        figure=fig,
    )

    ax_strict = fig.add_subplot(gs[0, 0])
    plot_two_set_venn(
        panel_data.union_size,
        panel_data.strict_size,
        panel_data.strict_overlap,
        labels=(label_values["panelA_left"], label_values["panelA_right"]),
        title=label_values["panelA_title"],
        title_fontsize=style.venn_title_fontsize,
        set_label_fontsizes=venn_label_sizes,
        ax=ax_strict,
        set_colors=style.venn_colors,
        alpha=style.venn_alpha,
    )

    ax_relaxed = fig.add_subplot(gs[0, 1])
    plot_two_set_venn(
        panel_data.union_size,
        panel_data.relaxed_size,
        panel_data.relaxed_overlap,
        labels=(label_values["panelB_left"], label_values["panelB_right"]),
        title=label_values["panelB_title"],
        title_fontsize=style.venn_title_fontsize,
        set_label_fontsizes=venn_label_sizes,
        ax=ax_relaxed,
        set_colors=style.venn_colors,
        alpha=style.venn_alpha,
    )

    fig.tight_layout()
    offset_strict = style.panel_label_offsets[0]
    offset_relaxed = style.panel_label_offsets[1]
    _place_panel_label(
        fig,
        ax_strict,
        style.panel_labels[0],
        x_pad=offset_strict[0],
        y_pad=offset_strict[1],
        fontsize=int(style.panel_label_fontsize),
    )
    _place_panel_label(
        fig,
        ax_relaxed,
        style.panel_labels[1],
        x_pad=offset_relaxed[0],
        y_pad=offset_relaxed[1],
        fontsize=int(style.panel_label_fontsize),
    )
    return fig, {"strict": ax_strict, "relaxed": ax_relaxed}


# --------------------------
# Directory helpers
# --------------------------

def _iter_sup_fig3_search_dirs(start: Path) -> List[Path]:
    search: List[Path] = []
    seen: Set[Path] = set()
    cur = start
    while True:
        if cur not in seen:
            search.append(cur)
            seen.add(cur)
        nested = cur / "Sup_Figure_3"
        if nested.exists() and nested not in seen:
            search.append(nested)
            seen.add(nested)
        if cur.parent == cur:
            break
        cur = cur.parent
    return search


def _locate_sup_fig3_root(start_dir: Optional[PathLike] = None) -> Path:
    start = _as_path(start_dir) if start_dir is not None else Path.cwd()
    start = start.resolve()
    for candidate in _iter_sup_fig3_search_dirs(start):
        if (candidate / "data_files").exists():
            return candidate.resolve()
    raise FileNotFoundError(
        f"Could not locate a Sup_Figure_3 directory with data_files relative to {start}"
    )


def resolve_sup_fig3_inputs(start_dir: Optional[PathLike] = None,
                            outputs_subdir: Optional[PathLike] = None) -> Dict[str, Path]:
    base_dir = _locate_sup_fig3_root(start_dir)
    data_dir = base_dir / "data_files"
    figure_dir = _ensure_dir(base_dir / "figure_files")
    figure_data_dir = _ensure_dir(base_dir / "figure_data")
    if outputs_subdir is None or str(outputs_subdir).strip() in {"", ".", "./"}:
        out_dir = figure_dir
    else:
        out_dir = _ensure_dir(figure_dir / str(outputs_subdir))

    required = {
        "ac_csv": data_dir / "ac_pfam_both_hits_final.csv",
        "bd_csv": data_dir / "bd_pfam_both_hits_final.csv",
        "zip_path": data_dir / "20251027_Mouse_Ovary_Manuscript_LCMS.zip",
    }
    for label, path in required.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing required file for {label}: {path}")

    paths: Dict[str, Path] = {
        "base_dir": base_dir,
        "data_dir": data_dir,
        "figure_dir": figure_dir,
        "figure_data_dir": figure_data_dir,
        "out_dir": out_dir,
    }
    paths.update(required)
    return paths


# --------------------------
# Loaders
# --------------------------

def load_pfam_transcripts(csv_path: PathLike, id_column: str = "transcript_id") -> Set[str]:
    csv_path = _as_path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(csv_path)
    df = pd.read_csv(csv_path)
    if id_column not in df.columns:
        for col in ["Transcript", "transcript", "tx_id", "tx", "id"]:
            if col in df.columns:
                id_column = col
                break
        else:
            raise KeyError(f"Transcript id column not found in {csv_path}")
    tx = df[id_column].dropna().astype(str).map(_normalize_transcript_id)
    return set(tx.tolist())


def _load_pfam_metadata_frame(csv_path: PathLike,
                              columns: Sequence[str]) -> pd.DataFrame:
    csv_path = _as_path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(csv_path)
    df = pd.read_csv(csv_path)
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise KeyError(f"Missing PFAM metadata columns {missing} in {csv_path}")
    subset = df.loc[:, list(columns)].copy()

    def _coerce_tx(value: object) -> Optional[str]:
        if pd.isna(value):
            return None
        return _normalize_transcript_id(str(value))

    subset.loc[:, "transcript_id"] = subset["transcript_id"].map(_coerce_tx)
    subset = subset.dropna(subset=["transcript_id"])
    subset.loc[:, "transcript_id"] = subset["transcript_id"].astype(str)
    return subset


def build_pfam_union_metadata_table(ac_csv: PathLike,
                                    bd_csv: PathLike,
                                    columns: Sequence[str] = (
                                        "transcript_id",
                                        "pfam_id",
                                        "target_name",
                                        "description",
                                    )) -> pd.DataFrame:
    """Return a transcript-level metadata table covering the PFAM union."""
    columns = tuple(columns)
    if "transcript_id" not in columns:
        columns = ("transcript_id",) + columns
    ordered_columns: Tuple[str, ...] = tuple(dict.fromkeys(columns))
    frames = [
        _load_pfam_metadata_frame(ac_csv, ordered_columns),
        _load_pfam_metadata_frame(bd_csv, ordered_columns),
    ]
    combined = pd.concat(frames, ignore_index=True)
    combined = combined.drop_duplicates(subset="transcript_id", keep="first")
    combined = combined.sort_values("transcript_id").reset_index(drop=True)
    return combined.loc[:, ordered_columns]


def export_pfam_union_metadata_table(ac_csv: PathLike,
                                     bd_csv: PathLike,
                                     out_path: PathLike,
                                     columns: Sequence[str] = (
                                         "transcript_id",
                                         "pfam_id",
                                         "target_name",
                                         "description",
                                     )) -> Path:
    """Write the PFAM union metadata table to a tab-separated file."""
    df = build_pfam_union_metadata_table(ac_csv, bd_csv, columns=columns)
    out_path = _as_path(out_path)
    _ensure_dir(out_path.parent)
    df.to_csv(out_path, sep="\t", index=False)
    return out_path


@dataclass
class MassSpecConfig:
    id_column: str = "T: Leading razor protein"
    intensity_cols_1: Tuple[str, str, str] = ("Intensity 1A", "Intensity 1B", "Intensity 1C")
    intensity_cols_2: Tuple[str, str, str] = ("Intensity 2A", "Intensity 2B", "Intensity 2C")

    def all_intensity_cols(self) -> Tuple[str, ...]:
        return tuple(self.intensity_cols_1) + tuple(self.intensity_cols_2)


def load_mass_spec_transcripts_from_txt(txt_path: PathLike,
                                        cfg: Optional[MassSpecConfig] = None,
                                        require_detected_in_any: Optional[str] = None) -> Set[str]:
    cfg = cfg or MassSpecConfig()
    txt_path = _as_path(txt_path)
    if not txt_path.exists():
        raise FileNotFoundError(txt_path)

    df = pd.read_csv(txt_path, sep="\t")

    if require_detected_in_any is not None:
        if require_detected_in_any == "1":
            cols = cfg.intensity_cols_1
        elif require_detected_in_any == "2":
            cols = cfg.intensity_cols_2
        elif require_detected_in_any == "any":
            cols = cfg.all_intensity_cols()
        else:
            raise ValueError("require_detected_in_any in {None, '1', '2', 'any'}")
        missing = [col for col in cols if col not in df.columns]
        if missing:
            raise KeyError(f"Missing intensity columns: {missing}")
        mask = False
        for col in cols:
            mask = mask | (pd.to_numeric(df[col], errors="coerce").fillna(0) > 0)
        df = df.loc[mask]

    id_col = cfg.id_column if cfg.id_column in df.columns else None
    if id_col is None:
        for col in [cfg.id_column, "Leading razor protein", "Proteins", "T: Proteins"]:
            if col in df.columns:
                id_col = col
                break
    if id_col is None:
        raise KeyError("No id column found for mass-spec file")

    tx = df[id_col].dropna().astype(str).map(_parse_transcript_from_leading_razor)
    tx = tx.dropna()
    return set(tx.tolist())


def load_mass_spec_transcripts_from_zip(zip_path: PathLike,
                                        member_name: str = "20251027_peptides_initial_and_uniqueProteins_filtered.txt",
                                        **kwargs) -> Set[str]:
    zip_path = _as_path(zip_path)
    if not zip_path.exists():
        raise FileNotFoundError(zip_path)
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        target = None
        if member_name in names:
            target = member_name
        else:
            low = member_name.lower()
            for name in names:
                if name.lower().endswith(low):
                    target = name
                    break
        if target is None:
            raise FileNotFoundError(f"{member_name} not found in zip")
        out_dir = _ensure_dir(zip_path.parent / "_extracted")
        out_path = out_dir / os.path.basename(target)
        with zf.open(target) as src, out_path.open("wb") as dst:
            dst.write(src.read())
    return load_mass_spec_transcripts_from_txt(out_path, **kwargs)


# --------------------------
# Stringency engine from peptides.txt
# --------------------------

def _get_intensity_cols(df: pd.DataFrame) -> List[str]:
    return [col for col in df.columns if col.startswith("Intensity ") and col != "Intensity"]


def build_ms_transcripts_from_peptides_raw(zip_path: PathLike,
                                           member_name: str = "peptides.txt",
                                           unique_only: bool = True,
                                           require_any_intensities: bool = True,
                                           pep_max: Optional[float] = None,
                                           source: str = "leading",
                                           min_peptides: Optional[int] = None,
                                           remove_reverse_contam: bool = True) -> Set[str]:
    zip_path = _as_path(zip_path)
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(member_name) as handle:
            df = pd.read_csv(io.BytesIO(handle.read()), sep='\t')

    if remove_reverse_contam:
        if "Reverse" in df.columns:
            df = df[df["Reverse"] != "+"]
        if "Potential contaminant" in df.columns:
            df = df[df["Potential contaminant"] != "+"]

    if unique_only and "Unique (Proteins)" in df.columns:
        df = df[df["Unique (Proteins)"].astype(str) == "yes"]

    if pep_max is not None and "PEP" in df.columns:
        df = df[pd.to_numeric(df["PEP"], errors="coerce").fillna(1.0) <= float(pep_max)]

    if require_any_intensities:
        cols = _get_intensity_cols(df)
        if cols:
            ints = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
            df = df[(ints > 0).any(axis=1)]

    if min_peptides is None:
        if source == "all" and "Proteins" in df.columns:
            txs: Set[str] = set()
            for value in df["Proteins"].astype(str):
                txs.update(_split_proteins_to_txs(value))
            return txs
        col = "Leading razor protein" if "Leading razor protein" in df.columns else None
        if col is None:
            raise KeyError("Leading razor protein column missing in peptides.txt")
        return set(
            df[col].astype(str).map(_parse_transcript_from_leading_razor).dropna().tolist()
        )

    seq_col = "Sequence" if "Sequence" in df.columns else None
    if seq_col is None:
        raise KeyError("Sequence column not found in peptides.txt")

    tx_to_peps: Dict[str, Set[str]] = {}
    if source == "all" and "Proteins" in df.columns:
        iterator = df.iterrows()
        for _, row in iterator:
            txs = _split_proteins_to_txs(row.get("Proteins", ""))
            seq = row.get(seq_col, "")
            for tx in txs:
                if not tx:
                    continue
                tx_to_peps.setdefault(tx, set()).add(seq)
    else:
        col = "Leading razor protein" if "Leading razor protein" in df.columns else None
        if col is None:
            raise KeyError("Leading razor protein column missing in peptides.txt")
        for _, row in df.iterrows():
            tx = _parse_transcript_from_leading_razor(row.get(col, ""))
            seq = row.get(seq_col, "")
            if tx:
                tx_to_peps.setdefault(tx, set()).add(seq)
    return {tx for tx, peptides in tx_to_peps.items() if len(peptides) >= int(min_peptides)}

# Panel C helper
def _build_ms_expand_all_any(zip_path: PathLike,
                             member_name: str = "peptides.txt") -> Set[str]:
    zip_path = _as_path(zip_path)
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(member_name) as handle:
            df = pd.read_csv(io.BytesIO(handle.read()), sep='\t')

    if "Reverse" in df.columns:
        df = df[df["Reverse"] != "+"]
    if "Potential contaminant" in df.columns:
        df = df[df["Potential contaminant"] != "+"]

    int_cols = [c for c in df.columns if c.startswith("Intensity ") and c != "Intensity"]
    if int_cols:
        ints = df[int_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        df = df[(ints > 0).any(axis=1)]

    if "Proteins" not in df.columns:
        raise KeyError("peptides.txt is missing the 'Proteins' column")

    txs: Set[str] = set()
    for value in df["Proteins"].astype(str):
        txs.update(_split_proteins_to_txs(value))
    return txs

# --------------------------
# Panel data preparation
# --------------------------

def prepare_sup_fig3_panel_data(ac_csv: PathLike,
                                bd_csv: PathLike,
                                zip_path: PathLike,
                                sensitivity_preset: str = "leading_2peps") -> SupFigure3PanelData:
    ac_csv = _as_path(ac_csv)
    bd_csv = _as_path(bd_csv)
    zip_path = _as_path(zip_path)

    pfam_A = load_pfam_transcripts(ac_csv)
    pfam_B = load_pfam_transcripts(bd_csv)
    pfam_union = pfam_A | pfam_B

    ms_strict = load_mass_spec_transcripts_from_zip(
        zip_path,
        member_name="20251027_peptides_initial_and_uniqueProteins_filtered.txt",
    )

    if sensitivity_preset == "leading_2peps":
        ms_relaxed = build_ms_transcripts_from_peptides_raw(
            zip_path,
            member_name="peptides.txt",
            unique_only=False,
            require_any_intensities=True,
            pep_max=None,
            source="leading",
            min_peptides=2,
            remove_reverse_contam=True,
        )
        relaxed_label = "MS (non-unique, leading, ≥2 peptides)"
        suffix = "relaxed_leading2peps"
        panel_title = "Panel B — Sensitivity (leading, ≥2 peptides)"
    elif sensitivity_preset == "expandAll_2peps":
        ms_relaxed = build_ms_transcripts_from_peptides_raw(
            zip_path,
            member_name="peptides.txt",
            unique_only=False,
            require_any_intensities=True,
            pep_max=None,
            source="all",
            min_peptides=2,
            remove_reverse_contam=True,
        )
        relaxed_label = "MS (non-unique, expand-all, ≥2 peptides)"
        suffix = "relaxed_expandAll2peps"
        panel_title = "Panel B — Sensitivity (expand-all, ≥2 peptides)"
    else:
        raise ValueError("Unknown sensitivity_preset; use 'leading_2peps' or 'expandAll_2peps'")

    ms_expand_any = _build_ms_expand_all_any(zip_path)
    strict_detected = pfam_union & ms_strict
    missed = pfam_union - strict_detected
    panelC_detected = missed & ms_expand_any
    panelC_not_detected = missed - panelC_detected

    return SupFigure3PanelData(
        pfam_union=pfam_union,
        ms_strict=ms_strict,
        ms_relaxed=ms_relaxed,
        relaxed_label=relaxed_label,
        relaxed_suffix=suffix,
        panel_titles=("Panel A — Strict", panel_title),
        sensitivity_preset=sensitivity_preset,
        panelC_detected=panelC_detected,
        panelC_not_detected=panelC_not_detected,
        panelC_variant_label="expand-all (≥1 peptide, non-unique)",
        panelC_suffix="panelC_expandAll_any",
    )

def export_sup_fig3_panel_tables(panel_data: SupFigure3PanelData,
                                 out_dir: PathLike,
                                 figure_data_dir: Optional[PathLike] = None) -> None:
    out_dir = _ensure_dir(out_dir)
    strict_overlap = sorted(panel_data.pfam_union & panel_data.ms_strict)
    relaxed_overlap = sorted(panel_data.pfam_union & panel_data.ms_relaxed)

    pd.Series(strict_overlap).to_csv(
        out_dir / "pfam_union_with_ms_STRICT.csv",
        index=False,
        header=["transcript_id"],
    )
    pd.Series(relaxed_overlap).to_csv(
        out_dir / f"pfam_union_with_ms_{panel_data.relaxed_suffix.upper()}.csv",
        index=False,
        header=["transcript_id"],
    )
    pd.Series(sorted(panel_data.panelC_detected)).to_csv(
        out_dir / "missed_detected_panelC_expandAll_any.csv",
        index=False,
        header=["transcript_id"],
    )
    pd.Series(sorted(panel_data.panelC_not_detected)).to_csv(
        out_dir / "missed_NOT_detected_panelC_expandAll_any.csv",
        index=False,
        header=["transcript_id"],
    )

    panelA_counts = {
        "pfam_union_size": panel_data.union_size,
        "ms_size": panel_data.strict_size,
        "overlap": panel_data.strict_overlap,
    }
    panelB_counts = {
        "pfam_union_size": panel_data.union_size,
        "ms_size": panel_data.relaxed_size,
        "overlap": panel_data.relaxed_overlap,
        "preset": panel_data.sensitivity_preset,
    }
    with open(out_dir / "panelA_counts.json", "w") as handle:
        json.dump(panelA_counts, handle, indent=2)
    with open(out_dir / "panelB_counts.json", "w") as handle:
        json.dump(panelB_counts, handle, indent=2)
    with open(out_dir / "panelC_expandAll_any_counts.json", "w") as handle:
        json.dump(
            {
                "pfam_union_size": panel_data.union_size,
                "panelA_overlap_strict": panel_data.strict_overlap,
                "panelA_missed_total": len(panel_data.missed_after_strict),
                "detected_at_all_under_variant": panel_data.panelC_detected_count,
                "not_detected_under_variant": panel_data.panelC_not_detected_count,
                "variant": "expand_all_any",
                "variant_label": panel_data.panelC_variant_label,
            },
            handle,
            indent=2,
        )
    if figure_data_dir is not None:
        figure_data_dir = _ensure_dir(figure_data_dir)
        summary_rows = [
            {
                "panel": "A_strict",
                "pfam_union_size": panel_data.union_size,
                "ms_size": panel_data.strict_size,
                "overlap": panel_data.strict_overlap,
                "variant_label": "strict",
            },
            {
                "panel": f"B_{panel_data.sensitivity_preset}",
                "pfam_union_size": panel_data.union_size,
                "ms_size": panel_data.relaxed_size,
                "overlap": panel_data.relaxed_overlap,
                "variant_label": panel_data.relaxed_label,
            },
            {
                "panel": "C_expandAll_any",
                "pfam_union_size": panel_data.union_size,
                "ms_size": len(panel_data.missed_after_strict),
                "overlap": panel_data.panelC_detected_count,
                "variant_label": panel_data.panelC_variant_label,
            },
        ]
        pd.DataFrame(summary_rows).to_csv(
            Path(figure_data_dir) / "Sup_Figure_3_panel_summary.csv",
            index=False,
        )

def _render_single_panel(panel_data: SupFigure3PanelData,
                         labels: Tuple[str, str],
                         title: str,
                         ms_size: int,
                         overlap: int,
                         out_dir: Path,
                         stem: str,
                         style: SupFigure3PanelStyle) -> Tuple[Path, Path, Path]:
    fig, ax = plt.subplots(figsize=style.single_panel_size)
    venn_label_sizes = _get_set_label_fontsizes(style)
    plot_two_set_venn(
        panel_data.union_size,
        ms_size,
        overlap,
        labels=labels,
        title=title,
        title_fontsize=style.venn_title_fontsize,
        set_label_fontsizes=venn_label_sizes,
        ax=ax,
        set_colors=style.venn_colors,
        alpha=style.venn_alpha,
    )
    fig.tight_layout()
    png = out_dir / f"{stem}.png"
    fig.savefig(png, dpi=style.dpi)
    plt.close(fig)
    return png

def generate_sup_fig3_panels(panel_data: SupFigure3PanelData,
                             out_dir: PathLike,
                             panel_style: Optional[SupFigure3PanelStyle] = None,
                             figure_data_dir: Optional[PathLike] = None,
                             panel_labels: Optional[SupFigure3PanelLabels] = None) -> Dict[str, str]:
    out_dir = _ensure_dir(out_dir)
    style = panel_style or SupFigure3PanelStyle()
    apply_sup_fig3_style(style)
    label_values = resolve_sup_fig3_panel_labels(panel_data, panel_labels)

    pA_png = _render_single_panel(
        panel_data,
        labels=(label_values["panelA_left"], label_values["panelA_right"]),
        title=label_values["panelA_title"],
        ms_size=panel_data.strict_size,
        overlap=panel_data.strict_overlap,
        out_dir=out_dir,
        stem="supp_fig3_panelA_strict",
        style=style,
    )
    pB_png = _render_single_panel(
        panel_data,
        labels=(label_values["panelB_left"], label_values["panelB_right"]),
        title=label_values["panelB_title"],
        ms_size=panel_data.relaxed_size,
        overlap=panel_data.relaxed_overlap,
        out_dir=out_dir,
        stem=f"supp_fig3_panelB_{panel_data.relaxed_suffix}",
        style=style,
    )

    fig_combo, _ = build_sup_fig3_gridspec(panel_data, style=style, panel_labels=panel_labels)
    combo_png = out_dir / f"supp_fig3_panels_{panel_data.relaxed_suffix}_gridspec.png"
    fig_combo.savefig(combo_png, dpi=style.dpi)
    plt.close(fig_combo)

    figC, axC = plt.subplots(figsize=style.panel_c_figsize)
    bars = [panel_data.panelC_detected_count, panel_data.panelC_not_detected_count]
    bar_labels = [label_values["panelC_detected"], label_values["panelC_not_detected"]]
    colors = style.panel_c_colors
    axC.barh([0, 1], bars, color=colors)
    axC.set_yticks([0, 1], bar_labels)
    axC.set_xlabel(label_values["panelC_xlabel"])
    axC.set_title(label_values["panelC_title"])
    for idx, value in enumerate(bars):
        axC.text(value, idx, f"  {value}", va="center", ha="left", fontsize=10)
    if style.panel_c_hide_box:
        for spine_name in ("top", "right"):
            spine = axC.spines.get(spine_name)
            if spine is not None:
                spine.set_visible(False)
        for spine_name in ("bottom", "left"):
            spine = axC.spines.get(spine_name)
            if spine is not None:
                spine.set_visible(True)
    figC.tight_layout()
    panelC_png = out_dir / f"supp_fig3_panelC_{panel_data.panelC_suffix}.png"
    figC.savefig(panelC_png, dpi=style.dpi)
    plt.close(figC)

    export_sup_fig3_panel_tables(panel_data, out_dir, figure_data_dir=figure_data_dir)

    return {
        "panelA_png": str(pA_png),
        "panelB_png": str(pB_png),
        "panels_gridspec_png": str(combo_png),
        "panelC_png": str(panelC_png),
        "composed_png": str(combo_png),
    }


# --------------------------
# Plotting helpers
# --------------------------

def plot_two_set_venn(U_size: int,
                      M_size: int,
                      UM_overlap: int,
                      labels: Tuple[str, str] = ("PFAM (A ∪ B)", "Mass spec (unique)"),
                      title: Optional[str] = None,
                      title_fontsize: Optional[float] = None,
                      set_label_fontsizes: Optional[Tuple[Optional[float], Optional[float]]] = None,
                      ax=None,
                      set_colors: Tuple[str, str] = ("#4C72B0", "#DD8452"),
                      alpha: float = 0.6):
    if venn2 is None:
        raise ImportError("matplotlib_venn is required; install it with `pip install matplotlib-venn`.")
    if ax is None:
        _, ax = plt.subplots(figsize=(5.5, 5.5))
    u_only = max(0, U_size - UM_overlap)
    m_only = max(0, M_size - UM_overlap)
    intersection = max(0, UM_overlap)
    venn = venn2(
        subsets=(u_only, m_only, intersection),
        set_labels=labels,
        set_colors=set_colors,
        ax=ax,
        alpha=alpha,
    )
    if title:
        if title_fontsize is not None:
            ax.set_title(title, fontsize=title_fontsize)
        else:
            ax.set_title(title)
    if set_label_fontsizes:
        set_labels = getattr(venn, "set_labels", None)
        if set_labels:
            for idx, size in enumerate(set_label_fontsizes):
                if size is None or idx >= len(set_labels):
                    continue
                label_text = set_labels[idx]
                if label_text is not None:
                    label_text.set_fontsize(size)
    return ax


# --------------------------
# Two-panel generator
# --------------------------

def make_supp_fig3_panels(ac_csv: PathLike,
                          bd_csv: PathLike,
                          zip_path: PathLike,
                          out_dir: PathLike,
                          sensitivity_preset: str = "leading_2peps",
                          dpi: int = 300,
                          panel_style: Optional[SupFigure3PanelStyle] = None,
                          figure_data_dir: Optional[PathLike] = None,
                          panel_labels: Optional[SupFigure3PanelLabels] = None) -> Dict[str, str]:
    style = panel_style or SupFigure3PanelStyle(dpi=dpi)
    panel_data = prepare_sup_fig3_panel_data(
        ac_csv=ac_csv,
        bd_csv=bd_csv,
        zip_path=zip_path,
        sensitivity_preset=sensitivity_preset,
    )
    if figure_data_dir is None:
        figure_data_dir = _ensure_dir(Path(out_dir).parent / "figure_data")
    return generate_sup_fig3_panels(
        panel_data,
        out_dir=out_dir,
        panel_style=style,
        figure_data_dir=figure_data_dir,
        panel_labels=panel_labels,
    )
