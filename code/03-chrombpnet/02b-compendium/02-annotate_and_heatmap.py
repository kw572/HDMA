#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import os
import pickle
import re
import sys
import types
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

HOCOMOCO14_MEME_URL = (
    "https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE/"
    "formatted_motifs/H14CORE_meme_format.meme"
)
HOCOMOCO14_ANNOTATION_URL = (
    "https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE/"
    "H14CORE_annotation.jsonl"
)
JASPAR2026_VERTEBRATES_MEME_URL = (
    "https://jaspar.elixir.no/download/data/2026/CORE/"
    "JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt"
)
JASPAR2026_CORE_METADATA_URL = (
    "https://mencius.uio.no/JASPAR/JASPAR_metadata/2026/"
    "ultimate_metadata_table_CORE.tsv"
)


def bootstrap_crested() -> Path | None:
    try:
        import crested  # noqa: F401
        return None
    except ImportError:
        pass

    candidate_roots = []
    env_repo = os.environ.get("CRESTED_REPO")
    if env_repo:
        candidate_roots.append(Path(env_repo))

    script_path = Path(__file__).resolve()
    candidate_roots.extend(
        [
            script_path.parents[4] / "CREsted",
            script_path.parents[5] / "CREsted",
            Path.cwd().parent / "CREsted",
        ]
    )

    for repo_root in candidate_roots:
        src_dir = repo_root / "src"
        if src_dir.exists():
            sys.path.insert(0, str(src_dir))
            if (src_dir / "crested").exists():
                return repo_root

    raise ImportError(
        "Unable to import `crested`. Activate an environment with CREsted installed "
        "or set CRESTED_REPO to a CREsted checkout."
    )


CRESTED_REPO_ROOT = bootstrap_crested()
try:
    import crested  # type: ignore  # noqa: E402
except ImportError:
    if CRESTED_REPO_ROOT is None:
        raise

    crested = types.ModuleType("crested")
    crested.__path__ = [str(CRESTED_REPO_ROOT / "src" / "crested")]
    sys.modules["crested"] = crested


def load_module_from_path(module_name: str, module_path: Path):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load module spec for {module_name} from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def ensure_loguru() -> None:
    if "loguru" in sys.modules:
        return

    loguru_module = types.ModuleType("loguru")
    logger = logging.getLogger("crested_modisco_fallback")
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(levelname)s %(message)s"))
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    loguru_module.logger = logger
    sys.modules["loguru"] = loguru_module


def ensure_motif_runtime_stubs() -> None:
    ensure_loguru()

    if "anndata" not in sys.modules:
        anndata_module = types.ModuleType("anndata")

        def _missing_read_h5ad(*args, **kwargs):
            raise ImportError("anndata is unavailable in this motif-only runtime")

        anndata_module.read_h5ad = _missing_read_h5ad
        sys.modules["anndata"] = anndata_module

    if "scanpy" not in sys.modules:
        scanpy_module = types.ModuleType("scanpy")
        sys.modules["scanpy"] = scanpy_module

    if "crested.utils" not in sys.modules:
        utils_pkg = types.ModuleType("crested.utils")
        utils_pkg.__path__ = [str(CRESTED_REPO_ROOT / "src" / "crested" / "utils")]
        sys.modules["crested.utils"] = utils_pkg


def load_crested_modisco_api():
    try:
        return crested.tl.modisco
    except Exception:
        if CRESTED_REPO_ROOT is None:
            raise

        ensure_motif_runtime_stubs()
        src_root = CRESTED_REPO_ROOT / "src"
        tl_root = src_root / "crested" / "tl"
        modisco_root = tl_root / "modisco"
        utils_root = src_root / "crested" / "utils"

        if "crested.tl" not in sys.modules:
            tl_pkg = types.ModuleType("crested.tl")
            tl_pkg.__path__ = [str(tl_root)]
            sys.modules["crested.tl"] = tl_pkg

        if "crested.tl.modisco" not in sys.modules:
            modisco_pkg = types.ModuleType("crested.tl.modisco")
            modisco_pkg.__path__ = [str(modisco_root)]
            sys.modules["crested.tl.modisco"] = modisco_pkg

        load_module_from_path("crested.utils._logging", utils_root / "_logging.py")
        load_module_from_path(
            "crested.tl.modisco._modisco_utils",
            modisco_root / "_modisco_utils.py",
        )
        return load_module_from_path(
            "crested.tl.modisco._tfmodisco",
            modisco_root / "_tfmodisco.py",
        )


def _read_html_table_with_bs4(source: str) -> pd.DataFrame:
    from bs4 import BeautifulSoup

    with open(source, encoding="utf-8") as handle:
        soup = BeautifulSoup(handle.read(), "html.parser")

    table = soup.find("table")
    if table is None:
        raise ValueError(f"No table found in HTML file: {source}")

    header_cells = table.find_all("th")
    headers = [cell.get_text(strip=True) for cell in header_cells]

    rows = []
    for row in table.find_all("tr"):
        cells = row.find_all("td")
        if not cells:
            continue
        values = [cell.get_text(strip=True) for cell in cells]
        rows.append(values)

    if not rows:
        raise ValueError(f"No data rows found in HTML file: {source}")

    width = max(len(headers), max(len(row) for row in rows))
    if not headers:
        headers = [f"col_{idx}" for idx in range(width)]
    elif len(headers) < width:
        headers.extend(f"extra_{idx}" for idx in range(len(headers), width))

    padded_rows = [row + [""] * (width - len(row)) for row in rows]
    df = pd.DataFrame(padded_rows, columns=headers[:width])
    numeric_prefixes = ("qval", "pval", "num_seqlets")
    for col in df.columns:
        if col.startswith(numeric_prefixes):
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def install_html_parser_fallback(modisco_api) -> None:
    try:
        import lxml  # noqa: F401
        return
    except ImportError:
        pass

    def _fallback_reader(source: str):
        try:
            return pd.read_html(source)[0]
        except ImportError:
            return _read_html_table_with_bs4(source)
        except ValueError as exc:
            return f"Error: {exc}"

    modisco_api.read_html_to_dataframe = _fallback_reader
    modisco_utils = sys.modules.get("crested.tl.modisco._modisco_utils")
    if modisco_utils is not None:
        modisco_utils.read_html_to_dataframe = _fallback_reader


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate CREsted-style compendium motifs and generate a clustermap."
    )
    parser.add_argument("--compendium-dir", required=True, type=Path)
    parser.add_argument("--modisco-root", required=True, type=Path)
    parser.add_argument("--annotation-dir", required=True, type=Path)
    parser.add_argument("--plots-dir", required=True, type=Path)
    parser.add_argument("--importance-threshold", default=0.0, type=float)
    parser.add_argument("--heatmap-width", default=25.0, type=float)
    parser.add_argument("--heatmap-height", default=8.0, type=float)
    parser.add_argument("--annotation-pval-threshold", default=0.05, type=float)
    parser.add_argument("--annotation-min-score", default=6.0, type=float)
    return parser.parse_args()


def load_json(path: Path):
    with path.open() as handle:
        return json.load(handle)


def load_compendium(compendium_dir: Path) -> tuple[dict, np.ndarray, list[str], dict[str, str], pd.DataFrame]:
    with (compendium_dir / "all_patterns.pkl").open("rb") as handle:
        all_patterns = pickle.load(handle)
    pattern_matrix = np.load(compendium_dir / "pattern_matrix.npy")
    classes = load_json(compendium_dir / "classes.json")
    matched_modisco_h5 = load_json(compendium_dir / "matched_modisco_h5.json")
    manifest = pd.read_csv(compendium_dir / "pattern_manifest.tsv", sep="\t")
    return all_patterns, pattern_matrix, classes, matched_modisco_h5, manifest


def find_report_htmls(matched_modisco_h5: dict[str, str]) -> dict[str, Path]:
    html_paths: dict[str, Path] = {}
    for class_name, h5_path in matched_modisco_h5.items():
        h5_file = Path(h5_path)
        candidates = [
            h5_file.parent / "counts_modisco_report" / "motifs.html",
            h5_file.parent / f"{h5_file.stem}_report" / "motifs.html",
            h5_file.parent / "report" / "motifs.html",
        ]
        for candidate in candidates:
            if candidate.exists():
                html_paths[class_name] = candidate
                break
    return html_paths


def extract_class_name_from_instance_id(instance_id: str) -> str:
    parts = instance_id.split("_")
    if len(parts) < 4:
        return instance_id
    return "_".join(parts[:-3])


def build_html_paths_for_patterns(all_patterns: dict, html_paths_by_class: dict[str, Path]) -> list[list[str]]:
    pattern_html_paths: list[list[str]] = []
    for pattern_idx in all_patterns:
        html_paths = []
        for instance_key in all_patterns[pattern_idx]["instances"]:
            instance_id = all_patterns[pattern_idx]["instances"][instance_key]["id"]
            class_name = extract_class_name_from_instance_id(instance_id)
            html_path = html_paths_by_class.get(class_name)
            html_paths.append(str(html_path) if html_path is not None else "")
        pattern_html_paths.append(html_paths)
    return pattern_html_paths


def download_if_missing(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    urllib.request.urlretrieve(url, dest)
    return dest


def load_hocomoco14_annotations(annotation_jsonl: Path) -> dict[str, dict]:
    mapping: dict[str, dict] = {}
    with annotation_jsonl.open() as handle:
        for line in handle:
            record = json.loads(line)
            motif_name = record["name"]
            tf_name = record.get("tf")
            masterlist_info = record.get("masterlist_info", {})
            species_info = masterlist_info.get("species", {}).get("HUMAN", {})
            human_gene_symbol = (
                species_info.get("gene_symbol")
            )
            synonyms = species_info.get("gene_synonyms", [])
            mapping[motif_name] = {
                "tf": tf_name,
                "gene_symbol": human_gene_symbol,
                "gene_synonyms": synonyms,
                "tfclass_family": masterlist_info.get("tfclass_family"),
                "tfclass_subfamily": masterlist_info.get("tfclass_subfamily"),
                "tfclass_class": masterlist_info.get("tfclass_class"),
                "tfclass_superclass": masterlist_info.get("tfclass_superclass"),
                "collection": record.get("collection"),
                "quality": record.get("quality"),
            }
    return mapping


def load_jaspar2026_metadata(metadata_tsv: Path) -> dict[str, dict]:
    table = pd.read_csv(metadata_tsv, sep="\t")
    mapping: dict[str, dict] = {}
    for _, row in table.iterrows():
        if str(row.get("collection", "")).strip().upper() != "CORE":
            continue
        matrix_id = str(row["matrix_id"]).strip()
        if not matrix_id:
            continue
        tf_name = str(row.get("name", "")).strip() or None
        family = str(row.get("family", "")).strip() or None
        tf_class = str(row.get("class", "")).strip() or None
        species = str(row.get("species", "")).strip() or None
        tax_group = str(row.get("tax_group", "")).strip() or None
        mapping[matrix_id] = {
            "tf": tf_name,
            "gene_symbol": tf_name,
            "gene_synonyms": [],
            "tfclass_family": family,
            "tfclass_subfamily": None,
            "tfclass_class": tf_class,
            "tfclass_superclass": None,
            "collection": row.get("collection"),
            "quality": None,
            "species": species,
            "tax_group": tax_group,
        }
    return mapping


def parse_meme_database(meme_path: Path) -> list[dict]:
    motifs: list[dict] = []
    current_id: str | None = None
    current_name: str | None = None
    current_ppm: list[list[float]] = []
    expected_rows = 0

    def finalize_current():
        nonlocal current_id, current_name, current_ppm, expected_rows
        if current_id is None or not current_ppm:
            return
        ppm = np.array(current_ppm, dtype=float)
        motifs.append(
            {
                "id": current_id,
                "name": current_name or current_id,
                "ppm": ppm,
            }
        )
        current_id = None
        current_name = None
        current_ppm = []
        expected_rows = 0

    with meme_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("MOTIF "):
                finalize_current()
                parts = line.split(maxsplit=2)
                current_id = parts[1] if len(parts) > 1 else None
                current_name = parts[2] if len(parts) > 2 else current_id
                continue
            if line.startswith("letter-probability matrix:"):
                width_match = re.search(r"\bw=\s*(\d+)", line)
                if width_match:
                    expected_rows = int(width_match.group(1))
                continue
            if current_id is not None and expected_rows > 0:
                if line.startswith("URL "):
                    finalize_current()
                    continue
                row = [float(x) for x in line.split()]
                if len(row) == 4:
                    current_ppm.append(row)
                    if len(current_ppm) == expected_rows:
                        finalize_current()

    finalize_current()
    return motifs


def score_hocomoco14_matches(
    all_patterns: dict,
    modisco_utils_module,
    hocomoco_patterns: list[dict],
    min_score: float,
    top_n: int = 5,
) -> dict[int, dict[str, list[str]]]:
    pattern_match_dict: dict[int, dict[str, list[str]]] = {}
    hocomoco_ids = [motif["id"] for motif in hocomoco_patterns]

    for pattern_idx, pattern_info in all_patterns.items():
        representative = pattern_info["pattern"]
        scores = modisco_utils_module.match_score_patterns(
            representative,
            hocomoco_patterns,
            use_ppm=True,
        ).reshape(-1)
        ranked = sorted(
            zip(hocomoco_ids, scores, strict=False),
            key=lambda pair: float(pair[1]),
            reverse=True,
        )
        filtered = [(motif_id, float(score)) for motif_id, score in ranked if float(score) >= min_score][:top_n]
        if filtered:
            pattern_match_dict[int(pattern_idx)] = {
                "matches": [motif_id for motif_id, _ in filtered],
                "scores": [score for _, score in filtered],
            }

    return pattern_match_dict


def build_annotation_table(
    manifest: pd.DataFrame,
    all_patterns: dict,
    matched_modisco_h5: dict[str, str],
    modisco_api,
    annotation_pval_threshold: float,
    annotation_min_score: float,
    annotation_dir: Path,
) -> tuple[pd.DataFrame, dict]:
    del matched_modisco_h5, annotation_pval_threshold
    annotation_db = os.environ.get("CRESTED_ANNOTATION_DB", "hocomoco14").strip().lower()
    metadata = {
        "annotation_mode": "metadata_only",
        "n_patterns_with_matches": 0,
    }

    annotations = manifest.copy()
    annotations["motif_matches"] = pd.NA
    annotations["motif_match_scores"] = pd.NA
    annotations["tf_family_candidates"] = pd.NA
    annotations["tf_candidates"] = pd.NA
    annotations["annotation_status"] = "metadata_only"
    annotations["annotation_note"] = "No HOCOMOCO v14 match exceeded the annotation score threshold."

    if annotation_db == "jaspar2026":
        cache_dir = annotation_dir / "jaspar2026_cache"
        meme_path = download_if_missing(
            os.environ.get("JASPAR2026_MEME_URL", JASPAR2026_VERTEBRATES_MEME_URL),
            cache_dir / "JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt",
        )
        metadata_path = download_if_missing(
            os.environ.get("JASPAR2026_METADATA_URL", JASPAR2026_CORE_METADATA_URL),
            cache_dir / "JASPAR2026_CORE_metadata.tsv",
        )
        annotation_meta = load_jaspar2026_metadata(metadata_path)
        annotation_patterns = parse_meme_database(meme_path)
        metadata["annotation_mode"] = "jaspar2026"
        metadata["jaspar2026_meme_url"] = os.environ.get("JASPAR2026_MEME_URL", JASPAR2026_VERTEBRATES_MEME_URL)
        metadata["jaspar2026_metadata_url"] = os.environ.get("JASPAR2026_METADATA_URL", JASPAR2026_CORE_METADATA_URL)
        metadata["annotation_min_score"] = annotation_min_score
        metadata["jaspar2026_meme_file"] = str(meme_path)
        metadata["jaspar2026_metadata_file"] = str(metadata_path)
        metadata["n_jaspar2026_motifs"] = len(annotation_patterns)
    else:
        cache_dir = annotation_dir / "hocomoco14_cache"
        meme_path = download_if_missing(
            os.environ.get("HOCOMOCO14_MEME_URL", HOCOMOCO14_MEME_URL),
            cache_dir / "H14CORE_meme_format.meme",
        )
        annotation_jsonl = download_if_missing(
            os.environ.get("HOCOMOCO14_ANNOTATION_URL", HOCOMOCO14_ANNOTATION_URL),
            cache_dir / "H14CORE_annotation.jsonl",
        )
        annotation_meta = load_hocomoco14_annotations(annotation_jsonl)
        annotation_patterns = parse_meme_database(meme_path)
        metadata["annotation_mode"] = "hocomoco14"
        metadata["hocomoco14_meme_url"] = os.environ.get("HOCOMOCO14_MEME_URL", HOCOMOCO14_MEME_URL)
        metadata["hocomoco14_annotation_url"] = os.environ.get("HOCOMOCO14_ANNOTATION_URL", HOCOMOCO14_ANNOTATION_URL)
        metadata["annotation_min_score"] = annotation_min_score
        metadata["hocomoco14_meme_file"] = str(meme_path)
        metadata["hocomoco14_annotation_file"] = str(annotation_jsonl)
        metadata["n_hocomoco14_motifs"] = len(annotation_patterns)

    modisco_utils_module = sys.modules["crested.tl.modisco._modisco_utils"]
    pattern_match_dict = score_hocomoco14_matches(
        all_patterns=all_patterns,
        modisco_utils_module=modisco_utils_module,
        hocomoco_patterns=annotation_patterns,
        min_score=annotation_min_score,
    )

    metadata["n_patterns_with_matches"] = len(pattern_match_dict)

    annotations["annotation_status"] = "no_match"
    annotations["annotation_note"] = "No motif-database match exceeded the annotation score threshold."

    for pattern_idx, match_info in pattern_match_dict.items():
        mask = annotations["pattern_idx"] == int(pattern_idx)
        matches = sorted(set(match_info.get("matches", [])))
        scores = match_info.get("scores", [])
        annotations.loc[mask, "motif_matches"] = ";".join(matches) if matches else pd.NA
        annotations.loc[mask, "motif_match_scores"] = (
            ";".join(f"{score:.3f}" for score in scores) if scores else pd.NA
        )

        tf_family_candidates: list[str] = []
        tf_candidates: list[str] = []
        for motif_name in matches:
            motif_meta = annotation_meta.get(motif_name, {})
            for family_key in [
                "tfclass_family",
                "tfclass_subfamily",
                "tfclass_class",
                "tfclass_superclass",
            ]:
                family_name = motif_meta.get(family_key)
                if family_name:
                    tf_family_candidates.append(str(family_name))
            for candidate in [
                motif_meta.get("gene_symbol"),
                motif_meta.get("tf"),
                *motif_meta.get("gene_synonyms", []),
            ]:
                if candidate:
                    tf_candidates.append(candidate)

        tf_family_candidates = sorted(dict.fromkeys(tf_family_candidates))
        tf_candidates = sorted(dict.fromkeys(tf_candidates))
        annotations.loc[mask, "tf_family_candidates"] = (
            ";".join(tf_family_candidates) if tf_family_candidates else pd.NA
        )
        annotations.loc[mask, "tf_candidates"] = (
            ";".join(tf_candidates) if tf_candidates else pd.NA
        )
        annotations.loc[mask, "annotation_status"] = (
            "hocomoco14_family_match"
            if tf_family_candidates
            else "hocomoco14_tf_match"
            if tf_candidates
            else "hocomoco14_motif_match"
        )
        annotations.loc[mask, "annotation_note"] = (
            "Matched representative motif directly against the configured motif database and expanded TF family and TF names from the official release metadata."
        )

        if annotation_db == "jaspar2026":
            annotations.loc[mask, "annotation_status"] = (
                "jaspar2026_family_match"
                if tf_family_candidates
                else "jaspar2026_tf_match"
                if tf_candidates
                else "jaspar2026_motif_match"
            )
            annotations.loc[mask, "annotation_note"] = (
                "Matched representative motif directly against JASPAR 2026 CORE vertebrate PFMs and expanded TF family and TF names from the official JASPAR 2026 CORE metadata table."
            )
        else:
            annotations.loc[mask, "annotation_status"] = (
                "hocomoco14_family_match"
                if tf_family_candidates
                else "hocomoco14_tf_match"
                if tf_candidates
                else "hocomoco14_motif_match"
            )
            annotations.loc[mask, "annotation_note"] = (
                "Matched representative motif directly against HOCOMOCO v14 H14CORE and expanded TF family and TF names from the official annotation JSONL."
            )

    return annotations, metadata


def build_pattern_label(row: pd.Series) -> str:
    pattern_idx = int(row["pattern_idx"])
    tf_family_candidates = row.get("tf_family_candidates")
    tf_candidates = row.get("tf_candidates")
    motif_matches = row.get("motif_matches")

    if pd.notna(tf_family_candidates) and str(tf_family_candidates).strip():
        primary = str(tf_family_candidates).split(";")[0]
        return f"{primary} [p{pattern_idx}]"

    if pd.notna(tf_candidates) and str(tf_candidates).strip():
        primary = str(tf_candidates).split(";")[0]
        return f"{primary} [p{pattern_idx}]"

    if pd.notna(motif_matches) and str(motif_matches).strip():
        primary = str(motif_matches).split(";")[0]
        return f"{primary} [p{pattern_idx}]"

    representative_id = str(row.get("representative_id", f"pattern_{pattern_idx}"))
    return f"{representative_id} [p{pattern_idx}]"


def build_labeled_matrix(
    pattern_matrix: np.ndarray,
    classes: list[str],
    annotations: pd.DataFrame,
    importance_threshold: float,
) -> tuple[pd.DataFrame, list[int]]:
    max_importance = np.max(np.abs(pattern_matrix), axis=0)
    keep = max_importance > importance_threshold
    filtered = pattern_matrix[:, keep]

    if filtered.size == 0:
        raise ValueError(
            "No patterns remain after importance-threshold filtering; lower CRESTED_HEATMAP_IMPORTANCE_THRESHOLD."
        )

    kept_pattern_indices = np.where(keep)[0]
    annotation_lookup = (
        annotations.set_index("pattern_idx", drop=False)
        .reindex(kept_pattern_indices)
        .reset_index(drop=True)
    )
    columns = [build_pattern_label(row) for _, row in annotation_lookup.iterrows()]
    return pd.DataFrame(filtered, index=classes, columns=columns), kept_pattern_indices.tolist()


def zscore_rows(data: pd.DataFrame) -> pd.DataFrame:
    values = data.to_numpy(dtype=float)
    means = values.mean(axis=1, keepdims=True)
    stds = values.std(axis=1, keepdims=True)
    stds[stds == 0] = 1.0
    zvalues = (values - means) / stds
    return pd.DataFrame(zvalues, index=data.index, columns=data.columns)


def ppm_to_ic_matrix(ppm: np.ndarray) -> np.ndarray:
    ppm = np.asarray(ppm, dtype=float)
    ppm = np.clip(ppm, 1e-9, 1.0)
    entropy = -(ppm * np.log2(ppm)).sum(axis=1, keepdims=True)
    info = 2.0 - entropy
    return ppm * info


def add_pwm_logos_to_clustermap(
    grid,
    ppms: list[np.ndarray],
    logo_height_fraction: float = 0.28,
    logo_y_padding: float = 0.18,
) -> None:
    import logomaker

    col_order = (
        grid.dendrogram_col.reordered_ind
        if getattr(grid, "dendrogram_col", None) is not None and grid.dendrogram_col is not None
        else list(range(len(ppms)))
    )
    ordered_ppms = [ppms[idx] for idx in col_order]

    fig = grid.fig
    heatmap_pos = grid.ax_heatmap.get_position()
    n = len(ordered_ppms)
    if n == 0:
        return

    logo_width = heatmap_pos.width / n
    logo_height = heatmap_pos.height * logo_height_fraction
    logo_y = heatmap_pos.y0 - logo_height * (1.0 + logo_y_padding)

    fig.set_size_inches(fig.get_size_inches()[0], fig.get_size_inches()[1] * (1 + logo_height_fraction))

    for i, ppm in enumerate(ordered_ppms):
        x0 = heatmap_pos.x0 + logo_width * i
        ax = fig.add_axes([x0, logo_y, logo_width, logo_height])
        ic_matrix = ppm_to_ic_matrix(ppm)
        logo_df = pd.DataFrame(ic_matrix, columns=["A", "C", "G", "T"])
        logomaker.Logo(logo_df, ax=ax)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)


def add_pwm_logos_to_clustermap_rows(
    grid,
    ppms: list[np.ndarray],
    logo_width_fraction: float = 0.22,
    logo_x_padding: float = 0.04,
) -> None:
    import logomaker

    row_order = (
        grid.dendrogram_row.reordered_ind
        if getattr(grid, "dendrogram_row", None) is not None and grid.dendrogram_row is not None
        else list(range(len(ppms)))
    )
    ordered_ppms = [ppms[idx] for idx in row_order]

    fig = grid.fig
    heatmap_pos = grid.ax_heatmap.get_position()
    n = len(ordered_ppms)
    if n == 0:
        return

    logo_height = heatmap_pos.height / n
    logo_width = heatmap_pos.width * logo_width_fraction
    logo_x = heatmap_pos.x0 - logo_width * (1.0 + logo_x_padding)

    fig.set_size_inches(fig.get_size_inches()[0] * (1 + logo_width_fraction), fig.get_size_inches()[1])

    for i, ppm in enumerate(ordered_ppms):
        y0 = heatmap_pos.y0 + heatmap_pos.height - logo_height * (i + 1)
        ax = fig.add_axes([logo_x, y0, logo_width, logo_height])
        ic_matrix = ppm_to_ic_matrix(ppm)
        logo_df = pd.DataFrame(ic_matrix, columns=["A", "C", "G", "T"])
        logomaker.Logo(logo_df, ax=ax)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)


def save_clustermap(
    data: pd.DataFrame,
    output_path: Path,
    width: float,
    height: float,
    title: str,
    colorbar_label: str,
    center: float = 0,
    pattern_ppms: list[np.ndarray] | None = None,
    logo_axis: str = "column",
) -> None:
    sns.set_theme(style="white")
    row_cluster = data.shape[0] > 1
    col_cluster = data.shape[1] > 1
    width = min(width, max(14.0, data.shape[1] * 1.2))
    grid = sns.clustermap(
        data,
        cmap="coolwarm",
        center=center,
        method="average",
        figsize=(width, height),
        yticklabels=data.index.tolist(),
        xticklabels=data.columns.tolist(),
        dendrogram_ratio=(0.05, 0.2),
        cbar_pos=(1.05, 0.4, 0.01, 0.3),
        row_cluster=row_cluster,
        col_cluster=col_cluster,
    )
    colorbar = grid.ax_heatmap.collections[0].colorbar
    colorbar.set_label(colorbar_label, rotation=270, labelpad=20)
    grid.ax_heatmap.set_xlabel("Class")
    grid.ax_heatmap.set_ylabel("Annotated motif")
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
    grid.ax_heatmap.set_xticklabels(
        grid.ax_heatmap.get_xticklabels(),
        rotation=90,
        ha="center",
        fontsize=11,
    )
    if pattern_ppms is not None:
        if logo_axis == "row":
            add_pwm_logos_to_clustermap_rows(grid, pattern_ppms)
        else:
            add_pwm_logos_to_clustermap(grid, pattern_ppms)
    grid.fig.suptitle(title, y=1.02)
    grid.savefig(output_path, bbox_inches="tight")
    plt.close(grid.fig)


def main(args: argparse.Namespace) -> None:
    all_patterns, pattern_matrix, classes, matched_modisco_h5, manifest = load_compendium(args.compendium_dir)
    args.annotation_dir.mkdir(parents=True, exist_ok=True)
    args.plots_dir.mkdir(parents=True, exist_ok=True)

    modisco_api = load_crested_modisco_api()
    install_html_parser_fallback(modisco_api)
    annotations, annotation_metadata = build_annotation_table(
        manifest=manifest,
        all_patterns=all_patterns,
        matched_modisco_h5=matched_modisco_h5,
        modisco_api=modisco_api,
        annotation_pval_threshold=args.annotation_pval_threshold,
        annotation_min_score=args.annotation_min_score,
        annotation_dir=args.annotation_dir,
    )
    annotation_metadata["compendium_dir"] = str(args.compendium_dir)
    annotation_metadata["plots_dir"] = str(args.plots_dir)
    annotation_metadata["annotation_dir"] = str(args.annotation_dir)

    annotations.to_csv(args.annotation_dir / "pattern_annotations.tsv", sep="\t", index=False)
    with (args.annotation_dir / "annotation_summary.json").open("w") as handle:
        json.dump(annotation_metadata, handle, indent=2)

    counts_df, kept_pattern_indices = build_labeled_matrix(
        pattern_matrix=pattern_matrix,
        classes=classes,
        annotations=annotations,
        importance_threshold=args.importance_threshold,
    )
    pattern_ppms = [
        np.asarray(all_patterns[str(pattern_idx)]["pattern"]["ppm"], dtype=float)
        for pattern_idx in kept_pattern_indices
    ]
    counts_df_t = counts_df.T
    counts_df_t.to_csv(args.plots_dir / "pattern_matrix_annotated_counts.tsv", sep="\t")

    zscore_df_t = zscore_rows(counts_df_t)
    zscore_df_t.to_csv(args.plots_dir / "pattern_matrix_annotated_znorm.tsv", sep="\t")

    save_clustermap(
        data=counts_df_t,
        output_path=args.plots_dir / "pattern_clustermap_counts_annotated.png",
        width=args.heatmap_width,
        height=max(args.heatmap_height, len(counts_df_t.index) * 0.45),
        title="CREsted-style motif compendium clustermap (counts, annotated, transposed)",
        colorbar_label="Motif count / importance",
        center=0,
        pattern_ppms=pattern_ppms,
        logo_axis="row",
    )
    save_clustermap(
        data=zscore_df_t,
        output_path=args.plots_dir / "pattern_clustermap_znorm_annotated.png",
        width=args.heatmap_width,
        height=max(args.heatmap_height, len(zscore_df_t.index) * 0.45),
        title="CREsted-style motif compendium clustermap (motif-row z-score, annotated, transposed)",
        colorbar_label="Motif-row z-score",
        center=0,
        pattern_ppms=pattern_ppms,
        logo_axis="row",
    )

    print(json.dumps(
        {
            "annotation_table": str(args.annotation_dir / "pattern_annotations.tsv"),
            "annotation_summary": str(args.annotation_dir / "annotation_summary.json"),
            "counts_tsv": str(args.plots_dir / "pattern_matrix_annotated_counts.tsv"),
            "znorm_tsv": str(args.plots_dir / "pattern_matrix_annotated_znorm.tsv"),
            "counts_clustermap_png": str(args.plots_dir / "pattern_clustermap_counts_annotated.png"),
            "znorm_clustermap_png": str(args.plots_dir / "pattern_clustermap_znorm_annotated.png"),
            "n_patterns": int(pattern_matrix.shape[1]),
            "n_classes": int(pattern_matrix.shape[0]),
            "annotation_mode": annotation_metadata["annotation_mode"],
        },
        indent=2,
    ))


if __name__ == "__main__":
    args = parse_args()
    main(args)
