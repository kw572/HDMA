#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import os
import pickle
import sys
import types
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


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


def find_motif_to_tf_file() -> Path | None:
    env_path = os.environ.get("CRESTED_MOTIF_TO_TF_FILE")
    if env_path and Path(env_path).exists():
        return Path(env_path)

    candidate_paths = [
        Path.home() / ".cache" / "crested" / "motif_db" / "motif_tf_collection.tsv",
        Path.home() / ".local" / "share" / "crested" / "motif_db" / "motif_tf_collection.tsv",
    ]
    for candidate in candidate_paths:
        if candidate.exists():
            return candidate
    return None


def choose_annotation_columns(motif_to_tf_df: pd.DataFrame) -> list[str]:
    preferred = [
        "Direct_annot",
        "Motif_similarity_annot",
        "Orthology_annot",
        "TF",
        "Gene",
    ]
    cols = [col for col in preferred if col in motif_to_tf_df.columns]
    if cols:
        return cols
    return [col for col in motif_to_tf_df.columns if "annot" in col.lower()]


def build_annotation_table(
    manifest: pd.DataFrame,
    all_patterns: dict,
    matched_modisco_h5: dict[str, str],
    modisco_api,
    annotation_pval_threshold: float,
) -> tuple[pd.DataFrame, dict]:
    html_paths_by_class = find_report_htmls(matched_modisco_h5)
    metadata = {
        "annotation_mode": "metadata_only",
        "n_report_htmls": len(html_paths_by_class),
        "n_patterns_with_matches": 0,
    }

    annotations = manifest.copy()
    annotations["motif_matches"] = pd.NA
    annotations["tf_candidates"] = pd.NA
    annotations["annotation_status"] = "metadata_only"
    annotations["annotation_note"] = "No motif match HTML reports were found for this compendium input set."

    if not html_paths_by_class:
        return annotations, metadata

    motif_to_tf_file = find_motif_to_tf_file()
    pattern_html_paths = build_html_paths_for_patterns(all_patterns, html_paths_by_class)
    if any(path == "" for paths in pattern_html_paths for path in paths):
        annotations["annotation_status"] = "partial_html_coverage"
        annotations["annotation_note"] = (
            "Some merged pattern instances map to datasets without motifs.html, so TF annotation was skipped."
        )
        metadata["annotation_mode"] = "partial_html_coverage"
        return annotations, metadata

    pattern_match_dict = modisco_api.find_pattern_matches(
        all_patterns,
        pattern_html_paths,
        p_val_thr=annotation_pval_threshold,
    )
    metadata["n_patterns_with_matches"] = len(pattern_match_dict)

    pattern_tf_dict = {}
    if motif_to_tf_file is not None:
        motif_to_tf_df = modisco_api.read_motif_to_tf_file(str(motif_to_tf_file))
        cols = choose_annotation_columns(motif_to_tf_df)
        if cols:
            pattern_tf_dict, _ = modisco_api.create_pattern_tf_dict(
                pattern_match_dict,
                motif_to_tf_df,
                all_patterns,
                cols,
            )
            metadata["annotation_mode"] = "motif_and_tf"
            metadata["motif_to_tf_file"] = str(motif_to_tf_file)
        else:
            metadata["annotation_mode"] = "motif_only"
            metadata["motif_to_tf_file"] = str(motif_to_tf_file)
    else:
        metadata["annotation_mode"] = "motif_only"

    annotations["annotation_status"] = "no_match"
    annotations["annotation_note"] = "No motif database match passed the p-value threshold."

    for pattern_idx, match_info in pattern_match_dict.items():
        mask = annotations["pattern_idx"] == int(pattern_idx)
        matches = sorted(set(match_info.get("matches", [])))
        annotations.loc[mask, "motif_matches"] = ";".join(matches) if matches else pd.NA
        annotations.loc[mask, "annotation_status"] = "motif_match"
        annotations.loc[mask, "annotation_note"] = "Matched motif database entries from modisco report HTML."

        if pattern_idx in pattern_tf_dict:
            tfs = sorted(set(pattern_tf_dict[pattern_idx].get("tfs", [])))
            annotations.loc[mask, "tf_candidates"] = ";".join(tfs) if tfs else pd.NA
            annotations.loc[mask, "annotation_status"] = "motif_and_tf_match"
            annotations.loc[mask, "annotation_note"] = "Matched motif database entries and mapped them to TF candidates."

    return annotations, metadata


def save_clustermap(
    pattern_matrix: np.ndarray,
    classes: list[str],
    output_path: Path,
    importance_threshold: float,
    width: float,
    height: float,
) -> None:
    max_importance = np.max(np.abs(pattern_matrix), axis=0)
    keep = max_importance > importance_threshold
    filtered = pattern_matrix[:, keep]

    if filtered.size == 0:
        raise ValueError(
            "No patterns remain after importance-threshold filtering; lower CRESTED_HEATMAP_IMPORTANCE_THRESHOLD."
        )

    data = pd.DataFrame(filtered, index=classes)
    sns.set_theme(style="white")
    row_cluster = data.shape[0] > 1
    col_cluster = data.shape[1] > 1
    grid = sns.clustermap(
        data,
        cmap="coolwarm",
        center=0,
        method="average",
        figsize=(width, height),
        yticklabels=classes,
        xticklabels=True,
        dendrogram_ratio=(0.05, 0.2),
        cbar_pos=(1.05, 0.4, 0.01, 0.3),
        row_cluster=row_cluster,
        col_cluster=col_cluster,
    )
    colorbar = grid.ax_heatmap.collections[0].colorbar
    colorbar.set_label("Motif importance", rotation=270, labelpad=20)
    grid.ax_heatmap.set_xlabel("Merged pattern index")
    grid.ax_heatmap.set_ylabel("Class")
    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
    grid.fig.suptitle("CREsted-style motif compendium clustermap", y=1.02)
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
    )
    annotation_metadata["compendium_dir"] = str(args.compendium_dir)
    annotation_metadata["plots_dir"] = str(args.plots_dir)
    annotation_metadata["annotation_dir"] = str(args.annotation_dir)

    annotations.to_csv(args.annotation_dir / "pattern_annotations.tsv", sep="\t", index=False)
    with (args.annotation_dir / "annotation_summary.json").open("w") as handle:
        json.dump(annotation_metadata, handle, indent=2)

    save_clustermap(
        pattern_matrix=pattern_matrix,
        classes=classes,
        output_path=args.plots_dir / "pattern_clustermap.png",
        importance_threshold=args.importance_threshold,
        width=args.heatmap_width,
        height=args.heatmap_height,
    )

    print(json.dumps(
        {
            "annotation_table": str(args.annotation_dir / "pattern_annotations.tsv"),
            "annotation_summary": str(args.annotation_dir / "annotation_summary.json"),
            "clustermap_png": str(args.plots_dir / "pattern_clustermap.png"),
            "n_patterns": int(pattern_matrix.shape[1]),
            "n_classes": int(pattern_matrix.shape[0]),
            "annotation_mode": annotation_metadata["annotation_mode"],
        },
        indent=2,
    ))


if __name__ == "__main__":
    args = parse_args()
    main(args)
