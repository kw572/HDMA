#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import pickle
import re
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd


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
            try:
                import crested  # noqa: F401
                return repo_root
            except ImportError:
                continue

    raise ImportError(
        "Unable to import `crested`. Activate an environment with CREsted installed "
        "or set CRESTED_REPO to a CREsted checkout."
    )


CRESTED_REPO_ROOT = bootstrap_crested()
import crested  # noqa: E402


def load_module_from_path(module_name: str, module_path: Path):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load module spec for {module_name} from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def load_crested_modisco_api():
    try:
        return crested.tl.modisco
    except Exception as exc:
        if CRESTED_REPO_ROOT is None:
            raise

        src_root = CRESTED_REPO_ROOT / "src"
        tl_root = src_root / "crested" / "tl"
        modisco_root = tl_root / "modisco"

        tl_pkg = sys.modules.get("crested.tl")
        if tl_pkg is None:
            tl_pkg = types.ModuleType("crested.tl")
            tl_pkg.__path__ = [str(tl_root)]
            sys.modules["crested.tl"] = tl_pkg

        modisco_pkg = sys.modules.get("crested.tl.modisco")
        if modisco_pkg is None:
            modisco_pkg = types.ModuleType("crested.tl.modisco")
            modisco_pkg.__path__ = [str(modisco_root)]
            sys.modules["crested.tl.modisco"] = modisco_pkg

        load_module_from_path(
            "crested.tl.modisco._modisco_utils",
            modisco_root / "_modisco_utils.py",
        )
        tfmodisco_module = load_module_from_path(
            "crested.tl.modisco._tfmodisco",
            modisco_root / "_tfmodisco.py",
        )
        return tfmodisco_module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a CREsted-style motif compendium from ChromBPNet TF-MoDISco outputs."
    )
    parser.add_argument("--modisco-root", required=True, type=Path)
    parser.add_argument("--keep-file", required=True, type=Path)
    parser.add_argument("--keep-file-fallback", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--dataset-filter-regex", default="", type=str)
    parser.add_argument("--sim-threshold", default=6.0, type=float)
    parser.add_argument("--trim-ic-threshold", default=0.025, type=float)
    parser.add_argument("--discard-ic-threshold", default=0.1, type=float)
    parser.add_argument("--pattern-parameter", default="seqlet_count", type=str)
    parser.add_argument("--normalize-pattern-matrix", default=0, type=int)
    parser.add_argument("--write-similarity-matrix", default=0, type=int)
    return parser.parse_args()


def read_keep_file(path: Path) -> list[str]:
    if not path.exists():
        return []
    datasets = []
    with path.open() as handle:
        for line in handle:
            fields = line.strip().split("\t")
            if fields and fields[0] and fields[0] != "NA":
                datasets.append(fields[0])
    return datasets


def discover_datasets(modisco_root: Path) -> list[str]:
    if not modisco_root.exists():
        return []
    datasets = []
    for child in sorted(modisco_root.iterdir()):
        if not child.is_dir():
            continue
        if (child / "counts_modisco_output.h5").exists():
            datasets.append(child.name)
    return datasets


def collect_datasets(args: argparse.Namespace) -> list[str]:
    datasets = read_keep_file(args.keep_file)
    if not datasets:
        datasets = read_keep_file(args.keep_file_fallback)
    if not datasets:
        datasets = discover_datasets(args.modisco_root)

    if args.dataset_filter_regex:
        pattern = re.compile(args.dataset_filter_regex)
        datasets = [dataset for dataset in datasets if pattern.search(dataset)]

    return sorted(dict.fromkeys(datasets))


def build_matched_files(
    datasets: list[str],
    modisco_root: Path,
) -> tuple[dict[str, str], list[dict[str, str]]]:
    matched_files: dict[str, str] = {}
    missing: list[dict[str, str]] = []

    for dataset in datasets:
        modisco_h5 = modisco_root / dataset / "counts_modisco_output.h5"
        if modisco_h5.exists():
            matched_files[dataset] = str(modisco_h5)
        else:
            missing.append({"dataset": dataset, "expected_path": str(modisco_h5)})

    return matched_files, missing


def pattern_indices(all_patterns: dict[str, dict]) -> list[str]:
    return sorted(all_patterns.keys(), key=lambda idx: int(idx))


def build_pattern_manifest(all_patterns: dict[str, dict]) -> pd.DataFrame:
    rows = []
    for pattern_idx in pattern_indices(all_patterns):
        pattern = all_patterns[pattern_idx]
        class_names = sorted(pattern["classes"].keys())
        class_seqlets = {
            class_name: int(pattern["classes"][class_name]["n_seqlets"])
            for class_name in class_names
        }
        rows.append(
            {
                "pattern_idx": int(pattern_idx),
                "representative_id": pattern["pattern"]["id"],
                "representative_ic": float(pattern["ic"]),
                "pos_pattern": bool(pattern.get("pos_pattern", False)),
                "n_classes": len(class_names),
                "classes": ",".join(class_names),
                "class_seqlets_json": json.dumps(class_seqlets, sort_keys=True),
                "n_instances": len(pattern["instances"]),
            }
        )
    return pd.DataFrame(rows)


def build_pattern_matrix_df(
    classes: list[str],
    pattern_matrix: np.ndarray,
    all_patterns: dict[str, dict],
) -> pd.DataFrame:
    ordered_indices = pattern_indices(all_patterns)
    if pattern_matrix.shape[1] != len(ordered_indices):
        raise ValueError(
            "Pattern matrix column count does not match the number of merged patterns."
        )
    columns = [f"pattern_{idx}" for idx in ordered_indices]
    return pd.DataFrame(pattern_matrix, index=classes, columns=columns)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    modisco_api = load_crested_modisco_api()

    datasets = collect_datasets(args)
    if not datasets:
        raise SystemExit("No datasets found from keep files or TF-MoDISco directories.")

    matched_files, missing = build_matched_files(datasets, args.modisco_root)
    if not matched_files:
        raise SystemExit("No TF-MoDISco H5 files were found for the requested datasets.")

    classes = list(matched_files.keys())

    all_patterns = modisco_api.process_patterns(
        matched_files=matched_files,
        sim_threshold=args.sim_threshold,
        trim_ic_threshold=args.trim_ic_threshold,
        discard_ic_threshold=args.discard_ic_threshold,
        verbose=True,
    )

    pattern_matrix = modisco_api.create_pattern_matrix(
        classes=classes,
        all_patterns=all_patterns,
        normalize=bool(args.normalize_pattern_matrix),
        pattern_parameter=args.pattern_parameter,
    )

    manifest_df = build_pattern_manifest(all_patterns)
    pattern_matrix_df = build_pattern_matrix_df(classes, pattern_matrix, all_patterns)

    (args.output_dir / "matched_modisco_h5.json").write_text(
        json.dumps(matched_files, indent=2) + "\n"
    )
    (args.output_dir / "missing_modisco_h5.json").write_text(
        json.dumps(missing, indent=2) + "\n"
    )
    (args.output_dir / "classes.json").write_text(json.dumps(classes, indent=2) + "\n")

    with (args.output_dir / "all_patterns.pkl").open("wb") as handle:
        pickle.dump(all_patterns, handle)

    np.save(args.output_dir / "pattern_matrix.npy", pattern_matrix)
    manifest_df.to_csv(args.output_dir / "pattern_manifest.tsv", sep="\t", index=False)
    pattern_matrix_df.to_csv(
        args.output_dir / "pattern_matrix.tsv",
        sep="\t",
        index_label="class",
    )

    summary = {
        "modisco_root": str(args.modisco_root),
        "output_dir": str(args.output_dir),
        "n_requested_datasets": len(datasets),
        "n_matched_datasets": len(matched_files),
        "n_missing_datasets": len(missing),
        "n_patterns": int(len(all_patterns)),
        "pattern_matrix_shape": list(pattern_matrix.shape),
        "sim_threshold": args.sim_threshold,
        "trim_ic_threshold": args.trim_ic_threshold,
        "discard_ic_threshold": args.discard_ic_threshold,
        "pattern_parameter": args.pattern_parameter,
        "normalize_pattern_matrix": bool(args.normalize_pattern_matrix),
    }

    if args.write_similarity_matrix:
        similarity_matrix, similarity_indices = modisco_api.calculate_similarity_matrix(
            all_patterns
        )
        np.save(args.output_dir / "pattern_similarity_matrix.npy", similarity_matrix)
        (args.output_dir / "pattern_similarity_indices.json").write_text(
            json.dumps(list(similarity_indices), indent=2) + "\n"
        )
        summary["pattern_similarity_shape"] = list(similarity_matrix.shape)

    (args.output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
