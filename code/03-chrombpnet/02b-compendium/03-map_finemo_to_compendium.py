#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
import pickle
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Attach 02b compendium motif-group annotations to Fi-NeMo hits."
    )
    parser.add_argument("--compendium-dir", required=True, type=Path)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--hits-tsv", required=True, type=Path)
    parser.add_argument("--output-tsv-gz", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    return parser.parse_args()


def first_semicolon_value(value: str | None, fallback: str) -> str:
    if value is None:
        return fallback
    value = str(value).strip()
    if not value:
        return fallback
    return value.split(";")[0]


def build_pattern_name(annotation_row: dict[str, str]) -> str:
    return first_semicolon_value(
        annotation_row.get("tf_candidates"),
        first_semicolon_value(
            annotation_row.get("motif_matches"),
            annotation_row.get("representative_id", "unknown_pattern"),
        ),
    )


def build_pattern_family(annotation_row: dict[str, str]) -> str:
    return first_semicolon_value(
        annotation_row.get("tf_family_candidates"),
        build_pattern_name(annotation_row),
    )


def parse_instance_id(instance_id: str) -> tuple[str, str, str]:
    if "_pos_patterns_" in instance_id:
        dataset, suffix = instance_id.rsplit("_pos_patterns_", 1)
        return dataset, "pos_patterns", suffix
    if "_neg_patterns_" in instance_id:
        dataset, suffix = instance_id.rsplit("_neg_patterns_", 1)
        return dataset, "neg_patterns", suffix
    raise ValueError(f"Unexpected instance id format: {instance_id}")


def load_annotation_rows(annotation_tsv: Path) -> dict[int, dict[str, str]]:
    with annotation_tsv.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = {}
        for row in reader:
            rows[int(row["pattern_idx"])] = row
    return rows


def build_alias_map(compendium_dir: Path) -> dict[tuple[str, str], dict[str, str]]:
    with (compendium_dir / "all_patterns.pkl").open("rb") as handle:
        all_patterns = pickle.load(handle)

    annotation_rows = load_annotation_rows(compendium_dir / "annotation" / "pattern_annotations.tsv")
    alias_map: dict[tuple[str, str], dict[str, str]] = {}

    for pattern_idx_str, pattern_info in all_patterns.items():
        pattern_idx = int(pattern_idx_str)
        annotation_row = annotation_rows[pattern_idx]
        pattern_name = build_pattern_name(annotation_row)
        pattern_family = build_pattern_family(annotation_row)
        pattern_label = f"{pattern_family} [p{pattern_idx}]"

        for instance_id in pattern_info["instances"]:
            dataset, pattern_group, pattern_number = parse_instance_id(instance_id)
            aliases = {
                f"{pattern_group}.pattern_{pattern_number}",
                f"{pattern_group}_pattern_{pattern_number}",
                f"{pattern_group}_{pattern_number}",
                f"{pattern_group}/pattern_{pattern_number}",
                f"pattern_{pattern_number}",
            }
            record = {
                "dataset": dataset,
                "motif_instance_id": instance_id,
                "compendium_pattern_idx": str(pattern_idx),
                "compendium_pattern_id": f"p{pattern_idx}",
                "compendium_pattern_name": pattern_name,
                "compendium_pattern_family": pattern_family,
                "compendium_pattern_label": pattern_label,
                "representative_id": annotation_row.get("representative_id", ""),
                "annotation_status": annotation_row.get("annotation_status", ""),
                "motif_matches": annotation_row.get("motif_matches", ""),
                "tf_candidates": annotation_row.get("tf_candidates", ""),
                "tf_family_candidates": annotation_row.get("tf_family_candidates", ""),
            }
            for alias in aliases:
                alias_map.setdefault((dataset, alias), record)

    return alias_map


def annotate_hits(
    hits_tsv: Path,
    dataset: str,
    alias_map: dict[tuple[str, str], dict[str, str]],
    output_tsv_gz: Path,
) -> dict[str, int]:
    output_tsv_gz.parent.mkdir(parents=True, exist_ok=True)
    matched = 0
    total = 0

    with hits_tsv.open() as in_handle, gzip.open(output_tsv_gz, "wt", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        extra_fields = [
            "motif_lookup_key",
            "motif_instance_id",
            "compendium_pattern_idx",
            "compendium_pattern_id",
            "compendium_pattern_name",
            "compendium_pattern_family",
            "compendium_pattern_label",
            "representative_id",
            "annotation_status",
            "motif_matches",
            "tf_candidates",
            "tf_family_candidates",
        ]
        writer = csv.DictWriter(
            out_handle,
            delimiter="\t",
            fieldnames=fieldnames + extra_fields,
            extrasaction="ignore",
        )
        writer.writeheader()

        for row in reader:
            total += 1
            motif_lookup = row.get("motif_name_orig") or row.get("motif_name") or ""
            row["motif_lookup_key"] = motif_lookup
            match = alias_map.get((dataset, motif_lookup))
            if match is not None:
                matched += 1
                row.update(match)
            else:
                for key in extra_fields[1:]:
                    row.setdefault(key, "")
            writer.writerow(row)

    return {
        "dataset": dataset,
        "n_hits": total,
        "n_hits_mapped_to_compendium": matched,
        "n_hits_unmapped": total - matched,
    }


def main() -> None:
    args = parse_args()
    alias_map = build_alias_map(args.compendium_dir)
    summary = annotate_hits(
        hits_tsv=args.hits_tsv,
        dataset=args.dataset,
        alias_map=alias_map,
        output_tsv_gz=args.output_tsv_gz,
    )
    with args.summary_json.open("w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
