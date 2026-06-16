#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path


def sanitize_variants(barcode: str):
    values = {barcode, barcode.replace("#", "_"), barcode.replace("_", "#")}
    if "#" in barcode:
        values.add(barcode.split("#")[-1])
    if "_" in barcode:
        values.add(barcode.split("_")[-1])
    return [value for value in values if value]


def load_reference_seqlevels(chromsizes_path: Path):
    seqlevels = []
    with chromsizes_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0]:
                seqlevels.append(parts[0])
    return seqlevels


def normalize_chrom_name(chrom: str, reference_seqlevels: set[str]):
    if not reference_seqlevels:
        return chrom
    if chrom in reference_seqlevels:
        return chrom
    alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    if alt in reference_seqlevels:
        return alt
    return None


def build_barcode_maps(mapping_path: Path):
    exact = {}
    bare = {}
    with mapping_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    for idx, row in enumerate(rows):
        barcode = row["cell_name"]
        for variant in sanitize_variants(barcode):
            exact.setdefault(variant, set()).add(idx)
        bare_key = barcode.split("#")[-1].split("_")[-1]
        bare.setdefault(bare_key, set()).add(idx)

    return rows, exact, {key: values for key, values in bare.items() if len(values) == 1}


def resolve_index(barcode: str, exact: dict[str, set[int]], bare: dict[str, set[int]]):
    candidates = set()
    hit = exact.get(barcode)
    if hit:
        candidates.update(hit)
    bare_hit = bare.get(barcode)
    if bare_hit:
        candidates.update(bare_hit)
    if len(candidates) == 1:
        return next(iter(candidates))
    return None


def open_fragment_input(path: Path):
    with path.open("rb") as handle:
        magic = handle.read(2)
    if path.suffix in {".gz", ".bgz"} or magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fragments", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--chromsizes", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--stats-out", required=True)
    args = parser.parse_args()

    fragments_path = Path(args.fragments).expanduser().resolve()
    mapping_path = Path(args.mapping).expanduser().resolve()
    chromsizes_path = Path(args.chromsizes).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    stats_out = Path(args.stats_out).expanduser().resolve()

    outdir.mkdir(parents=True, exist_ok=True)
    stats_out.parent.mkdir(parents=True, exist_ok=True)

    rows, exact_map, bare_map = build_barcode_maps(mapping_path)
    reference_seqlevels = set(load_reference_seqlevels(chromsizes_path))

    writers = {}
    counts = {
        "total": 0,
        "written": 0,
        "unmatched": 0,
        "malformed": 0,
        "dropped_chrom": 0,
        "normalized_chrom": 0,
    }

    def get_writer(stem: str):
        writer = writers.get(stem)
        if writer is None:
            writer = (outdir / f"{stem}.tsv").open("w", encoding="utf-8")
            writers[stem] = writer
        return writer

    try:
        with open_fragment_input(fragments_path) as handle:
            for line in handle:
                counts["total"] += 1
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    counts["malformed"] += 1
                    continue

                chrom = normalize_chrom_name(fields[0], reference_seqlevels)
                if chrom is None:
                    counts["dropped_chrom"] += 1
                    continue
                if chrom != fields[0]:
                    counts["normalized_chrom"] += 1
                fields[0] = chrom

                idx = resolve_index(fields[3], exact_map, bare_map)
                if idx is None:
                    counts["unmatched"] += 1
                    continue

                stem = rows[idx]["output_stem"]
                writer = get_writer(stem)
                writer.write("\t".join(fields[:4]) + "\n")
                counts["written"] += 1
    finally:
        for writer in writers.values():
            writer.close()

    with stats_out.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in counts.items():
            handle.write(f"{key}\t{value}\n")

    print(f"@ wrote stats: {stats_out}")
    for key, value in counts.items():
        print(f"@ {key}: {value}")


if __name__ == "__main__":
    main()
