#!/usr/bin/env python3

import argparse
import gzip
import subprocess
import tempfile
from pathlib import Path

import pyBigWig


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert sorted ATAC fragment BED/TSV files into unstranded cut-site bigWigs."
    )
    parser.add_argument("--fragments", required=True, help="Sorted fragments file.")
    parser.add_argument("--chromsizes", required=True, help="Two-column chrom sizes file.")
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output prefix; writes <prefix>_unstranded.bw.",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def load_chromsizes(path: Path):
    chroms = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[0]:
                continue
            chroms.append((parts[0], int(parts[1])))
    if not chroms:
        raise ValueError(f"No chromosome sizes found in {path}")
    return chroms


def flush_run(writer, chrom, start, end, value):
    if chrom is None or start is None or end is None or value <= 0:
        return
    writer.addEntries([chrom], [start], ends=[end], values=[float(value)])


def main():
    args = parse_args()

    fragments_path = Path(args.fragments).expanduser().resolve()
    chromsizes_path = Path(args.chromsizes).expanduser().resolve()
    output_prefix = Path(args.output_prefix).expanduser().resolve()
    output_path = Path(f"{output_prefix}_unstranded.bw")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    chromsizes = load_chromsizes(chromsizes_path)
    chromsize_map = dict(chromsizes)
    chrom_order = {chrom: idx for idx, (chrom, _) in enumerate(chromsizes)}

    with tempfile.TemporaryDirectory(prefix="chrombpnet_bw_") as tmpdir_name:
        tmpdir = Path(tmpdir_name)
        raw_cutsites = tmpdir / "cutsites.tsv"
        sorted_cutsites = tmpdir / "cutsites.sorted.tsv"

        with raw_cutsites.open("w", encoding="utf-8") as out_handle:
            with open_text(fragments_path) as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue

                    chrom = parts[0]
                    if chrom not in chromsize_map:
                        continue

                    start = int(parts[1])
                    end = int(parts[2])
                    chrom_end = chromsize_map[chrom]
                    if start < 0 or end <= start or start >= chrom_end:
                        continue

                    out_handle.write(f"{chrom_order[chrom]}\t{chrom}\t{start}\n")

                    right_cut = end - 1
                    if 0 <= right_cut < chrom_end:
                        out_handle.write(f"{chrom_order[chrom]}\t{chrom}\t{right_cut}\n")

        subprocess.run(
            ["sort", "-k1,1n", "-k3,3n", str(raw_cutsites), "-o", str(sorted_cutsites)],
            check=True,
        )

        current_chrom = None
        current_pos = None
        current_count = 0

        run_chrom = None
        run_start = None
        run_end = None
        run_value = 0

        with pyBigWig.open(str(output_path), "w") as bw:
            bw.addHeader(chromsizes)

            with sorted_cutsites.open("r", encoding="utf-8") as handle:
                for line in handle:
                    _, chrom, pos_text = line.rstrip("\n").split("\t")
                    pos = int(pos_text)

                    if current_chrom is None:
                        current_chrom = chrom
                        current_pos = pos
                        current_count = 1
                        continue

                    if chrom == current_chrom and pos == current_pos:
                        current_count += 1
                        continue

                    if run_chrom == current_chrom and run_end == current_pos and run_value == current_count:
                        run_end = current_pos + 1
                    else:
                        flush_run(bw, run_chrom, run_start, run_end, run_value)
                        run_chrom = current_chrom
                        run_start = current_pos
                        run_end = current_pos + 1
                        run_value = current_count

                    current_chrom = chrom
                    current_pos = pos
                    current_count = 1

            if current_chrom is not None:
                if run_chrom == current_chrom and run_end == current_pos and run_value == current_count:
                    run_end = current_pos + 1
                else:
                    flush_run(bw, run_chrom, run_start, run_end, run_value)
                    run_chrom = current_chrom
                    run_start = current_pos
                    run_end = current_pos + 1
                    run_value = current_count

            flush_run(bw, run_chrom, run_start, run_end, run_value)

    print(output_path)


if __name__ == "__main__":
    main()
