import csv
import hashlib
import os
import re
import subprocess
import sys
import tempfile
import threading
from io import BytesIO, StringIO
from glob import glob

import altair as alt
import logomaker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import torch
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
try:
    from scipy.stats import ttest_rel
except ImportError:
    ttest_rel = None
from pyfaidx import Fasta

# Resolve repo-local imports relative to this file so Streamlit works from Docker/systemd.
APP_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_ROOT = os.path.dirname(APP_DIR)
CHROMBPNET_DIR = os.path.join(CODE_ROOT, "03-chrombpnet")
sys.path.insert(0, CHROMBPNET_DIR)

from bpnetlite import BPNet
from bpnetlite.attribute import deep_lift_shap
from bpnetlite.bpnet import CountWrapper
from tangermeme_utils.wrappers import ChromBPNetWrapper


if os.environ.get("CHROMBPNET_EMBEDDED_APP") != "1":
    st.set_page_config(
        page_title="ChromBPNet Genome Browser",
        page_icon="DNA",
        layout="wide",
    )


# =========================
# DEVICE
# =========================
if torch.backends.mps.is_available():
    DEVICE = torch.device("mps")
elif torch.cuda.is_available():
    DEVICE = torch.device("cuda")
else:
    DEVICE = torch.device("cpu")

if DEVICE.type == "cuda":
    torch.backends.cudnn.benchmark = True


INPUTLEN = 2114
OUTPUTLEN = 1000
PROFILE_WINDOW = 1000
DEFAULT_WINDOW_HALF_WIDTH = 250
VALID_BASES = {"A", "C", "G", "T", "N"}
DEFAULT_LOGO_YLIM = 0.03
PROFILE_SHIFT = (INPUTLEN - OUTPUTLEN) // 2
DEFAULT_FIMO_ENABLED = False
DEFAULT_MOTIF_DATABASE = "JASPAR2022 vertebrates"
DEFAULT_FIMO_PVALUE = 1e-4
FIMO_BATCH_SIZE = 32

MOTIF_DATABASE_SPECS = {
    "JASPAR2022 vertebrates": {"filename": "JASPAR2022_vertebrates.pfm"},
    "JASPAR2020 vertebrates": {"filename": "JASPAR2020_vertebrates.pfm"},
    "HOCOMOCOv11 human": {"filename": "HOCOMOCOv11_HUMAN.pfm"},
    "HOCOMOCOv11 mouse": {"filename": "HOCOMOCOv11_MOUSE.pfm"},
}

MOTIF_TRACK_COLORS = [
    "#2563eb",
    "#f97316",
    "#16a34a",
    "#dc2626",
    "#7c3aed",
    "#0891b2",
]


# =========================
# PATHS
# =========================
DEFAULT_PROJ_IN = os.environ.get("CHROMBPNET_PROJECT_DIR", os.path.join(CHROMBPNET_DIR, "data", "local_run_chr", "work"))
DEFAULT_MODEL_RUN = "bias_1-col_aspn_ogna_thresh0.4"
DEFAULT_CELL_TYPE = "6-f13a1b_serpine1_gdf5_anxa2a_emilin3a_hi"
DEFAULT_FOLD_MODE = "1 fold"
GENOME_SPECS = {
    "GRCz11": {
        "label": "Zebrafish (GRCz11)",
        "default_coordinate": "chr7:8320256-8321214",
        "candidates": [
            os.environ.get("CHROMBPNET_GENOME_FASTA"),
            os.path.join(CHROMBPNET_DIR, "data", "local_run_chr", "refs", "Danio_rerio.GRCz11.chr.dna.toplevel.fa"),
            "/workspace/genome/ref/Danio_rerio.GRCz11.chr.dna.toplevel.fa",
            "/data/genome/ref/Danio_rerio.GRCz11.chr.dna.toplevel.fa",
        ],
        "glob_patterns": [
            "/workspace/genome/ref/Danio_rerio.GRCz11.chr.dna.toplevel.fa*",
            "/data/genome/ref/Danio_rerio.GRCz11.chr.dna.toplevel.fa*",
        ],
    },
    "mm10": {
        "label": "Mouse (mm10)",
        "default_coordinate": "chr13:55766039-55767038",
        "candidates": [
            os.path.join(CHROMBPNET_DIR, "data", "local_run_chr", "refs", "mm10", "mm10.fa"),
            "/data/projects/active/HDMA/code/03-chrombpnet/data/local_run_chr/refs/mm10/mm10.fa",
        ],
        "glob_patterns": [],
    },
    "mm39": {
        "label": "Mouse (mm39)",
        "default_coordinate": "chr13:55913852-55914851",
        "candidates": [
            os.path.join(CHROMBPNET_DIR, "data", "local_run_chr", "refs", "mm39", "mm39.fa"),
            "/data/projects/active/HDMA/code/03-chrombpnet/data/local_run_chr/refs/mm39/mm39.fa",
        ],
        "glob_patterns": [],
    },
}
DEFAULT_GENOME_KEY = "GRCz11"


def resolve_genome_fasta(genome_key):
    spec = GENOME_SPECS[genome_key]
    candidates = [candidate for candidate in spec["candidates"] if candidate]
    for pattern in spec["glob_patterns"]:
        candidates.extend(sorted(glob(pattern)))

    for candidate in candidates:
        if candidate and os.path.isfile(candidate):
            return candidate
    return candidates[0] if candidates else ""


def apply_genome_selection(genome_key):
    st.session_state.genome_key = genome_key
    st.session_state.genome_fasta = resolve_genome_fasta(genome_key)
    st.session_state.coordinate_region = GENOME_SPECS[genome_key]["default_coordinate"]
    st.session_state.coordinate_input_text = st.session_state.coordinate_region
    st.session_state.sequence_label = ""
    st.session_state.active_original_seq = ""
    st.session_state.active_seq = ""
    st.session_state.results = None
    st.session_state.zoom_center_genomic = None
    st.session_state.browser_seq = ""
    st.session_state.browser_start_1 = 1
    st.session_state.browser_end_1 = 1
    st.session_state.browser_active_start_1 = 1
    st.session_state.accumulated_matches = []
    st.session_state.pending_search_selection = []
    st.session_state.browser_highlight_intervals = []


# =========================
# LOAD MODELS / DATA
# =========================

@st.cache_resource
def load_models(model_paths):
    pred_models = []
    count_models = []
    for model_path in model_paths:
        pred_models.append(ChromBPNetWrapper(BPNet.from_chrombpnet(filename=model_path)).to(DEVICE).eval())
        count_models.append(CountWrapper(BPNet.from_chrombpnet(filename=model_path)).to(DEVICE).eval())
    return pred_models, count_models


@st.cache_resource
def load_genome(genome_path):
    return Fasta(genome_path)


@st.cache_data
def discover_models(project_path):
    model_root = os.path.join(project_path, "01-models")
    discovered = {}
    if not os.path.isdir(model_root):
        return discovered

    for run_name in sorted(os.listdir(model_root)):
        run_dir = os.path.join(model_root, run_name)
        if not os.path.isdir(run_dir):
            continue

        for cell_type in sorted(os.listdir(run_dir)):
            cell_dir = os.path.join(run_dir, cell_type)
            if not os.path.isdir(cell_dir):
                continue

            for fold in sorted(os.listdir(cell_dir)):
                model_path = os.path.join(cell_dir, fold, "models", "chrombpnet_nobias.h5")
                if os.path.isfile(model_path):
                    discovered.setdefault(run_name, {}).setdefault(cell_type, {})[fold] = model_path

    return discovered



@st.cache_data
def discover_project_paths():
    candidates = []
    roots = [
        os.path.join(CHROMBPNET_DIR, "data"),
        "/workspace/HDMA/03-chrombpnet/data",
        "/data/projects/active/HDMA/code/03-chrombpnet/data",
    ]

    seen_roots = set()
    for root in roots:
        real_root = os.path.realpath(root)
        if real_root in seen_roots or not os.path.isdir(real_root):
            continue
        seen_roots.add(real_root)

        for entry in sorted(os.listdir(real_root)):
            work_dir = os.path.join(real_root, entry, "work")
            if os.path.isdir(work_dir):
                display_label = "ucmab_14dpf" if entry == "local_run_chr" else entry
                candidates.append((display_label, os.path.realpath(work_dir)))

    deduped = []
    seen_labels = set()
    for label, work_dir in candidates:
        if label not in seen_labels:
            deduped.append((label, work_dir))
            seen_labels.add(label)
    return deduped


def get_project_option_labels(project_options):
    return [label for label, _ in project_options]


def resolve_project_path(project_options, selected_label, fallback_path):
    for label, work_dir in project_options:
        if label == selected_label:
            return work_dir
    return fallback_path


def get_fold_model_paths(model_catalog, run_name, cell_type, fold_mode):
    fold_map = model_catalog.get(run_name, {}).get(cell_type, {})
    folds = sorted(fold_map)
    if fold_mode == "1 fold":
        folds = folds[:1]
    return [fold_map[fold] for fold in folds], folds


@st.cache_data
def discover_motif_database_roots():
    patterns = [
        os.path.join(sys.prefix, "lib", "python*", "site-packages", "data", "motif_databases"),
        os.path.join(sys.prefix, "lib", "python*", "dist-packages", "data", "motif_databases"),
        "/opt/conda/lib/python*/site-packages/data/motif_databases",
        "/usr/local/lib/python*/site-packages/data/motif_databases",
        "/usr/local/lib/python*/dist-packages/data/motif_databases",
        "/usr/lib/python*/dist-packages/data/motif_databases",
    ]
    candidates = []
    for pattern in patterns:
        candidates.extend(glob(pattern))

    roots = []
    seen = set()
    for candidate in candidates:
        real_candidate = os.path.realpath(candidate)
        if os.path.isdir(real_candidate) and real_candidate not in seen:
            roots.append(real_candidate)
            seen.add(real_candidate)
    return roots


def resolve_motif_database_path(database_label):
    spec = MOTIF_DATABASE_SPECS.get(database_label)
    if spec is None:
        raise ValueError(f"Unknown motif database: {database_label}")

    for root in discover_motif_database_roots():
        candidate = os.path.join(root, spec["filename"])
        if os.path.isfile(candidate):
            return candidate

    raise FileNotFoundError(
        f"Could not find the motif database for '{database_label}'. Checked: {', '.join(discover_motif_database_roots()) or 'no known roots'}"
    )


def motif_label_from_raw_name(raw_name):
    if ".H11MO." in raw_name:
        return raw_name.split(".H11MO.", 1)[0]

    motif_label = raw_name
    if "_" in raw_name:
        motif_label = raw_name.split("_", 1)[1]
    if "." in motif_label:
        motif_label = motif_label.split(".", 1)[1]
    return motif_label


def sanitize_motif_identifier(text, fallback="motif"):
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", str(text or "")).strip("_")
    return cleaned or fallback


@st.cache_data
def discover_available_motif_databases():
    available = []
    for label in MOTIF_DATABASE_SPECS:
        try:
            resolve_motif_database_path(label)
        except FileNotFoundError:
            continue
        available.append(label)
    return available


@st.cache_data
def prepare_fimo_database(database_label):
    source_path = resolve_motif_database_path(database_label)
    if source_path.endswith(".meme"):
        return source_path

    cache_root = os.path.join(tempfile.gettempdir(), "chrombpnet_fimo_databases")
    os.makedirs(cache_root, exist_ok=True)

    fingerprint = hashlib.sha1(
        f"pfm2meme-v4:{source_path}:{os.path.getmtime(source_path)}:{os.path.getsize(source_path)}".encode("utf-8")
    ).hexdigest()[:16]
    stem = os.path.splitext(os.path.basename(source_path))[0]
    converted_path = os.path.join(cache_root, f"{stem}.{fingerprint}.meme")
    if os.path.isfile(converted_path):
        return converted_path

    motifs = []
    current_name = None
    current_rows = []
    with open(source_path, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if current_name and current_rows:
                    motifs.append((current_name, current_rows))
                current_name = line[1:].strip()
                current_rows = []
                continue

            values = [float(value) for value in line.replace("[", " ").replace("]", " ").split()]
            if len(values) != 4:
                raise ValueError(f"Unsupported motif row in {source_path}: {line}")
            current_rows.append(values)

    if current_name and current_rows:
        motifs.append((current_name, current_rows))
    if not motifs:
        raise ValueError(f"No motifs were parsed from {source_path}")

    with open(converted_path, "w", encoding="utf-8") as handle:
        handle.write("MEME version 4\n\n")
        handle.write("ALPHABET= ACGT\n\n")
        handle.write("strands: + -\n\n")
        handle.write("Background letter frequencies\n")
        handle.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        for motif_idx, (raw_name, rows) in enumerate(motifs, start=1):
            motif_label = sanitize_motif_identifier(motif_label_from_raw_name(raw_name))
            motif_id = sanitize_motif_identifier(raw_name, fallback=f"motif_{motif_idx}")
            nsites = max(1, int(round(max(sum(row) for row in rows))))
            handle.write(f"MOTIF {motif_id} {motif_label}\n")
            handle.write(f"letter-probability matrix: alength= 4 w= {len(rows)} nsites= {nsites} E= 0\n")
            for row in rows:
                total = float(sum(row))
                if total <= 0:
                    probs = [0.25, 0.25, 0.25, 0.25]
                else:
                    probs = [value / total for value in row]
                handle.write(" ".join(f"{prob:.8f}" for prob in probs) + "\n")
            handle.write("\n")
    return converted_path


@st.cache_data
def prepare_fimo_database_batches(database_label, batch_size=FIMO_BATCH_SIZE):
    meme_path = prepare_fimo_database(database_label)
    with open(meme_path, encoding="utf-8") as handle:
        lines = handle.readlines()

    motif_start_indexes = [idx for idx, line in enumerate(lines) if line.startswith("MOTIF ")]
    if not motif_start_indexes:
        return [meme_path]

    header = lines[: motif_start_indexes[0]]
    motif_start_indexes.append(len(lines))
    cache_root = os.path.join(tempfile.gettempdir(), "chrombpnet_fimo_batches")
    os.makedirs(cache_root, exist_ok=True)

    fingerprint = hashlib.sha1(
        f"{meme_path}:{os.path.getmtime(meme_path)}:{os.path.getsize(meme_path)}:{batch_size}".encode("utf-8")
    ).hexdigest()[:16]
    stem = os.path.splitext(os.path.basename(meme_path))[0]

    batch_paths = []
    for batch_number, start_idx in enumerate(range(0, len(motif_start_indexes) - 1, int(batch_size)), start=1):
        batch_path = os.path.join(cache_root, f"{stem}.{fingerprint}.batch{batch_number:03d}.meme")
        if not os.path.isfile(batch_path):
            chunk_lines = list(header)
            chunk_boundaries = motif_start_indexes[start_idx : start_idx + int(batch_size) + 1]
            for left, right in zip(chunk_boundaries[:-1], chunk_boundaries[1:]):
                chunk_lines.extend(lines[left:right])
            with open(batch_path, "w", encoding="utf-8") as handle:
                handle.writelines(chunk_lines)
        batch_paths.append(batch_path)

    return batch_paths


def summarize_motif_label(hit):
    label = hit.get("motif_alt_id") or hit.get("motif_id") or hit.get("pattern_name") or "motif"
    if label == ".":
        label = hit.get("motif_id") or hit.get("pattern_name") or "motif"
    label = str(label).replace("_", " ")
    return label if len(label) <= 24 else f"{label[:21]}..."


def parse_fimo_tsv(fimo_stdout):
    raw_lines = [line.strip() for line in fimo_stdout.splitlines() if line.strip()]
    header_line = None
    data_lines = []
    for line in raw_lines:
        if line.startswith("#pattern name"):
            header_line = line[1:]
            continue
        if line.startswith("#"):
            continue
        data_lines.append(line)

    if header_line is None:
        header_line = "pattern name\tsequence name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched sequence"

    rows = [header_line, *data_lines]
    if len(rows) == 1:
        return []

    reader = csv.DictReader(StringIO("\n".join(rows)), delimiter="\t")
    hits = []
    for row in reader:
        start_text = row.get("start")
        stop_text = row.get("stop") or row.get("end")
        score_text = row.get("score")
        pvalue_text = row.get("p-value") or row.get("pvalue")
        if not start_text or not stop_text or not score_text or not pvalue_text:
            continue

        start_idx = max(0, int(start_text) - 1)
        end_idx = max(start_idx + 1, int(stop_text))
        pvalue = float(pvalue_text)
        qvalue_text = str(row.get("q-value") or row.get("qvalue") or "").strip()
        qvalue = float(qvalue_text) if qvalue_text and qvalue_text.lower() != "nan" else None
        matched_sequence = row.get("matched_sequence") or row.get("matched sequence") or ""
        hit = {
            "pattern_name": row.get("pattern name", ""),
            "motif_id": row.get("motif_id") or row.get("pattern name", ""),
            "motif_alt_id": row.get("motif_alt_id", ""),
            "label": row.get("motif_alt_id") or row.get("motif_id") or row.get("pattern name") or "motif",
            "start_idx": start_idx,
            "end_idx": end_idx,
            "strand": row.get("strand", "."),
            "score": float(score_text),
            "pvalue": pvalue,
            "qvalue": qvalue,
            "matched_sequence": matched_sequence,
        }
        hit["display_label"] = summarize_motif_label(hit)
        hits.append(hit)

    hits.sort(key=lambda item: (item["pvalue"], -(item["end_idx"] - item["start_idx"]), item["start_idx"]))
    return hits


@st.cache_data(show_spinner=False)
def scan_sequence_with_fimo(seq, database_label, pvalue_threshold, max_hits=None):
    clean_seq = sanitize_sequence(seq)
    invalid = validate_sequence(clean_seq)
    if invalid:
        raise ValueError(f"FIMO scan received invalid bases: {', '.join(invalid)}")
    if not clean_seq:
        return {
            "database_label": database_label,
            "hits": [],
            "shown_hits": 0,
            "total_hits": 0,
            "pvalue_threshold": float(pvalue_threshold),
            "sequence_length": 0,
        }

    motif_db_paths = prepare_fimo_database_batches(database_label)
    with tempfile.TemporaryDirectory(prefix="chrombpnet_fimo_") as tmpdir:
        fasta_path = os.path.join(tmpdir, "query.fa")
        with open(fasta_path, "w", encoding="utf-8") as handle:
            handle.write(">query\n")
            handle.write(clean_seq)
            handle.write("\n")

        all_hits = []
        for motif_db_path in motif_db_paths:
            result = subprocess.run(
                ["fimo", "--text", "--thresh", f"{float(pvalue_threshold):.6g}", motif_db_path, fasta_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            if result.returncode != 0:
                stderr = (result.stderr or "").strip()
                raise RuntimeError(
                    f"FIMO failed on motif batch {os.path.basename(motif_db_path)}: {stderr or 'unknown error'}"
                )
            all_hits.extend(parse_fimo_tsv(result.stdout))

    all_hits.sort(key=lambda item: (item["pvalue"], -(item["end_idx"] - item["start_idx"]), item["start_idx"]))
    if max_hits is not None:
        all_hits = all_hits[: max(0, int(max_hits))]
    return {
        "database_label": database_label,
        "hits": all_hits,
        "shown_hits": len(all_hits),
        "total_hits": len(all_hits),
        "pvalue_threshold": float(pvalue_threshold),
        "sequence_length": len(clean_seq),
    }


# =========================
# COORDINATE / SEQUENCE
# =========================
def parse_coordinate_region(region):
    chrom, coords = region.split(":")
    start_1, end_1 = map(int, coords.split("-"))
    return chrom, start_1, end_1


def extract_sequence(region, genome_path):
    chrom, start_1, end_1 = parse_coordinate_region(region)
    genome = load_genome(genome_path)
    return str(genome[chrom][start_1 - 1 : end_1]).upper()


def sanitize_sequence(seq):
    return "".join(seq.upper().replace("U", "T").split())


def read_single_fasta(path):
    with open(path) as handle:
        return parse_fasta_text(handle.read())


def parse_fasta_text(text):
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("No sequence found.")

    header = "custom"
    seq_lines = []
    for line in lines:
        if line.startswith(">"):
            if seq_lines:
                raise ValueError("Please provide exactly one FASTA sequence.")
            header = line[1:].strip() or "custom"
        else:
            seq_lines.append(line)

    seq = sanitize_sequence("".join(seq_lines if seq_lines else lines))
    invalid = validate_sequence(seq)
    if invalid:
        raise ValueError(f"Sequence contains unsupported characters: {invalid}")
    if not seq:
        raise ValueError("No sequence found.")
    return header, seq


def validate_sequence(seq):
    return sorted(set(seq) - VALID_BASES)


def dinucleotide_style_padding(seq, length, rng):
    seq = "".join(base for base in seq if base in "ACGT")
    if length <= 0:
        return ""
    if len(seq) == 0:
        return "A" * length
    if len(seq) == 1:
        return seq * length

    transitions = {base: [] for base in "ACGT"}
    for left, right in zip(seq[:-1], seq[1:]):
        transitions[left].append(right)

    starts = list(seq[:-1])
    current = rng.choice(starts)
    out = []
    while len(out) < length:
        out.append(current)
        choices = transitions.get(current, [])
        current = rng.choice(choices) if choices else rng.choice(starts)

    return "".join(out[:length])


def center_resize(seq, padding_seed=0):
    if len(seq) > INPUTLEN:
        start = (len(seq) - INPUTLEN) // 2
        seq = seq[start : start + INPUTLEN]
    elif len(seq) < INPUTLEN:
        pad = INPUTLEN - len(seq)
        rng = np.random.default_rng(padding_seed)
        left_pad = dinucleotide_style_padding(seq, pad // 2, rng)
        right_pad = dinucleotide_style_padding(seq, pad - pad // 2, rng)
        seq = left_pad + seq + right_pad
    return seq


def one_hot_encode_dna(seq, randomize_n=True, seed=0):
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}

    rng = np.random.default_rng(seed)

    seq = seq.upper()
    arr = np.zeros((4, len(seq)), dtype=np.float32)

    for i, base in enumerate(seq):
        if base == "N":
            if randomize_n:
                base = rng.choice(["A", "C", "G", "T"])
            else:
                base = "A"

        if base in mapping:
            arr[mapping[base], i] = 1.0

    return torch.tensor(arr, dtype=torch.float32)


def normalize_attr_matrix(attr):
    attr = np.asarray(attr)
    if attr.ndim != 2:
        raise ValueError(f"Expected 2D attribution array, got shape {attr.shape}")
    if attr.shape[0] == 4 and attr.shape[1] != 4:
        return attr.T
    if attr.shape[1] == 4:
        return attr
    raise ValueError(f"Unexpected attribution shape {attr.shape}; expected 4 x L or L x 4.")


def extract_attr_window(attr, seq_len):
    attr = normalize_attr_matrix(attr)
    if seq_len >= INPUTLEN:
        return attr
    pad = INPUTLEN - seq_len
    left = pad // 2
    return attr[left : left + seq_len]


def high_z_attr_windows(attr, window_size=10, z_threshold=1.0):
    attr = normalize_attr_matrix(attr)
    if len(attr) == 0:
        return [], np.array([])

    window_size = max(1, min(int(window_size), len(attr)))
    per_base_score = np.abs(attr).sum(axis=1)
    window_scores = np.convolve(per_base_score, np.ones(window_size) / window_size, mode="valid")
    score_std = float(window_scores.std())
    if score_std == 0:
        zscores = np.zeros_like(window_scores)
    else:
        zscores = (window_scores - float(window_scores.mean())) / score_std

    high_windows = [(idx, idx + window_size) for idx, zscore in enumerate(zscores) if zscore > float(z_threshold)]
    return merge_intervals(high_windows), zscores


def intervals_overlap(left, right):
    return left[0] < right[1] and right[0] < left[1]


def filter_motif_scan_by_attr(motif_scan, attr, window_size=10, z_threshold=1.0):
    if motif_scan is None:
        return None

    high_windows, zscores = high_z_attr_windows(attr, window_size=window_size, z_threshold=z_threshold)
    filtered_hits = [
        hit
        for hit in motif_scan.get("hits", [])
        if any(intervals_overlap((int(hit["start_idx"]), int(hit["end_idx"])), window) for window in high_windows)
    ]

    filtered_scan = dict(motif_scan)
    filtered_scan["hits"] = filtered_hits
    filtered_scan["shown_hits"] = len(filtered_hits)
    filtered_scan["unfiltered_total_hits"] = int(motif_scan.get("total_hits", len(motif_scan.get("hits", []))))
    filtered_scan["total_hits"] = len(filtered_hits)
    filtered_scan["high_z_intervals"] = high_windows
    filtered_scan["high_z_window_count"] = len(high_windows)
    filtered_scan["z_threshold"] = float(z_threshold)
    filtered_scan["z_window_size"] = int(window_size)
    filtered_scan["max_zscore"] = float(zscores.max()) if len(zscores) else None
    return filtered_scan


def clip_int(value, lower, upper):
    return max(lower, min(int(value), upper))


def genomic_window(region_start_1, region_end_1, center_1, half_width):
    left = max(region_start_1, center_1 - half_width)
    right = min(region_end_1, center_1 + half_width)
    if right <= left:
        right = min(region_end_1, left + 1)
    return left, right


def genome_to_index(genomic_coord_1, region_start_1):
    return genomic_coord_1 - region_start_1


def index_to_genome(index_0, region_start_1):
    return region_start_1 + index_0


def make_position_ruler(start_1, end_1):
    return "".join(str(pos % 10) for pos in range(start_1, end_1 + 1))


def find_sequence_matches(seq, query):
    clean_query = sanitize_sequence(query)
    if not clean_query:
        return []

    matches = []
    query_len = len(clean_query)
    for start in range(0, len(seq) - query_len + 1):
        window = seq[start : start + query_len]
        if all(q == "N" or s == q for s, q in zip(window, clean_query)):
            matches.append((start, start + query_len))
    return matches


def merge_intervals(intervals):
    if not intervals:
        return []

    merged = []
    for start, end in sorted(intervals):
        if not merged or start >= merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [tuple(interval) for interval in merged]


def selected_intervals_from_rows(rows):
    return merge_intervals(
        [(row["start_idx"], row["end_idx"]) for row in rows if row.get("selected")]
    )


def extract_centered_genome_sequence(region, genome_path, target_len=INPUTLEN):
    chrom, start_1, end_1 = parse_coordinate_region(region)
    genome = load_genome(genome_path)
    chrom_len = len(genome[chrom])

    center_0 = ((start_1 - 1) + end_1) // 2
    seq_start_0 = center_0 - (target_len // 2)
    seq_end_0 = seq_start_0 + target_len

    if seq_start_0 < 0:
        seq_end_0 -= seq_start_0
        seq_start_0 = 0
    if seq_end_0 > chrom_len:
        shift = seq_end_0 - chrom_len
        seq_start_0 = max(0, seq_start_0 - shift)
        seq_end_0 = chrom_len

    seq = str(genome[chrom][seq_start_0:seq_end_0]).upper().replace("U", "T")
    return chrom, seq_start_0 + 1, seq_end_0, seq


def extend_deleted_sequence_from_genome(seq, chrom, start_1, end_1, genome_path, target_len=INPUTLEN):
    seq = sanitize_sequence(seq)
    if len(seq) >= target_len:
        return seq[:target_len], start_1, start_1 + target_len - 1

    genome = load_genome(genome_path)
    chrom_seq = genome[chrom]
    chrom_len = len(chrom_seq)

    missing = target_len - len(seq)
    take_left = missing // 2
    take_right = missing - take_left

    left_start_0 = max(0, start_1 - 1 - take_left)
    left_end_0 = start_1 - 1
    right_start_0 = end_1
    right_end_0 = min(chrom_len, end_1 + take_right)

    left_flank = str(chrom_seq[left_start_0:left_end_0]).upper().replace("U", "T")
    right_flank = str(chrom_seq[right_start_0:right_end_0]).upper().replace("U", "T")

    still_missing = target_len - (len(left_flank) + len(seq) + len(right_flank))
    if still_missing > 0:
        extra_left = min(still_missing, left_start_0)
        left_start_0 -= extra_left
        left_flank = str(chrom_seq[left_start_0:left_end_0]).upper().replace("U", "T")
        still_missing = target_len - (len(left_flank) + len(seq) + len(right_flank))
    if still_missing > 0:
        extra_right = min(still_missing, chrom_len - right_end_0)
        right_end_0 += extra_right
        right_flank = str(chrom_seq[right_start_0:right_end_0]).upper().replace("U", "T")

    extended = sanitize_sequence(left_flank + seq + right_flank)
    if len(extended) > target_len:
        extended = extended[:target_len]
        right_end_0 = left_start_0 + len(extended)

    return extended, left_start_0 + 1, left_start_0 + len(extended)


# =========================
# MODEL
# =========================

ATTRIBUTION_LOCK = threading.Lock()


def reset_attr_model_state(model):
    for module in model.modules():
        handles = getattr(module, "handles", None)
        if handles:
            for handle in handles:
                try:
                    handle.remove()
                except Exception:
                    pass
            module.handles = []

        for hook_name in ("_forward_hooks", "_forward_pre_hooks", "_backward_hooks"):
            hook_dict = getattr(module, hook_name, None)
            if hook_dict:
                hook_dict.clear()

        if hasattr(module, "_NON_LINEAR_OPS"):
            try:
                delattr(module, "_NON_LINEAR_OPS")
            except Exception:
                pass


def first_tensor_output(output):
    if isinstance(output, (tuple, list)):
        if not output:
            raise ValueError("Model returned an empty output tuple/list.")
        return output[0]
    return output


def run_model_batch(seqs, pred_models, count_models, batch_size=None):
    resized = [center_resize(seq) for seq in seqs]
    batch_size = max(1, int(batch_size or len(resized) or 1))

    profiles_by_fold = []
    counts_by_fold = []
    with torch.inference_mode():
        for pred_model, count_model in zip(pred_models, count_models):
            fold_profiles = []
            fold_counts = []
            for start in range(0, len(resized), batch_size):
                chunk = resized[start : start + batch_size]
                x = torch.stack([one_hot_encode_dna(seq) for seq in chunk], dim=0).to(DEVICE, non_blocking=True)
                y = first_tensor_output(pred_model(x))
                logcounts = first_tensor_output(count_model(x))
                fold_profiles.append(y.detach().cpu().numpy().reshape(len(chunk), -1))
                fold_counts.append(np.exp(logcounts.detach().cpu().numpy().reshape(len(chunk), -1)[:, 0]))
            profiles_by_fold.append(np.concatenate(fold_profiles, axis=0))
            counts_by_fold.append(np.concatenate(fold_counts, axis=0))

    profiles = np.stack(profiles_by_fold, axis=0)
    counts = np.stack(counts_by_fold, axis=0).astype(np.float64)
    return profiles, counts


def summarize_model_batch(profiles, counts, seq_index):
    seq_profiles = profiles[:, seq_index, :]
    seq_counts = counts[:, seq_index]
    return (
        seq_profiles.mean(axis=0),
        float(seq_counts.mean()),
        seq_profiles.std(axis=0),
        float(seq_counts.std(ddof=0)),
        seq_profiles,
        seq_counts,
    )


def run_model(seq, pred_models, count_models):
    profiles, counts = run_model_batch([seq], pred_models, count_models, batch_size=1)
    return summarize_model_batch(profiles, counts, 0)


def get_attr(seq, model_paths):
    seq = center_resize(seq)
    attr_device = DEVICE if DEVICE.type == "cuda" else torch.device("cpu")
    x = one_hot_encode_dna(seq).unsqueeze(0).to(attr_device)

    attrs = []
    with ATTRIBUTION_LOCK:
        for model_path in model_paths:
            attr_model = CountWrapper(BPNet.from_chrombpnet(filename=model_path)).to(attr_device).eval()
            reset_attr_model_state(attr_model)
            attr = deep_lift_shap(
                attr_model,
                x,
                random_state=0,
                device=str(attr_device),
            )
            attrs.append(normalize_attr_matrix(attr[0].detach().cpu().numpy()))
            del attr_model
            if attr_device.type == "cuda":
                torch.cuda.empty_cache()

    attrs = np.stack(attrs, axis=0)
    return attrs.mean(axis=0), attrs.std(axis=0)


def run_all_outputs(before_seq, after_seq, pred_models, count_models, model_paths):
    prediction_batch_size = max(1, len(model_paths))
    profile_batch, count_batch = run_model_batch(
        [before_seq, after_seq],
        pred_models,
        count_models,
        batch_size=prediction_batch_size,
    )
    prof_before, cnt_before, prof_before_std, cnt_before_std, prof_before_folds, cnt_before_folds = summarize_model_batch(profile_batch, count_batch, 0)
    prof_after, cnt_after, prof_after_std, cnt_after_std, prof_after_folds, cnt_after_folds = summarize_model_batch(profile_batch, count_batch, 1)
    attr_before, attr_before_std = get_attr(before_seq, model_paths)
    attr_after, attr_after_std = get_attr(after_seq, model_paths)
    attr_before = extract_attr_window(attr_before, len(before_seq))
    attr_after = extract_attr_window(attr_after, len(after_seq))
    attr_before_std = extract_attr_window(attr_before_std, len(before_seq))
    attr_after_std = extract_attr_window(attr_after_std, len(after_seq))
    return {
        "prof_before": prof_before,
        "prof_after": prof_after,
        "prof_before_std": prof_before_std,
        "prof_after_std": prof_after_std,
        "prof_before_folds": prof_before_folds,
        "prof_after_folds": prof_after_folds,
        "cnt_before": cnt_before,
        "cnt_after": cnt_after,
        "cnt_before_std": cnt_before_std,
        "cnt_after_std": cnt_after_std,
        "cnt_before_folds": cnt_before_folds,
        "cnt_after_folds": cnt_after_folds,
        "attr_before": attr_before,
        "attr_after": attr_after,
        "attr_before_std": attr_before_std,
        "attr_after_std": attr_after_std,
    }


# =========================
# EDIT
# =========================
def replace_region(seq, start_0, end_0_exclusive, new_seq):
    return seq[:start_0] + new_seq + seq[end_0_exclusive:]


def edit_intervals(seq, intervals, operation, replacement):
    updated = seq
    for start, end in sorted(intervals, reverse=True):
        new_seq = "" if operation == "Delete" else replacement
        updated = replace_region(updated, start, end, new_seq)
    return updated


def normalize_editor_text(value):
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return str(value)


def normalize_match_rows(rows):
    normalized_rows = []
    for row in rows:
        normalized = dict(row)
        normalized["selected"] = bool(normalized.get("selected", False))
        normalized["row_action"] = normalized.get("row_action") or "Keep"
        normalized["row_replacement"] = normalize_editor_text(normalized.get("row_replacement", ""))
        normalized_rows.append(normalized)
    return normalized_rows


def merge_editor_state(editor_df, editor_key):
    editor_state = st.session_state.get(editor_key, {})
    edited_rows = editor_state.get("edited_rows", {}) if isinstance(editor_state, dict) else {}
    if not edited_rows:
        return editor_df

    merged_df = editor_df.copy()
    for row_idx, updates in edited_rows.items():
        try:
            idx = int(row_idx)
        except (TypeError, ValueError):
            continue
        if idx < 0 or idx >= len(merged_df):
            continue
        for column_name, value in updates.items():
            if column_name not in merged_df.columns:
                continue
            if column_name == "row_replacement":
                value = normalize_editor_text(value)
            merged_df.iat[idx, merged_df.columns.get_loc(column_name)] = value
    return merged_df


def edit_rows(seq, rows):
    updated = seq
    for row in sorted(rows, key=lambda item: item["start_idx"], reverse=True):
        operation = row["row_action"]
        replacement = sanitize_sequence(normalize_editor_text(row.get("row_replacement", ""))) if operation == "Replace" else ""
        new_seq = "" if operation == "Delete" else replacement
        updated = replace_region(updated, row["start_idx"], row["end_idx"], new_seq)
    return updated


# =========================
# PLOTTING
# =========================
def plot_profile_overlay(prof_before, prof_after):
    start = max(0, (len(prof_before) - PROFILE_WINDOW) // 2)
    end = start + PROFILE_WINDOW
    x = np.arange(start, end)

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(x, prof_before[start:end], label="Before", linewidth=2)
    ax.plot(x, prof_after[start:end], label="After", linewidth=2)
    ax.set_title(f"Predicted profile, center {PROFILE_WINDOW} bp")
    ax.set_xlabel("Position in output profile")
    ax.set_ylabel("Signal")
    ax.legend()
    return fig


def compute_paired_count_test(before_counts, after_counts):
    before = np.asarray(before_counts, dtype=np.float64)
    after = np.asarray(after_counts, dtype=np.float64)
    n_folds = int(min(len(before), len(after)))

    if n_folds < 2:
        return {"pvalue": None, "statistic": None, "n_folds": n_folds, "message": "Need at least 2 folds for a paired t-test."}

    if ttest_rel is None:
        return {"pvalue": None, "statistic": None, "n_folds": n_folds, "message": "scipy is unavailable for the paired t-test."}

    result = ttest_rel(after[:n_folds], before[:n_folds], nan_policy="omit")
    return {
        "pvalue": float(result.pvalue) if result.pvalue is not None else None,
        "statistic": float(result.statistic) if result.statistic is not None else None,
        "n_folds": n_folds,
        "message": None,
    }


def plot_count_summary(before_counts, after_counts):
    before = np.asarray(before_counts, dtype=np.float64)
    after = np.asarray(after_counts, dtype=np.float64)
    n_folds = int(min(len(before), len(after)))

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    x = np.array([0, 1], dtype=np.float64)
    means = np.array([before.mean(), after.mean()], dtype=np.float64)
    stds = np.array([before.std(ddof=0), after.std(ddof=0)], dtype=np.float64)

    ax.bar(x, means, color=["#2563eb", "#f97316"], alpha=0.85, width=0.6)
    ax.errorbar(x, means, yerr=stds, fmt="none", ecolor="#334155", capsize=5, linewidth=1.4)

    for idx in range(n_folds):
        ax.plot(x, [before[idx], after[idx]], color="#94a3b8", alpha=0.55, linewidth=1.1, marker="o", markersize=3.5)

    ax.set_xticks(x)
    ax.set_xticklabels(["Before", "After"])
    ax.set_ylabel("Predicted count")
    ax.set_title("Count change across folds")
    return fig


def format_count_datapoints(before_counts, after_counts):
    before = np.asarray(before_counts, dtype=np.float64)
    after = np.asarray(after_counts, dtype=np.float64)
    lines = []
    for idx, (before_value, after_value) in enumerate(zip(before, after), start=1):
        lines.append(
            f"Fold {idx}: before={before_value:.4f}, after={after_value:.4f}, delta={after_value - before_value:.4f}"
        )
    return lines



def plot_profile_overlay_marked(prof_before, prof_after, highlight_windows=None, selected_windows=None, deleted_windows=None, prof_before_std=None, prof_after_std=None, prof_before_folds=None, prof_after_folds=None, window_offset=PROFILE_SHIFT):
    start = max(0, (len(prof_before) - PROFILE_WINDOW) // 2)
    end = start + PROFILE_WINDOW
    x = np.arange(start, end)

    fig, ax = plt.subplots(figsize=(12, 4))
    before_window = prof_before[start:end]
    after_window = prof_after[start:end]
    unedited_profile = np.allclose(before_window, after_window)

    before_mean_color = "#2563eb"
    before_fold_color = "#93c5fd"
    after_mean_color = "#f97316"
    after_fold_color = "#fdba74"

    if prof_before_folds is not None:
        for idx, fold_profile in enumerate(prof_before_folds):
            label = "Before folds" if idx == 0 else None
            ax.plot(x, fold_profile[start:end], linewidth=1.0, color=before_fold_color, alpha=0.55, label=label)

    if prof_after_folds is not None:
        same_as_before = prof_before_folds is not None and np.allclose(prof_before_folds[:, start:end], prof_after_folds[:, start:end])
        if not same_as_before:
            for idx, fold_profile in enumerate(prof_after_folds):
                label = "After folds" if idx == 0 else None
                ax.plot(x, fold_profile[start:end], linewidth=1.0, color=after_fold_color, alpha=0.55, label=label)

    ax.plot(x, before_window, label="Before mean", linewidth=2.2, color=before_mean_color)
    if unedited_profile:
        ax.plot(x, after_window, label="After mean", linewidth=2.2, color=before_mean_color, alpha=0.9)
    else:
        ax.plot(x, after_window, label="After mean", linewidth=2.2, color=after_mean_color)

    for window_start, window_end in selected_windows or []:
        out_start = max(start, window_start - window_offset)
        out_end = min(end, window_end - window_offset)
        if out_end > out_start:
            ax.axvspan(out_start, out_end, color="#d1d5db", alpha=0.5)

    for window_start, window_end in highlight_windows or []:
        out_start = max(start, window_start - window_offset)
        out_end = min(end, window_end - window_offset)
        if out_end > out_start:
            ax.axvspan(out_start, out_end, color="#dbeafe", alpha=0.45)

    ax.set_title(f"Predicted profile across folds, center {PROFILE_WINDOW} bp")
    ax.set_xlabel("Position in output profile")
    ax.set_ylabel("Signal")
    ax.legend()
    return fig


def add_window_labels(ax, selected_windows=None, deleted_windows=None):
    label_lines = []

    if selected_windows:
        label_lines.extend([f"Selected: {item}" for item in selected_windows[:4]])
        if len(selected_windows) > 4:
            label_lines.append(f"+ {len(selected_windows) - 4} more selected")

    if deleted_windows:
        label_lines.extend([f"Deleted: {item}" for item in deleted_windows[:4]])
        if len(deleted_windows) > 4:
            label_lines.append(f"+ {len(deleted_windows) - 4} more deleted")

    if not label_lines:
        return

    ax.text(
        0.01,
        0.98,
        "\n".join(label_lines),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "#cbd5e1"},
    )


def add_logo_highlights(ax, highlights=None, selected_highlights=None):
    for start, end in selected_highlights or []:
        if end > start:
            ax.axvspan(start - 0.5, end - 0.5, color="#d1d5db", alpha=0.5)
    for start, end in highlights or []:
        if end > start:
            ax.axvspan(start - 0.5, end - 0.5, color="#dbeafe", alpha=0.45)


def assign_motif_track_rows(motif_hits):
    placed = []
    row_end_positions = []
    for hit in sorted(motif_hits, key=lambda item: (item["start_idx"], item["end_idx"], item["pvalue"])):
        placed_hit = dict(hit)
        for row_idx, last_end in enumerate(row_end_positions):
            if placed_hit["start_idx"] >= last_end:
                row_end_positions[row_idx] = placed_hit["end_idx"]
                placed_hit["track_row"] = row_idx
                break
        else:
            row_end_positions.append(placed_hit["end_idx"])
            placed_hit["track_row"] = len(row_end_positions) - 1
        placed.append(placed_hit)
    return placed, max(1, len(row_end_positions))


def build_overlap_groups(motif_hits):
    if not motif_hits:
        return []

    groups = []
    current_group = [dict(sorted(motif_hits, key=lambda item: (item["start_idx"], item["end_idx"], item["pvalue"]))[0])]
    current_end = current_group[0]["end_idx"]

    for hit in sorted(motif_hits, key=lambda item: (item["start_idx"], item["end_idx"], item["pvalue"]))[1:]:
        hit_copy = dict(hit)
        if hit_copy["start_idx"] < current_end:
            current_group.append(hit_copy)
            current_end = max(current_end, hit_copy["end_idx"])
        else:
            groups.append(current_group)
            current_group = [hit_copy]
            current_end = hit_copy["end_idx"]

    groups.append(current_group)
    return groups


def summarize_overlap_groups(motif_hits):
    summaries = []
    for group_idx, group in enumerate(build_overlap_groups(motif_hits)):
        group_start = min(item["start_idx"] for item in group)
        group_end = max(item["end_idx"] for item in group)
        ordered_group = sorted(
            group,
            key=lambda item: (
                -(min(item["end_idx"], group_end) - max(item["start_idx"], group_start)),
                item["pvalue"],
                item["start_idx"],
                item["end_idx"],
            ),
        )
        representative = dict(ordered_group[0])
        representative_start = int(representative["start_idx"])
        representative_end = int(representative["end_idx"])
        representative["group_id"] = f"group_{group_idx}"
        representative["overlap_count"] = len(ordered_group)
        representative["group_hits"] = ordered_group
        representative["start_idx"] = group_start
        representative["end_idx"] = group_end
        representative["representative_overlap_bp"] = representative_end - representative_start
        summaries.append(representative)
    return summaries


def draw_motif_track(
    ax,
    motif_scan,
    seq_len,
    tick_positions,
    tick_labels,
    axis_label,
    highlights=None,
    selected_highlights=None,
):
    hits = summarize_overlap_groups(motif_scan.get("hits", []))
    add_logo_highlights(ax, highlights=highlights, selected_highlights=selected_highlights)
    ax.set_xlim(-0.5, seq_len - 0.5)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=0)
    ax.set_xlabel(axis_label)
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#94a3b8")

    if not hits:
        ax.set_ylim(-0.5, 0.5)
        ax.text(
            0.5,
            0.5,
            f"No FIMO hits in {motif_scan['database_label']} at p <= {motif_scan['pvalue_threshold']:.1e}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9,
            color="#475569",
        )
        ax.set_title("Motif hits", fontsize=10)
        return

    placed_hits, row_count = assign_motif_track_rows(hits)
    ax.set_ylim(-0.6, row_count - 0.4)
    ax.invert_yaxis()

    for row_idx in range(row_count):
        ax.hlines(row_idx, -0.5, seq_len - 0.5, color="#e2e8f0", linewidth=0.6, zorder=0)

    for hit_idx, hit in enumerate(placed_hits):
        width = max(0.8, hit["end_idx"] - hit["start_idx"])
        x0 = hit["start_idx"] - 0.5
        y0 = hit["track_row"] - 0.26
        color = MOTIF_TRACK_COLORS[hit_idx % len(MOTIF_TRACK_COLORS)]
        ax.add_patch(
            Rectangle(
                (x0, y0),
                width,
                0.52,
                facecolor=color,
                edgecolor="#0f172a",
                linewidth=0.6,
                alpha=0.8,
            )
        )
        ax.text(
            x0 + (width / 2.0),
            hit["track_row"],
            str(hit.get("overlap_count", 1)),
            ha="center",
            va="center",
            fontsize=8.5,
            fontweight="bold",
            color="white",
            clip_on=True,
        )

    shown = len(hits)
    total = motif_scan.get("total_hits", shown)
    suffix = f"{shown} overlap groups" if total == shown else f"{shown} overlap groups from {total} total hits"
    ax.set_title(
        f"Motif hits: {motif_scan['database_label']} ({suffix}, p <= {motif_scan['pvalue_threshold']:.1e})",
        fontsize=10,
    )


def plot_logo_matrix(
    attr,
    title,
    genomic_start_1,
    highlights=None,
    selected_highlights=None,
    y_limit=None,
):
    attr = normalize_attr_matrix(attr)
    df = pd.DataFrame(attr, columns=["A", "C", "G", "T"])
    fig_width = max(10, min(22, len(df) / 10))
    fig, ax = plt.subplots(figsize=(fig_width, 2.8))
    fig.subplots_adjust(left=0.055, right=0.995, bottom=0.22, top=0.88)
    logomaker.Logo(df, ax=ax)
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_xlim(-0.5, len(df) - 0.5)
    if y_limit is not None:
        ax.set_ylim(-float(y_limit), float(y_limit))
    ax.set_title(title)
    tick_count = min(6, len(df))
    tick_positions = np.linspace(0, len(df) - 1, num=tick_count, dtype=int)
    tick_labels = [str(genomic_start_1 + int(pos)) for pos in tick_positions]
    add_logo_highlights(ax, highlights=highlights, selected_highlights=selected_highlights)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=0)
    ax.set_xlabel("Genomic coordinate")
    return fig


def motif_hits_dataframe(motif_hits, coordinate_start_1):
    rows = []
    for hit in motif_hits:
        rows.append(
            {
                "motif": hit.get("display_label", hit.get("label", "motif")),
                "full_label": hit.get("label", hit.get("display_label", "motif")).replace("_", " "),
                "start": coordinate_start_1 + int(hit["start_idx"]),
                "end": coordinate_start_1 + int(hit["end_idx"]) - 1,
                "strand": hit.get("strand", "."),
                "pvalue": hit.get("pvalue"),
                "score": hit.get("score"),
                "matched_sequence": hit.get("matched_sequence", ""),
            }
        )
    return pd.DataFrame(rows)


def motif_group_chart_dataframe(motif_scan, coordinate_start_1):
    rows = []
    summarized_hits = summarize_overlap_groups(motif_scan.get("hits", []))
    placed_hits, row_count = assign_motif_track_rows(summarized_hits)
    sequence_length = int(motif_scan.get("sequence_length") or 0)
    track_start = int(coordinate_start_1)
    track_end = int(coordinate_start_1 + max(sequence_length, 1) - 1)
    for idx, hit in enumerate(placed_hits):
        group_hits = hit.get("group_hits", [])
        cluster_members = []
        for member in group_hits:
            member_overlap_bp = min(int(member["end_idx"]), int(hit["end_idx"])) - max(int(member["start_idx"]), int(hit["start_idx"]))
            cluster_members.append(
                f"{member.get('display_label', member.get('label', 'motif'))} "
                f"({coordinate_start_1 + int(member['start_idx'])}-{coordinate_start_1 + int(member['end_idx']) - 1}, "
                f"overlap={member_overlap_bp} bp, "
                f"p={float(member.get('pvalue', 1.0)):.2e})"
            )
        rows.append(
            {
                "group_id": hit["group_id"],
                "track_row": int(hit.get("track_row", 0)),
                "track_label": f"track_{int(hit.get('track_row', 0)) + 1}",
                "start": coordinate_start_1 + int(hit["start_idx"]),
                "end": coordinate_start_1 + int(hit["end_idx"]) - 1,
                "track_start": track_start,
                "track_end": track_end,
                "box_start": int(hit["start_idx"]),
                "box_end": int(hit["end_idx"]),
                "midpoint": (coordinate_start_1 + int(hit["start_idx"]) + coordinate_start_1 + int(hit["end_idx"]) - 1) / 2.0,
                "span_bp": max(1, int(hit["end_idx"]) - int(hit["start_idx"])),
                "motif": hit.get("display_label", hit.get("label", "motif")),
                "full_label": hit.get("label", hit.get("display_label", "motif")).replace("_", " "),
                "overlap_count": int(hit.get("overlap_count", 1)),
                "representative_overlap_bp": int(hit.get("representative_overlap_bp", hit["end_idx"] - hit["start_idx"])),
                "pvalue": float(hit.get("pvalue", 1.0)),
                "score": float(hit.get("score", 0.0)),
                "strand": hit.get("strand", "."),
                "matched_sequence": hit.get("matched_sequence", ""),
                "color": MOTIF_TRACK_COLORS[idx % len(MOTIF_TRACK_COLORS)],
                "cluster_members": " | ".join(cluster_members) if cluster_members else "None",
                "row_count": row_count,
            }
        )
    return pd.DataFrame(rows)


def build_interactive_motif_chart(chart_df, axis_label):
    row_count = max(1, int(chart_df["row_count"].max())) if not chart_df.empty else 1
    chart_height = max(140, min(280, row_count * 42 + 28))
    domain_start = int(chart_df["track_start"].min())
    domain_end = int(chart_df["track_end"].max())
    base = (
        alt.Chart(chart_df)
        .encode(
            x=alt.X(
                "start:Q",
                title=axis_label,
                scale=alt.Scale(domain=[domain_start, domain_end], zero=False),
            ),
            x2="end:Q",
            y=alt.Y(
                "track_label:N",
                sort=[f"track_{idx + 1}" for idx in range(row_count)],
                title=None,
                axis=None,
            ),
            tooltip=[
                alt.Tooltip("full_label:N", title="Representative motif"),
                alt.Tooltip("start:Q", title="Start"),
                alt.Tooltip("end:Q", title="End"),
                alt.Tooltip("pvalue:Q", title="Lowest p-value", format=".2e"),
                alt.Tooltip("overlap_count:Q", title="Overlapping motifs"),
                alt.Tooltip("cluster_members:N", title="Motifs in cluster"),
            ],
        )
    )

    boxes = base.mark_bar(
        cornerRadius=5,
        stroke="#0f172a",
        strokeWidth=1.0,
        size=26,
    ).encode(
        color=alt.Color("color:N", scale=None, legend=None),
        opacity=alt.value(0.9),
    )

    counts = base.mark_text(fontSize=13, fontWeight="bold", color="white").encode(
        x=alt.X("midpoint:Q"),
        text="overlap_count:Q",
    )

    labels = base.transform_filter(
        "datum.span_bp >= 22"
    ).mark_text(
        fontSize=10,
        fontWeight="bold",
        color="#111827",
        align="center",
        baseline="bottom",
        dy=-12,
    ).encode(
        x=alt.X("midpoint:Q"),
        text=alt.Text("motif:N"),
    )

    return (
        (boxes + counts + labels)
        .properties(height=chart_height)
        .configure_view(stroke="#d1d5db")
        .configure_axis(grid=False)
    )


def get_selected_group_id(selection_state, selection_name="motif_pick"):
    try:
        selected = selection_state["selection"][selection_name]
    except Exception:
        return None

    if not selected:
        return None
    if isinstance(selected, dict):
        value = selected.get("group_id")
        if isinstance(value, list):
            return value[0] if value else None
        return value
    if isinstance(selected, list) and selected:
        first = selected[0]
        if isinstance(first, dict):
            value = first.get("group_id")
            if isinstance(value, list):
                return value[0] if value else None
            return value
    return None


def render_interactive_motif_panel(
    motif_scan,
    coordinate_start_1,
    axis_label,
    chart_key,
):
    chart_df = motif_group_chart_dataframe(motif_scan, coordinate_start_1)
    if chart_df.empty:
        unfiltered_total = motif_scan.get("unfiltered_total_hits", motif_scan.get("total_hits", 0))
        if unfiltered_total:
            st.caption(
                f"No motif hits overlap {motif_scan.get('z_window_size', 10)} bp attribution windows "
                f"with z > {motif_scan.get('z_threshold', 1.0):.1f}."
            )
        else:
            st.caption(f"No FIMO hits in {motif_scan['database_label']} at p <= {motif_scan['pvalue_threshold']:.1e}")
        return

    st.altair_chart(
        build_interactive_motif_chart(chart_df, axis_label),
        width="stretch",
        key=chart_key,
    )
    st.caption(
        f"Motif hits: {len(chart_df)} overlap groups from {motif_scan['total_hits']} filtered hits "
        f"({motif_scan.get('unfiltered_total_hits', motif_scan['total_hits'])} FIMO hits before filtering). "
        f"Filter: {motif_scan.get('z_window_size', 10)} bp attribution windows with z > {motif_scan.get('z_threshold', 1.0):.1f}. "
        "Hover over a box to see the full motif list in that cluster."
    )


def make_text_page(title, lines, max_chars_per_line=120):
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis("off")
    wrapped_lines = []
    for line in lines:
        line = str(line)
        if len(line) <= max_chars_per_line:
            wrapped_lines.append(line)
            continue
        for i in range(0, len(line), max_chars_per_line):
            wrapped_lines.append(line[i : i + max_chars_per_line])
    fig.text(0.05, 0.97, title, fontsize=16, fontweight="bold", va="top")
    fig.text(0.05, 0.94, "\n".join(wrapped_lines), fontsize=9, va="top", family="monospace")
    return fig


def build_pdf_report(figures):
    pdf_buffer = BytesIO()
    with PdfPages(pdf_buffer) as pdf:
        for fig in figures:
            pdf.savefig(fig, bbox_inches="tight")
    pdf_buffer.seek(0)
    return pdf_buffer.getvalue()


# =========================
# APP STATE
# =========================
if "project_path" not in st.session_state:
    st.session_state.project_path = DEFAULT_PROJ_IN
if "model_run" not in st.session_state:
    st.session_state.model_run = DEFAULT_MODEL_RUN
if "cell_type" not in st.session_state:
    st.session_state.cell_type = DEFAULT_CELL_TYPE
if "fold_mode" not in st.session_state:
    st.session_state.fold_mode = DEFAULT_FOLD_MODE
if "project_option" not in st.session_state:
    st.session_state.project_option = os.path.basename(os.path.dirname(DEFAULT_PROJ_IN.rstrip("/")))
if "genome_key" not in st.session_state:
    st.session_state.genome_key = DEFAULT_GENOME_KEY
if "genome_fasta" not in st.session_state:
    st.session_state.genome_fasta = resolve_genome_fasta(st.session_state.genome_key)
if "coordinate_region" not in st.session_state:
    st.session_state.coordinate_region = GENOME_SPECS[st.session_state.genome_key]["default_coordinate"]
if "coordinate_input_text" not in st.session_state:
    st.session_state.coordinate_input_text = st.session_state.coordinate_region
if "sequence_label" not in st.session_state:
    st.session_state.sequence_label = ""
if "active_original_seq" not in st.session_state:
    st.session_state.active_original_seq = ""
if "active_seq" not in st.session_state:
    st.session_state.active_seq = ""
if "active_start_1" not in st.session_state:
    st.session_state.active_start_1 = 1
if "active_end_1" not in st.session_state:
    st.session_state.active_end_1 = 1
if "original_start_1" not in st.session_state:
    st.session_state.original_start_1 = 1
if "original_end_1" not in st.session_state:
    st.session_state.original_end_1 = 1
if "results" not in st.session_state:
    st.session_state.results = None
if "zoom_center_genomic" not in st.session_state:
    st.session_state.zoom_center_genomic = None
if "zoom_half_width" not in st.session_state:
    st.session_state.zoom_half_width = DEFAULT_WINDOW_HALF_WIDTH
if "selection_start_genomic" not in st.session_state:
    st.session_state.selection_start_genomic = None
if "selection_end_genomic" not in st.session_state:
    st.session_state.selection_end_genomic = None
if "match_query" not in st.session_state:
    st.session_state.match_query = ""
if "accumulated_matches" not in st.session_state:
    st.session_state.accumulated_matches = []
if "logo_y_limit" not in st.session_state:
    st.session_state.logo_y_limit = DEFAULT_LOGO_YLIM
if "fimo_enabled" not in st.session_state:
    st.session_state.fimo_enabled = DEFAULT_FIMO_ENABLED
if "motif_database" not in st.session_state:
    st.session_state.motif_database = DEFAULT_MOTIF_DATABASE
if "fimo_pvalue" not in st.session_state:
    st.session_state.fimo_pvalue = DEFAULT_FIMO_PVALUE
if "last_search_added" not in st.session_state:
    st.session_state.last_search_added = ""
if "pending_search_selection" not in st.session_state:
    st.session_state.pending_search_selection = []
if "deleted_windows" not in st.session_state:
    st.session_state.deleted_windows = []
if "output_selected_intervals" not in st.session_state:
    st.session_state.output_selected_intervals = []
if "output_selected_labels" not in st.session_state:
    st.session_state.output_selected_labels = []
if "browser_highlight_intervals" not in st.session_state:
    st.session_state.browser_highlight_intervals = []
if "browser_seq" not in st.session_state:
    st.session_state.browser_seq = ""
if "browser_start_1" not in st.session_state:
    st.session_state.browser_start_1 = 1
if "browser_end_1" not in st.session_state:
    st.session_state.browser_end_1 = 1
if "browser_active_start_1" not in st.session_state:
    st.session_state.browser_active_start_1 = 1


st.title("ChromBPNet Genome Browser")
st.caption("Browse a genomic window, accumulate motif matches across searches, and iteratively delete selected intervals while keeping averaged fold predictions up to date.")

with st.sidebar:
    st.subheader("Inputs")
    project_options = discover_project_paths()
    option_labels = get_project_option_labels(project_options)
    if not option_labels:
        option_labels = [os.path.basename(os.path.dirname(st.session_state.project_path.rstrip("/"))) or "default"]
    project_option_value = st.session_state.project_option if st.session_state.project_option in option_labels else option_labels[0]
    project_option_input = st.selectbox(
        "Project dataset",
        option_labels,
        index=option_labels.index(project_option_value),
        key="project_dataset_select",
    )
    project_path_input = resolve_project_path(project_options, project_option_input, st.session_state.project_path)
    project_model_catalog = discover_models(project_path_input.strip())
    project_run_options = sorted(project_model_catalog)
    if project_option_input != st.session_state.project_option:
        st.session_state.project_option = project_option_input
        st.session_state.project_path = project_path_input.strip()
        st.session_state.model_run = project_run_options[0] if project_run_options else ""
        project_cell_options = sorted(project_model_catalog.get(st.session_state.model_run, {})) if st.session_state.model_run else []
        st.session_state.cell_type = project_cell_options[0] if project_cell_options else ""
        st.rerun()

    model_catalog = project_model_catalog
    run_options = sorted(model_catalog)
    if run_options:
        model_run_value = st.session_state.model_run if st.session_state.model_run in run_options else run_options[0]
        model_run_input = st.selectbox(
            "Model run",
            run_options,
            index=run_options.index(model_run_value),
            key=f"model_run_select::{project_option_input}",
        )
        cell_options = sorted(model_catalog[model_run_input])
        cell_value = st.session_state.cell_type if st.session_state.cell_type in cell_options else cell_options[0]
        cell_type_input = st.selectbox(
            "Cell type",
            cell_options,
            index=cell_options.index(cell_value),
            key=f"cell_type_select::{project_option_input}::{model_run_input}",
        )
        fold_mode_options = ["1 fold", "All folds"]
        fold_mode_value = st.session_state.fold_mode if st.session_state.fold_mode in fold_mode_options else DEFAULT_FOLD_MODE
        fold_mode_input = st.selectbox(
            "Fold mode",
            fold_mode_options,
            index=fold_mode_options.index(fold_mode_value),
        )
    else:
        model_run_input = ""
        cell_type_input = ""
        fold_mode_input = DEFAULT_FOLD_MODE
        st.warning("No chrombpnet_nobias.h5 models found under project_path/01-models.")

    genome_options = list(GENOME_SPECS)
    genome_labels = {key: GENOME_SPECS[key]["label"] for key in genome_options}
    genome_key_value = st.session_state.genome_key if st.session_state.genome_key in genome_options else DEFAULT_GENOME_KEY
    genome_key_input = st.selectbox(
        "Genome assembly",
        genome_options,
        index=genome_options.index(genome_key_value),
        format_func=lambda key: genome_labels[key],
        key="genome_assembly_select",
        on_change=lambda: apply_genome_selection(st.session_state.genome_assembly_select),
    )
    if genome_key_input != st.session_state.genome_key:
        st.rerun()

    coordinate_input = st.text_input("Genome coordinate", key="coordinate_input_text")
    reload_clicked = st.button("Load inputs", use_container_width=True)

    if reload_clicked:
        st.session_state.project_option = project_option_input
        st.session_state.project_path = project_path_input.strip()
        st.session_state.model_run = model_run_input
        st.session_state.cell_type = cell_type_input
        st.session_state.fold_mode = fold_mode_input
        st.session_state.genome_key = genome_key_input
        st.session_state.genome_fasta = resolve_genome_fasta(genome_key_input)
        st.session_state.coordinate_region = coordinate_input.strip() or GENOME_SPECS[genome_key_input]["default_coordinate"]
        st.session_state.sequence_label = ""
        st.session_state.active_original_seq = ""
        st.session_state.active_seq = ""
        st.session_state.active_start_1 = 1
        st.session_state.active_end_1 = 1
        st.session_state.original_start_1 = 1
        st.session_state.original_end_1 = 1
        st.session_state.results = None
        st.session_state.zoom_center_genomic = None
        st.session_state.zoom_half_width = DEFAULT_WINDOW_HALF_WIDTH
        st.session_state.selection_start_genomic = None
        st.session_state.selection_end_genomic = None
        st.session_state.match_query = ""
        st.session_state.accumulated_matches = []
        st.session_state.last_search_added = ""
        st.session_state.pending_search_selection = []
        st.session_state.deleted_windows = []
        st.session_state.output_selected_intervals = []
        st.session_state.output_selected_labels = []
        st.session_state.browser_highlight_intervals = []
        st.session_state.browser_seq = ""
        st.session_state.browser_start_1 = 1
        st.session_state.browser_end_1 = 1
        st.session_state.browser_active_start_1 = 1
        st.rerun()


project_path = st.session_state.project_path
genome_fasta = st.session_state.genome_fasta or resolve_genome_fasta(st.session_state.genome_key)
coordinate_region = st.session_state.coordinate_region
model_catalog = discover_models(project_path)

model_files, model_folds = get_fold_model_paths(
    model_catalog,
    st.session_state.model_run,
    st.session_state.cell_type,
    st.session_state.fold_mode,
)
if not model_files:
    st.error(f"No model folds found for the selected project/run/cell type under `{project_path}`.")
    st.stop()

missing_model_files = [model_file for model_file in model_files if not os.path.exists(model_file)]
if missing_model_files:
    st.error("Model file(s) not found:\n" + "\n".join(missing_model_files))
    st.stop()

if not os.path.exists(genome_fasta):
    st.error(f"Genome FASTA not found:\n{genome_fasta}")
    st.stop()

pred_models, count_models = load_models(tuple(model_files))

if not st.session_state.active_original_seq:
    try:
        chrom, active_start_1, active_end_1, seq = extract_centered_genome_sequence(
            coordinate_region, genome_fasta, target_len=INPUTLEN
        )
        st.session_state.sequence_label = chrom
        st.session_state.active_start_1 = active_start_1
        st.session_state.active_end_1 = active_end_1
        st.session_state.original_start_1 = active_start_1
        st.session_state.original_end_1 = active_end_1
        st.session_state.active_original_seq = sanitize_sequence(seq)
    except Exception as exc:
        st.error(f"Could not load sequence.\n\n{exc}")
        st.stop()

if not st.session_state.active_seq:
    st.session_state.active_seq = st.session_state.active_original_seq

chrom = st.session_state.sequence_label or "custom"
active_original_seq = sanitize_sequence(st.session_state.active_original_seq)
active_seq = sanitize_sequence(st.session_state.active_seq)
active_start_1 = st.session_state.active_start_1
active_end_1 = st.session_state.active_end_1
region_center_1 = active_start_1 + (len(active_seq) // 2)

if st.session_state.zoom_center_genomic is None:
    st.session_state.zoom_center_genomic = region_center_1

zoom_center_1 = clip_int(st.session_state.zoom_center_genomic, active_start_1, active_end_1)
zoom_half_width = clip_int(
    st.session_state.zoom_half_width,
    25,
    max(25, min(DEFAULT_WINDOW_HALF_WIDTH * 4, (len(active_seq) // 2) or 25)),
)

window_start_1, window_end_1 = genomic_window(active_start_1, active_end_1, zoom_center_1, zoom_half_width)
window_start_idx = genome_to_index(window_start_1, active_start_1)
window_end_idx = genome_to_index(window_end_1, active_start_1) + 1
visible_seq = active_seq[window_start_idx:window_end_idx]

if not st.session_state.browser_seq:
    st.session_state.browser_seq = visible_seq
    st.session_state.browser_start_1 = window_start_1
    st.session_state.browser_end_1 = window_end_1
    st.session_state.browser_active_start_1 = active_start_1

browser_seq = st.session_state.browser_seq
browser_start_1 = st.session_state.browser_start_1
browser_end_1 = st.session_state.browser_end_1
browser_active_start_1 = st.session_state.browser_active_start_1
browser_start_idx = genome_to_index(browser_start_1, browser_active_start_1)

if st.session_state.selection_start_genomic is None:
    st.session_state.selection_start_genomic = max(window_start_1, zoom_center_1 - 10)
if st.session_state.selection_end_genomic is None:
    st.session_state.selection_end_genomic = min(window_end_1, st.session_state.selection_start_genomic + 9)

selection_start_1 = clip_int(st.session_state.selection_start_genomic, window_start_1, window_end_1)
selection_end_1 = clip_int(st.session_state.selection_end_genomic, selection_start_1, window_end_1)
st.session_state.selection_start_genomic = selection_start_1
st.session_state.selection_end_genomic = selection_end_1

if st.session_state.results is None:
    with st.spinner("Loading predictions..."):
        st.session_state.results = run_all_outputs(
            active_original_seq,
            active_seq,
            pred_models,
            count_models,
            tuple(model_files),
        )

results = st.session_state.results
fallback_highlight = (
    genome_to_index(selection_start_1, browser_start_1),
    genome_to_index(selection_end_1, browser_start_1) + 1,
)

with st.sidebar:
    st.subheader("Session")
    st.write(f"Genome assembly: `{GENOME_SPECS[st.session_state.genome_key]['label']}`")
    st.write(f"Project: `{project_path}`")
    st.write(f"Model run: `{st.session_state.model_run}`")
    st.write(f"Cell type: `{st.session_state.cell_type}`")
    st.write(f"Fold mode: `{st.session_state.fold_mode}`")
    st.write(f"Folds used: `{len(model_folds)}`")
    st.write(f"Chromosome: `{chrom}`")
    st.write(f"Active span: `{chrom}:{active_start_1}-{active_end_1}`")
    st.write(f"Visible window: `{chrom}:{window_start_1}-{window_end_1}`")
    st.write(f"Device: `{DEVICE}`")
    if len(active_seq) < INPUTLEN:
        st.warning("Chromosome boundaries prevented a full 2114 bp extraction.")
    st.caption("Genome coordinates shorter than 2114 bp are extended to a centered 2114 bp window when possible.")
    current_logo_y_limit = max(0.005, min(float(st.session_state.logo_y_limit), 0.75))
    logo_y_limit = st.slider(
        "Logo y-axis scale",
        min_value=0.005,
        max_value=0.75,
        value=current_logo_y_limit,
        step=0.005,
        format="%.3f",
    )
    st.session_state.logo_y_limit = logo_y_limit
    motif_database_options = discover_available_motif_databases()
    st.subheader("Motif scan")
    st.session_state.fimo_enabled = st.checkbox(
        "Run FIMO on displayed sequences",
        value=bool(st.session_state.fimo_enabled),
        help="Adds a motif-hit track below each logo plot using the selected motif database.",
    )
    if motif_database_options:
        current_motif_database = (
            st.session_state.motif_database
            if st.session_state.motif_database in motif_database_options
            else motif_database_options[0]
        )
        st.session_state.motif_database = st.selectbox(
            "Motif database",
            options=motif_database_options,
            index=motif_database_options.index(current_motif_database),
            disabled=not st.session_state.fimo_enabled,
        )
    else:
        st.warning("No JASPAR or HOCOMOCO motif databases were found for FIMO.")
    st.session_state.fimo_pvalue = st.number_input(
        "FIMO p-value threshold",
        min_value=1e-8,
        max_value=1.0,
        value=float(st.session_state.fimo_pvalue),
        format="%.1e",
        disabled=not st.session_state.fimo_enabled or not motif_database_options,
    )


st.subheader("1. Genome browser")
with st.form("zoom_form"):
    zoom_cols = st.columns(2)
    with zoom_cols[0]:
        zoom_center_input = st.number_input(
            "Move window center",
            min_value=active_start_1,
            max_value=active_end_1,
            value=zoom_center_1,
            step=1,
        )
    with zoom_cols[1]:
        zoom_half_width_input = st.number_input(
            "Zoom half-width (bp each side)",
            min_value=25,
            max_value=max(25, min(DEFAULT_WINDOW_HALF_WIDTH * 4, (len(active_seq) // 2) or 25)),
            value=zoom_half_width,
            step=25,
        )
    zoom_submit = st.form_submit_button("Apply zoom")

if zoom_submit:
    st.session_state.zoom_center_genomic = int(zoom_center_input)
    st.session_state.zoom_half_width = int(zoom_half_width_input)
    st.session_state.selection_start_genomic = None
    st.session_state.selection_end_genomic = None
    new_window_start_1, new_window_end_1 = genomic_window(
        active_start_1,
        active_end_1,
        int(zoom_center_input),
        int(zoom_half_width_input),
    )
    new_window_start_idx = genome_to_index(new_window_start_1, active_start_1)
    new_window_end_idx = genome_to_index(new_window_end_1, active_start_1) + 1
    st.session_state.browser_seq = active_seq[new_window_start_idx:new_window_end_idx]
    st.session_state.browser_start_1 = new_window_start_1
    st.session_state.browser_end_1 = new_window_end_1
    st.session_state.browser_active_start_1 = active_start_1
    st.session_state.browser_highlight_intervals = selected_intervals_from_rows(st.session_state.accumulated_matches)
    st.rerun()

match_query_input = st.text_input(
    "Search sequence in visible window",
    value=st.session_state.match_query,
    help="Type a DNA string such as ACGTG. Matches are searched in the currently visible edited sequence.",
)

clean_match_query = sanitize_sequence(match_query_input)
query_invalid = validate_sequence(clean_match_query) if clean_match_query else []

if match_query_input != st.session_state.match_query:
    st.session_state.match_query = match_query_input
    st.session_state.pending_search_selection = []

matches = []
search_options = []
if clean_match_query and not query_invalid:
    matches = find_sequence_matches(browser_seq, clean_match_query)
    search_options = [
        f"{idx + 1}. {chrom}:{browser_start_1 + start}-{browser_start_1 + end - 1} ({end - start} bp)"
        for idx, (start, end) in enumerate(matches)
    ]
    valid_pending = [item for item in st.session_state.pending_search_selection if item in search_options]
    if valid_pending != st.session_state.pending_search_selection:
        st.session_state.pending_search_selection = valid_pending
elif query_invalid:
    st.error(f"Search contains invalid bases: {', '.join(query_invalid)}")

if search_options:
    result_cols = st.columns([2.6, 1.0])
    with result_cols[0]:
        pending_selected = st.multiselect(
            "Search results",
            options=search_options,
            default=st.session_state.pending_search_selection,
        )
        st.session_state.pending_search_selection = pending_selected
    with result_cols[1]:
        add_matches_clicked = st.button("Add selected", use_container_width=True)
elif clean_match_query and not query_invalid:
    st.info("No visible matches found for the typed sequence.")
    add_matches_clicked = False
else:
    add_matches_clicked = False

new_match_rows = []
if add_matches_clicked:
    if not clean_match_query:
        st.warning("Enter a sequence before adding matches.")
    elif query_invalid:
        st.error(f"Search contains invalid bases: {', '.join(query_invalid)}")
    elif not st.session_state.pending_search_selection:
        st.warning("Choose at least one search result from the dropdown.")
    else:
        existing_keys = {
            (row["query"], row["start_idx"], row["end_idx"])
            for row in st.session_state.accumulated_matches
        }
        selected_indexes = [search_options.index(label) for label in st.session_state.pending_search_selection]
        for idx in selected_indexes:
            start, end = matches[idx]
            global_start = browser_start_idx + start
            global_end = browser_start_idx + end
            key = (clean_match_query, global_start, global_end)
            if key in existing_keys:
                continue
            new_match_rows.append(
                {
                    "selected": True,
                    "row_action": "Keep",
                    "row_replacement": "",
                    "query": clean_match_query,
                    "coordinate": f"{chrom}:{active_start_1 + global_start}-{active_start_1 + global_end - 1}",
                    "width_bp": global_end - global_start,
                    "sequence": active_seq[global_start:global_end],
                    "start_idx": global_start,
                    "end_idx": global_end,
                }
            )
        st.session_state.accumulated_matches.extend(new_match_rows)
        st.session_state.last_search_added = clean_match_query
        st.session_state.browser_highlight_intervals = selected_intervals_from_rows(st.session_state.accumulated_matches)

accumulated_df = pd.DataFrame(normalize_match_rows(st.session_state.accumulated_matches))
if not accumulated_df.empty:
    if "row_action" not in accumulated_df.columns:
        accumulated_df["row_action"] = "Keep"
    if "row_replacement" not in accumulated_df.columns:
        accumulated_df["row_replacement"] = ""

    header_cols = st.columns([0.9, 1.0, 1.2, 0.9, 1.8, 0.8, 1.0])
    header_labels = ["selected", "Edit", "New bases", "query", "coordinate", "width_bp", "sequence"]
    for col, label in zip(header_cols, header_labels):
        col.markdown(f"**{label}**")

    updated_rows = []
    for row in accumulated_df.to_dict("records"):
        row_id = f"{row['query']}_{row['start_idx']}_{row['end_idx']}"
        selected_key = f"match_selected_{row_id}"
        action_key = f"match_action_{row_id}"
        replacement_key = f"match_replacement_{row_id}"

        if selected_key not in st.session_state:
            st.session_state[selected_key] = bool(row.get("selected", False))
        if action_key not in st.session_state:
            st.session_state[action_key] = row.get("row_action", "Keep") or "Keep"
        if replacement_key not in st.session_state:
            st.session_state[replacement_key] = normalize_editor_text(row.get("row_replacement", ""))

        row_cols = st.columns([0.9, 1.0, 1.2, 0.9, 1.8, 0.8, 1.0])
        with row_cols[0]:
            selected_value = st.checkbox("Select", value=st.session_state[selected_key], key=selected_key, label_visibility="collapsed")
        with row_cols[1]:
            action_value = st.selectbox(
                "Edit",
                options=["Keep", "Delete", "Replace"],
                index=["Keep", "Delete", "Replace"].index(st.session_state[action_key] if st.session_state[action_key] in ["Keep", "Delete", "Replace"] else "Keep"),
                key=action_key,
                label_visibility="collapsed",
            )
        with row_cols[2]:
            replacement_value = st.text_input(
                "New bases",
                value=st.session_state[replacement_key],
                key=replacement_key,
                label_visibility="collapsed",
                disabled=action_value != "Replace",
            )
        with row_cols[3]:
            st.write(row["query"])
        with row_cols[4]:
            st.write(row["coordinate"])
        with row_cols[5]:
            st.write(row["width_bp"])
        with row_cols[6]:
            st.write(row["sequence"])

        updated_row = dict(row)
        updated_row["selected"] = bool(selected_value)
        updated_row["row_action"] = action_value
        updated_row["row_replacement"] = normalize_editor_text(replacement_value)
        updated_rows.append(updated_row)

    st.session_state.accumulated_matches = normalize_match_rows(updated_rows)

    remove_cols = st.columns([1.0, 1.0, 1.0, 2.0])
    with remove_cols[0]:
        remove_matches_clicked = st.button("Remove selected", use_container_width=True)
    with remove_cols[1]:
        apply_row_edits_clicked = st.button("Apply row edits", use_container_width=True)
    with remove_cols[2]:
        reset_clicked = st.button("Reset sequence", use_container_width=True)
    if remove_matches_clicked:
        st.session_state.accumulated_matches = [
            row for row in st.session_state.accumulated_matches if not row["selected"]
        ]
        st.session_state.browser_highlight_intervals = selected_intervals_from_rows(st.session_state.accumulated_matches)
        st.rerun()
    if reset_clicked:
        st.session_state.active_seq = st.session_state.active_original_seq
        st.session_state.active_start_1 = st.session_state.original_start_1
        st.session_state.active_end_1 = st.session_state.original_end_1
        st.session_state.zoom_center_genomic = st.session_state.original_start_1 + (len(st.session_state.active_original_seq) // 2)
        st.session_state.results = run_all_outputs(
            st.session_state.active_original_seq,
            st.session_state.active_original_seq,
            pred_models,
            count_models,
            tuple(model_files),
        )
        st.session_state.accumulated_matches = []
        st.session_state.deleted_windows = []
        st.session_state.output_selected_intervals = []
        st.session_state.output_selected_labels = []
        st.session_state.browser_highlight_intervals = []
        st.session_state.browser_seq = ""
        st.session_state.browser_start_1 = 1
        st.session_state.browser_end_1 = 1
        st.session_state.browser_active_start_1 = 1
        st.rerun()
    if apply_row_edits_clicked:
        rows_to_edit = [
            row for row in st.session_state.accumulated_matches if row.get("selected") and row.get("row_action") != "Keep"
        ]
        if not rows_to_edit:
            st.warning("Choose at least one selected row with Delete or Replace in the table.")
        else:
            invalid_row = None
            empty_replace_row = None
            for row in rows_to_edit:
                if row["row_action"] == "Replace":
                    clean_row_replacement = sanitize_sequence(normalize_editor_text(row.get("row_replacement", "")))
                    invalid_chars = validate_sequence(clean_row_replacement)
                    if not clean_row_replacement:
                        empty_replace_row = row["coordinate"]
                        break
                    if invalid_chars:
                        invalid_row = (row["coordinate"], invalid_chars)
                        break
                    row["row_replacement"] = clean_row_replacement
            if empty_replace_row:
                st.error(f"Enter replacement bases for {empty_replace_row}.")
            elif invalid_row:
                st.error(f"Invalid bases for {invalid_row[0]}: {', '.join(invalid_row[1])}")
            else:
                edited_intervals = merge_intervals([(row["start_idx"], row["end_idx"]) for row in rows_to_edit])
                edited_labels = [
                    f"{chrom}:{active_start_1 + start}-{active_start_1 + end - 1}"
                    for start, end in edited_intervals
                ]
                updated_seq = edit_rows(active_seq, rows_to_edit)
                new_start_1 = active_start_1
                new_end_1 = active_end_1
                if len(updated_seq) < INPUTLEN:
                    updated_seq, new_start_1, new_end_1 = extend_deleted_sequence_from_genome(
                        updated_seq,
                        chrom,
                        active_start_1,
                        active_end_1,
                        genome_fasta,
                        target_len=INPUTLEN,
                    )
                st.session_state.active_seq = updated_seq
                st.session_state.active_start_1 = new_start_1
                st.session_state.active_end_1 = new_end_1
                st.session_state.pending_search_selection = []
                st.session_state.deleted_windows = [
                    row["coordinate"] for row in rows_to_edit if row["row_action"] == "Delete"
                ]
                st.session_state.output_selected_intervals = edited_intervals
                st.session_state.output_selected_labels = edited_labels
                with st.spinner("Updating predictions and attributions..."):
                    st.session_state.results = run_all_outputs(
                        active_original_seq,
                        updated_seq,
                        pred_models,
                        count_models,
                        tuple(model_files),
                    )
                for row in st.session_state.accumulated_matches:
                    if row.get("selected"):
                        row["row_action"] = "Keep"
                        row["row_replacement"] = ""
                st.rerun()

selected_rows = [row for row in st.session_state.accumulated_matches if row["selected"]]
selected_intervals = [(row["start_idx"], row["end_idx"]) for row in selected_rows]
selected_intervals = merge_intervals(selected_intervals)
highlight_intervals = []
for start, end in st.session_state.browser_highlight_intervals:
    local_start = start - browser_start_idx
    local_end = end - browser_start_idx
    if local_end > 0 and local_start < len(browser_seq):
        highlight_intervals.append((max(0, local_start), min(len(browser_seq), local_end)))

if not highlight_intervals:
    highlight_intervals = [fallback_highlight]

browser_attr = extract_attr_window(get_attr(browser_seq, tuple(model_files))[0], len(browser_seq))
browser_motif_scan = None
if st.session_state.fimo_enabled and discover_available_motif_databases():
    browser_motif_scan = filter_motif_scan_by_attr(
        scan_sequence_with_fimo(
            browser_seq,
            st.session_state.motif_database,
            st.session_state.fimo_pvalue,
        ),
        browser_attr,
    )
browser_logo_fig = plot_logo_matrix(
    browser_attr,
    "Current logo",
    browser_start_1,
    highlights=highlight_intervals,
    y_limit=st.session_state.logo_y_limit,
)
st.pyplot(browser_logo_fig)
if browser_motif_scan is not None:
    render_interactive_motif_panel(
        browser_motif_scan,
        browser_start_1,
        "Genomic coordinate",
        "genome_browser_motif_track",
    )
st.code(
    f"{chrom}:{browser_start_1}-{browser_end_1}\n"
    f"{make_position_ruler(browser_start_1, browser_end_1)}\n"
    f"{browser_seq}"
)
browser_report_fig = make_text_page(
    "ChromBPNet Browser Snapshot",
    [
        "Sidebar",
        f"Genome assembly: {GENOME_SPECS[st.session_state.genome_key]['label']}",
        f"Project: {project_path}",
        f"Model run: {st.session_state.model_run}",
        f"Cell type: {st.session_state.cell_type}",
        f"Fold mode: {st.session_state.fold_mode}",
        f"Folds used: {len(model_folds)}",
        f"Logo y-axis scale: {st.session_state.logo_y_limit:.3f}",
        f"Device: {DEVICE}",
        "",
        f"Chromosome: {chrom}",
        f"Active span: {chrom}:{active_start_1}-{active_end_1}",
        f"Visible window: {chrom}:{browser_start_1}-{browser_end_1}",
        f"Search query: {clean_match_query or 'None'}",
        f"Count prediction before: {results['cnt_before']:.2f} +/- {results['cnt_before_std']:.2f}",
        f"Count prediction after: {results['cnt_after']:.2f} +/- {results['cnt_after_std']:.2f}",
        "",
        f"{chrom}:{browser_start_1}-{browser_end_1}",
        make_position_ruler(browser_start_1, browser_end_1),
        browser_seq,
    ],
)
can_edit_match = bool(selected_intervals)
if not can_edit_match:
    st.caption("Add matches from one or more searches, then select rows in the table for deletion or replacement.")


st.subheader("Outputs")
count_test = compute_paired_count_test(results["cnt_before_folds"], results["cnt_after_folds"])
count_datapoint_lines = format_count_datapoints(results["cnt_before_folds"], results["cnt_after_folds"])
output_cols = st.columns([1.0, 1.2])
with output_cols[0]:
    st.markdown(
        (
            "<div style='font-size:0.9rem; line-height:1.5;'>"
            f"<div><strong>Count prediction before</strong>: {results['cnt_before']:.2f} +/- {results['cnt_before_std']:.2f}</div>"
            f"<div><strong>Count prediction after</strong>: {results['cnt_after']:.2f} +/- {results['cnt_after_std']:.2f}</div>"
            f"<div><strong>Count delta</strong>: {results['cnt_after'] - results['cnt_before']:.2f}</div>"
            "</div>"
        ),
        unsafe_allow_html=True,
    )
    if count_test["pvalue"] is not None:
        st.markdown(
            (
                "<div style='font-size:0.82rem; color:#475569; line-height:1.45; margin-top:0.35rem;'>"
                f"<div><strong>Paired t-test p-value</strong>: {count_test['pvalue']:.3g}</div>"
                f"<div>t = {count_test['statistic']:.3f}, n = {count_test['n_folds']} folds</div>"
                "</div>"
            ),
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            f"<div style='font-size:0.82rem; color:#475569; line-height:1.45; margin-top:0.35rem;'>{count_test['message']}</div>",
            unsafe_allow_html=True,
        )
with output_cols[1]:
    count_fig = plot_count_summary(results["cnt_before_folds"], results["cnt_after_folds"])
    st.pyplot(count_fig)

output_selected_intervals = st.session_state.output_selected_intervals
output_selected_labels = st.session_state.output_selected_labels

profile_fig = plot_profile_overlay_marked(
    results["prof_before"],
    results["prof_after"],
    highlight_windows=output_selected_intervals,
    selected_windows=selected_intervals,
    deleted_windows=st.session_state.deleted_windows,
    prof_before_std=results["prof_before_std"],
    prof_after_std=results["prof_after_std"],
    prof_before_folds=results["prof_before_folds"],
    prof_after_folds=results["prof_after_folds"],
)
add_window_labels(
    profile_fig.axes[0],
    selected_windows=[row["coordinate"] for row in selected_rows] or output_selected_labels,
    deleted_windows=st.session_state.deleted_windows,
)
st.pyplot(profile_fig)
count_report_fig = make_text_page(
    "Count Summary",
    [
        f"Count prediction before: {results['cnt_before']:.4f} +/- {results['cnt_before_std']:.4f}",
        f"Count prediction after: {results['cnt_after']:.4f} +/- {results['cnt_after_std']:.4f}",
        f"Count delta: {results['cnt_after'] - results['cnt_before']:.4f}",
        f"Paired t-test p-value: {count_test['pvalue']:.6g}" if count_test["pvalue"] is not None else f"Paired t-test: {count_test['message']}",
        f"Paired t-test statistic: {count_test['statistic']:.6f}" if count_test["statistic"] is not None else "",
        f"Folds used: {count_test['n_folds']}",
        "",
        "Fold datapoints",
        *count_datapoint_lines,
    ],
)
browser_pdf_bytes = build_pdf_report([browser_report_fig, browser_logo_fig, count_report_fig, count_fig, profile_fig])
st.download_button(
    "Save browser as PDF",
    data=browser_pdf_bytes,
    file_name=f"{chrom}_{browser_start_1}_{browser_end_1}_browser.pdf",
    mime="application/pdf",
    use_container_width=True,
)
before_logo_fig = plot_logo_matrix(
    results["attr_before"],
    "Before logo",
    active_start_1,
    highlights=output_selected_intervals,
    selected_highlights=selected_intervals,
    y_limit=st.session_state.logo_y_limit,
)
st.pyplot(before_logo_fig)
before_motif_scan = (
    filter_motif_scan_by_attr(
        scan_sequence_with_fimo(
            active_original_seq,
            st.session_state.motif_database,
            st.session_state.fimo_pvalue,
        ),
        results["attr_before"],
    )
    if st.session_state.fimo_enabled and discover_available_motif_databases()
    else None
)
if before_motif_scan is not None:
    render_interactive_motif_panel(
        before_motif_scan,
        active_start_1,
        "Genomic coordinate",
        "genome_before_motif_track",
    )
after_logo_fig = plot_logo_matrix(
    results["attr_after"],
    "After logo",
    active_start_1,
    highlights=output_selected_intervals,
    y_limit=st.session_state.logo_y_limit,
)
st.pyplot(after_logo_fig)
after_motif_scan = (
    filter_motif_scan_by_attr(
        scan_sequence_with_fimo(
            active_seq,
            st.session_state.motif_database,
            st.session_state.fimo_pvalue,
        ),
        results["attr_after"],
    )
    if st.session_state.fimo_enabled and discover_available_motif_databases()
    else None
)
if after_motif_scan is not None:
    render_interactive_motif_panel(
        after_motif_scan,
        active_start_1,
        "Genomic coordinate",
        "genome_after_motif_track",
    )

report_summary_fig = make_text_page(
    "ChromBPNet Genome Browser Report",
    [
        "Sidebar",
        f"Genome assembly: {GENOME_SPECS[st.session_state.genome_key]['label']}",
        f"Project: {project_path}",
        f"Model run: {st.session_state.model_run}",
        f"Cell type: {st.session_state.cell_type}",
        f"Fold mode: {st.session_state.fold_mode}",
        f"Folds used: {len(model_folds)}",
        f"Chromosome: {chrom}",
        f"Active span: {chrom}:{active_start_1}-{active_end_1}",
        f"Visible window: {chrom}:{window_start_1}-{window_end_1}",
        f"Device: {DEVICE}",
        f"Logo y-axis scale: {st.session_state.logo_y_limit:.3f}",
        "",
        f"Count prediction before: {results['cnt_before']:.2f} +/- {results['cnt_before_std']:.2f}",
        f"Count prediction after: {results['cnt_after']:.2f} +/- {results['cnt_after_std']:.2f}",
        f"Count delta: {results['cnt_after'] - results['cnt_before']:.2f}",
        f"Paired t-test p-value: {count_test['pvalue']:.6g}" if count_test["pvalue"] is not None else f"Paired t-test: {count_test['message']}",
        f"Paired t-test statistic: {count_test['statistic']:.6f}" if count_test["statistic"] is not None else "",
        f"Deleted windows: {', '.join(st.session_state.deleted_windows) if st.session_state.deleted_windows else 'None'}",
        "",
        f"{chrom}:{browser_start_1}-{browser_end_1}",
        make_position_ruler(browser_start_1, browser_end_1),
        browser_seq,
    ],
)
count_detail_fig = make_text_page("Count Fold Datapoints", count_datapoint_lines)
pdf_bytes = build_pdf_report(
    [
        report_summary_fig,
        count_detail_fig,
        count_fig,
        browser_logo_fig,
        profile_fig,
        before_logo_fig,
        after_logo_fig,
    ]
)
st.download_button(
    "Save whole page as PDF",
    data=pdf_bytes,
    file_name=f"{chrom}_{active_start_1}_{active_end_1}_chrombpnet_report.pdf",
    mime="application/pdf",
    use_container_width=True,
)
