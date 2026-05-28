# 02b-Compendium Progress

## Scope

Set up and run a CREsted-style motif compendium pipeline at `code/03-chrombpnet/02b-compendium` that consumes `01-train_models` TF-MoDISco outputs on USC HPC.

## Current Status

- Added `02b-compendium` HPC wrapper and Python driver.
- Added `02b-compendium` output paths to `config.sh` and `config.ucmab.sh`.
- Added usage notes to `code/03-chrombpnet/README.md`.
- Local syntax checks passed for the new shell and Python files.
- Preparing HPC sync and first submit attempt.

## Files Added Or Updated

- `code/03-chrombpnet/02b-compendium/00-run_all_compendium.sbatch`
- `code/03-chrombpnet/02b-compendium/01-build_crested_compendium.sh`
- `code/03-chrombpnet/02b-compendium/01-build_crested_compendium.py`
- `code/03-chrombpnet/config.sh`
- `code/03-chrombpnet/config.ucmab.sh`
- `code/03-chrombpnet/README.md`

## Verification Log

### 2026-05-28 Local syntax

- `bash -n code/03-chrombpnet/02b-compendium/00-run_all_compendium.sbatch`
- `bash -n code/03-chrombpnet/02b-compendium/01-build_crested_compendium.sh`
- `python3 -m py_compile code/03-chrombpnet/02b-compendium/01-build_crested_compendium.py`
- Result: passed

## Errors And Fixes

### Error 1: local default `python3` missing scientific packages

- Symptom: `ModuleNotFoundError: No module named 'numpy'`
- Where: local smoke test of `01-build_crested_compendium.py`
- Resolution: treated as local environment mismatch, not a code bug. HPC wrapper was designed around an explicit env activation path instead of assuming system Python.

### Error 2: local `crested` env lacked backend and motif stack

- Symptom: CREsted import path required a backend, and the available local env also lacked `modiscolite`
- Resolution:
  - Added fallback bootstrap logic in `01-build_crested_compendium.py` to load the CREsted motif module directly from a checkout when `crested.tl` lazy import is blocked by backend setup.
  - Updated the wrapper default to `CRESTED_ENV_NAME=modiscolite`.
  - Documented that the runtime env must provide `numpy`, `pandas`, `modiscolite`, and `memelite`, with `CRESTED_REPO` available if CREsted is not installed into that same env.

## HPC Run Log

### 2026-05-28 pilot job `8997738`

- Checkout found at `/scratch1/kuangtse/HDMA`
- CREsted checkout found separately at `/scratch1/kuangtse/code/CREsted`
- HPC branch `mac-preprocess` could not fast-forward because the run checkout has diverged local commits and unrelated modified files
- To avoid overwriting that worktree, synced only `code/03-chrombpnet/02b-compendium/` into the existing checkout with `rsync`
- Submitted pilot job:
  - `sbatch --export=ALL,CRESTED_ENV_NAME=modiscolite,CRESTED_REPO=/scratch1/kuangtse/code/CREsted,CHROMBPNET_DATASET_FILTER_REGEX="^1_Jaw_Hyoid$" 00-run_all_compendium.sbatch`

## Next Actions

- Patch wrapper to auto-discover CREsted on HPC.
- Resync `02b-compendium` into the HPC checkout.
- Resubmit pilot job.
- Record the next runtime error or first successful motif-processing step.

## Additional Errors And Fixes

### Error 3: HPC checkout could not fast-forward to pushed branch

- Symptom: `fatal: Not possible to fast-forward, aborting.`
- Cause: `/scratch1/kuangtse/HDMA` is a dirty and diverged run checkout with local commits and many unrelated modified files.
- Resolution: avoided merge/rebase in-place and synced only the new `02b-compendium` folder into that run checkout.

### Error 4: pilot HPC job could not locate CREsted at runtime

- Job: `8997738`
- Symptom: `ImportError: Unable to import 'crested'. Activate an environment with CREsted installed or set CRESTED_REPO to a CREsted checkout.`
- Cause: even though a CREsted checkout exists at `/scratch1/kuangtse/code/CREsted`, the job runtime did not end up with a usable `CRESTED_REPO` path.
- Resolution: updated `01-build_crested_compendium.sh` to auto-detect common HPC CREsted checkout paths such as `/scratch1/$USER/code/CREsted` and export `CRESTED_REPO` before launching Python.

### Error 5: CREsted checkout path was available, but bootstrap still required a successful top-level package import

- Job: `8997741`
- Symptom: job log showed `crested_repo=/scratch1/kuangtse/code/CREsted`, followed by the same `Unable to import crested` failure.
- Cause: the Python bootstrap still tried `import crested` after adding the checkout `src` path, and the `modiscolite` environment is lighter than a full CREsted install.
- Resolution: changed `01-build_crested_compendium.py` so checkout discovery only requires `src/crested` to exist. If the top-level import still fails, the script now creates a minimal `crested` package namespace pointing at the checkout and lets the motif module load directly from source.

### Error 6: `modiscolite` HPC env is missing `loguru`

- Job: `8997745`
- Symptom: `ModuleNotFoundError: No module named 'loguru'`
- Cause: CREsted's motif helpers import `loguru` for logging, but the available HPC `modiscolite` env does not include it.
- Resolution: added a small fallback shim in `01-build_crested_compendium.py` that injects a minimal `loguru.logger` backed by Python's standard `logging` module before loading CREsted's motif source files.

### Error 7: remote job logs were accidentally deleted during folder resync

- Jobs affected: `8997745` log artifacts
- Symptom: `sacct` showed the job failed, but `/scratch1/kuangtse/HDMA/code/03-chrombpnet/02b-compendium/Log/` was empty when inspected afterward.
- Cause: the manual `rsync --delete` used to sync `02b-compendium/` from local to HPC removed the remote `Log/` directory because the local source tree does not keep job logs.
- Resolution: treat `Log/` as HPC runtime state and avoid syncing over it before collecting failure output. Future monitor cycles should inspect logs first, then sync only after the relevant traceback has been copied into `Progress.md`.

### Error 8: `modiscolite` HPC env is also missing `anndata` and `scanpy`

- Job: `8998068`
- Symptom: `ModuleNotFoundError: No module named 'anndata'`
- Cause: CREsted's `_tfmodisco.py` imports `anndata`, `scanpy`, and `crested.utils._logging` at module import time, even though the `process_patterns` path we are calling does not actually use the H5AD/scanpy functionality.
- Resolution: extended the bootstrap shim in `01-build_crested_compendium.py` to:
  - create minimal `anndata` and `scanpy` placeholder modules for this motif-only runtime
  - register `crested.utils` as a package namespace
  - load `crested.utils._logging` directly from the CREsted checkout before loading `_tfmodisco.py`

## Monitoring Notes

### 2026-05-28 13:13 PT

- Re-read the preserved `8998068` traceback and confirmed it is the pre-fix failure for missing `anndata`.
- Verified that the local and HPC copies of `01-build_crested_compendium.py` now match at 344 lines and include the `loguru`, `anndata`, and `scanpy` bootstrap shims.
- Conclusion: the next step is a fresh pilot resubmission; no additional code change is justified until a post-shim traceback or first successful processing milestone is captured.

### 2026-05-28 13:15 PT

- Submitted a fresh one-dataset pilot after confirming the synced runtime shims:
  - `sbatch --export=ALL,CRESTED_ENV_NAME=modiscolite,CHROMBPNET_DATASET_FILTER_REGEX="^1_Jaw_Hyoid$" 00-run_all_compendium.sbatch`
  - New job id: `8998081`
- Early monitoring result:
  - `8998081` entered `RUNNING` on `a01-04`
  - wrapper banner printed normally with the expected `modisco_root`, `output_dir`, and `crested_repo`
  - no immediate Python traceback appeared within the first ~40 seconds
- This is the first pilot that has cleanly progressed past the previous import-time failures for `loguru` and `anndata`.
- Next monitor step: wait for either the first motif-processing output or the next runtime traceback before changing code again.

### 2026-05-28 13:18 PT

- Final status for `8998081`:
  - `sacct` reports `COMPLETED` with exit code `0:0`
  - elapsed time: `00:00:43`
- Successful runtime milestones captured from `Log/chrombpnet-crested-compendium-8998081.out`:
  - opened `/scratch1/kuangtse/HDMA/code/03-chrombpnet/data/NCC_36hpf/work/01-models/modisco/bias_1_Jaw_Hyoid_thresh0.4/1_Jaw_Hyoid/counts_modisco_output.h5`
  - processed motif similarity matches for the `1_Jaw_Hyoid` pilot dataset
  - completed post-hoc merging with `32` patterns remaining after one iteration
- Output bundle confirmed at `/scratch1/kuangtse/HDMA/code/03-chrombpnet/data/NCC_36hpf/work/02b-compendium/crested_patterns`:
  - `all_patterns.pkl`
  - `classes.json`
  - `matched_modisco_h5.json`
  - `missing_modisco_h5.json`
  - `pattern_manifest.tsv`
  - `pattern_matrix.npy`
  - `pattern_matrix.tsv`
  - `run_summary.json`
- `run_summary.json` confirms:
  - `n_requested_datasets = 1`
  - `n_matched_datasets = 1`
  - `n_missing_datasets = 0`
  - `n_patterns = 32`
  - `pattern_matrix_shape = [1, 32]`
- Conclusion: the `02b-compendium` pilot now starts cleanly and completes successfully on USC HPC for the one-dataset validation case.

## Annotation And Heatmap Follow-up

### 2026-05-28 13:30 PT

- Added a new post-processing stage under `code/03-chrombpnet/02b-compendium`:
  - `02-annotate_and_heatmap.sh`
  - `02-annotate_and_heatmap.py`
- Updated `00-run_all_compendium.sbatch` so future `02b-compendium` runs automatically:
  - build the CREsted merged-pattern bundle
  - annotate merged motifs from available modisco report HTMLs
  - render a CREsted-style clustermap from the merged pattern matrix
- Pilot post-process outputs were generated successfully on HPC at:
  - `/scratch1/kuangtse/HDMA/code/03-chrombpnet/data/NCC_36hpf/work/02b-compendium/crested_patterns/annotation/pattern_annotations.tsv`
  - `/scratch1/kuangtse/HDMA/code/03-chrombpnet/data/NCC_36hpf/work/02b-compendium/crested_patterns/annotation/annotation_summary.json`
  - `/scratch1/kuangtse/HDMA/code/03-chrombpnet/data/NCC_36hpf/work/02b-compendium/crested_patterns/plots/pattern_clustermap.png`
- Pilot annotation result:
  - `annotation_mode = motif_only`
  - `n_report_htmls = 1`
  - `n_patterns_with_matches = 26`
- Reason for `motif_only` rather than TF-level annotation:
  - the run could recover motif-name matches from `motifs.html`
  - no cached `motif_tf_collection.tsv` was available in the current HPC env, so TF alias expansion was not added yet

### Error 9: annotation step failed because `pandas.read_html` required `lxml`

- Symptom: `ImportError: Missing optional dependency 'lxml'. Use pip or conda to install lxml.`
- Cause: the HPC `modiscolite` env lacked `lxml`, but `crested.tl.modisco.find_pattern_matches` relies on parsing `motifs.html` tables.
- Resolution: patched `02-annotate_and_heatmap.py` with a BeautifulSoup-based fallback HTML table parser and injected it into the CREsted motif runtime before calling `find_pattern_matches`.

### Error 10: fallback HTML parser left q-value columns as strings

- Symptom: `TypeError: '<' not supported between instances of 'str' and 'float'`
- Cause: q-value columns parsed from HTML needed explicit numeric coercion before CREsted compared them with the p-value threshold.
- Resolution: updated the fallback parser to coerce `qval*`, `pval*`, and `num_seqlets` columns to numeric types.

### Error 11: single-row pilot matrix could not be hierarchically clustered by seaborn

- Symptom: `ValueError: The number of observations cannot be determined on an empty distance matrix.`
- Cause: the one-dataset validation matrix has only one row, so row clustering is undefined.
- Resolution: updated the clustermap renderer to automatically disable row or column clustering whenever the filtered matrix dimension is `< 2`.

## Full Atlas Run

### 2026-05-28 14:28 PT

- Confirmed that local QC keep files are empty or absent for this workspace, so `02b-compendium` will use its fallback dataset discovery from the modisco directory.
- Confirmed the current full discovered dataset set under `bias_1_Jaw_Hyoid_thresh0.4`:
  - `1_Jaw_Hyoid`
  - `3_Ventral_oral_Joint_mandible`
  - `4_Frontonasal`
  - `6_Maxilla_I`
  - `7_Teeth`
  - `8_Maxilla_II`
  - `9_Intermediate_Jaw`
  - `Gills_Ventral_3-6Arch`
  - `Opercle`
  - `Ventral_2nd_Arch`
  - `Ventral_Intermediate`
- Submitted the full all-cluster `02b-compendium` run with no dataset filter:
  - `sbatch --export=ALL,CRESTED_ENV_NAME=modiscolite 00-run_all_compendium.sbatch`
  - Job id: `9002341`
- Early log check confirms the intended scope:
  - `dataset_filter=<all datasets>`
  - job state entered `RUNNING` on node `d05-37`

### 2026-05-28 14:40 PT

- Follow-up monitoring from the local workstation was blocked by a transient SSH connectivity failure to USC HPC.
- Exact symptom from a minimal probe:
  - `ssh: connect to host discovery.usc.edu port 22: Operation timed out`
  - verbose SSH showed consecutive timeouts to `10.72.0.13` and `10.72.0.14`
- Impact:
  - could not fetch the latest `sacct` status for job `9002341`
  - could not inspect the current `.out`/`.err` logs or output files during this monitor cycle
- Resolution:
  - no pipeline code change was indicated
  - leave the full-run monitor active and retry on the next heartbeat once cluster connectivity recovers

### 2026-05-28 14:45 PT

- SSH connectivity to USC HPC recovered on the next monitor cycle.
- Final status for the full all-cluster job `9002341`:
  - `sacct` reports `COMPLETED` with exit code `0:0`
  - elapsed time: `00:05:38`
- The full compendium run completed cleanly across all discovered datasets.
- `run_summary.json` confirms:
  - `n_requested_datasets = 11`
  - `n_matched_datasets = 11`
  - `n_missing_datasets = 0`
  - `n_patterns = 45`
  - `pattern_matrix_shape = [11, 45]`
- `annotation/annotation_summary.json` confirms:
  - `annotation_mode = motif_only`
  - `n_report_htmls = 11`
  - `n_patterns_with_matches = 39`
- The resulting `pattern_matrix.tsv` now contains 11 cluster rows:
  - `1_Jaw_Hyoid`
  - `3_Ventral_oral_Joint_mandible`
  - `4_Frontonasal`
  - `6_Maxilla_I`
  - `7_Teeth`
  - `8_Maxilla_II`
  - `9_Intermediate_Jaw`
  - `Gills_Ventral_3-6Arch`
  - `Opercle`
  - `Ventral_2nd_Arch`
  - `Ventral_Intermediate`
- Conclusion: `02b-compendium` now consolidates the full atlas ChromBPNet motif set in CREsted style and produces the expected multi-cluster heatmap input and motif-annotation outputs on HPC.

## HOCOMOCO14 Re-annotation

### 2026-05-28 14:54 PT

- Updated `02-annotate_and_heatmap.py` so annotation no longer depends on the older pre-rendered modisco report labels for TF naming.
- The annotation step now:
  - downloads the official HOCOMOCO v14 H14CORE MEME database from `https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE/formatted_motifs/H14CORE_meme_format.meme`
  - downloads the official HOCOMOCO v14 H14CORE metadata from `https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE/H14CORE_annotation.jsonl`
  - matches each merged CREsted representative motif directly against the HOCOMOCO14 motif database using CREsted's TOMTOM-style scoring
  - expands TF names and synonyms from the official HOCOMOCO14 JSONL
- Re-ran `02-annotate_and_heatmap.sh` on the completed 11-cluster compendium output.
- New annotation summary:
  - `annotation_mode = hocomoco14`
  - `n_hocomoco14_motifs = 1595`
  - `n_patterns_with_matches = 32`
  - `annotation_min_score = 6.0`
- Example successful HOCOMOCO14 TF annotations now present in `pattern_annotations.tsv`:
  - `CTCF/CTCFL`
  - `DLX1/DLX2/DLX4/DLX5/DLX6`
  - `ATF3/ATF6A/BATF2/BATF3/JUNB`
  - `RXRA/RXRB/RXRG/COT1`
- Updated outputs remain under:
  - `/scratch1/kuangtse/HDMA/code/03-chrombpnet/data/NCC_36hpf/work/02b-compendium/crested_patterns/annotation/`
  - `/scratch1/kuangtse/HDMA/code/03-chrombpnet/data/NCC_36hpf/work/02b-compendium/crested_patterns/plots/`

### Error 12: HOCOMOCO14 MEME parser initially failed on the `w=` header format

- Symptom: `ValueError: invalid literal for int() with base 10: ''`
- Cause: the HOCOMOCO14 MEME header uses `w= 10` formatting, while the first parser version assumed `w=10` without whitespace.
- Resolution: updated the parser to extract motif width with a regex that tolerates optional spaces after `w=`.
