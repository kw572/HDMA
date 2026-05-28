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
