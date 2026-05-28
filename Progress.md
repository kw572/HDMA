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

- Pending first login, sync, submit, and monitor cycle.

## Next Actions

- Commit and push the new pipeline files.
- Log into USC HPC.
- Update the HPC checkout from GitHub.
- Submit `02b-compendium/00-run_all_compendium.sbatch`.
- Record the first scheduler/runtime error and fix it.
- Keep iterating until the job starts cleanly.
