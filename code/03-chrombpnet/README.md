
Overview of the steps involved in training ChromBPNet models on the atlas
and analyzing outputs.

We show how to load and use trained ChromBPNet models in the tutorial at [`code/05-misc/04-ChromBPNet_use_cases.ipnyb`](https://github.com/GreenleafLab/HDMA/blob/main/code/05-misc/04-ChromBPNet_use_cases.ipynb) ([html](https://greenleaflab.github.io/HDMA/code/05-misc/04-ChromBPNet_use_cases.html)).

## Tools used

- ChromBPNet: https://github.com/kundajelab/chrombpnet
- TF-MoDISco (tfmodiscolite implementation): https://github.com/jmschrei/tfmodisco-lite
- gimmemotifs: https://gimmemotifs.readthedocs.io/en/master/
- Fi-NeMo: https://github.com/austintwang/finemo_gpu
- tangermeme: https://github.com/jmschrei/tangermeme


## `chrombpnet_utils`

Python utilities for downstream ChromBPNet related analysis, plotting, etc.


## `tangermeme_utils`

Python utilities for downstream ChromBPNet related analysis, particularly
for interfacing with the [tangermeme](https://tangermeme.readthedocs.io/en/latest/index.html) package.
These are mainly used for the *in silico* experimentation.


## 0. Preparing inputs, `00-inputs`

These scripts pertain to generating the necessary inputs for training ChromBPNet models.

- `01`: get fragments per cluster per sample from ArchR projects
- `01`: alternatively, `00-prepare_inputs/01-get_fragments_from_rds.R` can extract the same per-cluster fragment TSV layout directly from a Seurat `.rds` plus fragment file, which is useful when you want to bypass ArchR
- `02`: process fragments for calling peaks (generating pseudoreplicates)
- `03`: group fragments across samples per cell type, and sort fragments
- `04`: call peaks with macs2 on pseudoreplicates
- `04b`: call peaks on the fragments of each cluster (no pseudoreps), to generate p-value and fold-change signal bigwig tracks
- `05`: format peaks for ChromBPNet
- `06`: make bigwigs
- `07`: make chromosome splits for test/train/validation folds
- `08`: get GC-matched negative peaks
- `09`: tabix-index peaks for viewing in browser
- `10`: randomly subset backgrounds to 100 regions per cell type to use for downstream analysis

### Notes on local compatibility changes

The local edits used for running `00-inputs` on this Mac fall into two groups.

Mac-only / local shell compatibility:
- wrapper scripts were adjusted to tolerate paths with spaces, local `bash` behavior, and the absence of cluster-only assumptions such as `module load`
- several runners were changed to avoid fragile `parallel` usage for local execution
- local runs use the explicit conda env / local binary paths when needed

General portability / non-Mac-specific changes:
- `01` can now start from a Seurat `.rds` plus fragment file instead of requiring ArchR
- chromosome names are normalized so either `1` or `chr1`-style fragment inputs can be matched against the reference genome files
- downstream preprocessing no longer assumes human-only references; scripts now respect `chromsizes`, `ref_fasta`, and computed genome size from `config.sh`
- local ChromBPNet prep was patched to work with numeric chromosome names such as zebrafish `1`-`25`
- peak formatting now drops out-of-bounds summit-centered windows, which is useful for smaller genomes and is not specific to macOS
- empty blacklist files are handled cleanly, which is useful when no curated blacklist is available for the target genome

Zebrafish-specific configuration in the local test run:
- the current local config uses `Danio rerio` GRCz11 reference files and chromosome sizes
- `07` split generation was adapted to build folds from the chromosomes actually present in the configured genome, instead of relying on the original human chromosome lists




## 1. Training models, `01-train_models`

These scripts run the ChromBPNet workflow, which involves training a bias model,
using the bias model to train chromatin accessibility models in each cell type,
calculate contribution scores, generate model predictions, and perform motif discovery.

- `01`: train bias models
- `02`: train ChromBPNet models with a given bias model, initial interpretation with 30k subsample peaks
- `03`: perform model QC and make figures with performance metrics
- `04`: get contribution scores (_aka_ DeepLIFT or DeepSHAP scores) as h5 files per fold per cell type
- `05`: average contribs across folds, separately for profile and counts
- `06`: deep TF-MoDISco run with 1M seqlets on averaged contribs
- `07`: generate contribution score bigwigs, for viewing in browser
- `08`: predict accessibility in peak regions across folds, and average the predictions (bias-corrected) to generate bw
- `09`: prepare supplementary table with ChromBPNet performance and QC metrics

### Notes on interpretation outputs

- `02` training via `chrombpnet pipeline` already runs a built-in interpretation/QC pass.
- That training-time interpretation uses a peak subsample when there are more than `30,000` peaks, writing files under `auxiliary/interpret_subsample/`.
- `04` and beyond are still needed because they generate the final downstream interpretation products:
  - per-fold contribution exports on the provided peak set
  - averaged SHAPs across kept folds
  - final MoDISco runs on averaged fold outputs
  - final browser bigwigs and averaged predictions
- In other words, training-time SHAP/MoDISco is mainly for fold-level QC, while `04+` produces the dataset-level outputs used for downstream analysis.

### USC CARC notes

These are the current assumptions for running the ChromBPNet shell / sbatch workflow on USC CARC.

- Current config-driven local run paths use `code/03-chrombpnet/data/local_run_chr/work`, as defined in `config.sh`.
- GPU jobs now default to `--partition=gpu`; CPU helper jobs default to `--partition=main`.
- GPU request defaults use generic `--gres=gpu:1` rather than model-specific `gpu:a40:1`.
- The current default GPU module stack in the active training scripts is:
  - `PRE_MODULES='legacy/CentOS7 gcc/8.3.0'`
  - `CUDA_MODULE='cuda/11.2.0'`
  - `CUDNN_MODULE='cudnn/8.1.0.77-11.2'`
- On this cluster, the wrapper job itself may need a shorter submit-time override than the child jobs. A safe pattern has been:
  - `sbatch --time=12:00:00 00-run_all_training.sbatch`
- To limit a run to selected datasets, pass `CHROMBPNET_DATASET_FILTER_REGEX`, for example:
  - `sbatch --time=12:00:00 --export=ALL,CHROMBPNET_DATASET_FILTER_REGEX='^6-' 00-run_all_training.sbatch`
- The training/downstream submitters now avoid duplicate resubmission by checking for matching active SLURM job names.
- Incomplete fold reruns are cleaned automatically before resubmission if the fold is not already queued/running.
- The wrapper and downstream QC scripts now follow `config.sh` / `base_dir` instead of older hardcoded `data/local_run/work` paths.
- `03-model_QC.Rmd` now reconstructs cluster-level cell counts from unique barcodes in `cluster_fragments/fragments/*__sorted.tsv` when running the local workflow, so reruns do not depend on the broader manuscript metadata table.
- If `06` finished `counts_modisco_output.h5` but failed during `modisco report` or `modisco meme`, use `06-rerun_report_and_meme.sbatch` to repair those outputs without recomputing MoDISco motifs. It respects `BIAS_PARAMS`, `MEME_DB`, and optional `CHROMBPNET_DATASET_FILTER_REGEX`.




## 2. Assembly of motif compendium, `02-compendium`

These scripts pertain to assembling a motif compendium from across the atlas,
calling instances of each motif in each cell type, and then annotating those instances.

- `01`: convert MoDISco CWMs to PFMs to prepare the input for `gimme cluster`:
this is run _per organ_, and motifs from all cell types (models) within an organ are concatenated.
Done separately for positive and negative patterns.
- `02`: for each organ, run `gimme cluster` on PFMs of motifs from all cell types in that organ, to get broad clusterings.
Done separately for positive and negative patterns.
- `02b`: run `gimme cluster` again, using as input all the clusters from the _per organ_ runs of `gimme cluster` in Step 02.
Dine separately for positive and negative patterns. This produces "superclusters" of patterns from across cell types.
- `03`: merge similar MoDISco patterns within each supercluster together, using the `modiscolite.SimilarPatternCollapser` functionality.
This produces one MoDISco h5 object per supercluster, containing merged, non-redundant patterns.
- `04`: combine MoDISco h5 objects from the merging step into one h5 object, match to known patterns, assign unique index, and prepare as MEME format
- `05`: call hits in each cell type, using the unified motifs
- `06`: generate HTML reports to assess the quality of hit calling with unified motif set, filter & label hits using manual annotation
- `06b`: filter out low QC hits, reconcile overlapping hits, and get genomic annotations
- `07`: run NucleoATAC per cell type

### `02-compendium` versus `02b-compendium`

- `02-compendium` is the full atlas-wide motif unification workflow. It converts per-dataset MoDISco outputs into PFMs, clusters motifs within organs and then across organs, merges redundant patterns, annotates the unified compendium, calls motif hits across datasets, reconciles overlaps, and optionally runs NucleoATAC.
- `02b-compendium` is a narrower CREsted-style alternative. It builds a merged compendium directly from the per-dataset `counts_modisco_output.h5` files, annotates those merged patterns, renders a heatmap-style summary, then runs Fi-NeMo on each dataset and maps the resulting hits back onto the merged `02b` compendium IDs.
- Operationally, `02-compendium/00-run_all_compendium.sbatch` orchestrates the end-to-end atlas compendium plus downstream hit-calling and NucleoATAC stages. `02b-compendium/00-run_all_compendium.sbatch` currently runs only the compendium-build and annotation/heatmap phases; Fi-NeMo submission remains a separate `03-call_finemo_hits.sh` step.
- Use `02-compendium` when downstream analyses need the legacy unified motif catalog and reconciled hit set. Use `02b-compendium` when you want a CREsted / Fi-NeMo-centric view that preserves per-dataset motif instances while linking them back to merged motif families.


## 2b. CREsted-style motif compendium, `02b-compendium`

This is a lighter alternative to `02-compendium` that merges motifs across
datasets using CREsted's motif-processing utilities directly on the final
per-dataset TF-MoDISco H5 outputs from `01-train_models/06`.

- `00`: HPC wrapper for the CREsted-style compendium run
- `01`: collect kept datasets from `01-models/qc`, map them to `counts_modisco_output.h5`,
  run `crested.tl.modisco.process_patterns`, and write a compact merged-motif bundle
- `02`: read the merged-pattern bundle, add motif annotations when report HTMLs are available,
  and render a CREsted-style clustermap from the pattern matrix

Outputs are written under `base_dir/02b-compendium/crested_patterns` and include:
- `all_patterns.pkl`: CREsted merged-pattern object
- `pattern_manifest.tsv`: one row per merged pattern with representative motif metadata
- `pattern_matrix.tsv` / `pattern_matrix.npy`: class-by-pattern matrix from CREsted
- `matched_modisco_h5.json`: dataset-to-input H5 mapping
- `run_summary.json`: run parameters and counts
- `annotation/pattern_annotations.tsv`: merged-pattern annotations with motif/TF matches when available
- `annotation/annotation_summary.json`: notes on whether full motif annotation was possible
- `plots/pattern_clustermap.png`: CREsted-style clustered heatmap of merged motif importance

### USC CARC notes for `02b-compendium`

- The wrapper expects a Python environment with `numpy`, `pandas`, `modiscolite`, and `memelite` available.
- By default it activates `conda activate modiscolite`, assuming CREsted is available via `CRESTED_REPO` or already installed in that environment.
- If you prefer using an installed CREsted environment, make sure it also has motif dependencies and a backend-capable Keras stack available.
- You can override the Python environment with:
  - `CRESTED_ENV_NAME`
  - `CRESTED_VENV`
  - `CRESTED_PYTHON_BIN`
  - `CRESTED_REPO` if CREsted is available as a checkout rather than an installed package
- Optional motif-to-TF annotation can also use `CRESTED_MOTIF_TO_TF_FILE` if you want to point at a cached `motif_tf_collection.tsv` explicitly.
- Typical submit pattern:
  - `sbatch code/03-chrombpnet/02b-compendium/00-run_all_compendium.sbatch`
  - `sbatch --export=ALL,CHROMBPNET_DATASET_FILTER_REGEX='^(1_Jaw_Hyoid|4_Frontonasal)$' code/03-chrombpnet/02b-compendium/00-run_all_compendium.sbatch`


## 3. Analysis of TF binding site syntax, `03-syntax`

- `00`: compute the mean and summed importance score for all trimmed CWMs
- `01`: global overview of de novo motifs and predictive instances
- `02`: visualization of genomic tracks at example loci
- `03`: visualization of motif instances (hits) for a few examples
- `04`: in silico analysis of cooperativity
  - `04a`: assemble the set of composites to test
  - `04b`: run the in silico analysis
  - `04c`: summarize the results and generate plots
  - `04d`: make animations showing hard and soft synergy
- `05`: in silico analysis of cell type-specific cooperativity
  - `05a`: run the analysis
  - `05b`: visualize results
- `06`: *in silico* ablations of negative motifs
- `07`: compute motifs with enriched co-occurrence
- `08`: investigate overlaps in peaks and motif instance calls across clusters


_**NOTE**_: variant scoring using ChromBPNet models is done in the `code/06-variants` directory.
