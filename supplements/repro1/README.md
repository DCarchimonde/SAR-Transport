# SAR-Transport reproducibility supplement repro1

This independently versioned supplement accompanies **v1.0.0**. It adds the original training, qualification and export implementations, frozen configurations, environment records and prepared training inputs used for the reported study. The v1.0.0 tag, five release assets, task assignments, predictions and scientific results are unchanged. No new training, inference, numerical scoring or resampling was performed to assemble this supplement.

The underlying scientific result is `d5d3227b8a0269e2a896056ba41abf7d502d768f`; the manuscript/resource RC2 is `a21719a84edaab3176b65db1d8151ca335ac3a31`. Original source identities and SHA256 values are in `SOURCE_IDENTITY.tsv`. Files from the Ki implementation freeze `6a2b2b1c6f3d580eb0c26db34d3608cda038b2b2` retain their exact frozen bytes, including the older package initializer needed by that historical implementation.

## Files and use

1. Obtain both archives from [SAR-Transport v1.0.0](https://github.com/DCarchimonde/SAR-Transport/releases/tag/v1.0.0) and extract them into one resource directory.
2. Download this entire `supplements/repro1` directory, including `TRAINING_INPUTS.zip`. GitHub's repository archive or `git clone` includes it. Check `SHA256SUMS.txt` before use.
3. Use `REPRODUCIBILITY_ASSET_MAP.tsv` to find the implementation, configuration and inputs associated with each reported result. `RESTORE_MAP.tsv` points to unchanged files already in v1.0.0 instead of duplicating them.
4. If historical directory paths are needed, run the copy-only helper:

```text
python stage_layout.py --resource-root /path/to/v1-resource --destination /path/to/new-historical-layout
```

This helper does not load pickles, download raw activities, execute a study script or fit a model. `TRAINING_INPUTS.zip` restores the original `dataset.pkl`, `preprocessors.pkl`, identity universe and qualified/conflict observation tables. The two pickle hashes match the original input manifest. They are Python/scikit-learn artifacts, not portable task APIs; use the recorded Python and package versions. The uncompressed training inputs require about 518 MB of storage. The public JSON/TSV task API does not require these files.

## What is supported

| Reproduction level | Assets and limit |
| --- | --- |
| Load tasks, reconstruct task views and score supplied predictions | Existing v1.0.0 code, tables, labels and saved predictions. No historical model dependency is required. |
| Inspect the exact neural/RF generation implementation | Original `run_sar_transport.py`, `sar_transport_methods.py`, preprocessing code, protocol, freeze chain, method-adaptation log and environment records. |
| Refit the historical model families | Implementation, seeds, hyperparameters, prepared training inputs and the required ActFound backbone are supplied. Fits were not rerun or validated in this packaging pass. Numerical identity across environments is not promised. |
| Replay inference from the exact original trained checkpoints | Trained neural/RF checkpoints are **not included**. This layer is not supplied; do not describe the supplement as checkpoint replay. The original prediction tables and checkpoint hashes remain available for evaluation and provenance. Checkpoints are not required for scoring those tables or refitting from the supplied inputs. |
| Reconstruct original task derivation | IC50 qualification/features/export scripts, Ki identity precheck, frozen identity metadata/roles, numerical-qualification implementation and input maps are provided. Upstream raw activities remain at their pinned public source. No fresh ChEMBL rebuild or live-metadata replay was validated. |
| Reproduce P0 and Ki policy logic | Original P0 implementation and inputs are already in v1.0.0; the frozen Ki runner and its exact core/baseline modules are mapped here. P0 remains exploratory and Ki K1/K2 establish neither advantage nor equivalence. |

## Historical paths, guards and environments

The original IC50 environment is recorded in `historical/reports/sar_transport_v1/EXECUTION_ENVIRONMENT.json`: Python 3.12.14, PyTorch 2.5.1 CPU, RDKit 2025.3.5, NumPy 2.3.5, SciPy 1.17.0, scikit-learn 1.8.0, joblib 1.5.3 and matplotlib 3.10.8. Ki package versions are in the frozen `SAR_RESOURCE_KI_v1.json` supplied by v1.0.0 and restored by the map. P0 and the identity precheck additionally import pandas; its exact historical version is not recorded in the supplied freeze. This is a specific environment limitation, not an assertion of a locked complete environment.

The IC50 runner determines its root from its file location, creates `runs/sar_transport_v1`, and verifies its freeze records. Calibration precedes RF/neural fits; the frozen seeds are 17 and 43. The original commands are:

```text
python scripts/run_sar_transport.py --stage calibrate
python scripts/run_sar_transport.py --stage rf --seed 17
python scripts/run_sar_transport.py --stage rf --seed 43
python scripts/run_sar_transport.py --stage neural --model fscap --seed 17
```

The other recorded neural family names are `actfound`, `no_condition` and `conditional`, each with seeds 17 and 43. These are reproduction instructions for a researcher choosing to refit, not operations performed for this supplement. Historical result/CPU/stop guards are preserved and will need a separately documented output workspace for a deliberate rerun; never overwrite the archived results or describe exposed tasks as new independent validation. The Linux/POSIX resource monitor was not validated on Windows.

For raw IC50 input retrieval, `UPSTREAM_INPUTS.tsv` gives the exact destination, pinned URL and SHA256. The original qualification code verifies these bytes. No large upstream source files are duplicated here. Ki raw file paths and Git blob hashes are listed in the restored `reports/sar_resource_v1/KI_SOURCE_EXPECTED.json`; its runner's `fetch` stage uses the fixed Nelen repository commit. Source metadata snapshots already in v1.0.0 should be used for provenance instead of substituting current ChEMBL API responses.

Two original planning scripts retain their original filesystem assumptions: `provenance_precheck.py` sets `ROOT` to its historical scratch root; P0 `run_pilot.py` expects a sibling `sar_review_20261104/bundle` directory. The frozen outputs are supplied and are the inputs used by the Ki qualification/export. A fresh execution of either planning script requires setting up those layouts or recording a local path-only change. Neither script is advertised as a portable one-command raw-data rebuild. The precheck downloads live metadata, so rerunning it against today's API would not by itself reproduce the frozen metadata snapshot.

## Licensing and scope

The scope of v1.0.0 licenses remains unchanged: project-original code is MIT and ChEMBL-derived data retain CC BY-SA 3.0 attribution/share-alike terms. See the v1.0.0 `LICENSES/` directory. This supplement also contains the **GPL-3.0 ActFound backbone**, with its full original `LICENSE` beside the source. The MIT notice for project-original files does not relicense that third-party code; use involving the combined ActFound code must preserve the upstream GPL terms. No third-party pretrained weights are supplied. FS-CAP is represented by the recorded project task adaptation; unused upstream FS-CAP source files are not redistributed here.

The archive excludes manuscript prose, cover letters, author-approval documents, chats, C1/C2/R0/parallel-v2 experiments and unrelated private history. The stopped conditional-model/KG development, `RULE_UNAVAILABLE`, strict-scaffold budget block and previously exposed data remain unchanged. `STATIC_CHECKS.json` records syntax/hash checks only. No runtime or scientific verification is inferred from those checks.
