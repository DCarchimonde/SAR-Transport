# SAR-Transport reproducibility supplement repro1

This supplement accompanies [SAR-Transport v1.0.0](https://github.com/DCarchimonde/SAR-Transport/releases/tag/v1.0.0). It provides the historical training, qualification and export implementations, frozen configurations, environment records and prepared inputs for the reported reference results.

The scientific result commit is `d5d3227b8a0269e2a896056ba41abf7d502d768f`; the resource RC2 commit is `a21719a84edaab3176b65db1d8151ca335ac3a31`. `SOURCE_IDENTITY.tsv` records source identities and SHA-256 values. The Ki implementation is fixed at `6a2b2b1c6f3d580eb0c26db34d3608cda038b2b2`.

## Obtain and restore files

1. Download the code and task-data archives from v1.0.0 and extract them into one resource directory.
2. Download the entire `supplements/repro1` directory, including `TRAINING_INPUTS.zip`, and verify `SHA256SUMS.txt`.
3. Use `REPRODUCIBILITY_ASSET_MAP.tsv` to locate the implementation, configuration and inputs for each result. `RESTORE_MAP.tsv` locates dependencies supplied in v1.0.0.
4. Create the historical directory layout with:

```text
python stage_layout.py --resource-root /path/to/v1-resource --destination /path/to/new-historical-layout
```

`stage_layout.py` copies files into the requested layout. `TRAINING_INPUTS.zip` contains the original `dataset.pkl`, `preprocessors.pkl`, identity universe and qualified/conflict observation tables. It expands to approximately 518 MB. The pickle files require the recorded Python/scikit-learn environment; the public JSON/TSV task API can be used independently.

## Choose a workflow

| Workflow | Inputs and implementation |
| --- | --- |
| Load tasks and score saved predictions | Use the v1.0.0 evaluator, task tables, labels and saved predictions. Historical model dependencies are unnecessary for this workflow. |
| Inspect or refit neural/RF models | Use `run_sar_transport.py`, `sar_transport_methods.py`, preprocessing code, frozen protocols, seeds, prepared inputs and the ActFound backbone. Refitting requires the historical environment and a separate output workspace. |
| Reconstruct task derivation | Use the IC50 qualification/features/export scripts, Ki identity precheck, frozen metadata/roles, numerical-qualification code and input maps. Retrieve raw activities from their pinned sources. |
| Inspect or rerun reference policies | Use the P0 implementation and inputs in v1.0.0, and the frozen Ki runner and core/baseline modules identified by `RESTORE_MAP.tsv`. |

The distribution supplies prediction tables and refitting code/inputs, but not the original trained neural/RF checkpoints. Task loading and saved-prediction scoring were tested; refitting from this assembled supplement has not been validated.

## Environments and commands

The original IC50 environment is recorded in `historical/reports/sar_transport_v1/EXECUTION_ENVIRONMENT.json`: Python 3.12.14, PyTorch 2.5.1 CPU, RDKit 2025.3.5, NumPy 2.3.5, SciPy 1.17.0, scikit-learn 1.8.0, joblib 1.5.3 and matplotlib 3.10.8. Ki package versions are recorded in `SAR_RESOURCE_KI_v1.json`, supplied by v1.0.0 and restored by the map. P0 and the identity precheck also require pandas; its historical version was not recorded.

The IC50 runner resolves its root from its file location, creates `runs/sar_transport_v1`, and verifies freeze records. Run calibration before RF/neural fitting:

```text
python scripts/run_sar_transport.py --stage calibrate
python scripts/run_sar_transport.py --stage rf --seed 17
python scripts/run_sar_transport.py --stage rf --seed 43
python scripts/run_sar_transport.py --stage neural --model fscap --seed 17
```

Other neural family names are `actfound`, `no_condition` and `conditional`, each with seeds 17 and 43. Historical result, CPU and stop guards remain active. Use a separate output workspace for reruns. The resource monitor uses Linux/POSIX interfaces.

`UPSTREAM_INPUTS.tsv` specifies raw IC50 destinations, pinned URLs and SHA-256 values. The original qualification code verifies these bytes. Ki raw paths and Git blob hashes are in the restored `reports/sar_resource_v1/KI_SOURCE_EXPECTED.json`; the runner's `fetch` stage uses the pinned Nelen commit. Use the supplied metadata snapshots to reproduce the reported provenance.

Two scripts require historical layouts: `provenance_precheck.py` sets `ROOT` to its original scratch directory, and P0 `run_pilot.py` expects a sibling `sar_review_20261104/bundle` directory. Recreate these layouts or record a local path adjustment before execution. The precheck fetches live metadata, so use the supplied frozen outputs when reconstructing the reported Ki qualification.

## Licensing

Project-original code is MIT; ChEMBL-derived data retain CC BY-SA 3.0 attribution/share-alike terms. See v1.0.0 `LICENSES/`. The ActFound backbone is GPL-3.0 and includes its original license beside the source. Preserve those GPL terms when using the combined ActFound code. FS-CAP is supplied as the recorded task adaptation.

## Verification records

`MANIFEST.json` and `SHA256SUMS.txt` cover the distributed files. `STATIC_CHECKS.json` records syntax and source/input-hash checks; runtime coverage is described under the workflows above. Scientific result tables and historical source files retain their frozen identities.
