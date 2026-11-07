# Reproduction scope

| Level | Supplied assets and command | Current support |
|---|---|---|
| Load task inputs | Code wheel; `resource_data/CURRENT_TASK_MANIFEST.tsv`; JSON, structures, source values and metadata | Self-contained; no private repository needed. |
| Rebuild exported task views | `sar_transport_resource.rebuild.rebuild_exported_task(data_root, endpoint, task_id)` uses compressed TSVs and preserved manifests | Recreates task data, C8 order, query graph and weights from supplied tables; does not rerun curation. Fixed and acquisition example identity tested. |
| Score saved results | `score_predictions` or `score-saved`; saved predictions and separate truth | Scoring implementation is supplied. RC2 only scores one real task per track as software regression; full matrix command was not run. |
| Formal common4+4 integration | C4 sidecars, `FormalAcquisitionRunner`, example and schema | Pays saved C4 before callback; four remaining measurements; signature/compliance output. |
| Raw upstream eligibility/task derivation | Upstream fixed source, original SQL and identity/eligibility/group/role ledgers; original frozen project scripts | Source chain is documented, but the standalone distribution does not claim an independent end-to-end raw-ChEMBL recuration. Required historical runtime/input projections are not all bundled. |
| Reproduce historical neural/RF fitting | Original frozen private project contains training implementation/configurations/weights | Not bundled in this portable distribution. Saved outputs permit evaluation, not a claim that installing the wheel reproduces training. No pretrained weights or new inference included. |

The wheel only installs code. Data are delivered separately as `SAR_Transport_Data_v1.0.0.zip` (extract to `resource_data/`). The optional `SAR_Transport_Code_v1.0.0.zip` includes wheel, source, examples and tests (extract to `code/`). The combined resource distribution contains both. These are filenames, not claims of a public host or DOI.

### One task from portable tables

```python
from sar_transport_resource.rebuild import rebuild_exported_task
view = rebuild_exported_task('resource_data', 'IC50',
                            'fixed_c8_v1::CHEMBL1000620__CHEMBL1000621')
```

The table reconstruction reuses frozen identities and lists. Fixed C8 order comes from REFERENCE_LABELS, not sorted compound order. Qualification, fingerprints, max-min design and modelling are not run. `OBSERVABLE_COMPOUNDS` supplies only source values and structures; target query and pool truths are separate tables and not read by this routine.

### Full saved-method scoring (available, not executed in RC2)

```bash
python -m sar_transport_resource score-saved --data-root resource_data --manifest resource_data/IC50/TASK_MANIFEST.tsv --predictions-tsv resource_data/IC50/saved_predictions/tanimoto.tsv --split evaluation --out replayed_saved_tanimoto.json
```

The explicit historical manifest above supplies unchanged evaluation membership; users should inspect the current catalogue first. Both seed metrics must be calculated separately before averaging; do not silently ensemble predictions. This capability is not a fresh independent validation. All underlying outputs are exposed.

Historical P0 scripts and readmes are preserved as evidence, not advertised as standalone execution entry points. Their old project paths are not requirements for task loading or scoring. No third-party user or Windows validation is asserted.


## Reproducibility supplement (7 November 2026)

The independently versioned [repro1 supplement](supplements/repro1/README.md) adds the original historical training/qualification implementations, frozen configurations, prepared training inputs and a result-to-asset map. It documents the precise supported reproduction levels and the absence of original trained checkpoints. The v1.0.0 tag and release assets remain unchanged. No new scientific computations were run.
