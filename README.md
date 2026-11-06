# SAR-Transport 1.0.0

Resource repository: https://github.com/DCarchimonde/SAR-Transport

The repository provides unpacked source code and frozen code/data archives. Task data are supplied in the separate data archive. Download both versioned archives and extract them into the same directory before running the commands below. A wheel, a per-file manifest and SHA256SUMS identify the distribution.

SAR-Transport supplies paired-assay **tasks**, fixed analogue queries, permitted source measurements, target-reference budgets and dependence-group evaluation. It supports adding a predictor or a reference-selection callback without deriving these experimental-design decisions from upstream tables. The contribution is a derived task resource. The results do not establish neural superiority.

## Exact coverage

| Endpoint / track | Task definitions | Dependence groups | Saved prediction / trajectory coverage |
|---|---:|---:|---|
| IC50 fixed C8, all original partitions | 934: 686 development, 122 validation, 126 evaluation | 573 total | Reference predictions only for the 126 evaluation tasks / 92 groups. Each of 15 A–I/isotonic tables has 2,898 query rows. |
| IC50 acquisition P0 | 64 development tasks | 64 | Nine frozen policies per task; development evidence only. |
| Ki fixed C8 | 64 | 50 | A–D saved predictions; 39 targets, 7,873 query edges. |
| Ki acquisition | 64 | 50 | Nine saved policies per task; 39 targets, 3,369 query edges. |

`COVERAGE.tsv` separates task definitions, prediction tables, support-adaptation records and trajectories. The 252 ActFound task/seed adaptation records are not a 16th query-prediction table. Tracks overlap in groups and have different queries: do not pool them or sum their group counts as independent evidence.

## Install code, unpack data, and run

The wheel contains the Python core, **not the data**. The full resource distribution contains both `code/` and `resource_data/`. Alternatively, unpack `SAR_Transport_Code_v1.0.0.zip` and `SAR_Transport_Data_v1.0.0.zip` into the same working directory. Commands below run from the resource root and contain no personal or private-download paths.

```bash
python -m pip install --no-index code/dist/sar_transport_resource-1.0.0-py3-none-any.whl
python -m sar_transport_resource list --data-root resource_data --endpoint IC50 --split evaluation
python -m sar_transport_resource inspect --data-root resource_data --endpoint IC50 --task-id fixed_c8_v1::CHEMBL1000620__CHEMBL1000621
python -m unittest discover -s code/tests -v
```

The default catalogue is `resource_data/CURRENT_TASK_MANIFEST.tsv`; its `exposure_status` and `exposure_status_current` state current exposure. Original manifests and their former wording remain unchanged in `resource_data/{IC50,Ki}/TASK_MANIFEST.tsv` and the historical tables. No data recover blindness by being repackaged.

The standard-library core needs Python >=3.10. Optional baseline adapters need NumPy/SciPy/RDKit; loading, billing, task reconstruction and saved-prediction scoring do not. Clean-environment validation used Linux Python 3.12.14. Windows/macOS have not been tested.

## Add a predictor

```python
from sar_transport_resource import load_task_entry
entry = load_task_entry('resource_data', 'IC50',
                       'fixed_c8_v1::CHEMBL1000620__CHEMBL1000621')
print(entry.exposure_status)
view = entry.view
# Your predictor receives observable source values, structures and C8 labels.
# Return {query_id: finite_target_prediction}; query truth is supplied separately
# to score_predictions by the evaluator after predictions have been saved.
```

One regression command checks existing predictions on exactly one fixed and one acquisition task:

```bash
python code/examples/rc2_acceptance.py --data-root resource_data --out-dir my_software_check
```

Use a fresh output directory. This reads saved predictions; it does not train or infer a model, and its two example scores are not new scientific comparisons.

## Add a formal acquisition policy

`entry.initial_reference_contract` contains the exact saved `initial_reference_ids`, not a newly selected start. `FormalAcquisitionRunner` pays these four in the recorded order **before** invoking your callback. Your callback receives the four paid labels and can buy only four further distinct pool members. The total budget is eight. Repeats are free; query IDs, unknown IDs and requests beyond eight are refused and logged. An under-spent trajectory is marked `protocol_compliant: false`.

```python
from sar_transport_resource import FormalAcquisitionRunner, load_task_entry
entry = load_task_entry('resource_data', 'IC50',
                       'acquisition_sentinel_v1::CHEMBL1013099__CHEMBL1013112')
# A trusted evaluation driver supplies pool_truth_provider, not query labels.
# callback(paid_view, oracle) selects the remaining four and returns predictions.
# runner = FormalAcquisitionRunner(entry.view, entry.initial_reference_contract,
#                                 pool_truth_provider, 'new_trace.jsonl')
# result = runner.run(callback)
```

The runnable saved-policy example exposes the complete callback pattern:

```bash
python code/examples/formal_saved_policy.py --data-root resource_data --trace my_trace.jsonl --out my_formal_example.json
```

It replays a saved max-min path and returns saved predictions, without refitting. Output includes protocol/query signatures, measurement order, cost and compliance. There is no scientific advantage estimate from this software check.

The general `BudgetOracle` remains available for custom experiments. `code/examples/custom_policy.py` is the retained **non-formal** arbitrary-start demonstration; it explicitly reports `benchmark_comparable: false`. Such scores cannot be labelled as formal paper-track results. Signatures are cooperative research safeguards, not a security barrier against deliberately opening local truth files.

## Interpretation and provenance

Source/target order follows assay IDs, not measured chronology. Identity handling preserves stereochemistry, while Morgan fingerprints use radius 2, 2048 bits and `includeChirality=False`. Analogue edges are the upstream RASCAL edges, not universally single-cut edits. Acquisition measures the indirect value of bridging references for fixed unmeasured choices; it is not unrestricted lead optimization.

K1 and K2 established neither an advantage nor equivalence. The old 21.80% IC50 result is confined to its exposed tasks. A significant IC50 comparison and a non-significant Ki comparison do not establish a difference between endpoints. P0/KG and neural negative results remain intact. `REPRODUCTION_SCOPE.md`, `DATA_CARD.md` and `LICENSES/` specify what can be reproduced and redistributed. No public-development duration, community uptake, Windows validation or new scientific advantage is asserted.


## Distribution, citation and licences

Code and data are separate archives. Download both versioned assets, extract
them in one directory, and follow the commands above. `RELEASE_MANIFEST.json`
and `SHA256SUMS.txt` identify the exact files. The data archive also contains
the saved predictions, acquisition traces, qualification and exposure records,
and complete machine-readable scientific results. See `RESULTS.md`.

Project-original code is MIT, Copyright (c) 2026 Siyuan Tong; ChEMBL-derived
data and numerical results retain CC BY-SA 3.0. The upstream code notice is
preserved. See `LICENSE.md` and `LICENSES/` for the separate scopes.

Creators: Siyuan Tong (Independent Researcher, Wuhan, China) and Yuxuan Wang
(Faculty of Engineering, The University of New South Wales, Sydney, Australia).
Correspondence: z5653501@ad.unsw.edu.au. This research received no external
funding. The authors declare no competing interests. `CITATION.cff` provides
resource citation metadata. This resource distribution contains no paper
manuscript, author-authorization record, or complete private research history.
