# Data card

## Source and derived objects

Activities originate only from Jnelen/ChEMBL_MatchedPairsAnalysis commit `0a77cdd28abf40a1562aadce727d032cd668815c`: IC50_MinimalCuration and Ki_MinimalCuration_Unpruned. These activities were extracted from ChEMBL32 (doi:10.6019/CHEMBL.database.32). Official metadata subsequently obtained from ChEMBL37 are separately identified and must not be called the original activity-version metadata.

The resource adds directed source/target tasks, information restrictions, reference/query roles, original analogue graphs with degree weights, connected dependency groups, fixed splits, per-task measurement accounting, saved predictions and uniform decision evaluation. It does not generate new biological measurements. Author pairs use RASCAL similarity (Johnson >=0.6); they are **analogue pairs**, not necessarily single-cut matched molecular transformations.

## Identity and qualification

Structures preserve stereochemistry, charge and tautomer representation; atom maps are removed. Only carbon-free counterions are removed when exactly one organic component remains. Ambiguous multiple organic components are excluded. Canonical identity uses the smallest raw compound ID within a normalized isomeric-SMILES group. Original and normalized IDs are retained. The planning identity compatibility check is inherited, not rerun.

The upstream SQL enforces endpoint Ki or IC50, nM, exact `=` relation, non-null pChEMBL and no data-validity comment. Cache exports omit some SQL qualification fields; the source code/commit chain documents them. The new Ki pass checked Decimal numerical identities between Final and same-key cache records. Multiple distinct values for a canonical (assay,compound) key exclude that key and its edges. There is no averaging, potency-direction filter or replacement of censored/missing observations with numbers.

Ki numeric qualification read labels after the implementation freeze and before prediction. It removed 85 conflicting canonical keys. The fixed-C8 and acquisition applications each retain 64 tasks / 50 groups / 39 targets. Exact counting units, exclusions, original 52-group identity-only candidate stage and later 50-group numerical eligibility are kept separately.

## Fields

| Field | Definition |
|---|---|
| endpoint / potency_unit | IC50/pIC50 or Ki/pKi; p=-log10(value in mol/L), higher means more potent |
| source_assay / target_assay | Lexicographic ordering of assay IDs; not a temporal transition claim |
| target_chembl_id | Shared protein target; incompatible mapped targets excluded |
| compound_id / canonical_smiles | Canonical ID and retained stereochemical structure; alias ledger supplied |
| source_potency | Observed source value available for every predictor-visible compound |
| reference_pool_ids | Permitted target measurements; disjoint from query IDs |
| reference_labels | Fixed C8 labels, or the acquisition oracle's paid labels |
| query_ids / query_edges | Unique involved compounds and original unmodified decision graph |
| weight | Normalized inverse-degree edge weight, not an independent-sample count |
| group_id | Connected shared-assay/identity/document/DOI component; no group is split after numeric exclusions |
| split | Original IC50 development/validation/evaluation, P0 development, or frozen Ki application |
| query_signature | Hash of ordered query IDs and weighted graph, used to reject changed tasks |

Every exported potency is finite and qualified as a point measurement. Excluded/unknown/bounded records remain in eligibility ledgers, not input values. Missing metadata remain missing; no inferred substrate, timing, laboratory or technology is filled in. Raw assay IDs/documents are provenance, not predictor features. A researcher may load allowable assay descriptions from the accompanying metadata while respecting original training splits.

## Two tracks

Fixed C8 uses the original deterministic max-min design (Morgan radius2, 2048 bits, includeChirality=False; first smallest ID, lexical ties). IC50 C8 and queries are copied exactly. Ki C8 is generated once on qualified identities. Remaining original edges with both ends outside C8 define queries.

Acquisition uses a preassigned pool and fixed sentinel query graph. The same four diverse starting references are paid for before four additions. Random17/43/71 means diverse4+random4. Previously fixed roles only lose ineligible compounds/edges; no reallocation occurs. Original qualified compounds that are neither queried nor measured can remain in provenance without being predictor inputs.

## Scoring and uncertainty

For each edge, choose the larger predicted potency, with absolute predicted-gap tolerance 1e-8 resolved by smaller canonical ID. Loss is the observed better potency minus the observed chosen potency. Edge weights are proportional to 1/degree(a)+1/degree(b), normalized per task. Tasks are averaged within group, then groups equally. MAE uses unique query points; difference-MSE uses the weighted edges. Random seed metrics are averaged, not treated as independent observations.

The Ki K1/K2 family uses 5000 paired resamples of whole groups, seed20261104 and per-comparison 1.25/98.75 percentile limits. Both tracks share union-group draw multiplicities. The intervals condition on frozen choices. They do not include historical strategy/model selection uncertainty.

## Exposure, bias and use limits

All old IC50 evaluation results are exposed. P0 is development exploration. Ki identities, structures, provenance and published aggregate conclusions were known before freezing; programmatic value qualification is explicitly logged. Query scores were calculated only after all task-local predictions and acquisition traces were frozen. Ki results are now exposed too. Cross-task query/reference roles can overlap inside a group; task states do not share target information.

Ki groups touching exposed IC50 assays, raw/canonical compounds, normalized structures, DOI or document were removed whole. This is a same-database, source/compound-separated cross-endpoint application, not prospective replication, an independent database, all-new targets or strict new-scaffold validation.

Upstream data selection removes identical cross-assay values and differences exactly3 logarithmic units, and imposes assay/overlap/structure restrictions. Consequently the task population is selected and cannot estimate an unbiased ChEMBL-wide reversal frequency. Graph edges reuse compounds and physical observations; neither edge counts nor assay-pair counts are independent replication counts. A measured label can contain reporting error or unmodeled assay effects. Observed reversals are not demonstrated biological mechanisms. There is no held-out CYP confirmation in this resource.

Strict-scaffold sensitivity retains the historical status: eligibility passed; the frozen operation-budget preflight blocked fitting/scoring. It was not rerun using the resource budget.


## RC2 coverage and formal sidecars

934 refers to fixed-C8 IC50 task definitions (686/122/126 development/validation/evaluation). Saved A–I/isotonic predictions cover only the 126-task, 92-group evaluation partition: 15 tables, each 2,898 query rows. IC50 acquisition is only the 64-task P0 development pilot. Each Ki track contains 64 tasks / 50 groups / 39 targets. See `COVERAGE.tsv`; predictions were not generated to fill development/validation coverage.

The default `CURRENT_TASK_MANIFEST.tsv` gives current exposed status. `initial_reference_contract` links to a backward-compatible sidecar in `resource_data/contracts/{endpoint}/`. Task JSON and its targets, queries and weights are unchanged. `initial_reference_ids` are copied in order from the original P0 `initial4_ids` or Ki trace `initial4` fields; `protocol_signature` hashes the complete contract (excluding that signature), and `task_definition_sha256` binds the original JSON bytes. The formal runner pays common C4 before the callback and accepts four additional unique references. Other initial sets are custom, non-comparable tasks.

Sentinel roles are assigned by sorting the identity-qualified canonical IDs on `(SHA256(('sentinel_v0|' + ID).encode('utf-8')).hexdigest(), ID)`. The first floor(n/2) are sentinels, the rest references. Only original edges with two sentinel endpoints are queries; their unique participating IDs become the final query list. After assignment, stored lists are canonical-ID sorted. Numerical qualification only deletes; it does not reallocate roles. Assay-ID order defines direction, not time order. This track quantifies indirect bridging-reference value on fixed unmeasured analogue choices.

Ki intervals are nominal 95% family-wise intervals, constructed from the original fixed paired-group bootstrap. No finite-sample exact coverage or equivalence is claimed. The endpoint contrast is not a formal test of differences between IC50 and Ki effects.
