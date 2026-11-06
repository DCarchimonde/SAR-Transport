# Source and immutable evidence map

The only activity source is `Jnelen/ChEMBL_MatchedPairsAnalysis@0a77cdd28abf40a1562aadce727d032cd668815c`:
https://github.com/Jnelen/ChEMBL_MatchedPairsAnalysis/tree/0a77cdd28abf40a1562aadce727d032cd668815c

IC50_MinimalCuration and Ki_MinimalCuration_Unpruned activities derive from ChEMBL32 (DOI 10.6019/CHEMBL.database.32). Assay/document metadata were supplemented from ChEMBL37; these are not claimed as one identical database snapshot. Qualification retained finite exact measurements under the upstream equality/nM/source SQL chain. Multiple values at the same normalized key were excluded without averaging. Ki retains 85 excluded conflicting keys and the original whole-group overlap removals. Upstream activity-based exclusions are retained as selection biases, not a database-wide rate estimate.

* `resource_data/IC50/provenance/`: task/group/alias/assay maps, original C8, eligible/query edge lists and observation eligibility.
* `resource_data/Ki/provenance/`: source and alias maps, precheck groups and frozen sentinel/pool roles. No component is split after numerical exclusions.
* `resource_data/Ki/reference_results/`: immutable qualification/score/prediction/trajectory artefacts.
* `resource_data/IC50/P0_development_evidence/`: original development protocol, scripts, designs, trajectories and supplementary timing.
* `resource_data/contracts/`: new wrapper-only saved initial-reference contracts, with input hashes and unchanged query signatures.
* `SCIENTIFIC_ASSET_INVARIANCE.json` in the resource root: old predictions, scores, intervals, query files and weights remain byte-identical; new sidecars and entry-point overlays are listed separately.

Old SAR result: 0b22c0af92a5b9369da987e3e850a275e56949cf. Ki implementation freeze: 6a2b2b1c6f3d580eb0c26db34d3608cda038b2b2. Scientific resource result: d5d3227b8a0269e2a896056ba41abf7d502d768f. Original delivery: df4dd25a2f71dd656ba6d760defe61579df6477a. These are private provenance IDs, not public download endpoints. No reader needs private access to load the supplied exports.

Current exposure is carried by the default manifest. Historical manifests remain verbatim to preserve their temporal record; their former 'not scored yet' wording is not the current state. Original reports/archives are retained in the private project and are not overwritten by RC2.
