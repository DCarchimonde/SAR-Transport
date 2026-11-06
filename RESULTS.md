# Preserved scientific results and limits

No new fitting, model inference, matrix scoring or resampling was performed
for the version 1.0.0 release. The underlying scientific result is project
commit d5d3227b8a0269e2a896056ba41abf7d502d768f, following Ki implementation
freeze 6a2b2b1c6f3d580eb0c26db34d3608cda038b2b2.

| Frozen Ki primary contrast | Difference in choice loss | Nominal 95% family-wise interval |
|---|---:|---|
| K1: Tanimoto minus source-potency RBF residual calibration, fixed C8 | -0.0048797672 | [-0.0173834356, +0.0087088078] |
| K2: max-min minus random-seed metric mean, acquisition | +0.0020477307 | [-0.0047051949, +0.0091859601] |

Neither contrast established an advantage or equivalence. Ki fixed-track
Tanimoto MAE was lower than source-potency MAE, while choice loss was not
lower; these are descriptive outcomes, not extra formal comparisons.

The old 21.80% IC50 improvement is limited to its already exposed evaluation
tasks. A significant IC50 comparison and a nonsignificant Ki comparison do
not establish a difference between endpoint effects. IC50 acquisition P0 is
development evidence: KG was about 10.60% worse than max-min, and the KG
upgrade remains stopped. Conditional-neural superiority remains stopped.
RULE_UNAVAILABLE is preserved. Strict-scaffold eligibility passed, but the
frozen operation-budget preflight blocked fitting and scoring. No new-scaffold
result exists in this release. Earlier project failures and exposure chronology
are retained in the scoped scientific ledgers, without publishing private
communications or unrelated project archives.

Complete saved results are in `resource_data/IC50/reference_results/`,
`resource_data/IC50/P0_development_evidence/`, and
`resource_data/Ki/reference_results/`. Query truth, pool truth, predictions,
initial-reference contracts and exposure descriptions remain separate files.
`COVERAGE.tsv` identifies definitions, predictions and trajectories separately.
