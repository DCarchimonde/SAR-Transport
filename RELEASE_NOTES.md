# SAR-Transport v1.0.0

First versioned code/data resource distribution. Project-original code: MIT,
Copyright (c) 2026 Siyuan Tong. ChEMBL-derived data and numerical results:
CC BY-SA 3.0; upstream notices preserved.

Assets: code archive, separate task-data archive, installable wheel, full
RELEASE_MANIFEST.json and SHA256SUMS.txt. The data archive includes task
views, truth in separate files, saved predictions, acquisition traces,
coverage, eligibility and exposure records, and complete scientific results.
The code archive includes the evaluator, common-C4 runner, schemas, tests,
examples, source provenance and reproduction/limitations documentation.

934 IC50 fixed-C8 task definitions preserve 686/122/126 original partitions;
15 saved-prediction tables cover the 126-task evaluation partition only.
IC50 acquisition has 64 P0 development tasks. Each Ki track has 64 tasks,
50 groups and 39 targets. Formal acquisition pays the same saved four
references before acquiring at most four new references, for total cost eight.

K1/K2 established neither advantage nor equivalence. The exposed IC50 21.80%
result is not generalized to Ki. KG and conditional-neural stopping decisions,
RULE_UNAVAILABLE and the unfitted/unscored strict-scaffold budget block remain.
No new fit, inference, matrix scoring or resampling was run for this release.

The wheel and supplied code are not a complete historical neural/RF training
pipeline. Linux was checked in RC2; Windows/macOS are untested. No paper
manuscript, private author-authorization record or unrelated history is included.
No DOI or community-use statistics are claimed.
