# Implemented predictors and comparison scope

This is a new task adaptation, trained from scratch. It is not a reproduction of the published FS-CAP or ActFound benchmark. The fixed activity source and all task splits are shared by the methods. Source-assay potency is available for each query; target-assay potency is limited to the same eight references in every task.

## Common inputs and calibration

The affine map uses development-standardized potency, an unpenalized intercept and a slope penalized toward one. Negative slopes are allowed. Its regularization is selected once on validation group-equal selection loss. The same chosen map supplies residuals to the potency-RBF, Tanimoto and deep kernels, and supplies the RF baseline/context features. FS-CAP predicts a residual over this map. ActFound retains its native support-anchored relative prediction rather than replacing that prediction with an ordinary regression MLP.

Development-fitted TF-IDF/SVD encodes the original assay method descriptions. The condition vector also includes available categorical fields and missingness indicators. It does not invent temperature, incubation time, substrate or enzyme preparation. Common target-name encoding is available to every learning method. Deep kernels additionally share the same development target vocabulary and UNK handling. Assay IDs, document IDs and DOI are grouping keys, never predictor features.

## FS-CAP adaptation

The source is Rose-STL-Lab/FS-CAP at `5eb9df1542e465c7cd9b8a0f81a6a521c5352a8b`. The adapter retains separate context/query encoders, target-value multiplication of context fingerprints, mean context aggregation and a joint prediction head. Source potency, target name and assay descriptors are appended to both branches; target context values are also explicit context inputs. Encoder hidden width is 512, with five hidden layers and 128-dimensional embeddings. The predictor widths are 512/512/512/256/256/1. These dimensions were frozen before the valid fits.

One seed-17 attempt normalized features separately inside each task. That inadvertently removed task-constant metadata and did not match the official mixed-task batch normalization. The attempt is retained as an engineering failure. The correction concatenates the eight tasks before normalization, then aggregates contexts and computes a task-equal loss. It did not change the data, model widths, objective, learning rate, seed or stopping rule. This consumes one of the two permitted engineering restarts. See `FSCAP_ENGINEERING_RESTART.json` and the archived implementation.

## ActFound adaptation

The source is BFeng14/ActFound at `c96f0559ea7a2504e797b9d65b9f1af1aa899ddd`. The adapter imports the official two-layer, width-2048 fingerprint backbone, including its per-step batch normalization and leaky-ReLU. It preserves the support-relative prediction, original-Morgan Tanimoto neighbor selection, learned-embedding similarity weights, pairwise metaobjective, first-order/second-order schedule and five final-weight-only support updates. The initial per-step learning rates are 0.001 and are learned on development tasks. Query target values are passed only to the development metaobjective. Validation and final prediction use a separate interface without query target values. Each task starts from the saved initialization and restores its batch-statistic state afterward.

KNN-MAML, transfer mixing, inverse ensembles, repeated random supports and public pretrained weights are disabled. The wrapper is device-neutral; it does not invoke the original GPU-specific batch-statistic restoration. Its explicit restoration copies the same nontrainable state tensors.

There is an information-use limitation to this literal concatenation adapter: original per-task batch normalization centers away additive task-constant descriptors in the first layer. Source potency still varies within a task and remains available. Thus this implementation retains the original architecture and the requested appended inputs, but it is not evidence about every possible condition-aware redesign of ActFound. RF, FS-CAP and the full deep kernel have effective task-condition pathways; the direct I-versus-H comparison isolates the tested deep-kernel condition pathway. Any claim about superiority over all context-aware ActFound variants would exceed this experiment.

## Deep kernels

Both kernels use the same directed message-passing encoder (depth three, hidden size 128), source potency, common target information, a nonnegative learned Tanimoto coefficient and a 64-dimensional inner-product feature map. The full model uses a small assay-condition modulation layer. The ablation replaces assay-condition input with a learned task-invariant vector while retaining target information and the same nominal module layout. It therefore removes assay-dependent modulation without removing the common target information. The invariant vector has no assay identity lookup.

The pointwise target MSE is task-equal. The training kernel ridge coefficient is one. The four allowed inference ridge values are compared on validation predictions from the same saved weights; those solves are counted as calibration fits, not extra global neural training. The eight-reference solve is differentiable during training. No classification, cycle, correlation or pair-consistency objective is added.

## Evaluation scope

Validation checkpoint selection is task-equal absolute MAE. Final comparator and calibration selection use group-equal selection loss. These are distinct operations and both are recorded. Seeds are reported separately and averaged at the metric level. Fixed-model group bootstrap intervals do not include the uncertainty of retraining or validation selection.

The held-out units are connected assay/compound/document groups within one curated public resource. They are not prospective experiments, separate databases or necessarily new chemical scaffolds. Original source selection based partly on cross-assay activity differences prevents interpreting observed reversal rates as unbiased rates for all of ChEMBL.
