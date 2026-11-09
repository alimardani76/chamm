# ChaMM scientific review and finite completion plan

## Verdict

**The production run is credible and the current restricted positive verdict is supported relative to its tested controls. The scientific closeout still needs a small set of substantive robustness analyses. No demonstrated fatal target leakage or arithmetic error was found in the inspected implementation.**

Do not rebuild the production graph or full descriptor table. The next run should use the completed datasets and address conditional calibration, conservative validation groups, the nonmetal graph definition, and interpretation of disagreement magnitude.

## What was checked directly

The review inspected the scientific findings, methods, connectivity decision, model/calibration/analysis/sensitivity code, numerical inputs, and representative rendered figures.

- Reaggregated **825,210 structure-loss rows**, reproducing **90 all/metal/nonmetal model-summary rows**, including atom-weighted and structure-balanced MAE, denominators and improvements. Maximum MAE discrepancy was below 1e-15 e.
- Recomputed the three main 1,000-resample QMOF-level intervals from the supplied per-structure losses. They reproduce the reported improvements.
- Reaggregated **495,882 structure-level chemical-contrast rows**, reproducing the means and denominators of **108 contrast summaries**.
- Checked all **46,994 metal predictions for each of three pairs**: unique atom keys, correct fold assignments, and observed differences matching the recovered source charges. Recomputed the metal-domain main improvements and their intervals.
- Independently checked the saved **34-structure / 5,198-atom** accepted graphs from the raw lattice and fractional coordinates: **12,258 directed edges**, correct reverse edges, no duplicate physical edge records, and correct distances. Four extra indexed-neighbor entries are distinct periodic images.
- Confirmed that the displayed five-fold split keeps each QMOF ID together. Code inspection shows charge controls fitted on outer-training structures, training-only preprocessing, no outer-test early stopping, and no charge-derived quantities in the structural feature list.

Main held-out structure-balanced results reproduced from saved losses:

| Pair | Piecewise baseline MAE, e | Structural two-hop MAE, e | Improvement, e |
|---|---:|---:|---:|
| DDEC6–CM5 | 0.026285 | 0.010917 | 0.015368 |
| DDEC6–Bader | 0.092968 | 0.041689 | 0.051279 |
| CM5–Bader | 0.074354 | 0.043808 | 0.030547 |

The metal-only follow-up improvements also reproduce: approximately 0.009639, 0.006397 and 0.010307 e against its element-residual controls.

These checks do not regenerate the full global models: full atom-level global predictions, fitted models, training sample keys and the 240 MB feature table are omitted from the compact review ZIP. Global numerical verification therefore uses the supplied structure-level loss records. Metal numerical verification additionally uses the saved atom-level predictions and independent recovered charge table. Graph bookkeeping does not establish independent chemical bond truth.

## Four substantive closeout tasks

### 1. Test added chemical information against a stronger conditional charge-only model

The existing 12-bin per-element calibration is a useful bounded control. It does not by itself establish that remaining predictability requires environment information rather than a more flexible mapping of the input charge. The current structural model also adds central-element encodings, so the contribution of local neighbors is not isolated by the baseline-to-first-shell comparison alone.

Run a small, fixed, nested comparison on identical folds, training samples and weights: a center-only residual correction; a flexible correction using central identity and q_A only; that same model with first-shell descriptors; and that same model with two-hop descriptors. Predict the same signed difference with the same piecewise offset, so the charge-only and enriched models have a directly comparable target. The last comparison adds structural information conditional on q_A. Keep these explanatory conditional models separate from the existing charge-blind structural models.

Recompute the key signed and absolute chemical contrasts after the stronger charge-only correction. Do not claim that any finite control exhausts every possible charge transformation. If gains shrink or disappear, report that outcome rather than searching for a winning alternative.

### 2. Use metadata that actually exists for a stronger validation challenge

The methods currently say verified family/topology grouping is unavailable. A fresh inspection of `qmof_database/qmof.json` in the already pinned archive found verified source fields that the production run did not use:

| Source field | Available cohort records | Unique values | Records in groups spanning current folds |
|---|---:|---:|---:|
| Exact MOFkey | 5,408 | 5,114 | 441 |
| Exact MOFid | 5,276 | 5,082 | 287 |
| Topology field | 5,484 | 179 | 5,405 |
| Source DOI | 7,964 | 3,460 | 4,941 |
| Reduced formula | 9,171 | 8,468 | 1,030 |

This is **not evidence that the reported record-held-out evaluation is invalid**. That limited evaluation was disclosed correctly. It does mean that a conservative blocking check is feasible without inventing families or collecting new data.

One bounded proposal is to form connected components linking records with the same non-null exact MOFkey or the same exact source reduced formula. On the present cohort this produces **8,347 blocks**, largest size **13**. There are **1,414 records in nonsingleton blocks**, with **1,189 records in blocks that cross the original folds**. Formula equality is deliberately conservative blocking, not a claim that structures are identical. Exact MOFkeys also do not establish all broader chemical-family relationships.

Run the compact conditional comparison under these blocked folds, report overall and repeated-block results, and use uncertainty resampling consistent with those blocks. Correct the metadata-availability statement. Do not call this topology-held-out or universal family transfer. A topology-extrapolation campaign is outside this closeout.

### 3. Test the graph component on which the leading ligand claim depends

The existing +0.10 Å taper changes metal-centered CrystalNN selection. The original-cutoff stress graph uses the **same nonmetal cutoff** as the accepted graph. Stability of carbon contrasts under those alternatives consequently does not test the nonmetal bond rule itself; the current report acknowledges that limitation.

On the existing 96-structure subset, retain accepted metal-involving edges and perturb only the nonmetal radius multiplier to 1.20 and 1.30. Recompute the affected features and matched chemical contrasts. Fixed-model inference can test prediction sensitivity, but do not describe it as retrained alternative-model performance. Report both changes in environment membership and changes in the corresponding effects. Keep only physically interpretable alternatives; a wider cutoff that creates documented secondary contacts is a stress test, not an equally valid graph.

This needs no new whole-cohort neighbor computation or algorithm search.

### 4. Finish the chemical atlas and claim-to-evidence link

The present positive story focuses on **signed residual shifts**. These demonstrate environment-associated direction of charge-mapping departure. They do not automatically identify environments with greater absolute disagreement or greater variability. Absolute-residual contrasts are already saved, but they need a clear interpretation alongside the signed results.

Provide a compact supported-environment atlas with raw signed/absolute scheme differences, calibrated signed departures, calibrated absolute departures, support counts, structure-level uncertainty and graph-sensitivity status. Use the stronger conditional control where supported. Distinguish recurring directional shifts from high disagreement magnitude, and avoid naming signed shifts as generic increases in disagreement. Unsupported classes remain visibly unsupported; small effects are not promoted because their intervals exclude zero.

Update the main figures and findings from the final tables. The work can close with a restricted, positive, mixed or negative result; closure does not require a broader claim than the data support.

## Work that does not need repeating

No new DFT, experimental validation, deep learning, complete neighbor-method sweep, full descriptor rebuild, or expanded hyperparameter search is required for these within-cohort questions. Keep the current cohort, operational graph, production descriptors and original results as the frozen starting point. The metal-only follow-up remains a disclosed response to a global negative result; it should not be portrayed as independent confirmation.

## Final acceptance boundary

Scientific closeout is justified when the conditional-calibration comparison, conservative group challenge, nonmetal sensitivity and signed/magnitude atlas are completed and reflected honestly in the figures and claims. If a test weakens the story, narrow the conclusion and finish. Broader family/topology extrapolation, independent observables and causal electronic mechanisms remain outside the demonstrated scope, rather than becoming an endless list of mandatory new experiments.

The companion Codex prompt gives an executable assignment. Supporting JSON files retain the numerical checks and metadata counts. `verified_source_metadata.csv` is a review-derived join of the pinned public metadata to the 9,171 cohort IDs; no family or topology labels were invented.
