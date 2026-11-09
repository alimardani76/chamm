# ChaMM final scientific findings

**Final verdict: Restricted positive: local chemical information survives tested charge conditioning and identifier/composition blocking; central-metal and magnitude claims remain qualified.** The finite closeout is complete. The production graph, cohort and descriptor table remain frozen. No new whole-cohort graph build, DFT, topology campaign or model search was performed.

## Conditional comparison: chemistry is not explained by the stronger scalar control

Every competitor uses the same fitting atoms, weights and outer-test atoms. The common offset is the original training-fitted 12-bin charge calibration. The center comparator is an element residual mean; the other three competitors use the same fixed tree with identity + q_A, then first shell, then first shell + two-hop paths. q_A appears only in this explicitly conditioned follow-up. No scheme is reference truth.

Structure-balanced held-out MAE in e, all **9,171 structures / 1,420,532 atoms**:

| split | pair | piecewise | center | charge_only | conditional_first | conditional_twohop | conditional_gain | ci_low | ci_high | folds_improved |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| original | ddec6_cm5 | 0.02628 | 0.02635 | 0.02443 | 0.01315 | 0.01022 | 0.01420 | 0.01407 | 0.01434 | 5 |
| original | ddec6_bader | 0.09297 | 0.09289 | 0.08457 | 0.04184 | 0.03889 | 0.04568 | 0.04511 | 0.04628 | 5 |
| original | cm5_bader | 0.07435 | 0.07445 | 0.06738 | 0.03969 | 0.03899 | 0.02838 | 0.02792 | 0.02886 | 5 |
| blocked | ddec6_cm5 | 0.02628 | 0.02635 | 0.02434 | 0.01309 | 0.01018 | 0.01416 | 0.01400 | 0.01432 | 5 |
| blocked | ddec6_bader | 0.09295 | 0.09286 | 0.08454 | 0.04197 | 0.03885 | 0.04569 | 0.04495 | 0.04640 | 5 |
| blocked | cm5_bader | 0.07436 | 0.07445 | 0.06750 | 0.03981 | 0.03914 | 0.02837 | 0.02782 | 0.02893 | 5 |

The main increment is conditional two-hop versus flexible charge-only. It remains positive in all five folds for every pair under both splits, with paired intervals above zero. Central-element residual recentering alone does not explain the result. All atom-weighted metrics, common-element support, unsupported calibrations, and per-fold increments are retained in `results/closeout_v1/atlas/`. The center-to-first comparison jointly adds q_A and local features; only charge-only-to-first isolates first-shell information conditional on q_A. Path increments are separately reported in `paired_increments.csv`; they are not the old charge-blind path increments.

## Conservative grouping challenge

The prior statement that grouping metadata was unavailable is corrected. The pinned `qmof_database/qmof.json` contains source reduced formulas and exact MOFkeys. All cohort IDs and those grouping fields match the supplied reviewer metadata. There are 5,408 non-null MOFkey records / 5,114 distinct keys, and 8,468 exact reduced formulas.

Linking equal non-null MOFkey **or** equal exact source reduced formula produces **8,347 components**, largest 13, with **1,414 nonsingleton records** and 1,189 records in components spanning the old folds. The deterministic new folds contain [1835, 1834, 1834, 1834, 1834] structures. This is identifier/composition-blocked validation; formula equality is not structural identity, and this does not prove universal family or topology separation.

Blocked conditional gain over charge-only, in e:

| pair | population | n_structures | n_components | gain | ci_low | ci_high |
| --- | --- | --- | --- | --- | --- | --- |
| ddec6_cm5 | all | 9171 | 8347 | 0.01416 | 0.01400 | 0.01432 |
| ddec6_cm5 | has_mofkey | 5408 | 4842 | 0.01441 | 0.01420 | 0.01462 |
| ddec6_cm5 | nonsingleton | 1414 | 590 | 0.01306 | 0.01252 | 0.01357 |
| ddec6_bader | all | 9171 | 8347 | 0.04569 | 0.04495 | 0.04640 |
| ddec6_bader | has_mofkey | 5408 | 4842 | 0.04370 | 0.04281 | 0.04473 |
| ddec6_bader | nonsingleton | 1414 | 590 | 0.04664 | 0.04386 | 0.04979 |
| cm5_bader | all | 9171 | 8347 | 0.02837 | 0.02782 | 0.02893 |
| cm5_bader | has_mofkey | 5408 | 4842 | 0.02646 | 0.02572 | 0.02725 |
| cm5_bader | nonsingleton | 1414 | 590 | 0.02687 | 0.02449 | 0.02924 |

The nonsingleton challenge population and MOFkey-present population both retain the gain; the singleton majority therefore does not conceal failure on those challenge sets. Blocked intervals resample components using sums of per-QMOF paired losses divided by resampled QMOF counts. They retain a structure-balanced estimand, rather than weighting every component mean equally.

## Metal conclusions are narrower

For the global conditional models, blocked metal results over **46,994 atoms / 9,165 structures** are:

| pair | charge_only | conditional_twohop | gain | ci_low | ci_high | folds_improved |
| --- | --- | --- | --- | --- | --- | --- |
| ddec6_cm5 | 0.02694 | 0.02202 | 0.00493 | 0.00451 | 0.00533 | 5 |
| ddec6_bader | 0.02044 | 0.02851 | -0.00807 | -0.00857 | -0.00749 | 0 |
| cm5_bader | 0.02959 | 0.02868 | 0.00091 | 0.00049 | 0.00137 | 5 |

DDEC6−Bader worsens consistently for central metals; the CM5−Bader gain is small. A broad three-pair central-metal benefit is not retained. A failure of this global population model does not prove absence of metal chemistry. The earlier all-metal follow-up remains separately disclosed in `results/metal_control/`; it was motivated by the previous global failure, was not an independent replication, and does not substitute for the present conditional test. No second domain-specific ladder was added in this closeout.

## Signed shifts, magnitude and spread are different results

The atlas reports raw signed differences, raw absolute differences, signed departures after flexible charge-only correction, absolute calibrated departures, and within-structure/environment residual SD. All are empirical quantities in e. Neither absolute departure nor residual spread is calibrated physical uncertainty or charge-scheme error.

- **Directional departure without greater calibrated magnitude:** O-attached versus strictly C/H-only carbon, DDEC6−Bader: signed contrast -0.07210 [-0.07491, -0.06946] e, but absolute-departure contrast -0.00112 [-0.00342, +0.00126] e, over 7,601 paired structures. The magnitude contrast is unsupported; the signed shift must not be called increased disagreement.
- **Greater absolute calibrated departure without a recurring signed departure:** N-attached/O-absent versus C/H-only carbon, CM5−Bader: signed contrast +0.00068 [-0.00378, +0.00540] e, while absolute departure increases by +0.10216 [+0.09904, +0.10558] e, over 6,689 paired structures. Opposite directional departures can cancel in their mean.
- **A weakened earlier headline:** N-bound versus C-bound H, DDEC6−CM5, changes from the earlier piecewise-residual estimate -0.06788 e to -0.02187 [-0.02266, -0.02099] e after the stronger charge-only control on the original folds. A residual contrast remains, but the former magnitude is not retained.

These examples also distinguish raw disagreement from calibrated departure (blocked, paired structure means, e):

| pair | contrast | n_structures | raw_signed | raw_absolute | calibrated_signed | calibrated_absolute | residual_sd |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ddec6_bader | C:O-attached vs C:C/H-only | 7601 | -0.77166 | 0.73594 | -0.07210 | -0.00112 | -0.03488 |
| cm5_bader | C:N-attached vs C:C/H-only | 6689 | -0.36672 | 0.34561 | 0.00068 | 0.10216 | 0.03716 |

The full signed/magnitude/spread contrasts, weak and null effects, ambiguity fractions, atom counts, structure counts and component support are in the atlas. A compact atlas uses at least 30 structures and 30 represented validation components; unsupported rare classes remain in the full table. No physical error tolerance was invented. The three pair differences are algebraically dependent and do not constitute three independent discoveries.

## Actual nonmetal-cutoff test

On the existing **96 structures / 16,165 atoms**, all accepted metal-involving edges were retained exactly, with no new metal neighbor calculation. Nonmetal cutoff 1.20 changes 0 first shells; 1.30 adds 12 unoriented physical edges, changes 18 first shells and 38 two-hop descriptors, but changes **0 headline environment memberships**.

The added contacts surround six H sites: H–O distances 1.2125–1.2200 Å, nearly linear near-midpoint O–H–O pairs. The accepted graph leaves those H sites isolated; the wider graph gives H degree two. Shared-proton/secondary-contact interpretation remains ambiguous, and 1.30 is a stress configuration, not equally established bond truth. These flagged proton-sharing environments do not support a universal proton-bonding claim.

Every headline contrast was recomputed using the **same stronger original-fold charge-only residuals**, on each graph's available pairs and on the intersection across all three graphs. The maximum common-pair effect change across all reported quantities is 0 e. Available/common atom and structure denominators and paired changes with QMOF intervals are saved. The unchanged effects establish stability only for these tested categories and settings on this subset. Some H/O/N contrasts have limited paired subset support; the evidence-status table labels that limitation. No retrained alternate-graph performance is claimed.

## Boundary of the final claim

The supported contribution is reproducible, interpretable **within-cohort local chemical structure in charge-definition disagreement conditional on the tested finite models**, robust to this identifier/composition challenge. It is strongest for common nonmetal environments and has no universal direction across schemes. Model capacity and global population weighting remain limitations, especially for metals. The conditional comparison tests these fixed regressors, not every conceivable scalar transformation.

Exact MOFkeys/formulas leave broader related-record dependence unresolved. Bootstrap intervals condition on fitted models and these validation units. The graph is an operational representation with documented ambiguity, not independent bond truth. Charge accuracy, redox, universal family/topology transfer, physical uncertainty calibration and causal electronic mechanisms remain outside the evidence. All four finite tests are finished; these outside-scope topics are not unfinished requirements for this closeout.

See `SCIENTIFIC_CLOSEOUT.md`, `RESULT_TABLES.md`, `METHODS.md`, `FIGURE_CAPTIONS.md`, and `results/closeout_v1/claim_to_evidence.csv` for the final evidence record.
