# Final numerical tables

All values in e unless counts. Full atom and structure MAEs, fold increments, component-resampling inputs and support are in results/closeout_v1/.

## Nested comparison, all structures

| split | pair | piecewise | center | charge_only | conditional_first | conditional_twohop | conditional_gain | ci_low | ci_high | folds_improved |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| original | ddec6_cm5 | 0.02628 | 0.02635 | 0.02443 | 0.01315 | 0.01022 | 0.01420 | 0.01407 | 0.01434 | 5 |
| original | ddec6_bader | 0.09297 | 0.09289 | 0.08457 | 0.04184 | 0.03889 | 0.04568 | 0.04511 | 0.04628 | 5 |
| original | cm5_bader | 0.07435 | 0.07445 | 0.06738 | 0.03969 | 0.03899 | 0.02838 | 0.02792 | 0.02886 | 5 |
| blocked | ddec6_cm5 | 0.02628 | 0.02635 | 0.02434 | 0.01309 | 0.01018 | 0.01416 | 0.01400 | 0.01432 | 5 |
| blocked | ddec6_bader | 0.09295 | 0.09286 | 0.08454 | 0.04197 | 0.03885 | 0.04569 | 0.04495 | 0.04640 | 5 |
| blocked | cm5_bader | 0.07436 | 0.07445 | 0.06750 | 0.03981 | 0.03914 | 0.02837 | 0.02782 | 0.02893 | 5 |

## Blocked challenge support

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

## Global central-metal result

| pair | charge_only | conditional_twohop | gain | ci_low | ci_high | folds_improved |
| --- | --- | --- | --- | --- | --- | --- |
| ddec6_cm5 | 0.02694 | 0.02202 | 0.00493 | 0.00451 | 0.00533 | 5 |
| ddec6_bader | 0.02044 | 0.02851 | -0.00807 | -0.00857 | -0.00749 | 0 |
| cm5_bader | 0.02959 | 0.02868 | 0.00091 | 0.00049 | 0.00137 | 5 |

## Signed and magnitude headline contrasts, blocked folds

| pair | left | right | quantity | n_structures | blocked_effect | blocked_ci_low | blocked_ci_high | validation_status | graph_status | graph_common_structures |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ddec6_cm5 | C:O-attached | C:C/H-only | calibrated_signed | 7601 | 0.00374 | 0.00309 | 0.00440 | positive in both validation schemes | unchanged in tested subset | 74 |
| ddec6_cm5 | C:O-attached | C:C/H-only | calibrated_absolute | 7601 | -0.00251 | -0.00299 | -0.00202 | negative in both validation schemes | unchanged in tested subset | 74 |
| ddec6_cm5 | C:N-attached | C:C/H-only | calibrated_signed | 6689 | -0.03540 | -0.03629 | -0.03448 | negative in both validation schemes | unchanged in tested subset | 70 |
| ddec6_cm5 | C:N-attached | C:C/H-only | calibrated_absolute | 6689 | 0.01264 | 0.01192 | 0.01335 | positive in both validation schemes | unchanged in tested subset | 70 |
| ddec6_cm5 | H:O-bound | H:C-bound | calibrated_signed | 3676 | 0.00371 | 0.00307 | 0.00432 | positive in both validation schemes | unchanged in tested subset | 57 |
| ddec6_cm5 | H:O-bound | H:C-bound | calibrated_absolute | 3676 | 0.00644 | 0.00613 | 0.00674 | positive in both validation schemes | unchanged in tested subset | 57 |
| ddec6_cm5 | H:N-bound | H:C-bound | calibrated_signed | 1745 | -0.02058 | -0.02144 | -0.01963 | negative in both validation schemes | unchanged in tested subset (limited paired support) | 14 |
| ddec6_cm5 | H:N-bound | H:C-bound | calibrated_absolute | 1745 | 0.01312 | 0.01237 | 0.01389 | positive in both validation schemes | unchanged in tested subset (limited paired support) | 14 |
| ddec6_cm5 | O:metal-adjacent | O:no-metal | calibrated_signed | 4622 | -0.00829 | -0.01006 | -0.00673 | negative in both validation schemes | unchanged in tested subset | 52 |
| ddec6_cm5 | O:metal-adjacent | O:no-metal | calibrated_absolute | 4622 | 0.00020 | -0.00083 | 0.00120 | weak/null or inconsistent; no recurring direction established | unchanged in tested subset | 52 |
| ddec6_cm5 | N:metal-adjacent | N:no-metal | calibrated_signed | 3216 | 0.01739 | 0.01520 | 0.01948 | positive in both validation schemes | unchanged in tested subset (limited paired support) | 25 |
| ddec6_cm5 | N:metal-adjacent | N:no-metal | calibrated_absolute | 3216 | -0.01361 | -0.01499 | -0.01216 | negative in both validation schemes | unchanged in tested subset (limited paired support) | 25 |
| ddec6_bader | C:O-attached | C:C/H-only | calibrated_signed | 7601 | -0.07210 | -0.07491 | -0.06946 | negative in both validation schemes | unchanged in tested subset | 74 |
| ddec6_bader | C:O-attached | C:C/H-only | calibrated_absolute | 7601 | -0.00112 | -0.00342 | 0.00126 | weak/null or inconsistent; no recurring direction established | unchanged in tested subset | 74 |
| ddec6_bader | C:N-attached | C:C/H-only | calibrated_signed | 6689 | -0.20842 | -0.21282 | -0.20383 | negative in both validation schemes | unchanged in tested subset | 70 |
| ddec6_bader | C:N-attached | C:C/H-only | calibrated_absolute | 6689 | 0.12735 | 0.12413 | 0.13032 | positive in both validation schemes | unchanged in tested subset | 70 |
| ddec6_bader | H:O-bound | H:C-bound | calibrated_signed | 3676 | 0.00627 | 0.00503 | 0.00747 | positive in both validation schemes | unchanged in tested subset | 57 |
| ddec6_bader | H:O-bound | H:C-bound | calibrated_absolute | 3676 | -0.00107 | -0.00180 | -0.00026 | weak/null or inconsistent; no recurring direction established | unchanged in tested subset | 57 |
| ddec6_bader | H:N-bound | H:C-bound | calibrated_signed | 1745 | -0.02474 | -0.02731 | -0.02211 | negative in both validation schemes | unchanged in tested subset (limited paired support) | 14 |
| ddec6_bader | H:N-bound | H:C-bound | calibrated_absolute | 1745 | 0.01857 | 0.01710 | 0.02025 | positive in both validation schemes | unchanged in tested subset (limited paired support) | 14 |
| ddec6_bader | O:metal-adjacent | O:no-metal | calibrated_signed | 4622 | -0.03348 | -0.03698 | -0.03014 | negative in both validation schemes | unchanged in tested subset | 52 |
| ddec6_bader | O:metal-adjacent | O:no-metal | calibrated_absolute | 4622 | -0.04916 | -0.05173 | -0.04639 | negative in both validation schemes | unchanged in tested subset | 52 |
| ddec6_bader | N:metal-adjacent | N:no-metal | calibrated_signed | 3216 | 0.02939 | 0.02020 | 0.03817 | positive in both validation schemes | unchanged in tested subset (limited paired support) | 25 |
| ddec6_bader | N:metal-adjacent | N:no-metal | calibrated_absolute | 3216 | -0.02139 | -0.02826 | -0.01459 | negative in both validation schemes | unchanged in tested subset (limited paired support) | 25 |
| cm5_bader | C:O-attached | C:C/H-only | calibrated_signed | 7601 | -0.11184 | -0.11450 | -0.10904 | negative in both validation schemes | unchanged in tested subset | 74 |
| cm5_bader | C:O-attached | C:C/H-only | calibrated_absolute | 7601 | 0.05603 | 0.05409 | 0.05811 | positive in both validation schemes | unchanged in tested subset | 74 |
| cm5_bader | C:N-attached | C:C/H-only | calibrated_signed | 6689 | 0.00068 | -0.00378 | 0.00540 | weak/null or inconsistent; no recurring direction established | unchanged in tested subset | 70 |
| cm5_bader | C:N-attached | C:C/H-only | calibrated_absolute | 6689 | 0.10216 | 0.09904 | 0.10558 | positive in both validation schemes | unchanged in tested subset | 70 |
| cm5_bader | H:O-bound | H:C-bound | calibrated_signed | 3676 | -0.03877 | -0.04003 | -0.03752 | negative in both validation schemes | unchanged in tested subset | 57 |
| cm5_bader | H:O-bound | H:C-bound | calibrated_absolute | 3676 | 0.02712 | 0.02623 | 0.02806 | positive in both validation schemes | unchanged in tested subset | 57 |
| cm5_bader | H:N-bound | H:C-bound | calibrated_signed | 1745 | 0.09484 | 0.09275 | 0.09678 | positive in both validation schemes | unchanged in tested subset (limited paired support) | 14 |
| cm5_bader | H:N-bound | H:C-bound | calibrated_absolute | 1745 | 0.07306 | 0.07084 | 0.07508 | positive in both validation schemes | unchanged in tested subset (limited paired support) | 14 |
| cm5_bader | O:metal-adjacent | O:no-metal | calibrated_signed | 4622 | -0.01483 | -0.01842 | -0.01167 | negative in both validation schemes | unchanged in tested subset | 52 |
| cm5_bader | O:metal-adjacent | O:no-metal | calibrated_absolute | 4622 | -0.03288 | -0.03561 | -0.03022 | negative in both validation schemes | unchanged in tested subset | 52 |
| cm5_bader | N:metal-adjacent | N:no-metal | calibrated_signed | 3216 | -0.06642 | -0.07375 | -0.05926 | negative in both validation schemes | unchanged in tested subset (limited paired support) | 25 |
| cm5_bader | N:metal-adjacent | N:no-metal | calibrated_absolute | 3216 | -0.05782 | -0.06281 | -0.05234 | negative in both validation schemes | unchanged in tested subset (limited paired support) | 25 |
