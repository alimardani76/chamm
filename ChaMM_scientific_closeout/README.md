# ChaMM scientific closeout

1. Read SCIENTIFIC_REVIEW.md for the checked findings, limits and remaining scientific work.
2. Open local Codex in the existing completed ChaMM project folder and paste LOCAL_CODEX_SCIENTIFIC_CLOSEOUT_PROMPT.md.
3. Keep the completed production datasets and caches. The prompt executes a finite follow-up using those outputs; it does not restart the graph build.

Supporting review evidence:
- numerical_review_evidence.json: recomputed summary metrics/intervals and independent graph bookkeeping.
- verified_source_metadata.csv: all 9,171 cohort IDs joined to fields actually present in the pinned QMOF metadata member.
- metadata_summary.json: exact metadata coverage and cross-fold repeated groups.
- proposed_blocking_summary.json: counts for the proposed identifier/composition grouping rule.

The metadata was read from qmof_database/qmof.json in the already checksum-verified QMOF Figshare v18 archive, file 59573735, MD5 0d89aaf66f2c306e86e47fd91cb1346e. The exact metadata values are source data; the proposed grouping rule is a reviewer-defined stress test, not an independently established family label.

The uploaded production review ZIP was inspected without modification. These follow-up files are new. No full model retraining or descriptor reconstruction was performed during this review.
