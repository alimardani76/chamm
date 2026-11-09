# Full production reproduction

The release provides every located input needed by the completed project.
The full fresh fit is intentionally not an advertised passed release check.
Use the original runtime records in `requirements-observed.txt` and
`reports/runtime_environment.json` as provenance. Binary checkpoint replay
is tested in that existing Windows Python 3.13.13 runtime. A portable full
environment lock has not been validated.

Retrieve and verify the release deposit first (`python reproduce/retrieve.py`).
All file hashes, including the actual source archive, are in
`production_files.json`; the archive MD5 is also checked by stage 14 against
publisher metadata. The original geometry cache excludes electronic data.

The original production sequence is preserved below. Execute only in a copy.
Stages can write/replace outputs, and existing parts or checkpoints can be
reused. Move the outputs being regenerated to an archive in that copy first
if you require a fresh computation. Do not erase the supplied scientific
record. Retain the original `results/models/training_samples.parquet`,
`results/models/model_configuration.json`, `results/grouped_splits.csv`,
metadata join, and fixed `results/sensitivity/subset_manifest.csv`.

```sh
python code/13_recover_charge_data.py --archive work/qmof_database.zip --metadata inputs/figshare_article_v18.json --benchmark inputs/expanded_benchmark_records.json --outdir recovered_charges
python code/14_prepare_source.py
python code/16_charge_analysis.py
python code/16b_calibration_check.py
python code/17_build_descriptors.py --workers 2
python code/18_grouped_models.py
python code/19_results.py
python code/24_conditional_closeout.py
python code/25_closeout_atlas.py
python code/26_nonmetal_closeout.py
python code/25_closeout_atlas.py --verify-compact
```

Create `recovered_charges/` before stage 13 if necessary. Stage 13 recovers
source charges; stage 16 serializes charge targets. Stages 14/17 consume
the cached geometry/accepted graphs. To independently rerun contact selection,
archive the accepted graph caches in the reproduction copy first; this
expensive whole-cohort calculation was not repeated for the release.
The accepted graph configuration and original decisions are in
`reports/CONNECTIVITY_DECISION.md` and `data/descriptor_manifest.json`.

Stages 15/15b/15c document the prior finite connectivity decision; stages
20/23 and the earlier metal graph outputs remain distinct from final
conditional models. Stage 26 changes nonmetal cutoffs only on the frozen
96-structure subset, holding accepted metal edges and charge-only residuals
fixed. No broad analysis campaign is needed to reproduce the manuscript
figures or compact numerical estimates.

Checks actually executed for this release, with their practical limits,
are recorded in `RELEASE_REPORT.md`. The full input replay script reads
pickled model checkpoints only from the checksum-verified release deposit.
