# ChaMM release validation

Release identifier: v1.0.0, https://github.com/alimardani76/chamm/releases/tag/v1.0.0 .
The exact source commit is the commit referenced by this release tag; resolve
it with `git rev-parse v1.0.0^{commit}`. The delivered local release report
records the actual pushed SHA and final CI outcome after publication.

The original project and supplied development package remain unchanged.
The main Introduction, Abstract and Conclusions have been rewritten, and
Methods and Results revised around conditional chemical information,
quantitative carbon contrasts and the negative central-metal result.
All five main and five SI figures have consistent typography, exact saved
inputs and PDF/SVG/400-dpi PNG exports. Main and SI tables are generated.

Checks completed before initial publication:

- All 329 scientific artifact-index entries match the original local SHA256.
- Compact reconstruction: 1,320 summary rows, 1,584 increments, 180 headline
  contrasts, 5,910 atlas rows, principal intervals, split separation and
  exact original training samples reconcile to better than `1e-12`.
- The compact figure/source validator passes 99 checks.
- Main and SI compile through pdfLaTeX/BibTeX with no undefined references
  or overfull boxes. Rendered main/SI pages were inspected; a crowded main
  figure model-axis label was corrected.
- The actual full production inputs are deposited separately with an archive
  SHA256 and 28,458 file-level hashes. The regenerable numeric matrix is
  excluded. No required scientific input is known missing.
- Text-source credential-pattern scan found no credentials. Curated files
  contain only this project's manuscript, code, numerical evidence and
  reproduction metadata. No unrelated repository remote was changed.
- Primary-reference metadata and claim-specific sources were reviewed;
  PACMAN pagination was corrected and relevant correction notices added.

Figure and compact numerical reproduction are distinct from full refitting.
The full-input audit replays three accepted-descriptor structures and one
saved model; it does not retrain models. An end-to-end refit, a whole-cohort
graph rebuild, and a portable full-model installation lock are not claimed.
Original observed dependencies remain distinct from the clean figure
environment. Final clean-checkout, retrieval, replay and CI results are
recorded after the published commit is verified.

Remaining author decisions: journal and journal-specific format; confirmation
of supplied affiliations, correspondence, equal contribution and conflicts
statement; original software/manuscript licensing; funding/acknowledgements
and contribution declarations; coauthor approval and submission. None has
been invented or treated as a missing scientific experiment.
