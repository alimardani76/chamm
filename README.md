# ChaMM

**Local Chemical Environments Explain Charge-Definition Disagreement in Metal–Organic Frameworks**

Hosein Alimardani, Shayan Abaei, and Mehrdad Asgari.

Local chemistry adds predictive information beyond the tested element and
source-charge controls in 9,171 QMOF structures / 1,420,532 atoms. The result
is strongest for nonmetals. The global conditional DDEC6–Bader central-metal
result is negative and is retained. Signed departure, absolute departure, and
within-environment residual spread are distinct empirical quantities; none
is a physical charge-error estimate. See `reports/SCIENTIFIC_CLOSEOUT.md`,
`reports/SCIENTIFIC_FINDINGS.md`, and `reports/CLAIM_TO_EVIDENCE.md`.

Release: https://github.com/alimardani76/chamm/releases/tag/v1.0.0 .
The fixed cohort, graph `621045db53a892f1`, ordered targets, splits, and
scientific conclusions are unchanged. No new model search is part of this release.

## Repository contents

| Path | Contents |
|---|---|
| `paper/` | Main/SI LaTeX, bibliography, vector figure masters, editable exports and compiled PDFs |
| `figures/` | Exact plotting CSVs, SHA256 provenance and working Python generators |
| `code/` | Original production and finite closeout analysis source (original path retained for imports) |
| `results/` | Compact numerical evidence, losses, contrasts, manifests, calibration and training samples |
| `data/` | Cohort and descriptor manifests; full atom tables in release deposit |
| `reports/`, `evidence/` | Scientific decisions, numerical claims, primary-reference checks and provenance |
| `reproduce/` | Retrieval, input verification, build commands and dependency records |

## Reproduction levels

Use a clean clone, Python 3.13, and a virtual environment. Commands below are
from the repository root. On Windows use `.venv/Scripts/python`; on Unix use
`.venv/bin/python`. `python` below means that environment's interpreter.

### 1. Figures and manuscript

```sh
python -m venv .venv
```

Activate it with `.\.venv\Scripts\Activate.ps1` in PowerShell or
`. .venv/bin/activate` in Bash. Then run:

```sh
python -m pip install -r reproduce/requirements-figures.txt
python figures/source/make_figures.py
python figures/source/make_tables.py
python figures/source/validate_package.py
python reproduce/build_paper.py
```

The complete installed figure dependency set is recorded in
`reproduce/requirements-figures-lock.txt`. These commands
rebuild all five main and five SI figures as vector PDF, editable SVG text,
and 400-dpi PNG previews from frozen saved estimates. No model fit or new
uncertainty calculation is involved. Main/SI tables are generated from the
same saved numerical results. The final command additionally needs pdfLaTeX
and BibTeX (TeX Live or MiKTeX). `paper/COMPILATION.md` covers Overleaf.

### 2. Compact numerical reconstruction

Using the same figure environment, run:

```sh
python code/25_closeout_atlas.py --verify-compact
```

This reaggregates 1,320 model-summary rows, 1,584 paired-increment rows,
180 headline contrasts and 5,910 atlas rows. It recalculates the six principal
paired bootstrap intervals, checks training/test separation and component
fold separation, and verifies the original training samples exactly.
It does not refit models. Numerical agreement is checked to `1e-12`.

### 3. Full model and descriptor inputs

The actual completed production project supplies a separate 2.01-GB ZIP
release asset, `ChaMM_Production_Inputs_v1.zip`; Git does not contain this
asset. It contains the full frozen atom tables, pinned QMOF archive, all
9,171 accepted graph/structure caches and descriptor parts, legacy and final
saved predictions, trained model checkpoints, folds and training inputs.
The regenerable closeout numeric matrix is excluded. No required scientific
input is known missing. Extract into a separate reproduction clone if you
intend to rerun stages that write scientific outputs.

```sh
python reproduce/retrieve.py
python -m pip install -r reproduce/requirements-observed.txt
python reproduce/check_full_inputs.py
```

`retrieve.py` downloads the real release asset, checks its size and SHA256,
and checks every extracted file against `production_files.json`. If already
downloaded, pass `--archive /path/to/ChaMM_Production_Inputs_v1.zip`.
Interrupted downloads retain a partial file and resume where the server
supports HTTP Range. `python reproduce/verify_remote.py` checks public asset
metadata, the server digest when available, and the live download prefix.
The audit checks all atom keys and three signed targets, replays accepted
descriptors for three declared structures, and replays one saved blocked
model on its full held-out fold. These are input and replay checks, not an
end-to-end training rerun. The original observed environment is provenance,
not a claim that those pins install on every platform. The independently
installed figure environment is recorded separately; its NumPy major
version differs from the original model runtime.

For source recovery and descriptor regeneration, the preserved execution
order is stages 13, 14, 16, 17, 18, then the finite closeout stages 24–26.
Read `reproduce/FULL_REPRODUCTION.md` before running these writing stages.
An end-to-end refit and whole-cohort graph rebuild were not performed during
release preparation. Do not interpret existing checkpoint reuse as a fresh fit.

## Attribution and author decisions

QMOF data retain CC BY 4.0 attribution; see `THIRD_PARTY.md` and the manuscript
bibliography. Original code/manuscript licensing, journal selection, author
metadata confirmation, funding and contribution declarations require author
decisions. No missing funding or roles have been invented.
