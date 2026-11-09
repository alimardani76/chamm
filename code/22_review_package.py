"""Assemble completed numerical results, artifact index and compact review ZIP."""
from pathlib import Path
import argparse,json,hashlib,zipfile,platform
import pandas as pd,pyarrow.parquet as pq
ROOT=Path(__file__).resolve().parents[1]

def markdown(df):
    names=df.columns.tolist();rows=['| '+' | '.join(names)+' |','| '+' | '.join(['---']*len(names))+' |']
    for row in df.itertuples(index=False,name=None):
        rows.append('| '+' | '.join(f'{v:.4f}' if isinstance(v,float) else str(v) for v in row)+' |')
    return '\n'.join(rows)

def main():
    if (ROOT/'results/closeout_v1').exists():
        from importlib import import_module
        return import_module('27_closeout_review').generate()
    p=argparse.ArgumentParser();p.add_argument('--verdict',required=True);a=p.parse_args()
    manifest=json.loads((ROOT/'data/descriptor_manifest.json').read_text());cfg=json.loads((ROOT/'results/models/model_configuration.json').read_text())
    assert manifest['structures']==9171 and manifest['atoms']==1420532
    assert len(list((ROOT/'results/models').glob('oof_fold*.parquet')))==15
    assert len(list((ROOT/'figures').glob('0*.png')))==5
    metrics=pd.read_csv(ROOT/'results/science/grouped_metrics.csv');effects=pd.read_csv(ROOT/'results/science/paired_environment_contrasts.csv')
    sens=pd.read_csv(ROOT/'results/sensitivity/heldout_graph_metrics.csv')
    main=metrics[metrics.stratum.isin(['all','metal'])&metrics.model.isin(['offset','affine','nonlinear','affine_ridge_first','affine_tree_first','affine_tree_second','nonlinear_tree_first','nonlinear_tree_second'])]
    columns=['pair','stratum','model','baseline','n_atoms','n_structures','mae_atom','mae_structure','improvement_structure','improvement_ci_low','improvement_ci_high','n_folds_improved']
    strong=effects[(effects.control=='nonlinear')&(effects.quantity=='signed')&effects.supported_denominator]
    metal=pd.read_csv(ROOT/'results/metal_control/grouped_metrics.csv')
    text='# Numerical results for review\n\nAll charge units are e. Full tables retain other models, elements, flags and negative results. Improvement uses the named charge-only control for residual models. Confidence intervals resample held-out QMOF IDs.\n\n'+markdown(main[columns])+'\n\n## Bounded metal-only population control\n\n'+markdown(metal[metal.stratum.eq('all_metals')][['pair','model','n_atoms','n_structures','mae_atom','mae_structure','improvement_structure','ci_low','ci_high','n_folds_improved']])+'\n\n## Prespecified within-structure signed chemical contrasts after piecewise calibration\n\n'+markdown(strong[['pair','left','right','n_paired_structures','effect_e','ci_low','ci_high']])+'\n\n## Fixed-model graph sensitivity (global models; dedicated metal sensitivity in results/metal_control/graph_sensitivity.csv)\n\n'+markdown(sens[(sens.model=='nonlinear_tree_second')&sens.stratum.isin(['all','metal'])][['variant','pair','stratum','n_structures','accepted_mae_structure','alternative_mae_structure','baseline_mae_structure','alternative_minus_accepted','delta_ci_low','delta_ci_high']])+'\n'
    (ROOT/'reports/RESULT_TABLES.md').write_text(text,encoding='utf-8')
    # Coverage flags count participation, not validated errors.
    coverage=pd.DataFrame([dict(flag=k,n_atoms=v,fraction_of_common_atoms=v/manifest['atoms']) for k,v in manifest['flags'].items()])
    coverage.to_csv(ROOT/'results/science/connectivity_coverage.csv',index=False)
    cols=[];predictors=set(cfg['first_features']+cfg['second_features']);stored_structural=set(pq.read_schema(ROOT/'data/structural_features.parquet').names)
    for block,path in [('structural','data/structural_features.parquet'),('target','data/charge_targets.parquet')]:
        for field in pq.read_schema(ROOT/path):
            cols.append(dict(block=block,column=field.name,type=str(field.type),structural_predictor=block=='structural' and field.name in predictors,
                meaning=('two-hop non-backtracking path count, not unique atoms' if field.name.startswith('twohop') else 'review diagnostic, not bond truth' if field.name in manifest['flags'] else 'see METHODS.md and code/chamm_graph.py')))
    for name in sorted(predictors-stored_structural):
        cols.append(dict(block='derived structural',column=name,type='float32',structural_predictor=True,meaning='central element one-hot' if name.startswith('central_') else 'count parsed from neighbor_metal_formula'))
    pd.DataFrame(cols).to_csv(ROOT/'reports/data_dictionary.csv',index=False)
    from chamm_graph import info
    from pymatgen.analysis.local_env import CovalentRadius
    from pymatgen.core import Element
    ed=pd.read_csv(ROOT/'results/charge/element_distributions.csv')
    elrows=[]
    for e in sorted(ed.element.unique()):
        row=dict(element=e,**info(e));row['ase_covalent_radius']=row.pop('radius');row['crystal_default_radius']=CovalentRadius.radius[e] if e in CovalentRadius.radius else Element(e).atomic_radius;elrows.append(row)
    pd.DataFrame(elrows).to_csv(ROOT/'reports/element_classification_and_radii.csv',index=False)
    captions='''# Figure captions

1. **Scheme differences and near-zero signs.** Left: atom-weighted mean signed differences by central element; structure denominators are shown, and full atom/structure denominators and distribution quantiles are in the input table. These three differences are algebraically dependent. Right: sign opposition at four magnitude thresholds. The strong reduction after excluding near-zero values limits interpretations based on raw numerical sign.

2. **Affine adequacy and residual calibration.** Held-out mean signed residuals in charge deciles for C and O, comparing training-fitted affine and fixed 12-bin piecewise controls. Deciles are for description only; every calibration is trained without its test structures. Curvature motivates retaining the bounded nonlinear control. Scheme A is the x variable, without asserting physical reference status.

3. **Interpretable chemical contrasts beyond scheme mapping.** Prespecified within-structure, same-element differences in signed residuals after affine or piecewise calibration. Each QMOF contributes one paired difference; bars show 95% intervals from 1,000 structure resamples. Denominators are paired structures. O/N attachment and neighboring-metal categories describe graph environments; they do not establish reaction mechanisms or oxidation states. Absolute residual contrasts and unsupported small comparisons remain in the tables.

4. **Structure-held-out model ladder.** Structure-balanced MAE with 95% QMOF-resampled intervals. Top: all atoms, including element/affine/piecewise controls, a linear baseline and fixed tree ladders. Bottom: metals, comparing the global models against one metal-only follow-up with an extra training-fitted element residual offset. Blue is the global population objective; green is metal-only weighting. Global fitting is capped; the small metal follow-up fits all training metals. Every held-out atom in the relevant population is evaluated. Two-hop means path composition. Global metal failures are preserved in numerical tables. Correlated pair targets and related QMOF structures limit the independence of confirmations.

5. **Connectivity interpretation and sensitivity.** Left: the accepted Cu2 O6 environment in qmof-8a95c27, projected onto two local principal axes; labels are actual three-dimensional distances. The longer O pair is nearly opposite at 179.999 degrees. This is consistent with extended coordination and makes no oxidation-state claim. Right: matched change in structure-balanced held-out MAE for fixed piecewise-residual two-hop models on the outcome-blind 96-structure subset: global fits for all atoms, dedicated fits for central metals. The +0.10 Å taper is a bounded alternate configuration; the original cutoff is a stress comparison with documented chemical failures. Intervals resample QMOF IDs. Subset coverage is limited and not a population estimate.

Every figure has PNG and SVG versions; exact plotting inputs are in `figures/inputs/`, and `code/21_figures.py` generates the figures.
'''
    (ROOT/'reports/FIGURE_CAPTIONS.md').write_text(captions,encoding='utf-8')
    readme='''# ChaMM: charge-definition disagreement in QMOF

This review package examines whether DDEC6/CM5/Bader differences retain local chemical structure after ordinary element offsets and within-element charge mapping. No scheme is ground truth and no redox assignment is made.

Start with [scientific findings](reports/SCIENTIFIC_FINDINGS.md), [numerical tables](reports/RESULT_TABLES.md), [connectivity decision](reports/CONNECTIVITY_DECISION.md), [methods and limitations](reports/METHODS.md), and [figure captions](reports/FIGURE_CAPTIONS.md). The compact [review ZIP](ChaMM_review.zip) excludes the full-cohort atom tables, global atom-level predictions, source archive and large caches; those remain in the local project and are listed in the artifact index. Small metal-only prediction files and structure-level uncertainty inputs are included. The original charge CSV can be rebuilt from the pinned source when reproducing from the small ZIP alone.

## Frozen data and outputs

- `data/chamm_atom_charges_recovered.csv.gz`: original recovered charge-only table; `data/chamm_common_cohort.csv`: fixed 9,171-structure cohort.
- `data/charge_targets.parquet`: 1,420,532 atom-keyed charges, signed/absolute differences, range, robust sign indicators and fixed folds. Targets are separate from structural features.
- `data/structural_features.parquet`: all common-cohort atoms; geometry/chemistry descriptors and explicit validity/ambiguity flags. `data/descriptor_manifest.json` records exact coverage and configuration identity. `reports/data_dictionary.csv` describes blocks and model inclusion.
- `results/grouped_splits.csv`: saved five-fold structure assignments.
- `results/charge/`: descriptive distributions, training calibration parameters/knots and full charge-only out-of-fold predictions.
- `results/models/oof_fold{0..4}_{pair}.parquet`: full valid-graph held-out predictions for the fixed model ladder. `training_samples.parquet` saves capped fitting keys/weights. `model_configuration.json` gives exact features/configuration.
- `results/science/`: grouped/element/fold metrics, QMOF-level losses, signed/absolute environment contrasts, representative atom keys and coverage.
- `results/metal_control/`: bounded metal-only follow-up, saved OOF predictions, element residual controls, and its graph sensitivity. Global negative metal comparisons are retained separately.
- `results/connectivity/`: all 34 benchmark graphs, decisive contact ledger/distances, changed-site records, graph checks and mapped supercell evidence.
- `results/sensitivity/`: outcome-blind 96-structure manifest, alternate descriptors/predictions and matched graph sensitivity results.
- `figures/`: five PNG/SVG figures and plotting inputs. `reports/`: collaborator review documents and environment/artifact records.
- `work/`: pinned source archive, geometry cache, configuration-specific graph/descriptor parts, fitted models and resumable logs. These are regenerable, not part of the small ZIP.

## Reproduction and resuming

Run from the project root in PowerShell, with Python and the scientific packages in `requirements-observed.txt`. `reports/runtime_environment.json` records the actual Windows runtime; observed versions are a record, not a claim that every binary wheel is portable to every Python/platform.

```powershell
python -u code/14_prepare_source.py
# Needed only when starting from the small review ZIP without the recovered charge CSV:
if (-not (Test-Path 'data/chamm_atom_charges_recovered.csv.gz')) {
    python -u code/13_recover_charge_data.py --archive work/qmof_database.zip --metadata inputs/figshare_article_v18.json --benchmark inputs/expanded_benchmark_records.json --outdir work/recovered_charge_input
    if ($LASTEXITCODE -ne 0) { throw 'Charge recovery failed; inspect its report before continuing' }
    Copy-Item -LiteralPath work/recovered_charge_input/chamm_atom_charges_recovered.csv.gz -Destination data/chamm_atom_charges_recovered.csv.gz
}
python -u code/15_connectivity_panel.py
python -u code/15b_geometry_checks.py
python -u code/15c_freeze_decision.py
python -u code/16_charge_analysis.py
python -u code/16b_calibration_check.py
python -u code/17_build_descriptors.py --workers 2
python -u code/20_graph_sensitivity.py --build-only
python -u code/18_grouped_models.py
python -u code/19_results.py
python -u code/20_graph_sensitivity.py --evaluate-only
python -u code/23_metal_control.py
python -u code/21_figures.py
python -u code/22_review_package.py --verdict "VERDICT_FROM_SCIENTIFIC_FINDINGS"
```

Stage 14 reuses the pinned archive and completed geometry cache. The download resumes a `.partial` file and verifies publisher size/MD5. It never substitutes a newer release. Stage 17 processes only missing per-structure descriptor parts; graph caches include source hashes and configuration identity. Re-run the same command after an interruption. Inspect `work/production_events.jsonl` and `work/production_build.log` for failures; invalid graph atoms remain traceable and do not alter charge-only membership. Stage 18 reuses completed pair/fold partitions only with the saved matching configuration. It uses at most 160,000 fitting atoms per fold, with every available structure-element cell represented; test coverage remains full. Retain `work/models/` to reproduce fixed-model alternate-graph inference without refitting.

The accepted connectivity decision is a reviewed scientific input, not an algorithmic pass label. Its physical definition is frozen in `reports/CONNECTIVITY_DECISION.md`. The original recovery instructions and source metadata remain in `inputs/RECOVERY_README.md` and `inputs/figshare_article_v18.json`. Recovery and initial saved implementation packages are preserved, not repeatedly re-audited.

## Sources

QMOF [Figshare version 18](https://figshare.com/articles/dataset/QMOF_Database/13147324/18), [official repository](https://github.com/Andrew-S-Rosen/QMOF). The archive is CC BY 4.0; cite the two QMOF papers linked in METHODS.md. ChaMM changes here are cohort extraction, operational graphs, descriptors and statistical analysis; no new DFT/experiments.
'''
    readme=readme.replace('VERDICT_FROM_SCIENTIFIC_FINDINGS',a.verdict)
    (ROOT/'README.md').write_text(readme,encoding='utf-8')
    status=f'''# Completed production/review run — 9 November 2026

Scientific verdict: **{a.verdict}**. Read `reports/SCIENTIFIC_FINDINGS.md` for the demonstrated scope and limitations.

- Pinned source size/MD5 verified; geometry cache: 9,171 structures / 1,420,532 atoms.
- Connectivity decision, 34-structure comparison, decisive ledger, perturbations and 1,132 mapped supercell checks completed.
- Full structural table finalized: {manifest['atoms']:,} atoms / {manifest['structures']:,} structures; {manifest['valid_graph_atoms']:,} graph-valid atoms. Charge alignment passed for every key/element.
- Graph configuration: `{manifest['configuration_id']}`. Review flags are preserved; see `results/science/connectivity_coverage.csv`.
- Full charge-only results and training-only offset/affine/nonlinear calibration controls completed.
- Five structure-held-out folds, 15 pair/fold OOF partitions, fixed Ridge/tree ladders and explicit deterministic training caps completed.
- Bounded metal-only population control completed with all 46,994 held-out metal atoms and fixed unchanged tree settings; global negative metal results retained.
- QMOF-resampled metrics, prespecified signed/absolute chemical contrasts and fixed-model graph sensitivity on 96 structures completed.
- Five figures with numerical inputs/captions, scientific review documents and compact review ZIP delivered.

No required expensive job remains interrupted. Source/graph/feature caches and fitted models are retained for resuming or reproducing. Re-run commands are in README.md. Review readiness does not resolve independent chemical bond truth, related-structure leakage, new-topology transfer, or causal electronic mechanisms.
'''
    (ROOT/'RUN_STATUS.md').write_text(status,encoding='utf-8')
    # Artifact identities for local frozen tables/predictions. Caches intentionally not walked.
    paths=[ROOT/'data/charge_targets.parquet',ROOT/'data/structural_features.parquet',ROOT/'data/descriptor_manifest.json',ROOT/'data/chamm_common_cohort.csv']
    paths+=sorted((ROOT/'results').rglob('*.csv'))+sorted((ROOT/'results/models').glob('oof_fold*.parquet'))
    paths+=sorted((ROOT/'results/charge').glob('*oof.parquet'))+[ROOT/'results/models/training_samples.parquet']
    paths+=sorted((ROOT/'results/metal_control').glob('*oof.parquet'))
    items=[]
    for path in paths:
        h=hashlib.sha256()
        with path.open('rb') as f:
            for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
        items.append(dict(path=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=h.hexdigest()))
    (ROOT/'reports/artifact_index.json').write_text(json.dumps(items,indent=2))
    # Review bundle includes all numerical tables, structure-level uncertainty inputs, figures and code.
    include=[ROOT/'README.md',ROOT/'RUN_STATUS.md',ROOT/'requirements-observed.txt']
    for folder in ('code','reports','figures'):
        include += [f for f in (ROOT/folder).rglob('*') if f.is_file() and '__pycache__' not in f.parts]
    include += [f for f in (ROOT/'results').rglob('*') if f.is_file() and f.suffix in ('.csv','.json')]
    include += [ROOT/'results/science/structure_losses.parquet',ROOT/'results/science/structure_contrasts.parquet']
    include += sorted((ROOT/'results/metal_control').glob('*oof.parquet'))
    include += [ROOT/'inputs'/f for f in ('expanded_benchmark_records.json','selection.json','figshare_article_v18.json','RECOVERY_README.md')]
    include += [ROOT/'data/chamm_common_cohort.csv',ROOT/'data/descriptor_manifest.json',ROOT/'data/recovery_report.json',ROOT/'data/serialized_table_validation.json']
    with zipfile.ZipFile(ROOT/'ChaMM_review.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for path in sorted(set(include)):z.write(path,path.relative_to(ROOT).as_posix())
    print('Review ZIP MB',(ROOT/'ChaMM_review.zip').stat().st_size/1e6,flush=True)

if __name__=='__main__':main()
