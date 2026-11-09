"""Outcome-blind 96-structure subset; taper perturbation and cutoff stress comparison."""
import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,json,hashlib,pickle,time,gc
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np,pandas as pd
from importlib import import_module
from chamm_graph import *
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/sensitivity'

def one(q):
    r=json.loads((ROOT/'work/structures'/(q+'.json')).read_text());s=structure(r);cg=cutoff(s)
    zg,meta=nearest(s,indices=[i for i,x in enumerate(s) if x.specie.is_metal],distance_shift=0.1)
    h=hybrid(s,cg,zg);validate(s,h)
    a=pd.DataFrame(features(r,s,h,meta,cg));a.to_parquet(OUT/'parts_shift'/(q+'.parquet'),index=False)
    b=pd.DataFrame(features(r,s,cg,{'errors':{},'uncertainty':{}},cg));b.to_parquet(OUT/'parts_cutoff'/(q+'.parquet'),index=False)
    return q,len(s)

def build():
    OUT.mkdir(parents=True,exist_ok=True)
    cfg={'accepted_configuration_id':identity()[0],'source_md5':'0d89aaf66f2c306e86e47fd91cb1346e','shift':0.1,'cutoff_scale':1.25,'selection_seed':'ChaMM-sensitivity-20261109'}
    cp=OUT/'configuration.json'
    if cp.exists():assert json.loads(cp.read_text())==cfg,'Sensitivity cache configuration mismatch'
    else:cp.write_text(json.dumps(cfg,indent=2))
    for v in ('shift','cutoff'):(OUT/('parts_'+v)).mkdir(exist_ok=True)
    cohort=pd.read_csv(ROOT/'data/chamm_common_cohort.csv');refs=json.loads((ROOT/'inputs/expanded_benchmark_records.json').read_text())
    seeds={r['qmof_id'] for r in refs};other=sorted(set(cohort.qmof_id)-seeds,key=lambda q:hashlib.sha256(('ChaMM-sensitivity-20261109:'+q).encode()).hexdigest())[:62]
    ids=sorted(seeds|set(other));assert len(ids)==96
    manifest=cohort.set_index('qmof_id').loc[ids].reset_index();manifest['selection']=manifest.qmof_id.map(lambda q:'retained diagnostic panel' if q in seeds else 'deterministic ID hash')
    manifest.to_csv(OUT/'subset_manifest.csv',index=False)
    pending=[q for q in ids if not (OUT/'parts_shift'/(q+'.parquet')).exists() or not (OUT/'parts_cutoff'/(q+'.parquet')).exists()]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for k,f in enumerate(as_completed([pool.submit(one,q) for q in pending])):
            q,n=f.result()
            if k%10==0:print('sensitivity built',k+1,len(pending),q,n,flush=True)
    for v in ('shift','cutoff'):
        pd.concat([pd.read_parquet(OUT/('parts_'+v)/(q+'.parquet')) for q in ids],ignore_index=True).to_parquet(OUT/('features_'+v+'.parquet'),index=False)
    (OUT/'selection_and_methods.json').write_text(json.dumps(dict(structures=96,selection='34 retained diagnostic structures + 62 smallest SHA256(ChaMM-sensitivity-20261109:qmof_id); no charge outcomes consulted',
        alternative='accepted hybrid, both CrystalNN taper endpoints shifted +0.10 angstrom; tested on panel',
        stress='original cutoff 1.25; known shell failures, not asserted equally chemically valid',
        evaluation='fixed outer-fold models, replacing structural features only; no retuning',limitations='small selected subset, not a population prevalence sample; no topology grouping'),indent=2))

def evaluate():
    models=import_module('18_grouped_models');science=import_module('19_results')
    ids=pd.read_csv(OUT/'subset_manifest.csv').qmof_id.tolist()
    accepted=pd.read_parquet(ROOT/'data/structural_features.parquet',filters=[('qmof_id','in',ids)])
    counts=[];sitechanges=[];metrics=[];contrasts=[];summaries=[];predictions=[]
    for variant in ('shift','cutoff'):
        af,names,nfirst=models.feature_frames(OUT/('features_'+variant+'.parquet'))
        af['first_formula']=pd.read_parquet(OUT/('features_'+variant+'.parquet'),columns=['first_formula']).first_formula.to_numpy()
        assert af.graph_valid.all(),'Alternate graph failure: retained feature flags must be reviewed before sensitivity inference'
        compare=accepted[['qmof_id','atom_index','element','cn','first_formula','ambiguous']].merge(af[['qmof_id','atom_index','cn','first_formula','ambiguous']],on=['qmof_id','atom_index'],suffixes=('_accepted','_alternative'),validate='one_to_one')
        compare['variant']=variant;compare['changed']=compare.first_formula_accepted.ne(compare.first_formula_alternative)|compare.cn_accepted.ne(compare.cn_alternative)
        # Formula/CN changes are a lower bound on physical edge changes; compare full cached edge identities below for shift.
        sitechanges.append(compare)
        for el,g in compare.groupby('element'):counts.append(dict(variant=variant,element=el,n_atoms=len(g),n_structures=g.qmof_id.nunique(),changed_cn_or_formula=int(g.changed.sum())))
        alt_env=af[['qmof_id','atom_index']].copy();alt_env['environment']=science.environment(af)
        for pair in science.PAIRS:
            dd=pd.concat([pd.read_parquet(p,filters=[('qmof_id','in',ids)]) for p in sorted((ROOT/'results/models').glob(f'oof_fold*_{pair}.parquet'))],ignore_index=True)
            alt=dd.copy()
            for fold in range(5):
                mask=dd.fold.eq(fold);keys=dd.loc[mask,['qmof_id','atom_index']]
                sub=keys.merge(af,on=['qmof_id','atom_index'],validate='one_to_one')
                for control in ('offset','affine','nonlinear'):
                    for level in ('first','second'):
                        label=control+'_tree_'+level
                        with (ROOT/'work/models'/f'{pair}_fold{fold}_{label}.pkl').open('rb') as stream:artifact=pickle.load(stream)
                        wanted=artifact['names']
                        xx=sub.reindex(columns=wanted,fill_value=0).to_numpy(dtype=np.float32)
                        alt.loc[mask,label]=dd.loc[mask,control].to_numpy()+artifact['model'].predict(xx)
            for stratum,mask in [('all',np.ones(len(dd),bool)),('metal',dd.is_metal.to_numpy(bool)),('nonmetal',~dd.is_metal.to_numpy(bool))]:
                for control in ('offset','affine','nonlinear'):
                    for level in ('first','second'):
                        label=control+'_tree_'+level
                        e0=np.abs(dd.observed-dd[label]);e1=np.abs(alt.observed-alt[label]);eb=np.abs(dd.observed-dd[control])
                        agg=pd.DataFrame(dict(qmof_id=dd.loc[mask,'qmof_id'],accepted=e0[mask],alternative=e1[mask],baseline=eb[mask])).groupby('qmof_id').mean()
                        delta=agg.alternative-agg.accepted;ci=science.bootstrap(delta.to_numpy())
                        metrics.append(dict(variant=variant,pair=pair,stratum=stratum,model=label,n_atoms=int(mask.sum()),n_structures=len(agg),
                            accepted_mae_structure=agg.accepted.mean(),alternative_mae_structure=agg.alternative.mean(),baseline_mae_structure=agg.baseline.mean(),
                            alternative_minus_accepted=delta.mean(),delta_ci_low=ci[0],delta_ci_high=ci[1]))
            alt=alt.merge(alt_env,on=['qmof_id','atom_index'],validate='one_to_one')
            ss,cc,_=science.contrasts(alt,pair)
            for r in ss:r['variant']=variant
            for r in cc:r['variant']=variant
            summaries+=ss;contrasts+=cc
            keep=['qmof_id','atom_index','fold','observed','offset','affine','nonlinear']+[c+'_tree_'+l for c in ('offset','affine','nonlinear') for l in ('first','second')]
            pp=alt[keep].copy();pp['pair']=pair;pp['variant']=variant;predictions.append(pp)
        print('sensitivity evaluated',variant,flush=True);del af;gc.collect()
    pd.DataFrame(counts).to_csv(OUT/'changed_shell_counts.csv',index=False)
    pd.concat(sitechanges).to_parquet(OUT/'site_changes.parquet',index=False)
    pd.DataFrame(metrics).to_csv(OUT/'heldout_graph_metrics.csv',index=False)
    pd.DataFrame(contrasts).to_csv(OUT/'paired_environment_contrasts.csv',index=False)
    pd.DataFrame(summaries).to_csv(OUT/'environment_summaries.csv',index=False)
    pd.concat(predictions).to_parquet(OUT/'alternative_oof_predictions.parquet',index=False)
    compare_environment_contrasts()

def compare_environment_contrasts():
    """Matched structure-level chemical contrasts, same subset and charge-only residuals."""
    science=import_module('19_results');ids=pd.read_csv(OUT/'subset_manifest.csv').qmof_id.tolist()
    obs=pd.read_parquet(OUT/'alternative_oof_predictions.parquet');obs=obs[obs.variant.eq('shift')]
    accepted=pd.read_parquet(ROOT/'data/structural_features.parquet',filters=[('qmof_id','in',ids)])
    keys=['pair','control','left','right','quantity'];baseparts=[];basesumm=[]
    for pair,g in obs.groupby('pair'):
        features=accepted[['qmof_id','atom_index','element','is_metal','ambiguous']].copy();features['environment']=science.environment(accepted)
        dd=g[['qmof_id','atom_index','observed','offset','affine','nonlinear']].merge(features,on=['qmof_id','atom_index'],validate='one_to_one')
        _,cc,pp=science.contrasts(dd,pair);basesumm+=cc;baseparts+=pp
    base=pd.DataFrame(baseparts);pd.DataFrame(basesumm).to_csv(OUT/'accepted_subset_environment_contrasts.csv',index=False);rows=[]
    for variant in ('shift','cutoff'):
        af=pd.read_parquet(OUT/('features_'+variant+'.parquet'));features=af[['qmof_id','atom_index','element','is_metal','ambiguous']].copy();features['environment']=science.environment(af)
        altparts=[]
        for pair,g in obs.groupby('pair'):
            dd=g[['qmof_id','atom_index','observed','offset','affine','nonlinear']].merge(features,on=['qmof_id','atom_index'],validate='one_to_one')
            _,_,pp=science.contrasts(dd,pair);altparts+=pp
        matched=base.merge(pd.DataFrame(altparts),on=['qmof_id']+keys,suffixes=('_accepted','_alternative'),validate='one_to_one')
        for group,g in matched.groupby(keys):
            delta=g.effect_e_alternative-g.effect_e_accepted;ci=science.bootstrap(delta.to_numpy())
            rows.append(dict(zip(keys,group),variant=variant,n_shared_paired_structures=len(g),accepted_contrast_e=g.effect_e_accepted.mean(),alternative_contrast_e=g.effect_e_alternative.mean(),change_e=delta.mean(),ci_low=ci[0],ci_high=ci[1]))
    pd.DataFrame(rows).to_csv(OUT/'matched_environment_sensitivity.csv',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--build-only',action='store_true');p.add_argument('--evaluate-only',action='store_true');a=p.parse_args()
    if not a.evaluate_only:build()
    if not a.build_only:evaluate()
