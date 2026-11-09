"""Finite original-fold and identifier/composition-blocked conditional comparison."""
import os
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2'
import argparse,gc,hashlib,json,pickle,time,zipfile
from pathlib import Path
from importlib import import_module
import numpy as np,pandas as pd
from xgboost import XGBRegressor
from chamm_closeout import *

def prepare_metadata():
    p=OUT/'metadata';p.mkdir(parents=True,exist_ok=True)
    cohort=pd.read_csv(ROOT/'data/chamm_common_cohort.csv');old=pd.read_csv(ROOT/'results/grouped_splits.csv')
    reviewer=pd.read_csv(ROOT/'ChaMM_scientific_closeout/verified_source_metadata.csv',keep_default_na=False).set_index('qmof_id')
    with zipfile.ZipFile(ROOT/'work/qmof_database.zip') as z:records=json.loads(z.read('qmof_database/qmof.json'))
    source={r['qmof_id']:r for r in records if r['qmof_id'] in set(cohort.qmof_id)}
    assert set(source)==set(cohort.qmof_id)==set(reviewer.index)
    for q,r in source.items():
        info=r['info'];m=info.get('mofid') or {}
        assert reviewer.loc[q,'formula_reduced']==info['formula_reduced']
        assert reviewer.loc[q,'mofkey']==(m.get('mofkey') or '')
    meta=reviewer.loc[cohort.qmof_id].reset_index();assert np.array_equal(meta.fold,meta.qmof_id.map(old.set_index('qmof_id').fold))
    ids=meta.qmof_id.tolist();parent={q:q for q in ids}
    def find(q):
        while parent[q]!=q:parent[q]=parent[parent[q]];q=parent[q]
        return q
    for col in ('mofkey','formula_reduced'):
        for val,g in meta[meta[col].ne('')].groupby(col):
            root=find(g.qmof_id.iloc[0])
            for q in g.qmof_id.iloc[1:]:parent[find(q)]=root
    buckets={}
    for q in ids:buckets.setdefault(find(q),[]).append(q)
    comps={q:'cc_'+min(v) for v in buckets.values() for q in v};meta['component']=meta.qmof_id.map(comps)
    sizes=meta.groupby('component').size();meta['component_size']=meta.component.map(sizes);meta['has_mofkey']=meta.mofkey.ne('');meta['old_fold']=meta.pop('fold')
    oldspan=meta.groupby('component').old_fold.nunique();meta['component_crossed_old_folds']=meta.component.map(oldspan)>1
    order=sorted(sizes.index,key=lambda c:(-int(sizes[c]),hashlib.sha256(('ChaMM-closeout-blocks:'+c).encode()).hexdigest()))
    counts=np.zeros(5,dtype=int);mapping={}
    for c in order:f=int(np.argmin(counts));mapping[c]=f;counts[f]+=int(sizes[c])
    meta['blocked_fold']=meta.component.map(mapping)
    summary=dict(source='qmof_database/qmof.json in pinned existing archive; reviewer join exact formulas and MOFkeys matched for all cohort IDs',
        records=len(meta),mofkey_records=int(meta.has_mofkey.sum()),distinct_mofkeys=meta.loc[meta.has_mofkey,'mofkey'].nunique(),
        formula_records=int(meta.formula_reduced.ne('').sum()),distinct_formulas=meta.formula_reduced.nunique(),components=len(sizes),largest_component=int(sizes.max()),
        nonsingleton_records=int((meta.component_size>1).sum()),records_components_crossing_old_folds=int(meta.component_crossed_old_folds.sum()),fold_counts=counts.tolist(),
        grouping='connected components: identical nonempty exact MOFkey OR identical exact source reduced formula; no topology inference',
        fold_rule='descending component size; SHA256(ChaMM-closeout-blocks:component) tie order; assign least populated fold, fold-index tie break')
    assert (summary['mofkey_records'],summary['distinct_mofkeys'],summary['distinct_formulas'],summary['components'],summary['largest_component'],summary['nonsingleton_records'],summary['records_components_crossing_old_folds'])==(5408,5114,8468,8347,13,1414,1189)
    meta.to_csv(p/'component_mapping.csv',index=False);meta[['qmof_id','component','component_size','blocked_fold','has_mofkey']].to_csv(p/'blocked_folds.csv',index=False)
    reviewer.reset_index().to_csv(p/'verified_source_metadata.csv',index=False);(p/'grouping_summary.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)

def reconcile_legacy(t):
    p=OUT/'legacy_global_reconciliation.json'
    if p.exists():return
    cfg=import_module('19_results');summary=pd.read_csv(ROOT/'results/science/grouped_metrics.csv');maxerr=0.;targetmax=0.;rows=0
    for pair in PAIRS:
        blocks=[]
        for fold in range(5):
            d=pd.read_parquet(ROOT/'results/models'/f'oof_fold{fold}_{pair}.parquet',columns=['qmof_id','atom_index','fold','observed','is_metal']+cfg.MODELS)
            ref=t[t.fold==fold];assert np.array_equal(d.qmof_id,ref.qmof_id) and np.array_equal(d.atom_index,ref.atom_index) and d.fold.eq(fold).all()
            targetmax=max(targetmax,float(np.max(np.abs(d.observed.to_numpy()-ref['delta_'+pair].to_numpy()))));blocks.append(d)
        d=pd.concat(blocks,ignore_index=True);assert len(d)==1420532
        for stratum,dd in [('all',d),('metal',d[d.is_metal]),('nonmetal',d[~d.is_metal])]:
            loss=pd.DataFrame({'qmof_id':dd.qmof_id})
            for m in cfg.MODELS:loss[m]=np.abs(dd.observed-dd[m])
            sm=loss.groupby('qmof_id')[cfg.MODELS].mean().mean();am=loss[cfg.MODELS].mean()
            for m in cfg.MODELS:
                r=summary[(summary.pair==pair)&(summary.stratum==stratum)&(summary.model==m)].iloc[0]
                maxerr=max(maxerr,abs(sm[m]-r.mae_structure),abs(am[m]-r.mae_atom));rows+=1
        del d,blocks;gc.collect()
    assert maxerr<1e-12 and targetmax<1e-12
    p.write_text(json.dumps(dict(status='PASS',aggregate_rows=rows,atom_targets_checked=3*1420532,max_target_error=targetmax,max_metric_error=maxerr,note='one local reconciliation of global atom OOF targets/folds and saved aggregate MAEs; reviewer bootstrap/graph audits not repeated'),indent=2));print('Legacy atom OOF reconciliation PASS',flush=True)

def run(splits):
    t=load_targets();reconcile_legacy(t);X,names,nfirst,center,cfg=feature_matrix();t['graph_valid']=True
    meta=pd.read_csv(OUT/'metadata/component_mapping.csv',keep_default_na=False).set_index('qmof_id')
    oldsamples=pd.read_parquet(ROOT/'results/models/training_samples.parquet')
    atomkeys=pd.MultiIndex.from_frame(t[['qmof_id','atom_index']])
    for split in splits:
        dest=OUT/split;dest.mkdir(exist_ok=True);modeldir=WORK/split;modeldir.mkdir(exist_ok=True)
        folds=t.fold.to_numpy() if split=='original' else t.qmof_id.map(meta.blocked_fold).to_numpy()
        config=dict(schema='conditional-closeout-v1',split=split,tree=cfg['tree'],cap=160000,per_structure_element_cap=12,
            target='delta_AB about common training-fitted 12-bin piecewise offset; fixed pair direction',
            models=['sample-weighted center-element residual mean','identity + q_A','identity + q_A + first shell','identity + q_A + first shell + twohop paths'],
            training_weights='identical among competitors; original saved sample/weights or same rule regenerated within blocked training structures',
            descriptor_sha256=cfg['descriptor_table_sha256'],first_features=names[:nfirst],second_features=names[nfirst:],center_features=[names[j] for j in center],
            split_sha256=hashlib.sha256(pd.DataFrame({'qmof_id':t.qmof_id,'fold':folds}).drop_duplicates().to_csv(index=False).encode()).hexdigest())
        cp=dest/'configuration.json'
        if cp.exists():assert json.loads(cp.read_text())==config
        else:cp.write_text(json.dumps(config,indent=2))
        samples=[]
        for fold in range(5):
            if split=='original':
                s=oldsamples[oldsamples.outer_fold==fold];trainids=atomkeys.get_indexer(pd.MultiIndex.from_frame(s[['qmof_id','atom_index']]));w=s.fit_weight.to_numpy();assert (trainids>=0).all()
            else:
                tt=t[['qmof_id','element','structure_weight','graph_valid']].copy();tt['fold']=folds
                trainids,w=import_module('18_grouped_models').training_sample(tt,fold);s=t.loc[trainids,['qmof_id','atom_index']].copy();s['outer_fold']=fold;s['fit_weight']=w
            assert np.all(folds[trainids]!=fold);samples.append(s);testids=np.flatnonzero(folds==fold)
            for pair in PAIRS:
                op=dest/f'oof_fold{fold}_{pair}.parquet'
                if op.exists():print('reuse',split,op.name,flush=True);continue
                start=time.monotonic();base,supported,pars,knots=calibration(t,folds,fold,pair);y=t['delta_'+pair].to_numpy();qa=t[pair.split('_')[0]].to_numpy()
                if split=='original':
                    old=pd.read_parquet(ROOT/'results/models'/f'oof_fold{fold}_{pair}.parquet',columns=['nonlinear']);assert np.allclose(base[testids],old.nonlinear,atol=1e-9,rtol=0)
                pd.DataFrame(pars).to_csv(dest/f'calibration_fold{fold}_{pair}.csv',index=False);pd.DataFrame(knots).to_csv(dest/f'knots_fold{fold}_{pair}.csv',index=False)
                residual=y[trainids]-base[trainids];bias={}
                trainels=t.loc[trainids,'element'].to_numpy()
                for el in np.unique(trainels):mask=trainels==el;bias[el]=float(np.average(residual[mask],weights=w[mask]))
                out=t.loc[testids,['qmof_id','atom_index','element','is_metal']].copy();out['fold']=fold;out['observed']=y[testids];out['q_A']=qa[testids];out['piecewise']=base[testids];out['center']=base[testids]+t.loc[testids,'element'].map(bias).fillna(0).to_numpy();out['calibration_supported']=supported[testids];out['center_seen']=t.loc[testids,'element'].isin(bias).to_numpy()
                for label,columns in [('charge_only',center),('conditional_first',list(range(nfirst))),('conditional_twohop',list(range(len(names))))]:
                    model=XGBRegressor(**cfg['tree']);model.fit(model_input(X,trainids,columns,qa),residual,sample_weight=w)
                    out[label]=base[testids]+predict(model,X,testids,columns,qa)
                    with (modeldir/f'{pair}_fold{fold}_{label}.pkl').open('wb') as stream:pickle.dump(dict(model=model,columns=columns,structural_feature_names=[names[j] for j in columns],q_A_column=pair.split('_')[0],center_bias=bias,configuration=config['schema']),stream,protocol=5)
                    del model;gc.collect()
                out['residual_charge_only']=out.observed-out.charge_only;out['absolute_residual_charge_only']=np.abs(out.residual_charge_only)
                assert np.isfinite(out[list(MODELS)+['observed']].to_numpy()).all()
                temp=op.with_suffix('.parquet.partial');out.to_parquet(temp,index=False,compression='zstd');temp.replace(op)
                print(f'{split} fold {fold} {pair}: fit {len(trainids)}, test {len(testids)}, {time.monotonic()-start:.1f}s',flush=True)
        pd.concat(samples,ignore_index=True).to_parquet(dest/'training_samples.parquet',index=False,compression='zstd')
        (dest/'COMPLETE.json').write_text(json.dumps(dict(partitions=15,atoms_per_pair=1420532,structures=9171),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--split',choices=['original','blocked','both'],default='both');p.add_argument('--prepare-only',action='store_true');a=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);WORK.mkdir(parents=True,exist_ok=True)
    prepare_metadata()
    if not a.prepare_only:run(['original','blocked'] if a.split=='both' else [a.split])
