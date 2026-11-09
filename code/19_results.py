"""Structure-resampled paired OOF metrics and prespecified chemical contrasts."""
import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,pickle
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/science'
PAIRS=('ddec6_cm5','ddec6_bader','cm5_bader');SEED=20261109
MODELS=['offset','affine','nonlinear','affine_ridge_first']+[c+'_tree_'+l for c in ('offset','affine','nonlinear') for l in ('first','second')]

def bootstrap(values,seed=SEED,reps=1000):
    v=np.asarray(values);rng=np.random.default_rng(seed);out=[]
    for k in range(reps):out.append(np.mean(v[rng.integers(0,len(v),len(v))]))
    return np.quantile(out,[.025,.975])

def metrics(d,pair):
    records=[];grouped=[];folds=[]
    strata={'all':np.ones(len(d),bool),'metal':d.is_metal.to_numpy(bool),'nonmetal':~d.is_metal.to_numpy(bool),'unambiguous':~d.ambiguous.to_numpy(bool)}
    for el,n in d.groupby('element').qmof_id.nunique().items():
        if n>=50:strata[el]=d.element.eq(el).to_numpy()
    for label,mask in strata.items():
        dd=d.loc[mask];loss=pd.DataFrame({'qmof_id':dd.qmof_id.to_numpy(),'fold':dd.fold.to_numpy()})
        for model in MODELS:loss[model]=np.abs(dd.observed.to_numpy()-dd[model].to_numpy())
        st=loss.groupby('qmof_id')[MODELS].mean();n=loss.groupby('qmof_id').size();fo=loss.groupby('qmof_id').fold.first()
        st.index.name='qmof_id'
        for model in MODELS:
            baseline=model.split('_')[0] if '_' in model else ('offset' if model!='offset' else 'offset')
            values=st[model].to_numpy();improvement=st[baseline].to_numpy()-values
            ci=bootstrap(values);di=bootstrap(improvement)
            records.append(dict(pair=pair,stratum=label,model=model,baseline=baseline,n_atoms=len(dd),n_structures=len(st),
                mae_atom=loss[model].mean(),mae_structure=values.mean(),mae_structure_ci_low=ci[0],mae_structure_ci_high=ci[1],
                improvement_structure=improvement.mean(),improvement_ci_low=di[0],improvement_ci_high=di[1],
                improvement_atom=loss[baseline].mean()-loss[model].mean(),n_folds_improved=sum(improvement[fo.to_numpy()==f].mean()>0 for f in range(5))))
            if label in ('all','metal','nonmetal'):
                for q,v,delta in zip(st.index,values,improvement):grouped.append(dict(qmof_id=q,pair=pair,stratum=label,model=model,mae=v,improvement=delta,n_atoms=int(n[q])))
            for fold in range(5):
                vv=values[fo.to_numpy()==fold];imp=improvement[fo.to_numpy()==fold]
                folds.append(dict(pair=pair,stratum=label,model=model,fold=fold,n_structures=len(vv),mae_structure=vv.mean(),improvement_structure=imp.mean()))
    return records,grouped,folds

def environment(f):
    env=pd.Series('other',index=f.index)
    c=f.element.eq('C');env.loc[c&(f.first_n_O>0)]='C:O-attached';env.loc[c&(f.first_n_O==0)&(f.first_n_N>0)]='C:N-attached'
    env.loc[c&((f.first_n_C+f.first_n_H)>0)&(f.cn==(f.first_n_C+f.first_n_H))]='C:C/H-only'
    h=f.element.eq('H')&(f.cn==1)
    for el in ('C','N','O'):env.loc[h&(f['first_n_'+el]==1)]='H:'+el+'-bound'
    for el in ('O','N'):
        env.loc[f.element.eq(el)&(f.first_n_metal>0)]=el+':metal-adjacent'
        env.loc[f.element.eq(el)&(f.first_n_metal==0)]=el+':no-metal'
    m=f.is_metal.astype(bool)
    env.loc[m&(f.first_n_O>0)&(f.first_n_N==0)]='metal:O-no-N'
    env.loc[m&(f.first_n_N>0)&(f.first_n_O==0)]='metal:N-no-O'
    env.loc[m&(f.first_n_N>0)&(f.first_n_O>0)]='metal:mixed-ON'
    return env

CONTRASTS=[('C:O-attached','C:C/H-only'),('C:N-attached','C:C/H-only'),('H:O-bound','H:C-bound'),('H:N-bound','H:C-bound'),('O:metal-adjacent','O:no-metal'),('N:metal-adjacent','N:no-metal')]

def contrasts(d,pair):
    summaries=[];effects=[];perstructure=[]
    for control in ('offset','affine','nonlinear'):
        d['residual']=d.observed-d[control];d['absolute_residual']=np.abs(d.residual)
        g=d.groupby(['qmof_id','element','environment']).agg(signed=('residual','mean'),absolute=('absolute_residual','mean'),n_atoms=('observed','size'),ambiguous_fraction=('ambiguous','mean')).reset_index()
        for (el,env),gg in g.groupby(['element','environment']):
            summaries.append(dict(pair=pair,control=control,element=el,environment=env,n_structures=len(gg),n_atoms=gg.n_atoms.sum(),
                                  mean_signed_structure=gg.signed.mean(),mean_absolute_structure=gg.absolute.mean(),ambiguous_fraction_structure=gg.ambiguous_fraction.mean()))
        for left,right in CONTRASTS:
            a=g[g.environment==left].set_index('qmof_id');b=g[g.environment==right].set_index('qmof_id');ids=a.index.intersection(b.index)
            if not len(ids):continue
            for value in ('signed','absolute'):
                v=a.loc[ids,value].to_numpy()-b.loc[ids,value].to_numpy();ci=bootstrap(v)
                effects.append(dict(pair=pair,control=control,left=left,right=right,quantity=value,n_paired_structures=len(ids),n_left_atoms=int(a.loc[ids,'n_atoms'].sum()),n_right_atoms=int(b.loc[ids,'n_atoms'].sum()),
                                    effect_e=v.mean(),ci_low=ci[0],ci_high=ci[1],supported_denominator=len(ids)>=30))
                for q,x in zip(ids,v):perstructure.append(dict(qmof_id=q,pair=pair,control=control,left=left,right=right,quantity=value,effect_e=x))
    return summaries,effects,perstructure

def metal_contrasts(d,pair):
    """Prespecified same-element, between-environment associations, with QMOF resampling."""
    rows=[]
    for control in ('offset','affine','nonlinear'):
        dd=d[d.is_metal].copy();dd['r']=dd.observed-dd[control];dd['absr']=np.abs(dd.r)
        g=dd.groupby(['qmof_id','element','environment']).agg(signed=('r','mean'),absolute=('absr','mean'),n=('r','size')).reset_index()
        for el in ('Cu','Zn','Ni','Co','Fe','Mn','Zr','Al'):
            for left,right in [('metal:O-no-N','metal:N-no-O'),('metal:mixed-ON','metal:O-no-N')]:
                a=g[(g.element==el)&(g.environment==left)].set_index('qmof_id');b=g[(g.element==el)&(g.environment==right)].set_index('qmof_id')
                if len(a)<2 or len(b)<2:continue
                ids=a.index.union(b.index)
                for quantity in ('signed','absolute'):
                    av=a[quantity].reindex(ids).to_numpy();bv=b[quantity].reindex(ids).to_numpy();rng=np.random.default_rng(SEED);bs=[]
                    for k in range(1000):
                        ii=rng.integers(0,len(ids),len(ids))
                        if np.isfinite(av[ii]).any() and np.isfinite(bv[ii]).any():bs.append(np.nanmean(av[ii])-np.nanmean(bv[ii]))
                    ci=np.quantile(bs,[.025,.975]);rows.append(dict(pair=pair,control=control,element=el,left=left,right=right,quantity=quantity,
                        n_left_structures=len(a),n_right_structures=len(b),n_overlapping_structures=len(a.index.intersection(b.index)),
                        n_left_atoms=int(a.n.sum()),n_right_atoms=int(b.n.sum()),effect_e=np.nanmean(av)-np.nanmean(bv),ci_low=ci[0],ci_high=ci[1],
                        supported_denominator=min(len(a),len(b))>=30,interpretation='same element; observational between environments; other structural differences not controlled'))
    return rows

def twohop_increment():
    """Paired added-path improvement, using saved held-out QMOF losses only."""
    d=pd.read_parquet(OUT/'structure_losses.parquet');folds=pd.read_csv(ROOT/'results/grouped_splits.csv').set_index('qmof_id').fold;rows=[]
    for (pair,stratum),g in d.groupby(['pair','stratum']):
        pp=g.pivot(index='qmof_id',columns='model',values='mae')
        for control in ('offset','affine','nonlinear'):
            a=pp[control+'_tree_first'];b=pp[control+'_tree_second'];delta=a-b;ci=bootstrap(delta.to_numpy());fo=folds.reindex(pp.index).to_numpy()
            rows.append(dict(population_model='global',pair=pair,stratum=stratum,control=control,n_structures=len(pp),first_mae=a.mean(),second_mae=b.mean(),improvement_e=delta.mean(),ci_low=ci[0],ci_high=ci[1],n_folds_improved=sum(delta.to_numpy()[fo==f].mean()>0 for f in range(5))))
    for pair in PAIRS:
        p=ROOT/'results/metal_control'/(pair+'_oof.parquet')
        if not p.exists():continue
        g=pd.read_parquet(p);a=np.abs(g.observed-g.metal_tree_first);b=np.abs(g.observed-g.metal_tree_second)
        pp=pd.DataFrame(dict(qmof_id=g.qmof_id,first=a,second=b)).groupby('qmof_id').mean();delta=pp['first']-pp['second'];ci=bootstrap(delta.to_numpy());fo=folds.reindex(pp.index).to_numpy()
        rows.append(dict(population_model='metal_only',pair=pair,stratum='metal',control='piecewise_plus_element_residual',n_structures=len(pp),first_mae=pp['first'].mean(),second_mae=pp['second'].mean(),improvement_e=delta.mean(),ci_low=ci[0],ci_high=ci[1],n_folds_improved=sum(delta.to_numpy()[fo==f].mean()>0 for f in range(5))))
    pd.DataFrame(rows).to_csv(OUT/'twohop_increment.csv',index=False)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    flags=['ambiguous','zero_neighbors','metal_c_contact','metal_h_contact','metal_metal_contact','extended_metal_contact','cutoff_disagreement','radius_sensitivity']
    f=pd.read_parquet(ROOT/'data/structural_features.parquet',columns=['qmof_id','atom_index','element','cn','is_metal','first_n_O','first_n_N','first_n_C','first_n_H','first_n_metal','first_formula','neighbor_metal_formula']+flags)
    coverage=[]
    groups=[('all',f),('metal',f[f.is_metal.astype(bool)]),('nonmetal',f[~f.is_metal.astype(bool)])]+list(f.groupby('element'))
    for label,g in groups:
        for flag in flags:coverage.append(dict(stratum=label,flag=flag,n_atoms=len(g),n_structures=g.qmof_id.nunique(),n_flagged_atoms=int(g[flag].sum()),fraction_atoms=g[flag].mean(),n_flagged_structures=g.loc[g[flag].astype(bool),'qmof_id'].nunique()))
    pd.DataFrame(coverage).to_csv(OUT/'connectivity_coverage_by_stratum.csv',index=False)
    del groups
    f['environment']=environment(f);env=f[['qmof_id','atom_index','environment']]
    allmetrics=[];allgroups=[];allfolds=[];allsumm=[];alleffects=[];allcontrasts=[];allmetal=[]
    for pair in PAIRS:
        paths=sorted((ROOT/'results/models').glob(f'oof_fold*_{pair}.parquet'));assert len(paths)==5
        d=pd.concat([pd.read_parquet(p) for p in paths],ignore_index=True);assert not d.duplicated(['qmof_id','atom_index']).any()
        expected=json.loads((ROOT/'data/descriptor_manifest.json').read_text())['valid_graph_atoms']
        assert len(d)==expected and np.isfinite(d[['observed']+MODELS].to_numpy()).all(),'Incomplete/nonfinite held-out predictions'
        m,g,ff=metrics(d,pair);allmetrics+=m;allgroups+=g;allfolds+=ff
        d=d.merge(env,on=['qmof_id','atom_index'],validate='one_to_one');ss,ee,cc=contrasts(d,pair);allsumm+=ss;alleffects+=ee;allcontrasts+=cc
        allmetal+=metal_contrasts(d,pair)
        print(pair,'metrics and contrasts complete',flush=True)
    pd.DataFrame(allmetrics).to_csv(OUT/'grouped_metrics.csv',index=False)
    pd.DataFrame(allgroups).to_parquet(OUT/'structure_losses.parquet',index=False,compression='zstd')
    pd.DataFrame(allfolds).to_csv(OUT/'fold_metrics.csv',index=False)
    pd.DataFrame(allsumm).to_csv(OUT/'environment_summaries.csv',index=False)
    pd.DataFrame(alleffects).to_csv(OUT/'paired_environment_contrasts.csv',index=False)
    pd.DataFrame(allmetal).to_csv(OUT/'metal_environment_contrasts.csv',index=False)
    pd.DataFrame(allcontrasts).to_parquet(OUT/'structure_contrasts.parquet',index=False)
    # Traceable examples selected by median environment atom count per structure, without outcomes.
    t=pd.read_parquet(ROOT/'data/charge_targets.parquet',columns=['qmof_id','atom_index','ddec6','cm5','bader']);examples=[]
    for label,g in f.groupby('environment'):
        if label=='other':continue
        counts=g.groupby('qmof_id').size();q=counts.index[np.argmin(abs(counts-counts.median()))]
        examples.append(g[g.qmof_id==q].iloc[[0]])
    pd.concat(examples).merge(t,on=['qmof_id','atom_index']).to_csv(OUT/'representative_sites.csv',index=False)
    coefficients=[]
    for pair in PAIRS:
        for fold in range(5):
            with (ROOT/'work/models'/f'{pair}_fold{fold}_affine_ridge_first.pkl').open('rb') as stream:artifact=pickle.load(stream)
            pipe=artifact['model'];coef=pipe.named_steps['ridge'].coef_;scale=pipe.named_steps['standardscaler'].scale_
            for feature,c,s in zip(artifact['names'],coef,scale):coefficients.append(dict(pair=pair,fold=fold,feature=feature,coefficient_per_training_sd=c,coefficient_per_original_unit=c/s))
    cf=pd.DataFrame(coefficients);cf.to_csv(OUT/'ridge_coefficients.csv',index=False)
    cf.groupby(['pair','feature']).agg(mean_standardized=('coefficient_per_training_sd','mean'),min_standardized=('coefficient_per_training_sd','min'),max_standardized=('coefficient_per_training_sd','max'),mean_original_unit=('coefficient_per_original_unit','mean')).reset_index().to_csv(OUT/'ridge_coefficient_summary.csv',index=False)
    twohop_increment()
    print(pd.DataFrame(allmetrics).query("stratum in ['all','metal'] and model in ['offset','affine','nonlinear','affine_tree_first','affine_tree_second','nonlinear_tree_second']").to_string(index=False),flush=True)

if __name__=='__main__':main()
