"""Full descriptive cohort, fixed grouped folds, training-only offset/affine controls."""
import os
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2'
import json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from chamm_graph import info
ROOT=Path(__file__).resolve().parents[1]
PAIRS=[('ddec6','cm5'),('ddec6','bader'),('cm5','bader')]
SEED=20261109

def fit_control(train,a,b):
    rows=[]
    for el,g in train.groupby('element',observed=True):
        x=g[a].to_numpy();y=g[b].to_numpy();w=g['structure_weight'].to_numpy();w=w/w.sum()
        xm=w@x;ym=w@y;var=w@((x-xm)**2);cov=w@((x-xm)*(y-ym))
        nstruct=g.qmof_id.nunique();slope=cov/var if var>1e-6 and nstruct>=20 else 1.
        rows.append(dict(element=el,pair=a+'_'+b,offset=w@(x-y),intercept=ym-slope*xm,slope=slope,
                         n_atoms=len(g),n_structures=nstruct,x_variance=var,affine_estimable=bool(var>1e-6 and nstruct>=20)))
    return pd.DataFrame(rows).set_index('element')

def descriptive_distributions(t=None):
    if t is None:t=pd.read_parquet(ROOT/'data/charge_targets.parquet')
    observables=['ddec6','cm5','bader','charge_range']+['delta_'+a+'_'+b for a,b in PAIRS]+['abs_'+a+'_'+b for a,b in PAIRS]
    rows=[]
    groups=[('all',t),('metal',t[t.is_metal]),('nonmetal',t[~t.is_metal])]+list(t.groupby('element'))
    for label,g in groups:
        perstructure=g.groupby('qmof_id')[observables].mean()
        for col in observables:
            v=g[col].to_numpy();sv=perstructure[col].to_numpy()
            rows.append(dict(stratum=label,observable=col,n_atoms=len(g),n_structures=len(perstructure),mean_atom=v.mean(),median_atom=np.median(v),p05_atom=np.quantile(v,.05),p95_atom=np.quantile(v,.95),mean_structure=sv.mean(),median_structure_mean=np.median(sv),p05_structure_mean=np.quantile(sv,.05),p95_structure_mean=np.quantile(sv,.95)))
    pd.DataFrame(rows).to_csv(ROOT/'results/charge/full_descriptive_distributions.csv',index=False)
    # Element distributions use equal QMOF means within each element, not element fractions per whole structure.
    p=ROOT/'results/charge/element_distributions.csv'
    if p.exists():
        ed=pd.read_csv(p)
        for el,g in t.groupby('element'):
            for a,b in PAIRS:
                pair=a+'_'+b;sm=g.groupby('qmof_id')[['delta_'+pair,'abs_'+pair]].mean();mask=ed.element.eq(el)&ed.pair.eq(pair)
                ed.loc[mask,'structure_balanced_mean_signed']=sm['delta_'+pair].mean();ed.loc[mask,'structure_balanced_mean_absolute']=sm['abs_'+pair].mean()
        ed.to_csv(p,index=False)

def main():
    (ROOT/'results/charge').mkdir(parents=True,exist_ok=True)
    t=pd.read_csv(ROOT/'data/chamm_atom_charges_recovered.csv.gz')
    cohort=pd.read_csv(ROOT/'data/chamm_common_cohort.csv');assert len(t)==1420532 and t.qmof_id.nunique()==9171
    assert np.isfinite(t[['ddec6','cm5','bader']].to_numpy()).all()
    assert t.groupby('qmof_id').size().sort_index().equals(cohort.set_index('qmof_id').n_atoms.sort_index())
    assert not t.duplicated(['qmof_id','atom_index']).any()
    refs=json.loads((ROOT/'inputs/expanded_benchmark_records.json').read_text());q=refs[0]['qmof_id']
    r=refs[0];g=t[t.qmof_id==q].sort_values('atom_index');assert len(g)==len(r['structure']['sites'])
    direct_sites=[0,len(g)//2,len(g)-1]
    for i in direct_sites:
        s=r['structure']['sites'][i];assert g.iloc[i].element==s['species'][0]['element']
        assert np.allclose(g.iloc[i][['ddec6','cm5','bader']].to_numpy(dtype=float),[s['properties'][k] for k in ('pbe_ddec_charge','pbe_cm5_charge','pbe_bader_charge')],atol=1e-12)
    t['is_metal']=t.element.map(lambda e:info(e)['metal'])
    ids=sorted(t.qmof_id.unique());rng=np.random.default_rng(SEED);rng.shuffle(ids)
    folds=pd.DataFrame({'qmof_id':ids,'fold':np.arange(len(ids))%5}).sort_values('qmof_id')
    folds.to_csv(ROOT/'results/grouped_splits.csv',index=False);t['fold']=t.qmof_id.map(folds.set_index('qmof_id').fold)
    sizes=t.groupby('qmof_id').size();t['structure_weight']=1./t.qmof_id.map(sizes)
    a=t[['ddec6','cm5','bader']].to_numpy();t['charge_range']=np.ptp(a,axis=1)
    t['near_zero_any_005']=(np.abs(a)<.05).any(axis=1)
    for eps in (0.,.01,.05,.10):t[f'sign_opposition_{eps:.2f}']=(a.min(axis=1)<-eps)&(a.max(axis=1)>eps)
    for x,y in PAIRS:
        k=x+'_'+y;t['delta_'+k]=t[x]-t[y];t['abs_'+k]=np.abs(t['delta_'+k])
    summaries=[];strata=[]
    for el,g in t.groupby('element',observed=True):
        w=g.structure_weight.to_numpy();w=w/w.sum()
        for x,y in PAIRS:
            d=g['delta_'+x+'_'+y].to_numpy()
            summaries.append(dict(element=el,is_metal=info(el)['metal'],pair=x+'_'+y,n_atoms=len(g),n_structures=g.qmof_id.nunique(),
                mean_signed=d.mean(),median_signed=np.median(d),mean_absolute=np.abs(d).mean(),p05=np.quantile(d,.05),p95=np.quantile(d,.95),
                structure_balanced_mean_signed=g.groupby('qmof_id')['delta_'+x+'_'+y].mean().mean(),structure_balanced_mean_absolute=g.groupby('qmof_id')['abs_'+x+'_'+y].mean().mean(),
                sign_opposition_005=g['sign_opposition_0.05'].mean(),near_zero_any_005=g.near_zero_any_005.mean()))
    for label,g in [('all',t),('metal',t[t.is_metal]),('nonmetal',t[~t.is_metal])]:
        row=dict(stratum=label,n_atoms=len(g),n_structures=g.qmof_id.nunique(),charge_range_mean=g.charge_range.mean(),charge_range_median=g.charge_range.median())
        for eps in (0.,.01,.05,.10):row[f'sign_opposition_{eps:.2f}']=g[f'sign_opposition_{eps:.2f}'].mean()
        strata.append(row)
    pd.DataFrame(summaries).to_csv(ROOT/'results/charge/element_distributions.csv',index=False)
    pd.DataFrame(strata).to_csv(ROOT/'results/charge/strata_summary.csv',index=False)
    descriptive_distributions(t)
    baseline=[];fitrows=[]
    for fold in range(5):
        train=t[t.fold!=fold];test=t[t.fold==fold];p=test[['qmof_id','atom_index']].copy()
        for x,y in PAIRS:
            k=x+'_'+y;c=fit_control(train,x,y);c['fold']=fold;fitrows.append(c.reset_index())
            fallback=np.average(train[x]-train[y],weights=train.structure_weight)
            p['offset_'+k]=test.element.map(c.offset).fillna(fallback).to_numpy()
            slope=test.element.map(c.slope).fillna(1).to_numpy();intercept=test.element.map(c.intercept).fillna(-fallback).to_numpy()
            p['affine_'+k]=(1-slope)*test[x].to_numpy()-intercept
            p['baseline_supported_'+k]=test.element.map(c.n_structures).fillna(0).to_numpy()>=20
        baseline.append(p);print('charge baseline fold',fold,flush=True)
    pd.concat(baseline).to_parquet(ROOT/'results/charge/baseline_oof.parquet',index=False)
    pd.concat(fitrows).to_csv(ROOT/'results/charge/training_control_parameters.csv',index=False)
    t.to_parquet(ROOT/'data/charge_targets.parquet',index=False)
    identity=dict(structures=9171,atoms=1420532,finite=True,unique_keys=True,direct_check=dict(qmof_id=q,sites=direct_sites),
        seed=SEED,folds=5,weights='each structure total weight one; element controls normalize this weight within element',
        sign_definition='opposite charges beyond +/-0.05 e, with 0,0.01,0.10 sensitivity; near-zero is not opposition',
        algebra='delta_DC + delta_CB = delta_DB; range=max pairwise absolute difference')
    (ROOT/'results/charge/cohort_confirmation.json').write_text(json.dumps(identity,indent=2))
    print(pd.DataFrame(strata).to_string(index=False),flush=True)

if __name__=='__main__':main()
