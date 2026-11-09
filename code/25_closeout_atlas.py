"""Paired conditional metrics, component-bootstrap challenge and signed/magnitude atlas."""
import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
import gc,json
import numpy as np,pandas as pd
from chamm_closeout import *
DEST=OUT/'atlas'
COMPARISONS=[('piecewise','center','element_recentering'),('center','charge_only','added_q_A'),('center','conditional_first','joint_q_A_and_first'),('charge_only','conditional_first','first_given_q_A'),('conditional_first','conditional_twohop','added_twohop_given_q_A_first'),('charge_only','conditional_twohop','structure_given_q_A')]

def bootstrap_matrix(values,n_atoms,blocks):
    """One joint draw sequence for all competitors and paired increments."""
    values=np.asarray(values,float);n_atoms=np.asarray(n_atoms,float);labels,inv=np.unique(np.asarray(blocks),return_inverse=True);k=len(labels)
    ss=np.zeros((k,values.shape[1]));aa=ss.copy();ns=np.bincount(inv).astype(float);na=np.bincount(inv,weights=n_atoms)
    np.add.at(ss,inv,values);np.add.at(aa,inv,values*n_atoms[:,None])
    rng=np.random.default_rng(SEED);bs=np.empty((1000,values.shape[1]));ba=bs.copy()
    for j in range(1000):
        ix=rng.integers(0,k,k);bs[j]=ss[ix].sum(axis=0)/ns[ix].sum();ba[j]=aa[ix].sum(axis=0)/na[ix].sum()
    return bs,ba,k

def metrics_for_units(u,split,pair):
    records=[];increments=[];folds=[];fi=[]
    selections={'all':u.population.eq('all'),'metal':u.population.eq('metal'),'nonmetal':u.population.eq('nonmetal'),
                'has_mofkey':u.population.eq('all')&u.has_mofkey,'nonsingleton':u.population.eq('all')&(u.component_size>1),
                'has_mofkey_metal':u.population.eq('metal')&u.has_mofkey,'nonsingleton_metal':u.population.eq('metal')&(u.component_size>1),
                'has_mofkey_nonmetal':u.population.eq('nonmetal')&u.has_mofkey,'nonsingleton_nonmetal':u.population.eq('nonmetal')&(u.component_size>1)}
    for el in sorted(set(u.population)-{'all','metal','nonmetal'}):selections[el]=u.population.eq(el)
    for stratum,mask in selections.items():
        g=u[mask]
        if g.empty:continue
        vals=g[list(MODELS)].to_numpy();blocks=g.qmof_id if split=='original' else g.component
        bs,ba,nb=bootstrap_matrix(vals,g.n_atoms,blocks);ci=np.quantile(bs,[.025,.975],axis=0);aci=np.quantile(ba,[.025,.975],axis=0)
        for j,m in enumerate(MODELS):
            records.append(dict(split=split,pair=pair,stratum=stratum,model=m,n_atoms=int(g.n_atoms.sum()),n_structures=len(g),n_validation_blocks=nb,
                mae_atom=np.average(vals[:,j],weights=g.n_atoms),mae_structure=vals[:,j].mean(),ci_low=ci[0,j],ci_high=ci[1,j],atom_ci_low=aci[0,j],atom_ci_high=aci[1,j],
                calibration_unsupported_atoms=int(g.unsupported.sum()),center_unseen_atoms=int(g.unseen.sum())))
            for fold in range(5):
                h=g[g.fold==fold]
                if len(h):folds.append(dict(split=split,pair=pair,stratum=stratum,model=m,fold=fold,n_atoms=int(h.n_atoms.sum()),n_structures=len(h),mae_structure=h[m].mean(),mae_atom=np.average(h[m],weights=h.n_atoms)))
        for baseline,model,label in COMPARISONS:
            ib=MODELS.index(baseline);im=MODELS.index(model);v=vals[:,ib]-vals[:,im];bc=np.quantile(bs[:,ib]-bs[:,im],[.025,.975]);ac=np.quantile(ba[:,ib]-ba[:,im],[.025,.975])
            increments.append(dict(split=split,pair=pair,stratum=stratum,comparison=label,baseline=baseline,model=model,n_atoms=int(g.n_atoms.sum()),n_structures=len(g),n_validation_blocks=nb,
                improvement_structure=v.mean(),ci_low=bc[0],ci_high=bc[1],improvement_atom=np.average(v,weights=g.n_atoms),atom_ci_low=ac[0],atom_ci_high=ac[1],n_folds_improved=sum(g.loc[g.fold==f,baseline].mean()>g.loc[g.fold==f,model].mean() for f in range(5)),
                interpretation='joint q_A and neighborhood addition, not isolated chemistry' if label=='joint_q_A_and_first' else 'nested fixed-model comparison'))
            for fold in range(5):
                h=g[g.fold==fold]
                if len(h):fi.append(dict(split=split,pair=pair,stratum=stratum,comparison=label,baseline=baseline,model=model,fold=fold,n_structures=len(h),improvement_structure=(h[baseline]-h[model]).mean(),improvement_atom=np.average(h[baseline]-h[model],weights=h.n_atoms)))
    return records,increments,folds,fi

def loss_units(d,meta):
    df=d[['qmof_id','element','is_metal']].copy()
    for m in MODELS:df[m]=np.abs(d.observed-d[m]).to_numpy()
    df['n_atoms']=1;df['unsupported']=(~d.calibration_supported).astype(int).to_numpy();df['unseen']=(~d.center_seen).astype(int).to_numpy()
    cells=df.groupby(['qmof_id','element','is_metal'],observed=True)[list(MODELS)+['n_atoms','unsupported','unseen']].sum().reset_index();del df
    subgroups=[('all',cells),('metal',cells[cells.is_metal]),('nonmetal',cells[~cells.is_metal])]
    subgroups += [(el,g) for el,g in cells.groupby('element') if g.qmof_id.nunique()>=50]
    parts=[]
    for label,g in subgroups:
        u=g.groupby('qmof_id')[list(MODELS)+['n_atoms','unsupported','unseen']].sum()
        for m in MODELS:u[m]/=u.n_atoms
        u['population']=label;parts.append(u.reset_index())
    u=pd.concat(parts,ignore_index=True).merge(meta[['component','component_size','has_mofkey','old_fold','blocked_fold']],left_on='qmof_id',right_index=True,validate='many_to_one')
    return u

def atlas_summaries(units,vectors,split,pair,meta):
    rows=[];effects=[]
    for (el,env),g in units.groupby(['element','environment']):
        blocks=g.qmof_id if split=='original' else g.qmof_id.map(meta.component)
        bs,ba,nb=bootstrap_matrix(g[list(QUANTITIES)],g.n_atoms,blocks);ci=np.quantile(bs,[.025,.975],axis=0)
        for j,q in enumerate(QUANTITIES):
            rows.append(dict(split=split,pair=pair,element=el,environment=env,quantity=q,n_atoms=int(g.n_atoms.sum()),n_structures=len(g),n_validation_blocks=nb,
                mean_structure=g[q].mean(),mean_atom=np.average(g[q],weights=g.n_atoms),ci_low=ci[0,j],ci_high=ci[1,j],ambiguity_fraction_structure=g.ambiguous_fraction.mean(),ambiguity_fraction_atom=np.average(g.ambiguous_fraction,weights=g.n_atoms),adequate_support=len(g)>=30,
                meaning={'raw_signed':'signed scheme difference','raw_absolute':'raw disagreement magnitude','calibrated_signed':'directional departure from stronger charge-only model','calibrated_absolute':'absolute departure from stronger charge-only model','residual_sd':'within-structure/environment residual spread; not physical uncertainty'}[q]))
    for (left,right,q),g in vectors.groupby(['left','right','quantity']):
        blocks=g.qmof_id if split=='original' else g.qmof_id.map(meta.component);ci=block_ci(g.effect,blocks)
        effects.append(dict(split=split,pair=pair,left=left,right=right,quantity=q,n_structures=len(g),n_validation_blocks=pd.Series(blocks).nunique(),left_atoms=int(g.left_atoms.sum()),right_atoms=int(g.right_atoms.sum()),effect=g.effect.mean(),ci_low=ci[0],ci_high=ci[1],left_ambiguity=g.left_ambiguity.mean(),right_ambiguity=g.right_ambiguity.mean(),adequate_support=len(g)>=30))
    return rows,effects

def grouping_challenge(losses):
    rows=[]
    for pair in PAIRS:
        p=losses[losses.pair.eq(pair)&losses.population.eq('all')]
        merged=p[p.split.eq('original')].merge(p[p.split.eq('blocked')],on='qmof_id',suffixes=('_original','_blocked'),validate='one_to_one')
        for stratum,mask in [('all',np.ones(len(merged),bool)),('has_mofkey',merged.has_mofkey_blocked),('nonsingleton',merged.component_size_blocked>1)]:
            g=merged[mask]
            for label in ('charge_only','conditional_first','conditional_twohop','structure_given_q_A'):
                if label=='structure_given_q_A':
                    original=g.charge_only_original-g.conditional_twohop_original;blocked=g.charge_only_blocked-g.conditional_twohop_blocked;meaning='blocked minus original conditional gain'
                else:original=g[label+'_original'];blocked=g[label+'_blocked'];meaning='blocked minus original MAE'
                delta=blocked-original;ci=block_ci(delta,g.component_blocked)
                rows.append(dict(pair=pair,stratum=stratum,quantity=label,n_structures=len(g),n_validation_blocks=g.component_blocked.nunique(),original=original.mean(),blocked=blocked.mean(),change=delta.mean(),ci_low=ci[0],ci_high=ci[1],meaning=meaning))
    pd.DataFrame(rows).to_csv(DEST/'grouping_challenge.csv',index=False)

def main():
    DEST.mkdir(parents=True,exist_ok=True);meta=pd.read_csv(OUT/'metadata/component_mapping.csv',keep_default_na=False).set_index('qmof_id')
    fcols=['qmof_id','atom_index','element','cn','is_metal','ambiguous','first_n_O','first_n_N','first_n_C','first_n_H','first_n_metal']
    features=pd.read_parquet(ROOT/'data/structural_features.parquet',columns=fcols)
    metrics=[];increments=[];folds=[];foldincrements=[];losses=[];atlas=[];contrasts=[];envunits=[];vectors=[]
    for split in ('original','blocked'):
        assert (OUT/split/'COMPLETE.json').exists()
        for pair in PAIRS:
            d=pd.concat([pd.read_parquet(p) for p in sorted((OUT/split).glob(f'oof_fold*_{pair}.parquet'))],ignore_index=True)
            assert len(d)==1420532 and not d.duplicated(['qmof_id','atom_index']).any() and np.isfinite(d[list(MODELS)].to_numpy()).all()
            u=loss_units(d,meta);u['fold']=u.old_fold if split=='original' else u.blocked_fold;u['split']=split;u['pair']=pair
            a,b,c,e=metrics_for_units(u,split,pair);metrics+=a;increments+=b;folds+=c;foldincrements+=e;losses.append(u)
            eu=environment_units(features,d);v=paired_vectors(eu);ar,cr=atlas_summaries(eu,v,split,pair,meta);atlas+=ar;contrasts+=cr
            eu['split']=split;eu['pair']=pair;v['split']=split;v['pair']=pair;envunits.append(eu);vectors.append(v)
            print('closeout metrics/atlas',split,pair,flush=True);del d;gc.collect()
    pd.DataFrame(metrics).to_csv(DEST/'grouped_metrics.csv',index=False);pd.DataFrame(increments).to_csv(DEST/'paired_increments.csv',index=False);pd.DataFrame(folds).to_csv(DEST/'fold_metrics.csv',index=False);pd.DataFrame(foldincrements).to_csv(DEST/'fold_increments.csv',index=False)
    lu=pd.concat(losses,ignore_index=True);lu.to_parquet(DEST/'structure_losses.parquet',index=False,compression='zstd');grouping_challenge(lu)
    pd.DataFrame(atlas).to_csv(DEST/'environment_atlas.csv',index=False);pd.DataFrame(contrasts).to_csv(DEST/'headline_contrasts.csv',index=False)
    pd.concat(envunits,ignore_index=True).to_parquet(DEST/'environment_units.parquet',index=False,compression='zstd');pd.concat(vectors,ignore_index=True).to_parquet(DEST/'contrast_vectors.parquet',index=False,compression='zstd')
    (DEST/'COMPLETE.json').write_text(json.dumps(dict(bootstrap_replicates=1000,original_unit='QMOF',blocked_unit='identifier/composition component; ratio of QMOF sums/counts',conditional_training='same sampled atoms/weights for competitors',all_atoms_each_pair=1420532),indent=2));print('Atlas and conditional comparisons complete',flush=True)

def verify_compact():
    """Reaggregate saved numerical review inputs; does not regenerate models or graphs."""
    lu=pd.read_parquet(DEST/'structure_losses.parquet');metrics=pd.read_csv(DEST/'grouped_metrics.csv');increments=pd.read_csv(DEST/'paired_increments.csv');err=0.;nmetric=0
    groups={(s,p):g for (s,p),g in lu.groupby(['split','pair'])}
    def select(split,pair,stratum):
        g=groups[(split,pair)]
        if stratum.startswith('has_mofkey'):
            pop='metal' if stratum.endswith('_metal') else 'nonmetal' if stratum.endswith('_nonmetal') else 'all';return g[g.population.eq(pop)&g.has_mofkey]
        if stratum.startswith('nonsingleton'):
            pop='metal' if stratum.endswith('_metal') else 'nonmetal' if stratum.endswith('_nonmetal') else 'all';return g[g.population.eq(pop)&(g.component_size>1)]
        return g[g.population.eq(stratum)]
    for r in metrics.itertuples():
        g=select(r.split,r.pair,r.stratum);err=max(err,abs(g[r.model].mean()-r.mae_structure),abs(np.average(g[r.model],weights=g.n_atoms)-r.mae_atom));assert len(g)==r.n_structures and int(g.n_atoms.sum())==r.n_atoms;nmetric+=1
    for r in increments.itertuples():
        g=select(r.split,r.pair,r.stratum);delta=g[r.baseline]-g[r.model];err=max(err,abs(delta.mean()-r.improvement_structure),abs(np.average(delta,weights=g.n_atoms)-r.improvement_atom))
    booterr=0.
    for split in ['original','blocked']:
        for pair in PAIRS:
            g=select(split,pair,'all');delta=g.charge_only-g.conditional_twohop;ci=block_ci(delta,g.qmof_id if split=='original' else g.component)
            r=increments[(increments.split==split)&increments.pair.eq(pair)&increments.stratum.eq('all')&increments.comparison.eq('structure_given_q_A')].iloc[0];booterr=max(booterr,abs(ci[0]-r.ci_low),abs(ci[1]-r.ci_high))
    v=pd.read_parquet(DEST/'contrast_vectors.parquet');c=pd.read_csv(DEST/'headline_contrasts.csv');sm=v.groupby(['split','pair','left','right','quantity']).effect.agg(['mean','size']).reset_index();check=c.merge(sm,on=['split','pair','left','right','quantity'],validate='one_to_one');assert len(check)==len(c) and np.array_equal(check.n_structures,check['size']);contrast_err=float(np.max(np.abs(check.effect-check['mean'])))
    eu=pd.read_parquet(DEST/'environment_units.parquet');a=pd.read_csv(DEST/'environment_atlas.csv');means=eu.groupby(['split','pair','element','environment'])[list(QUANTITIES)].mean().reset_index().melt(id_vars=['split','pair','element','environment'],var_name='quantity',value_name='reaggregated');check=a.merge(means,on=['split','pair','element','environment','quantity'],validate='one_to_one');assert len(check)==len(a);atlaserr=float(np.max(np.abs(check.mean_structure-check.reaggregated)))
    meta=pd.read_csv(OUT/'metadata/component_mapping.csv').set_index('qmof_id');assert meta.groupby('component').blocked_fold.nunique().max()==1
    for split in ['original','blocked']:
        samples=pd.read_parquet(OUT/split/'training_samples.parquet');expected=samples.qmof_id.map(meta.old_fold if split=='original' else meta.blocked_fold);assert not (expected.to_numpy()==samples.outer_fold.to_numpy()).any();assert not samples.duplicated(['qmof_id','atom_index','outer_fold']).any()
    legacy=ROOT/'results/models/training_samples.parquet'
    oldmatch='not checked: original full-project sample file absent'
    if legacy.exists():assert pd.read_parquet(OUT/'original/training_samples.parquet').equals(pd.read_parquet(legacy));oldmatch=True
    assert max(err,booterr,contrast_err,atlaserr)<1e-12
    report=dict(status='PASS',model_summary_rows=nmetric,increment_rows=len(increments),headline_contrast_rows=len(c),atlas_rows=len(a),max_metric_error=err,max_principal_bootstrap_error=booterr,max_contrast_mean_error=contrast_err,max_atlas_mean_error=atlaserr,blocked_component_leakage=0,training_test_overlap=0,original_training_samples_exact=oldmatch,note='reaggregates compact numerical inputs; full model regeneration/target checks remain documented separately')
    (OUT/'numeric_reconciliation.json').write_text(json.dumps(report,indent=2));print(report,flush=True)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--verify-compact',action='store_true');a=p.parse_args()
    verify_compact() if a.verify_compact else main()
