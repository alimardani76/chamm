"""Bounded domain control: equal metal-structure weights, fixed tree, no tuning.
Uses current-outer-fold saved calibration coefficients/knots on training atoms;
never uses cross-fitted baselines from a different outer fold for training targets.
"""
import os
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2'
import json,pickle,gc
from pathlib import Path
from importlib import import_module
import numpy as np,pandas as pd
from xgboost import XGBRegressor
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/metal_control'
models=import_module('18_grouped_models');science=import_module('19_results')
PAIRS=science.PAIRS

def calibration(t,pair,fold,params,knots):
    a,b=pair.split('_');out=np.empty(len(t));ps=params[(params.pair==pair)&(params.fold==fold)].set_index('element')
    kk=knots[(knots.pair==pair)&(knots.fold==fold)]
    for el,indices in t.groupby('element').groups.items():
        ix=np.asarray(indices);x=t.loc[ix,a].to_numpy();k=kk[kk.element==el].sort_values('x')
        if len(k)>1:
            xx=k.x.to_numpy();yy=k.y.to_numpy();yp=np.interp(x,xx,yy);lo=x<xx[0];hi=x>xx[-1]
            yp[lo]=yy[0]+(x[lo]-xx[0])*(yy[1]-yy[0])/(xx[1]-xx[0])
            yp[hi]=yy[-1]+(x[hi]-xx[-1])*(yy[-1]-yy[-2])/(xx[-1]-xx[-2])
            out[ix]=x-yp
        elif el in ps.index:out[ix]=(1-ps.loc[el,'slope'])*x-ps.loc[el,'intercept']
        else:out[ix]=np.nan # only unsupported held-out elements; replaced by their saved test-fold control
    return out

def main():
    OUT.mkdir(parents=True,exist_ok=True);(ROOT/'work/metal_models').mkdir(exist_ok=True)
    if not (ROOT/'work/metal_features.parquet').exists():
        pd.read_parquet(ROOT/'data/structural_features.parquet',filters=[('is_metal','==',1)]).to_parquet(ROOT/'work/metal_features.parquet',index=False)
    f,names,nfirst=models.feature_frames(ROOT/'work/metal_features.parquet')
    t=pd.read_parquet(ROOT/'data/charge_targets.parquet',filters=[('is_metal','==',True)]).reset_index(drop=True)
    assert len(t)==46994 and np.array_equal(t.qmof_id,f.qmof_id) and np.array_equal(t.atom_index,f.atom_index)
    t['ambiguous']=f.ambiguous.to_numpy();X=f[names].to_numpy(dtype=np.float32);del f
    t['metal_weight']=1./t.qmof_id.map(t.groupby('qmof_id').size())
    cfg=json.loads((ROOT/'results/models/model_configuration.json').read_text());tree=cfg['tree']
    params=pd.read_csv(ROOT/'results/charge/training_control_parameters.csv');knots=pd.read_csv(ROOT/'results/charge/nonlinear_training_knots.csv')
    config=dict(purpose='Follow-up control for nonmetal-dominated global objective, not a model search',folds='same frozen five QMOF folds',tree=tree,
        fitting='all outer-training metal atoms; one total weight per structure; no sampling',
        residual_baseline='current training-fold piecewise calibration + weighted training metal-element residual mean',
        first_features=names[:nfirst],second_features=names[nfirst:],atoms=len(t),structures=t.qmof_id.nunique(),descriptor_configuration=cfg['descriptor_configuration_id'])
    (OUT/'configuration.json').write_text(json.dumps(config,indent=2));allmetrics=[];sensitivity=[];coefficients=[]
    for pair in PAIRS:
        blocks=[]
        for fold in range(5):
            train=np.flatnonzero(t.fold.to_numpy()!=fold);test=np.flatnonzero(t.fold.to_numpy()==fold)
            y=t['delta_'+pair].to_numpy();base=calibration(t,pair,fold,params,knots)
            held=pd.read_parquet(ROOT/'results/models'/f'oof_fold{fold}_{pair}.parquet',filters=[('is_metal','==',True)],columns=['qmof_id','atom_index','nonlinear'])
            held=t.loc[test,['qmof_id','atom_index']].merge(held,on=['qmof_id','atom_index'],validate='one_to_one')
            assert np.allclose(base[test][np.isfinite(base[test])],held.nonlinear.to_numpy()[np.isfinite(base[test])],atol=1e-9)
            base[test]=held.nonlinear.to_numpy()
            w=t.loc[train,'metal_weight'].to_numpy();w=w/w.mean();bias={}
            for el in t.loc[train,'element'].unique():
                ii=train[t.loc[train,'element'].to_numpy()==el];bias[el]=np.average(y[ii]-base[ii],weights=t.loc[ii,'metal_weight'])
                coefficients.append(dict(pair=pair,fold=fold,element=el,residual_element_offset=bias[el],n_training_structures=t.loc[ii,'qmof_id'].nunique()))
            centered=base+t.element.map(bias).fillna(0).to_numpy();residual=y[train]-centered[train]
            pred=t.loc[test,['qmof_id','atom_index','element','fold','ambiguous']].copy();pred['observed']=y[test];pred['nonlinear']=base[test];pred['element_residual']=centered[test]
            for level,ncols in [('first',nfirst),('second',len(names))]:
                model=XGBRegressor(**tree);model.fit(X[train,:ncols],residual,sample_weight=w)
                pred['metal_tree_'+level]=centered[test]+model.predict(X[test,:ncols])
                with (ROOT/'work/metal_models'/f'{pair}_fold{fold}_{level}.pkl').open('wb') as stream:pickle.dump(dict(model=model,names=names[:ncols]),stream,protocol=5)
            blocks.append(pred);print('metal control',pair,fold,'fit atoms',len(train),flush=True)
        d=pd.concat(blocks,ignore_index=True);assert len(d)==46994 and not d.duplicated(['qmof_id','atom_index']).any()
        d.to_parquet(OUT/(pair+'_oof.parquet'),index=False)
        for stratum in ['all_metals']+[el for el,n in d.groupby('element').qmof_id.nunique().items() if n>=30]:
            dd=d if stratum=='all_metals' else d[d.element==stratum]
            loss=pd.DataFrame({'qmof_id':dd.qmof_id,'fold':dd.fold})
            for label in ['nonlinear','element_residual','metal_tree_first','metal_tree_second']:loss[label]=np.abs(dd.observed-dd[label])
            st=loss.groupby('qmof_id')[['nonlinear','element_residual','metal_tree_first','metal_tree_second']].mean();fo=loss.groupby('qmof_id').fold.first().to_numpy()
            for label in ['nonlinear','element_residual','metal_tree_first','metal_tree_second']:
                v=st[label].to_numpy();im=st.element_residual.to_numpy()-v;ci=science.bootstrap(im);errci=science.bootstrap(v)
                allmetrics.append(dict(pair=pair,stratum=stratum,model=label,n_atoms=len(dd),n_structures=len(st),mae_atom=loss[label].mean(),mae_structure=v.mean(),mae_ci_low=errci[0],mae_ci_high=errci[1],improvement_structure=im.mean(),ci_low=ci[0],ci_high=ci[1],n_folds_improved=sum(im[fo==fold].mean()>0 for fold in range(5))))
        # Fixed metal-domain models evaluated on exactly the same outcome-blind sensitivity subset.
        for variant in ('shift','cutoff'):
            af,_,_=models.feature_frames(ROOT/'results/sensitivity'/('features_'+variant+'.parquet'));af=af[af.is_metal.astype(bool)]
            dd=d.merge(af,on=['qmof_id','atom_index'],suffixes=('','_alt'),validate='one_to_one');altpred=dd[['qmof_id','atom_index','fold','observed','element_residual']].copy()
            for level in ('first','second'):
                altpred['metal_tree_'+level]=np.nan
                for fold in range(5):
                    mask=dd.fold.eq(fold)
                    with (ROOT/'work/metal_models'/f'{pair}_fold{fold}_{level}.pkl').open('rb') as stream:artifact=pickle.load(stream)
                    xx=dd.reindex(columns=artifact['names'],fill_value=0).loc[mask].to_numpy(dtype=np.float32)
                    altpred.loc[mask,'metal_tree_'+level]=dd.loc[mask,'element_residual'].to_numpy()+artifact['model'].predict(xx)
                e0=np.abs(dd.observed-dd['metal_tree_'+level]);e1=np.abs(altpred.observed-altpred['metal_tree_'+level]);eb=np.abs(dd.observed-dd.element_residual)
                agg=pd.DataFrame(dict(qmof_id=dd.qmof_id,accepted=e0,alternative=e1,baseline=eb)).groupby('qmof_id').mean()
                delta=agg.alternative-agg.accepted;ci=science.bootstrap(delta.to_numpy())
                sensitivity.append(dict(pair=pair,variant=variant,model='metal_tree_'+level,n_atoms=len(dd),n_structures=len(agg),accepted_mae_structure=agg.accepted.mean(),alternative_mae_structure=agg.alternative.mean(),baseline_mae_structure=agg.baseline.mean(),alternative_minus_accepted=delta.mean(),ci_low=ci[0],ci_high=ci[1]))
            altpred.to_parquet(OUT/f'{pair}_{variant}_sensitivity_oof.parquet',index=False)
    pd.DataFrame(allmetrics).to_csv(OUT/'grouped_metrics.csv',index=False);pd.DataFrame(sensitivity).to_csv(OUT/'graph_sensitivity.csv',index=False)
    pd.DataFrame(coefficients).to_csv(OUT/'training_residual_element_offsets.csv',index=False)
    science.twohop_increment()
    print(pd.DataFrame(allmetrics).query("stratum=='all_metals'").to_string(index=False),flush=True)

if __name__=='__main__':main()
