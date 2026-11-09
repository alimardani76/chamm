"""Fixed five-fold structure evaluation; compact, capped training, full OOF evaluation.
All structural matrices are charge-blind. Affine controls are offsets only.
"""
import os
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['MKL_NUM_THREADS']='2'
import gc,json,pickle,time,argparse,re
from pathlib import Path
import numpy as np,pandas as pd,pyarrow.parquet as pq
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from xgboost import XGBRegressor
from chamm_graph import info
from importlib import import_module
charge=import_module('16_charge_analysis');calibration=import_module('16b_calibration_check');PAIRS=charge.PAIRS;SEED=charge.SEED
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/models'
METALS=('Zn','Cu','Co','Ni','Fe','Mn','Zr','Al','Mg','Ca','Cd','Ag','La','Gd','Dy','Pb','Pt','Ba')

def feature_frames(path):
    schema=pq.read_schema(path).names
    central=['atomic_number','group','period','electronegativity','covalent_radius','is_metal']
    first=['cn']+[c for c in schema if c.startswith('first_n_')]+[c for c in schema if c.startswith(('distance_','normalized_distance_','neighbor_X_','neighbor_abs_'))]+['mean_neighbor_metal_Z']
    second=['twohop_paths']+[c for c in schema if c.startswith('twohop_n_')]
    use=central+first+second
    df=pd.read_parquet(path,columns=['qmof_id','atom_index','element','graph_valid','ambiguous','neighbor_metal_formula']+use)
    for col in use:df[col]=df[col].astype('float32')
    # Exact common neighboring metal counts; less common identities are retained in the frozen table.
    cache={}
    for v in df.neighbor_metal_formula.unique():
        counts={el:int(n) for el,n in re.findall(r'([A-Z][a-z]*)(\d+)',v)}
        cache[v]=[counts.get(el,0) for el in METALS]
    arr=np.asarray([cache[v] for v in df.neighbor_metal_formula],dtype=np.float32)
    for k,el in enumerate(METALS):df['neighbor_metal_n_'+el]=arr[:,k];first.append('neighbor_metal_n_'+el)
    # One-hot central identity avoids treating atomic number as the only elemental axis.
    for el in sorted(df.element.unique()):df['central_'+el]=(df.element==el).astype('float32');central.append('central_'+el)
    names=central+first+second
    return df,names,len(central+first)

def training_sample(t,fold):
    tr=t[(t.fold!=fold)&t.graph_valid]
    # Up to 12 atoms per observed structure-element cell, then a global deterministic cap.
    rng=np.random.default_rng(SEED+fold);order=rng.permutation(len(tr));ranked=tr.iloc[order]
    mask=ranked.groupby(['qmof_id','element'],sort=False).cumcount().to_numpy()<12
    candidates=ranked.index.to_numpy()[mask]
    cap=160000
    if len(candidates)>cap:
        # Reserve one observation per available structure-element cell before filling the cap.
        reserved=ranked[mask].drop_duplicates(['qmof_id','element']).index.to_numpy()
        rest=np.setdiff1d(candidates,reserved,assume_unique=True);rng.shuffle(rest)
        candidates=np.concatenate([reserved,rest[:max(0,cap-len(reserved))]])
    chosen=np.sort(candidates)
    cell_total=tr.groupby(['qmof_id','element']).size()
    sample=t.loc[chosen];cell_sample=sample.groupby(['qmof_id','element']).size()
    idx=pd.MultiIndex.from_frame(sample[['qmof_id','element']])
    # Correct capped sampling back to full-training atom proportions, with each structure weight one.
    w=cell_total.reindex(idx).to_numpy()/cell_sample.reindex(idx).to_numpy()*sample.structure_weight.to_numpy()
    w=w/w.mean()
    return chosen,w.astype('float32')

def predict_batches(model,X,ids):
    pred=np.empty(len(ids),dtype=np.float32)
    for k in range(0,len(ids),80000):pred[k:k+80000]=model.predict(X[ids[k:k+80000]])
    return pred

def main():
    OUT.mkdir(parents=True,exist_ok=True);(ROOT/'work/models').mkdir(exist_ok=True)
    cols=['qmof_id','atom_index','element','fold','structure_weight','is_metal','ddec6','cm5','bader']+['delta_'+a+'_'+b for a,b in PAIRS]
    t=pd.read_parquet(ROOT/'data/charge_targets.parquet',columns=cols)
    f,names,nfirst=feature_frames(ROOT/'data/structural_features.parquet')
    assert np.array_equal(t.qmof_id,f.qmof_id) and np.array_equal(t.atom_index,f.atom_index) and np.array_equal(t.element,f.element)
    t['graph_valid']=f.graph_valid.to_numpy();t['ambiguous']=f.ambiguous.to_numpy()
    X=f[names].to_numpy(dtype=np.float32);del f;gc.collect()
    config=dict(seed=SEED,folds=5,training_cap=160000,per_structure_element_cap=12,weighting='capped sample corrected to inverse structure size',
        tree=dict(n_estimators=100,max_depth=4,learning_rate=.08,min_child_weight=20,subsample=1.0,colsample_bytree=1.0,
                  reg_lambda=5.0,tree_method='hist',max_bin=128,n_jobs=2,random_state=SEED,objective='reg:squarederror'),
        linear=dict(model='median imputation, standardized Ridge',alpha=100.),first_features=names[:nfirst],second_features=names[nfirst:],
        no_early_stopping=True,outer_fold_used_for_fit=False,graph_valid_atoms=int(t.graph_valid.sum()),all_atoms=len(t))
    config['descriptor_configuration_id']=json.loads((ROOT/'data/descriptor_manifest.json').read_text())['configuration_id']
    config['descriptor_table_sha256']=json.loads((ROOT/'data/descriptor_manifest.json').read_text())['table_sha256']
    cp=OUT/'model_configuration.json'
    if cp.exists():assert json.loads(cp.read_text())==config,'Refuse to reuse incompatible model partitions'
    else:cp.write_text(json.dumps(config,indent=2))
    gc.collect()
    ladder=['affine_ridge_first']+[c+'_tree_'+l for c in ('offset','affine','nonlinear') for l in ('first','second')]
    (OUT/'evaluated_ladder.json').write_text(json.dumps(dict(models=ladder,resource_choice='One first-shell affine-residual Ridge baseline; tree ladders for all three controls. Compatible first block retained; unused extra Ridge fits not included in scientific comparisons.'),indent=2))
    samples=[]
    for fold in range(5):
        trainids,w=training_sample(t,fold);testids=np.flatnonzero((t.fold.to_numpy()==fold)&t.graph_valid.to_numpy())
        samples.append(pd.DataFrame(dict(qmof_id=t.loc[trainids,'qmof_id'].to_numpy(),atom_index=t.loc[trainids,'atom_index'].to_numpy(),outer_fold=fold,fit_weight=w)))
        for a,bb in PAIRS:
            pair=a+'_'+bb;y=t['delta_'+pair].to_numpy();start=time.monotonic()
            op=OUT/f'oof_fold{fold}_{pair}.parquet'
            if op.exists():print('reuse',op.name,flush=True);continue
            # Refit offset/scale coefficients on this training fold only.
            c=charge.fit_control(t[t.fold!=fold],a,bb)
            off=t.element.map(c.offset).to_numpy();slope=t.element.map(c.slope).to_numpy();inter=t.element.map(c.intercept).to_numpy()
            fallback=np.average(t.loc[t.fold!=fold,'delta_'+pair],weights=t.loc[t.fold!=fold,'structure_weight'])
            off=np.nan_to_num(off,nan=fallback);slope=np.nan_to_num(slope,nan=1.);inter=np.nan_to_num(inter,nan=-fallback)
            affine=(1-slope)*t[a].to_numpy()-inter
            nonlinear,_=calibration.nonlinear_predictions(t[t.fold!=fold],t.reset_index(drop=True),a,bb,c)
            pred=t.loc[testids,['qmof_id','atom_index','element','is_metal','ambiguous']].copy();pred['fold']=fold;pred['observed']=y[testids]
            pred['offset']=off[testids];pred['affine']=affine[testids];pred['nonlinear']=nonlinear[testids]
            pred['baseline_supported']=t.loc[testids,'element'].map(c.n_structures).fillna(0).to_numpy()>=20
            for control,base in [('offset',off),('affine',affine),('nonlinear',nonlinear)]:
                residual=y[trainids]-base[trainids]
                for level,ncols in [('first',nfirst),('second',len(names))]:
                    for kind in (('tree',) if control=='nonlinear' else ('ridge','tree')):
                        if kind=='ridge' and (control!='affine' or level!='first'):continue
                        # The full ladder is small and fixed; no search or outer-test stopping.
                        label=control+'_'+kind+'_'+level
                        model=make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),Ridge(alpha=100.)) if kind=='ridge' else XGBRegressor(**config['tree'])
                        xx=X[trainids,:ncols]
                        if kind=='ridge':model.fit(xx,residual,ridge__sample_weight=w)
                        else:model.fit(xx,residual,sample_weight=w)
                        pp=predict_batches(model,X[:,:ncols],testids);pred[label]=base[testids]+pp
                        mp=ROOT/'work/models'/f'{pair}_fold{fold}_{label}.pkl'
                        with mp.open('wb') as stream:pickle.dump(dict(model=model,names=names[:ncols]),stream,protocol=5)
                        del model,xx;gc.collect()
                print(f'fold {fold} {pair} {control}; {len(trainids)} fitting atoms; {(time.monotonic()-start):.1f}s',flush=True)
            tmp=op.with_suffix('.parquet.partial');pred.to_parquet(tmp,index=False,compression='zstd');tmp.replace(op)
    pd.concat(samples,ignore_index=True).to_parquet(OUT/'training_samples.parquet',index=False)
    print('All 15 pair/fold OOF partitions completed',flush=True)

if __name__=='__main__':main()
