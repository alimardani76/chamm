"""Read-only input audit and representative saved-model/descriptor replay, no refit."""
from pathlib import Path
import sys,json,hashlib,pickle,importlib
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'code'))
import chamm_graph as graph
def main():
 t=pd.read_parquet(R/'data/charge_targets.parquet');f=pd.read_parquet(R/'data/structural_features.parquet')
 assert len(t)==len(f)==1420532 and t.qmof_id.nunique()==9171
 assert t[['qmof_id','atom_index','element']].equals(f[['qmof_id','atom_index','element']])
 for pair in ['ddec6_cm5','ddec6_bader','cm5_bader']:
  a,b=pair.split('_');assert np.max(np.abs(t['delta_'+pair]-(t[a]-t[b])))<1e-12
 assert t.is_metal.sum()==46994
 cid=json.loads((R/'data/descriptor_manifest.json').read_text())['configuration_id'];assert cid=='621045db53a892f1'
 # Replay the accepted descriptors on three preselected diagnostic structures.
 ids=pd.read_csv(R/'results/sensitivity/subset_manifest.csv').qmof_id.iloc[:3]
 max_descriptor=0.
 for q in ids:
  r=json.loads((R/'work/structures'/f'{q}.json').read_text());s=graph.structure(r)
  g=json.loads((R/f'work/graphs_{cid}'/f'{q}.json').read_text());graph.validate(s,g['graph'])
  rebuilt=pd.DataFrame(graph.features(r,s,g['graph'],g['metadata'],graph.cutoff(s)))
  saved=pd.read_parquet(R/f'work/features_{cid}'/f'{q}.parquet')
  for c in rebuilt:
   if rebuilt[c].dtype.kind in 'biuf' and saved[c].dtype.kind in 'biuf':
    assert np.array_equal(np.isnan(rebuilt[c].to_numpy(dtype=float)),np.isnan(saved[c].to_numpy(dtype=float))),c
    err=np.nanmax(np.abs(rebuilt[c].to_numpy(dtype=float)-saved[c].to_numpy(dtype=float))) if rebuilt[c].notna().any() else 0.
    max_descriptor=max(max_descriptor,float(err));assert np.allclose(rebuilt[c],saved[c],atol=1e-10,rtol=1e-10,equal_nan=True),c
   else:assert rebuilt[c].astype(str).equals(saved[c].astype(str)),c
 # Replay one already fitted blocked model on its complete held-out fold.
 close=importlib.import_module('chamm_closeout');X,names,nfirst,center,cfg=close.feature_matrix()
 meta=pd.read_csv(R/'results/closeout_v1/metadata/component_mapping.csv').set_index('qmof_id')
 folds=t.qmof_id.map(meta.blocked_fold).to_numpy();ix=np.flatnonzero(folds==0)
 pair='ddec6_bader';base,_,_,_=close.calibration(t,folds,0,pair)
 with (R/'work/closeout_v1/blocked/ddec6_bader_fold0_conditional_twohop.pkl').open('rb') as src:checkpoint=pickle.load(src)
 pred=base[ix]+close.predict(checkpoint['model'],X,ix,checkpoint['columns'],t.ddec6.to_numpy())
 oof=pd.read_parquet(R/'results/closeout_v1/blocked/oof_fold0_ddec6_bader.parquet')
 assert np.array_equal(oof.qmof_id,t.iloc[ix].qmof_id)
 err=float(np.max(np.abs(pred-oof.conditional_twohop)));assert err<1e-8,err
 result=dict(status='PASS',atoms=1420532,structures=9171,target_pairs=3,descriptor_structures=list(ids),max_descriptor_error=max_descriptor,model='blocked/ddec6_bader/fold0/conditional_twohop',model_atoms=len(ix),max_saved_prediction_error=err,scope='Input/target audit; three accepted-descriptor replays; one saved-model prediction replay. No model training and no full-cohort descriptor rebuild.')
 (R/'reproduce/full_input_check.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
if __name__=='__main__':main()
