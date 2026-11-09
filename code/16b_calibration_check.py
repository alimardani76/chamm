"""Bounded check of affine adequacy: 12-bin training-only piecewise-linear control."""
import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np,pandas as pd
from importlib import import_module
base=import_module('16_charge_analysis');ROOT=Path(__file__).resolve().parents[1]

def nonlinear_predictions(train,test,a,b,control):
    out=np.empty(len(test));knots=[]
    # Test dataframe uses a fresh contiguous index. Quantile binning is fit only on training charges.
    for el,te in test.groupby('element',observed=True):
        tr=train[train.element==el];c=control.loc[el] if el in control.index else None
        x=tr[a].to_numpy();y=tr[b].to_numpy();w=tr.structure_weight.to_numpy();xt=te[a].to_numpy()
        if c is None:yp=xt-np.average(train[a]-train[b],weights=train.structure_weight)
        elif tr.qmof_id.nunique()<20 or len(tr)<100 or len(np.unique(x))<20:yp=c.intercept+c.slope*xt
        else:
            bounds=np.unique(np.quantile(x,np.linspace(0,1,13)))
            bins=np.searchsorted(bounds[1:-1],x,side='right');xx=[];yy=[]
            for k in range(len(bounds)-1):
                mask=bins==k
                if mask.any():xx.append(np.average(x[mask],weights=w[mask]));yy.append(np.average(y[mask],weights=w[mask]))
            xx=np.array(xx);yy=np.array(yy);yp=np.interp(xt,xx,yy)
            if len(xx)>1:
                low=xt<xx[0];high=xt>xx[-1]
                yp[low]=yy[0]+(xt[low]-xx[0])*(yy[1]-yy[0])/(xx[1]-xx[0])
                yp[high]=yy[-1]+(xt[high]-xx[-1])*(yy[-1]-yy[-2])/(xx[-1]-xx[-2])
            knots.extend(dict(element=el,pair=a+'_'+b,x=float(x0),y=float(y0)) for x0,y0 in zip(xx,yy))
        out[te.index]=xt-yp
    return out,knots

def main():
    cols=['qmof_id','atom_index','element','fold','structure_weight','ddec6','cm5','bader']
    t=pd.read_parquet(ROOT/'data/charge_targets.parquet',columns=cols);pred=[];knots=[]
    for fold in range(5):
        train=t[t.fold!=fold];test=t[t.fold==fold].reset_index(drop=True);p=test[['qmof_id','atom_index']].copy()
        for a,b in base.PAIRS:
            control=base.fit_control(train,a,b);v,k=nonlinear_predictions(train,test,a,b,control)
            p['nonlinear_'+a+'_'+b]=v
            for row in k:row['fold']=fold
            knots.extend(k)
        pred.append(p);print('nonlinear fold',fold,flush=True)
    p=pd.concat(pred);p.to_parquet(ROOT/'results/charge/nonlinear_oof.parquet',index=False)
    pd.DataFrame(knots).to_csv(ROOT/'results/charge/nonlinear_training_knots.csv',index=False)
    d=t.merge(pd.read_parquet(ROOT/'results/charge/baseline_oof.parquet'),on=['qmof_id','atom_index']).merge(p,on=['qmof_id','atom_index'])
    rows=[]
    for el,g in d.groupby('element'):
        for a,b in base.PAIRS:
            pair=a+'_'+b;bins=pd.qcut(g[a],10,duplicates='drop')
            for interval,gg in g.groupby(bins,observed=True):
                delta=gg[a]-gg[b]
                rows.append(dict(element=el,pair=pair,qA_mean=gg[a].mean(),n_atoms=len(gg),n_structures=gg.qmof_id.nunique(),
                    affine_signed_residual=(delta-gg['affine_'+pair]).mean(),nonlinear_signed_residual=(delta-gg['nonlinear_'+pair]).mean()))
    pd.DataFrame(rows).to_csv(ROOT/'results/charge/calibration_residual_deciles.csv',index=False)
    print('control MAE',[(a+'_'+b,float(np.abs(d[a]-d[b]-d['nonlinear_'+a+'_'+b]).mean())) for a,b in base.PAIRS],flush=True)

if __name__=='__main__':main()
