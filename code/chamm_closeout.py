"""Finite closeout helpers. Frozen geometry; explicitly charge-conditioned models only."""
import os
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['MKL_NUM_THREADS']='2'
import gc,hashlib,json,warnings
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/closeout_v1';WORK=ROOT/'work/closeout_v1'
PAIRS=('ddec6_cm5','ddec6_bader','cm5_bader');MODELS=('piecewise','center','charge_only','conditional_first','conditional_twohop')
SEED=20261109
CONTRASTS=[('C:O-attached','C:C/H-only'),('C:N-attached','C:C/H-only'),('H:O-bound','H:C-bound'),('H:N-bound','H:C-bound'),('O:metal-adjacent','O:no-metal'),('N:metal-adjacent','N:no-metal')]

def load_targets():
    cols=['qmof_id','atom_index','element','ddec6','cm5','bader','is_metal','fold','structure_weight']+['delta_'+p for p in PAIRS]
    return pd.read_parquet(ROOT/'data/charge_targets.parquet',columns=cols)

def boot_ratio(sums,counts=None,reps=1000):
    """Resample validation blocks; ratio of structure sums/counts preserves QMOF estimand."""
    sums=np.asarray(sums,float);counts=np.ones(len(sums)) if counts is None else np.asarray(counts,float)
    rng=np.random.default_rng(SEED);values=[]
    for _ in range(reps):
        ii=rng.integers(0,len(sums),len(sums));values.append(sums[ii].sum()/counts[ii].sum())
    return np.quantile(values,[.025,.975])

def block_ci(values,components):
    g=pd.DataFrame({'v':np.asarray(values),'component':np.asarray(components)}).groupby('component').v.agg(['sum','size'])
    return boot_ratio(g['sum'],g['size'])

def calibration(t,folds,fold,pair):
    """Exactly the original weighted affine/12-bin rules; fit only current training structures.
    Grouped array implementation avoids repeated scans of the full atom table.
    """
    a,b=pair.split('_');xall=t[a].to_numpy();yall=t[b].to_numpy();wall=t.structure_weight.to_numpy();train=folds!=fold
    fallback=np.average(xall[train]-yall[train],weights=wall[train]);pred=np.empty(len(t));supported=np.zeros(len(t),bool);pars=[];knots=[]
    for el,ix0 in t.groupby('element',sort=True).indices.items():
        ix=np.asarray(ix0);tr=ix[train[ix]];xt=xall[ix]
        if not len(tr):
            pred[ix]=fallback;pars.append(dict(element=el,n_atoms=0,n_structures=0,offset=fallback,intercept=-fallback,slope=1.,supported=False));continue
        x=xall[tr];y=yall[tr];w=wall[tr];w=w/w.sum();xm=w@x;ym=w@y;var=w@((x-xm)**2);cov=w@((x-xm)*(y-ym));ns=t.iloc[tr].qmof_id.nunique()
        slope=cov/var if var>1e-6 and ns>=20 else 1.;inter=ym-slope*xm
        pars.append(dict(element=el,n_atoms=len(tr),n_structures=ns,offset=float(w@(x-y)),intercept=float(inter),slope=float(slope),supported=bool(ns>=20)))
        supported[ix]=ns>=20
        if ns<20 or len(tr)<100 or len(np.unique(x))<20:yp=inter+slope*xt
        else:
            bounds=np.unique(np.quantile(x,np.linspace(0,1,13)));bins=np.searchsorted(bounds[1:-1],x,side='right');xx=[];yy=[]
            for k in range(len(bounds)-1):
                mask=bins==k
                if mask.any():xx.append(np.average(x[mask],weights=w[mask]));yy.append(np.average(y[mask],weights=w[mask]))
            xx=np.asarray(xx);yy=np.asarray(yy);yp=np.interp(xt,xx,yy)
            if len(xx)>1:
                lo=xt<xx[0];hi=xt>xx[-1]
                yp[lo]=yy[0]+(xt[lo]-xx[0])*(yy[1]-yy[0])/(xx[1]-xx[0]);yp[hi]=yy[-1]+(xt[hi]-xx[-1])*(yy[-1]-yy[-2])/(xx[-1]-xx[-2])
            knots.extend(dict(element=el,x=float(x0),y=float(y0)) for x0,y0 in zip(xx,yy))
        pred[ix]=xt-yp
    return pred,supported,pars,knots

def feature_matrix():
    """One reusable numeric matrix; original feature vocabulary/column order frozen."""
    from importlib import import_module
    cfg=json.loads((ROOT/'results/models/model_configuration.json').read_text());names=cfg['first_features']+cfg['second_features']
    p=WORK/'structural_X.npy';mp=WORK/'matrix_identity.json';ident={'descriptor_sha256':cfg['descriptor_table_sha256'],'names':names}
    if p.exists() and mp.exists():assert json.loads(mp.read_text())==ident
    else:
        legacy=import_module('18_grouped_models')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',pd.errors.PerformanceWarning);f,actual,nfirst=legacy.feature_frames(ROOT/'data/structural_features.parquet')
        assert actual==names
        x=np.lib.format.open_memmap(p,mode='w+',dtype='float32',shape=(len(f),len(names)))
        for j,c in enumerate(names):x[:,j]=f[c].to_numpy(dtype=np.float32)
        x.flush();del x,f;gc.collect();mp.write_text(json.dumps(ident,indent=2))
    center=[names.index(c) for c in names if c.startswith('central_') or c in ['atomic_number','group','period','electronegativity','covalent_radius','is_metal']]
    return np.load(p,mmap_mode='r'),names,len(cfg['first_features']),center,cfg

def physical_edge_set(g):return {(i,*e[:4]) for i,a in enumerate(g) for e in a}

def model_input(X,ids,columns,qA):
    return np.column_stack((X[np.ix_(ids,columns)],np.asarray(qA)[ids])).astype(np.float32)

def predict(model,X,ids,columns,qA):
    out=np.empty(len(ids),np.float32)
    for k in range(0,len(ids),60000):out[k:k+60000]=model.predict(model_input(X,ids[k:k+60000],columns,qA))
    return out

QUANTITIES=('raw_signed','raw_absolute','calibrated_signed','calibrated_absolute','residual_sd')

def environment_units(features,predictions):
    from importlib import import_module
    f=features[['qmof_id','atom_index','element','ambiguous']].copy();f['environment']=import_module('19_results').environment(features)
    d=predictions[['qmof_id','atom_index','observed','residual_charge_only']].merge(f,on=['qmof_id','atom_index'],validate='one_to_one')
    d['raw_signed']=d.observed;d['raw_absolute']=np.abs(d.observed);d['calibrated_signed']=d.residual_charge_only;d['calibrated_absolute']=np.abs(d.residual_charge_only);d['r2']=d.residual_charge_only**2
    g=d.groupby(['qmof_id','element','environment'],observed=True).agg(raw_signed=('raw_signed','mean'),raw_absolute=('raw_absolute','mean'),calibrated_signed=('calibrated_signed','mean'),calibrated_absolute=('calibrated_absolute','mean'),r2=('r2','mean'),n_atoms=('atom_index','size'),ambiguous_fraction=('ambiguous','mean')).reset_index()
    g['residual_sd']=np.sqrt(np.maximum(0,g.r2-g.calibrated_signed**2));return g.drop(columns='r2')

def paired_vectors(units):
    pieces=[]
    for left,right in CONTRASTS:
        a=units[units.environment.eq(left)].set_index('qmof_id');b=units[units.environment.eq(right)].set_index('qmof_id');ids=a.index.intersection(b.index)
        for q in QUANTITIES:
            v=pd.DataFrame(dict(qmof_id=ids,left=left,right=right,quantity=q,effect=a.loc[ids,q].to_numpy()-b.loc[ids,q].to_numpy(),left_atoms=a.loc[ids,'n_atoms'].to_numpy(),right_atoms=b.loc[ids,'n_atoms'].to_numpy(),left_ambiguity=a.loc[ids,'ambiguous_fraction'].to_numpy(),right_ambiguity=b.loc[ids,'ambiguous_fraction'].to_numpy()))
            pieces.append(v)
    return pd.concat(pieces,ignore_index=True)
