"""Charge-blind periodic operational graph and compact chemical descriptors."""
import hashlib, json, warnings
from collections import Counter
from functools import lru_cache
import numpy as np
from ase.data import covalent_radii, atomic_numbers
from pymatgen.core import Structure, Element
from pymatgen.analysis.local_env import CrystalNN, IsayevNN

CRYSTAL=dict(weighted_cn=False,cation_anion=False,distance_cutoffs=(0.5,1.0),x_diff_weight=0.0,
             porous_adjustment=True,search_cutoff=7.0,fingerprint_length=None)
CONFIG=dict(version='operational-v1',nonmetal_scale=1.25,min_distance=0.4,metal_rule='metal_center_union',
            metal_classification='pymatgen.core.Element.is_metal',radii='ase.data.covalent_radii',crystal=CRYSTAL,
            metal_metal='selected by either metal center; retain and flag; cutoff proximity separate')
KEYELS=('H','B','C','N','O','F','Si','P','S','Cl','Se','Br','I')

@lru_cache(None)
def info(el):
    e=Element(el)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore'); x=float(e.X)
    r=float(covalent_radii[atomic_numbers[el]])
    assert np.isfinite(r) and r>0
    return dict(Z=e.Z,group=e.group,period=e.row,X=x,radius=r,metal=bool(e.is_metal))

def structure(record):
    raw=record['structure']; assert raw['sites']
    assert all(len(x['species'])==1 and x['species'][0]['occu']==1 for x in raw['sites'])
    # No charges, oxidation states, bond-order or spin metadata are supplied to NN.
    s=Structure(raw['lattice']['matrix'],[x['species'][0]['element'] for x in raw['sites']],
                [x['abc'] for x in raw['sites']],coords_are_cartesian=False)
    assert np.isfinite(s.cart_coords).all() and abs(s.lattice.volume)>1e-6
    return s

def key(e):return tuple(e[:4])
def reverse(i,e):return (i,-e[1],-e[2],-e[3])
def directed_set(g):return {(i,*e[:4]) for i,a in enumerate(g) for e in a}
def formula(s,edges):return ';'.join(f'{e}{n}' for e,n in sorted(Counter(s[x[0]].specie.symbol for x in edges).items(),key=lambda kv:atomic_numbers[kv[0]]))

def cutoff(s,scale=1.25):
    radii=np.array([info(x.specie.symbol)['radius'] for x in s]); g=[[] for _ in s]
    ii,jj,ims,dd=s.get_neighbor_list(r=2*scale*max(radii)+0.1)
    for i,j,im,d in zip(ii,jj,ims,dd):
        if 0.4<=d<=scale*(radii[i]+radii[j]):g[int(i)].append([int(j),*np.rint(im).astype(int).tolist(),float(d)])
    return g

def nearest(s,method='zero',indices=None,distance_shift=0.0):
    indices=list(range(len(s))) if indices is None else list(indices)
    cfg={**CRYSTAL,'distance_cutoffs':tuple(v+distance_shift for v in CRYSTAL['distance_cutoffs'])}
    if method=='default':cfg['x_diff_weight']=3.0
    obj=IsayevNN(tol=0.25,cutoff=13.,allow_pathological=False,extra_nn_info=False,compute_adj_neighbors=False) if method=='isayev' else CrystalNN(**cfg)
    g=[[] for _ in s]; errors={}; uncertainty={}; warning_counts=Counter()
    for i in indices:
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter('always')
            try:
                if method=='isayev':ns=obj.get_nn_info(s,i)
                else:
                    data=obj.get_nn_data(s,i)
                    cn=max(data.cn_weights,key=data.cn_weights.get)
                    ns=data.cn_nninfo[cn]
                    probs=sorted(data.cn_weights.values(),reverse=True)
                    uncertainty[i]={'cn_probability':probs[0],'cn_probability_gap':probs[0]-(probs[1] if len(probs)>1 else 0)}
                for n in ns:
                    j=int(n['site_index']); im=np.asarray(n['image']); assert np.allclose(im,np.rint(im),atol=1e-8)
                    im=np.rint(im).astype(int)
                    v=(s[j].frac_coords+im-s[i].frac_coords)@s.lattice.matrix
                    assert np.allclose(s[j].coords+im@s.lattice.matrix,n['site'].coords,atol=1e-6)
                    d=float(np.linalg.norm(v))
                    if d>=0.4:g[i].append([j,*im.tolist(),d])
            except Exception as exc:errors[i]=f'{type(exc).__name__}: {exc}'
            warning_counts.update(str(w.message) for w in ws)
    return g,dict(errors=errors,uncertainty=uncertainty,warnings=dict(warning_counts))

def hybrid(s,cg,zg,rule='metal_center_union'):
    metal=[info(x.specie.symbol)['metal'] for x in s]; sets=[{key(e) for e in a} for a in zg]
    out=[{} for _ in s]
    for i,a in enumerate(cg):
        for e in a:
            if not metal[i] and not metal[e[0]]:out[i][key(e)]=e
    for i,a in enumerate(zg):
        for e in a:
            j=e[0]
            if not (metal[i] or metal[j]):continue
            if rule=='both' and reverse(i,e) not in sets[j]:continue
            if rule=='metal_center_union' and not metal[i]:continue
            out[i][key(e)]=e
            rev=[i,-e[1],-e[2],-e[3],e[4]];out[j][key(rev)]=rev
    return [sorted(a.values(),key=key) for a in out]

def validate(s,g,reciprocal=True):
    sets=[{key(e) for e in a} for a in g]; missing=0; mult=0
    for i,a in enumerate(g):
        assert len(a)==len(sets[i]),'Duplicate physical edge'
        mult+=len(a)-len({e[0] for e in a})
        for e in a:
            assert not (i==e[0] and e[1:4]==[0,0,0])
            d=np.linalg.norm((s[e[0]].frac_coords+np.array(e[1:4])-s[i].frac_coords)@s.lattice.matrix)
            assert abs(d-e[4])<1e-6
            missing+=reverse(i,e) not in sets[e[0]]
    if reciprocal:assert missing==0,'Missing reverse'
    return dict(directed_edges=sum(map(len,g)),missing_reverse=missing,extra_images=mult)

def counts(counter,prefix):
    d={prefix+e:counter.get(e,0) for e in KEYELS}
    d[prefix+'metal']=sum(n for e,n in counter.items() if info(e)['metal'])
    d[prefix+'other_nonmetal']=sum(n for e,n in counter.items() if not info(e)['metal'] and e not in KEYELS)
    return d

def features(record,s,g,meta,cg):
    sy=[x.specie.symbol for x in s]; props=[info(e) for e in sy]
    first=[Counter(sy[e[0]] for e in a) for a in g]; output=[]
    errors=meta.get('errors',{}); uc=meta.get('uncertainty',{})
    cheap_low=cutoff(s,1.20); cheap_high=cutoff(s,1.30)
    for i,a in enumerate(g):
        p=props[i]; c=first[i]; second=Counter()
        for e in a:
            cc=first[e[0]].copy();cc[sy[i]]-=1; assert cc[sy[i]]>=0
            second.update(+cc)
        mc=Counter({el:n for el,n in c.items() if info(el)['metal']})
        d=dict(qmof_id=record['qmof_id'],atom_index=i,element=sy[i],atomic_number=p['Z'],group=p['group'],period=p['period'],
               electronegativity=p['X'],covalent_radius=p['radius'],is_metal=int(p['metal']),cn=len(a),
               first_formula=formula(s,a),neighbor_metal_formula=';'.join(f'{el}{n}' for el,n in sorted(mc.items())),
               twohop_paths=sum(second.values()),graph_valid=not bool(errors),zero_neighbors=not bool(a),
               graph_failure=bool(errors),metal_c_contact=bool((p['metal'] and c.get('C',0)) or (sy[i]=='C' and mc)),
               metal_h_contact=bool((p['metal'] and c.get('H',0)) or (sy[i]=='H' and mc)),
               metal_metal_contact=bool(p['metal'] and mc),
               cutoff_disagreement={key(e) for e in a}!={key(e) for e in cg[i]},
               radius_sensitivity=({key(e) for e in cheap_low[i] if not p['metal'] and not props[e[0]]['metal']}!={key(e) for e in cheap_high[i] if not p['metal'] and not props[e[0]]['metal']}),
               extended_metal_contact=any((p['metal'] or props[e[0]]['metal']) and e[4]/(p['radius']+props[e[0]]['radius'])>1.20 for e in a),
               cn_probability=uc.get(i,uc.get(str(i),{})).get('cn_probability',np.nan),
               cn_probability_gap=uc.get(i,uc.get(str(i),{})).get('cn_probability_gap',np.nan),
               geometric_metal_proximity=sum(props[e[0]]['metal'] for e in cg[i]))
        d.update(counts(c,'first_n_'));d.update(counts(second,'twohop_n_'))
        d['mean_neighbor_metal_Z']=np.mean([info(el)['Z'] for el in mc.elements()]) if mc else np.nan
        for label,vals in [('distance',[e[4] for e in a]),('normalized_distance',[e[4]/(p['radius']+props[e[0]]['radius']) for e in a])]:
            v=np.array(vals); d.update({label+'_'+k:float(f(v)) if len(v) else np.nan for k,f in [('mean',np.mean),('std',np.std),('min',np.min),('max',np.max)]})
        d['distance_cv']=d['distance_std']/d['distance_mean'] if a else np.nan
        x=np.array([props[e[0]]['X'] for e in a]);x=x[np.isfinite(x)]
        d['neighbor_X_mean']=float(x.mean()) if len(x) else np.nan
        d['neighbor_X_std']=float(x.std()) if len(x) else np.nan
        d['neighbor_abs_delta_X_mean']=float(np.abs(x-p['X']).mean()) if len(x) and np.isfinite(p['X']) else np.nan
        d['ambiguous']=bool(d['graph_failure'] or d['zero_neighbors'] or d['metal_c_contact'] or d['metal_h_contact'] or d['metal_metal_contact'] or d['extended_metal_contact'] or (np.isfinite(d['cn_probability_gap']) and d['cn_probability_gap']<0.1))
        output.append(d)
    return output

def identity():
    import importlib.metadata
    versions={m:importlib.metadata.version(m) for m in ('pymatgen','ase','numpy','scipy')}
    blob=json.dumps(dict(config=CONFIG,versions=versions,module_sha256=hashlib.sha256(__import__('pathlib').Path(__file__).read_bytes()).hexdigest()),sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:16],json.loads(blob)
