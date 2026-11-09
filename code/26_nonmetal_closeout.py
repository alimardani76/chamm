"""Only the 96 retained structures: accepted metal edges fixed; nonmetal scales 1.20/1.30."""
import argparse,json,hashlib,time
import numpy as np,pandas as pd
from chamm_closeout import *
from chamm_graph import structure,cutoff,features,validate,key,info
DEST=OUT/'nonmetal'

def undirected(g):
    records={}
    for i,a in enumerate(g):
        for e in a:
            k=(i,*e[:4]);rev=(e[0],i,-e[1],-e[2],-e[3]);records[min(k,rev)]=e[4]
    return records

def build():
    DEST.mkdir(parents=True,exist_ok=True);(DEST/'graphs').mkdir(exist_ok=True)
    manifest=pd.read_csv(ROOT/'results/sensitivity/subset_manifest.csv');cid=json.loads((ROOT/'data/descriptor_manifest.json').read_text())['configuration_id']
    config=dict(source_accepted_configuration=cid,subset_sha256=hashlib.sha256((ROOT/'results/sensitivity/subset_manifest.csv').read_bytes()).hexdigest(),nonmetal_scales=[1.2,1.25,1.3],minimum_distance=.4,metal_edges='accepted image-aware edges retained exactly; no CrystalNN calls',classification='pymatgen Element.is_metal')
    cp=DEST/'configuration.json'
    if cp.exists():assert json.loads(cp.read_text())==config
    else:cp.write_text(json.dumps(config,indent=2))
    manifest.to_csv(DEST/'subset_manifest.csv',index=False);frames={s:[] for s in (1.2,1.25,1.3)};edges=[];members=[];progress=[]
    legacyenv=__import__('importlib').import_module('19_results').environment
    for k,q in enumerate(manifest.qmof_id):
        r=json.loads((ROOT/'work/structures'/(q+'.json')).read_text());s=structure(r);cached=json.loads((ROOT/'work'/('graphs_'+cid)/(q+'.json')).read_text());assert cached['configid']==cid
        base=cached['graph'];baseedges=undirected(base);metal=[x.specie.is_metal for x in s];cg=cutoff(s)
        bf=pd.read_parquet(ROOT/'work'/('features_'+cid)/(q+'.parquet'));frames[1.25].append(bf);be=legacyenv(bf)
        for scale in (1.2,1.3):
            g=[{key(e):e for e in a if metal[i] or metal[e[0]]} for i,a in enumerate(base)]
            for i,a in enumerate(cutoff(s,scale)):
                for e in a:
                    if not metal[i] and not metal[e[0]]:g[i][key(e)]=e
            g=[sorted(a.values(),key=key) for a in g];validate(s,g)
            assert {(i,*e[:4]) for i,a in enumerate(g) for e in a if metal[i] or metal[e[0]]}=={(i,*e[:4]) for i,a in enumerate(base) for e in a if metal[i] or metal[e[0]]}
            af=pd.DataFrame(features(r,s,g,cached['metadata'],cg));frames[scale].append(af);ae=legacyenv(af)
            altedges=undirected(g)
            for status,ks,distances in [('added',altedges.keys()-baseedges.keys(),altedges),('removed',baseedges.keys()-altedges.keys(),baseedges)]:
                for ek in sorted(ks):
                    i,j,tx,ty,tz=ek;el1=s[i].specie.symbol;el2=s[j].specie.symbol;note='boundary contact; no independent bond label';anchor=np.nan;angle=np.nan
                    if status=='added' and {el1,el2}=={'H','O'}:
                        h=i if el1=='H' else j;a=base[h]
                        if a:
                            anchor=min(e[4] for e in a);note='additional H...O contact alongside accepted bond; secondary/shared-proton candidate; 1.30 is stress configuration'
                            old=min(a,key=lambda e:e[4]);new=next(e for e in g[h] if key(e) not in {key(v) for v in base[h]} and s[e[0]].specie.symbol=='O')
                            v=(s[old[0]].frac_coords+old[1:4]-s[h].frac_coords)@s.lattice.matrix;w=(s[new[0]].frac_coords+new[1:4]-s[h].frac_coords)@s.lattice.matrix;angle=float(np.degrees(np.arccos(np.clip(v@w/np.linalg.norm(v)/np.linalg.norm(w),-1,1))))
                    edges.append(dict(qmof_id=q,scale=scale,status=status,i=i,j=j,tx=tx,ty=ty,tz=tz,element_i=el1,element_j=el2,distance=distances[ek],accepted_H_anchor_distance=anchor,H_angle_degree=angle,interpretation=note))
            for i in range(len(s)):
                firstchanged={key(e) for e in base[i]}!={key(e) for e in g[i]}
                pathchanged=any(bf.iloc[i][c]!=af.iloc[i][c] for c in bf if c.startswith('twohop_'))
                members.append(dict(qmof_id=q,atom_index=i,element=s[i].specie.symbol,scale=scale,accepted_environment=be.iloc[i],alternative_environment=ae.iloc[i],accepted_cn=len(base[i]),alternative_cn=len(g[i]),changed_physical_first_shell=firstchanged,changed_twohop=pathchanged,changed_environment=be.iloc[i]!=ae.iloc[i]))
            (DEST/'graphs'/f'{q}_s{scale:.2f}.json').write_text(json.dumps(dict(qmof_id=q,source_accepted_configuration=cid,scale=scale,graph=g),separators=(',',':')))
        if (k+1)%20==0:print('nonmetal subset built',k+1,96,flush=True)
    for scale,fs in frames.items():pd.concat(fs,ignore_index=True).to_parquet(DEST/f'features_{scale:.2f}.parquet',index=False,compression='zstd')
    pd.DataFrame(edges).to_csv(DEST/'physical_edge_changes.csv',index=False);mm=pd.DataFrame(members);mm.to_parquet(DEST/'membership_changes.parquet',index=False)
    mm.groupby(['scale','element']).agg(n_atoms=('atom_index','size'),n_structures=('qmof_id','nunique'),first_shell_changed=('changed_physical_first_shell','sum'),twohop_changed=('changed_twohop','sum'),environment_changed=('changed_environment','sum')).reset_index().to_csv(DEST/'change_summary.csv',index=False)
    review_boundary_hydrogens()
    print('No metal NN run; physical edge changes',len(edges),flush=True)

def review_boundary_hydrogens():
    """Geometry ledger for the actually added H...O contacts, including shared-proton cases."""
    p=DEST/'physical_edge_changes.csv';d=pd.read_csv(p);cid=json.loads((ROOT/'data/descriptor_manifest.json').read_text())['configuration_id']
    for row in d.itertuples():
        if row.status!='added' or {row.element_i,row.element_j}!={'H','O'}:continue
        h=row.i if row.element_i=='H' else row.j;r=json.loads((ROOT/'work/structures'/(row.qmof_id+'.json')).read_text());s=structure(r)
        base=json.loads((ROOT/'work'/('graphs_'+cid)/(row.qmof_id+'.json')).read_text())['graph'];g=json.loads((DEST/'graphs'/f'{row.qmof_id}_s{row.scale:.2f}.json').read_text())['graph']
        d.loc[row.Index,'accepted_H_cn']=len(base[h]);d.loc[row.Index,'alternative_H_cn']=len(g[h])
        oks=[e for e in g[h] if s[e[0]].specie.symbol=='O']
        if not base[h] and len(oks)==2:
            v=[(s[e[0]].frac_coords+e[1:4]-s[h].frac_coords)@s.lattice.matrix for e in oks]
            angle=float(np.degrees(np.arccos(np.clip(v[0]@v[1]/np.linalg.norm(v[0])/np.linalg.norm(v[1]),-1,1))))
            d.loc[row.Index,'H_angle_degree']=angle;d.loc[row.Index,'O_O_distance']=float(np.linalg.norm(v[0]-v[1]));d.loc[row.Index,'interpretation']='near-midpoint O...H...O pair; accepted H isolated, wider graph gives H degree 2; shared-proton/secondary-contact ambiguity; 1.30 stress configuration'
        elif base[h]:
            d.loc[row.Index,'interpretation']='added H...O alongside accepted contact; secondary/shared-proton candidate; 1.30 stress configuration'
    d.to_csv(p,index=False)

def analyze():
    assert (OUT/'original/COMPLETE.json').exists();ids=pd.read_csv(DEST/'subset_manifest.csv').qmof_id.tolist();avail=[];common=[];unitsall=[];subset=[];members=[]
    for pair in PAIRS:
        preds=pd.concat([pd.read_parquet(p,filters=[('qmof_id','in',ids)],columns=['qmof_id','atom_index','observed','residual_charge_only','fold']) for p in sorted((OUT/'original').glob(f'oof_fold*_{pair}.parquet'))],ignore_index=True)
        preds['pair']=pair;subset.append(preds);pv={}
        for scale in (1.2,1.25,1.3):
            f=pd.read_parquet(DEST/f'features_{scale:.2f}.parquet');u=environment_units(f,preds);u['pair']=pair;u['scale']=scale;unitsall.append(u)
            pv[scale]=paired_vectors(u)
            for (left,right,quantity),g in pv[scale].groupby(['left','right','quantity']):
                ci=boot_ratio(g.effect.to_numpy());avail.append(dict(pair=pair,scale=scale,left=left,right=right,quantity=quantity,n_structures=len(g),left_atoms=int(g.left_atoms.sum()),right_atoms=int(g.right_atoms.sum()),effect=g.effect.mean(),ci_low=ci[0],ci_high=ci[1],adequate_support=len(g)>=30))
        for left,right in CONTRASTS:
            sets={s:set(df.loc[df.left.eq(left)&df.quantity.eq('raw_signed'),'qmof_id']) for s,df in pv.items()};commonids=set.intersection(*sets.values())
            for q in sorted(set.union(*sets.values())):members.append(dict(pair=pair,left=left,right=right,qmof_id=q,pair_at_120=q in sets[1.2],pair_at_125=q in sets[1.25],pair_at_130=q in sets[1.3],in_common=q in commonids))
            for quantity in QUANTITIES:
                base=pv[1.25];base=base[base.left.eq(left)&base.quantity.eq(quantity)&base.qmof_id.isin(commonids)].set_index('qmof_id').sort_index()
                if base.empty:continue
                for scale in (1.2,1.25,1.3):
                    a=pv[scale];a=a[a.left.eq(left)&a.quantity.eq(quantity)&a.qmof_id.isin(commonids)].set_index('qmof_id').reindex(base.index)
                    delta=a.effect-base.effect;ci=boot_ratio(delta.to_numpy());ec=boot_ratio(a.effect.to_numpy())
                    common.append(dict(pair=pair,scale=scale,left=left,right=right,quantity=quantity,n_common_structures=len(a),left_atoms=int(a.left_atoms.sum()),right_atoms=int(a.right_atoms.sum()),accepted_effect=base.effect.mean(),effect=a.effect.mean(),effect_ci_low=ec[0],effect_ci_high=ec[1],change=delta.mean(),change_ci_low=ci[0],change_ci_high=ci[1],adequate_support=len(a)>=30))
    pd.concat(subset).to_parquet(DEST/'subset_stronger_residuals.parquet',index=False);pd.concat(unitsall).to_parquet(DEST/'environment_units.parquet',index=False)
    pd.DataFrame(avail).to_csv(DEST/'available_pair_contrasts.csv',index=False);pd.DataFrame(common).to_csv(DEST/'common_pair_contrasts.csv',index=False);pd.DataFrame(members).to_csv(DEST/'pair_membership.csv',index=False)
    (DEST/'COMPLETE.json').write_text(json.dumps(dict(structures=96,source_residuals='original-fold flexible charge-only OOF',scales=[1.2,1.25,1.3],new_metal_NN_runs=0),indent=2));print('Nonmetal contrast sensitivity complete',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--build-only',action='store_true');p.add_argument('--analyze-only',action='store_true');a=p.parse_args()
    if not a.analyze_only:build()
    if not a.build_only:analyze()
