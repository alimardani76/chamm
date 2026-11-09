"""Bounded 34-structure comparison; archived baseline results reused for 18 seeds."""
import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time,argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import pandas as pd
from chamm_graph import *
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/connectivity';OUT.mkdir(parents=True,exist_ok=True)
OLD=ROOT/'work/retained/ChaMM_restart/ChaMM_connectivity_review/ChaMM_connectivity_review/results/pilot'

def oldgraph(d,name,n):
    rows=d['methods'][name]['sites'];return [[[e['j'],*e['image'],e['distance']] for e in rows[str(i)]['edges']] for i in range(n)]

def run(r):
    p=OUT/(r['qmof_id']+'.json')
    if p.exists():return r['qmof_id'],'cached'
    s=structure(r);t=time.monotonic(); methods={};timings={};metas={}
    cg=cutoff(s);methods['cutoff']=cg
    old=OLD/(r['qmof_id']+'.json')
    for name,oname in [('default','crystal_default'),('isayev','isayev_0.25')]:
        t0=time.monotonic()
        if old.exists():
            methods[name]=oldgraph(json.loads(old.read_text()),oname,len(s));timings[name]='reused retained 18-seed result'
        else:
            methods[name],metas[name]=nearest(s,name);timings[name]=time.monotonic()-t0
    t0=time.monotonic();zg,meta=nearest(s);methods['zero']=zg;timings['zero']=time.monotonic()-t0;metas['zero']=meta
    for rule in ('both','either','metal_center_union'):methods['hybrid_'+rule]=hybrid(s,cg,zg,rule)
    metals=[i for i,x in enumerate(s) if x.specie.is_metal]
    for shift in (-0.1,0.1):
        t0=time.monotonic();gg,mm=nearest(s,indices=metals,distance_shift=shift)
        name=f'hybrid_shift{shift:+.1f}';methods[name]=hybrid(s,cg,gg);metas[name]=mm;timings[name]=time.monotonic()-t0
    for scale in (1.20,1.30):methods[f'hybrid_radius{scale:.2f}']=hybrid(s,cutoff(s,scale),zg)
    checks={m:validate(s,g,reciprocal=m.startswith('hybrid') or m=='cutoff') for m,g in methods.items()}
    data=dict(qmof_id=r['qmof_id'],atoms=len(s),methods=methods,metadata=metas,timings=timings,checks=checks,seconds=time.monotonic()-t)
    p.write_text(json.dumps(data));return r['qmof_id'],data['seconds']

def summarize(records):
    rows=[];contacts=[];sensitivity=[]
    cases={'qmof-85a8986':[0,12,44],'qmof-8a95c27':[2],'qmof-9938276':[28],'qmof-8eb67f4':[0], 'qmof-be32c2a':[4],'qmof-d689f74':[61]}
    for r in records:
        d=json.loads((OUT/(r['qmof_id']+'.json')).read_text());s=structure(r)
        for m,g in d['methods'].items():
            for i,a in enumerate(g):
                rows.append(dict(qmof_id=r['qmof_id'],atom_index=i,element=s[i].specie.symbol,method=m,cn=len(a),formula=formula(s,a)))
                if i in cases.get(r['qmof_id'],[]):
                    for e in a:contacts.append(dict(qmof_id=r['qmof_id'],atom_index=i,element=s[i].specie.symbol,method=m,neighbor=e[0],neighbor_element=s[e[0]].specie.symbol,tx=e[1],ty=e[2],tz=e[3],distance=e[4]))
        base=d['methods']['hybrid_metal_center_union']
        for m in ('hybrid_both','hybrid_either','hybrid_shift-0.1','hybrid_shift+0.1','hybrid_radius1.20','hybrid_radius1.30'):
            for i,(a,b) in enumerate(zip(base,d['methods'][m])):
                if {key(e) for e in a}!={key(e) for e in b}:sensitivity.append(dict(qmof_id=r['qmof_id'],atom_index=i,element=s[i].specie.symbol,variant=m,base_cn=len(a),variant_cn=len(b),base_formula=formula(s,a),variant_formula=formula(s,b)))
    pd.DataFrame(rows).to_csv(OUT/'panel_sites.csv',index=False)
    pd.DataFrame(contacts).to_csv(OUT/'decisive_contacts.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(OUT/'changed_sites.csv',index=False)
    print(pd.DataFrame(rows).query("qmof_id in @cases").merge(pd.DataFrame([(q,i) for q,ids in cases.items() for i in ids],columns=['qmof_id','atom_index'])).to_string(index=False),flush=True)

if __name__=='__main__':
    records=json.loads((ROOT/'inputs/expanded_benchmark_records.json').read_text());t=time.monotonic()
    with ProcessPoolExecutor(max_workers=2) as pool:
        for f in as_completed([pool.submit(run,r) for r in records]):print(f.result(),flush=True)
    summarize(records);print('panel wall seconds',time.monotonic()-t,flush=True)
