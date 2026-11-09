"""Mapped supercells and decisive local vectors, not a new algorithm search."""
import json
from pathlib import Path
import numpy as np,pandas as pd
from chamm_graph import *
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/connectivity'

def vectors(s,g,i,mapping=None):
    return sorted((int(mapping[e[0]]) if mapping is not None else e[0],*np.round((s[e[0]].frac_coords+np.array(e[1:4])-s[i].frac_coords)@s.lattice.matrix,5)) for e in g[i])

def main():
    rs={r['qmof_id']:r for r in json.loads((ROOT/'inputs/expanded_benchmark_records.json').read_text())};checks=[];geometry=[]
    for q,i in [('qmof-6247dea',0),('qmof-8eb67f4',0),('qmof-be32c2a',4)]:
        s=structure(rs[q]);d=json.loads((OUT/(q+'.json')).read_text());g=d['methods']['hybrid_metal_center_union']
        ss=s.copy();ss.make_supercell([2,1,1]);prim=ss.cart_coords@np.linalg.inv(s.lattice.matrix)
        mapping=[]
        for x in prim:
            diff=x-s.frac_coords;dist=np.linalg.norm((diff-np.rint(diff))@s.lattice.matrix,axis=1);j=int(dist.argmin());assert dist[j]<1e-6;mapping.append(j)
        zg,meta=nearest(ss,indices=[j for j,x in enumerate(ss) if x.specie.is_metal]);sg=hybrid(ss,cutoff(ss),zg);validate(ss,sg)
        for k,j in enumerate(mapping):
            ok=vectors(s,g,j)==vectors(ss,sg,k,mapping);assert ok,(q,j,k)
            path1=sum(len(g[e[0]])-1 for e in g[j]);path2=sum(len(sg[e[0]])-1 for e in sg[k]);assert path1==path2
            checks.append(dict(qmof_id=q,primitive_site=j,supercell_site=k,repeat='2,1,1',neighbor_vectors_match=ok,twohop_paths_match=True))
    q='qmof-8a95c27';s=structure(rs[q]);g=json.loads((OUT/(q+'.json')).read_text())['methods']['hybrid_metal_center_union'];a=g[2]
    long=[e for e in a if e[4]>2.3];v=[(s[e[0]].frac_coords+e[1:4]-s[2].frac_coords)@s.lattice.matrix for e in long]
    angle=np.degrees(np.arccos(np.clip(np.dot(*v)/(np.linalg.norm(v[0])*np.linalg.norm(v[1])),-1,1)))
    geometry.append(dict(case='Cu2 longer O pair',value=float(angle),units='degree',interpretation='nearly opposite axial contacts; retained plausible coordination; no oxidation assignment'))
    for q,i in [('qmof-85a8986',0),('qmof-be32c2a',4),('qmof-8eb67f4',0)]:
        s=structure(rs[q]);d=json.loads((OUT/(q+'.json')).read_text());cg=d['methods']['cutoff'];a=cg[i]
        for e in a:
            geometry.append(dict(case=f'{q}/{i} -> {e[0]} image {e[1:4]} {s[e[0]].specie.symbol}',value=e[4],units='angstrom',
                interpretation='neighbor shell '+formula(s,cg[e[0]])))
    pd.DataFrame(checks).to_csv(OUT/'mapped_supercell_checks.csv',index=False)
    pd.DataFrame(geometry).to_csv(OUT/'decisive_geometry.csv',index=False)
    print('mapped checks',len(checks),'Cu axial angle',angle,flush=True)

if __name__=='__main__':main()
