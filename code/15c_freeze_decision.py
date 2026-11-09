"""Persist the bounded contact adjudication and execution identity."""
import json,time
from pathlib import Path
import numpy as np,pandas as pd
from chamm_graph import *
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/connectivity'

def main():
    rs={r['qmof_id']:r for r in json.loads((ROOT/'inputs/expanded_benchmark_records.json').read_text())};replay=[]
    for q,i in [('qmof-85a8986',44),('qmof-85a8986',0),('qmof-8a95c27',2),('qmof-9938276',28),('qmof-d689f74',61)]:
        s=structure(rs[q]);old=json.loads((OUT/(q+'.json')).read_text())
        for name in ('default','isayev'):
            gg,meta=nearest(s,name,indices=[i]);assert not meta['errors'];assert {key(e) for e in gg[i]}=={key(e) for e in old['methods'][name][i]}
            replay.append(dict(qmof_id=q,atom_index=i,method=name,retained_replay_match=True))
    pd.DataFrame(replay).to_csv(OUT/'retained_baseline_replay.csv',index=False)
    ledger=[
      ('qmof-85a8986',44,'C-C','1.3944, 1.3957','clear covalent','retained','short aromatic-like ligand C-C; default CrystalNN misses both; cutoff and zero-weight retain'),
      ('qmof-85a8986',12,'O-H / H...O','1.0025 / 1.6211','covalent / secondary hydrogen contact','short retained; longer rejected','nonmetal cutoff separates short O-H from hydrogen contact without using charge'),
      ('qmof-8a95c27',2,'Cu-O','1.9817-2.0185 + 2.3985 pair','plausible extended coordination','all six retained; longer pair interpretation flagged in ledger','long pair 179.99935 degrees; nearly opposite; no oxidation-state inference'),
      ('qmof-85a8986',0,'Dy-O/N vs C/H','O 2.3296-2.5118; N 2.6585,2.7144; C 3.1616-3.3346; H 2.7614','O/N plausible coordination; distant C/H secondary or unresolved','O7/N2 retained; C4/H1 rejected','distant C sites are ligand carbons adjacent to donor-bearing ligand skeleton; distances and screening support secondary interpretation; not a universal metal-C veto'),
      ('qmof-9938276',28,'Pb-O/N','2.6506-2.7696','plausible coordination; shell boundary uncertain','O4/N4 retained','coherent eight-contact shell narrowly outside cutoff; CN8 operational, not independent chemical truth'),
      ('qmof-8eb67f4',0,'Zr-O vs Zr-Zr','O 2.0705-2.2735; Zr 3.5340-3.5379','O coordination / cluster metal proximity','O8 retained; Zr4 proximity recorded separately','O bridges common to cluster; selected graph excludes these Zr links, not all metal-metal links'),
      ('qmof-be32c2a',4,'Zr-C/P','C 2.4040-2.5677; P 2.6154-2.6945','plausible organometallic atomic contacts; ligand-group CN differs','C11/P2 retained and metal-C flag set','ten C atoms form two five-membered cycles in induced carbon graph plus one additional C; atomic CN13 is not ligand-group CN'),
      ('qmof-d689f74',61,'Pt-C/N','C 1.9843-1.9875; N 2.2344-2.2383','clear plausible organometallic coordination','C2/N2 retained and metal-C flag set','short C contacts preserved; no blanket carbon deletion'),
    ]
    pd.DataFrame(ledger,columns=['qmof_id','atom_index','contact','distance_angstrom','interpretation','decision','reason']).to_csv(OUT/'contact_review_ledger.csv',index=False)
    import networkx as nx
    r=rs['qmof-be32c2a'];s=structure(r);d=json.loads((OUT/(r['qmof_id']+'.json')).read_text());g=d['methods']['cutoff'];cs={e[0] for e in g[4] if s[e[0]].specie.symbol=='C'}
    ng=nx.Graph();ng.add_nodes_from(cs)
    for i in cs:
        for e in g[i]:
            if e[0] in cs:ng.add_edge(i,e[0])
    cycles=nx.cycle_basis(ng);assert sorted(map(len,cycles))==[5,5]
    summaries=[]
    for p in OUT.glob('qmof-*.json'):
        d=json.loads(p.read_text());summaries.append(dict(qmof_id=d['qmof_id'],atoms=d['atoms'],zero_seconds=d['timings']['zero'],
                    metal_shift_seconds=d['timings']['hybrid_shift+0.1'],metal_sites=len(d['metadata']['hybrid_shift+0.1']['uncertainty']),
                    operational_edges=d['checks']['hybrid_metal_center_union']['directed_edges'],
                    missing_reverse=d['checks']['hybrid_metal_center_union']['missing_reverse'],extra_images=d['checks']['hybrid_metal_center_union']['extra_images']))
    pd.DataFrame(summaries).to_csv(OUT/'runtime_and_graph_summary.csv',index=False)
    cid,cfg=identity();elements=sorted({x['species'][0]['element'] for r in rs.values() for x in r['structure']['sites']})
    metadata=dict(configuration_id=cid,configuration=cfg,metal_elements_in_panel=[e for e in elements if info(e)['metal']],
        panel_structures=34,panel_atoms=5198,baseline_origin='18 seed baselines reused; 16 additions recomputed; 10 decisive seed-site calculations replayed with installed versions',
        carbon_cycles=cycles,decision='accept metal-center union hybrid for production with ambiguity flags')
    (OUT/'decision_metadata.json').write_text(json.dumps(metadata,indent=2))
    print('Decision frozen',cid,'carbon cycles',cycles,flush=True)

if __name__=='__main__':main()
