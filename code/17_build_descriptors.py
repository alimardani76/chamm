"""Resumable Windows-safe, two-worker production build. No electronic predictors."""
import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,json,time,hashlib,traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import pandas as pd
import pyarrow as pa,pyarrow.parquet as pq
from chamm_graph import *
ROOT=Path(__file__).resolve().parents[1]

def process(args):
    q,configid=args;cache=ROOT/'work'/('graphs_'+configid);parts=ROOT/'work'/('features_'+configid)
    dest=parts/(q+'.parquet');sp=ROOT/'work/structures'/(q+'.json')
    try:
        source=sp.read_bytes();sha=hashlib.sha256(source).hexdigest()
        gp=cache/(q+'.json');r=json.loads(source);s=structure(r);cg=cutoff(s);start=time.monotonic()
        if gp.exists():
            data=json.loads(gp.read_text());assert data['source_sha256']==sha and data['configid']==configid
            g=data['graph'];meta=data['metadata']
        else:
            metals=[i for i,x in enumerate(s) if x.specie.is_metal]
            zg,meta=nearest(s,indices=metals);g=hybrid(s,cg,zg)
            validate(s,g)
            data=dict(qmof_id=q,source_sha256=sha,configid=configid,graph=g,metadata=meta,seconds=time.monotonic()-start)
            tmp=gp.with_suffix('.tmp');tmp.write_text(json.dumps(data,separators=(',',':')));tmp.replace(gp)
        rows=features(r,s,g,meta,cg);df=pd.DataFrame(rows)
        # Stable schema even in parts where every distance is missing.
        for col in df:
            if col in ('qmof_id','element','first_formula','neighbor_metal_formula'):df[col]=df[col].astype(str)
            elif df[col].dtype.kind not in 'biu':df[col]=df[col].astype('float64')
        tmp=dest.with_suffix('.tmp');df.to_parquet(tmp,index=False,compression='zstd');tmp.replace(dest)
        return dict(qmof_id=q,atoms=len(s),valid=not bool(meta['errors']),errors=meta['errors'],seconds=time.monotonic()-start)
    except Exception as e:return dict(qmof_id=q,failed=repr(e),traceback=traceback.format_exc())

def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=2);p.add_argument('--limit',type=int);a=p.parse_args()
    assert (ROOT/'reports/CONNECTIVITY_DECISION.md').exists(),'Connectivity must be decided before production'
    cid,config=identity();cache=ROOT/'work'/('graphs_'+cid);parts=ROOT/'work'/('features_'+cid)
    cache.mkdir(exist_ok=True);parts.mkdir(exist_ok=True)
    (cache/'configuration.json').write_text(json.dumps(config,indent=2))
    cohort=pd.read_csv(ROOT/'data/chamm_common_cohort.csv');ids=cohort.qmof_id.tolist()
    if a.limit:ids=ids[:a.limit]
    pending=[q for q in ids if not (parts/(q+'.parquet')).exists()];start=time.monotonic();done=len(ids)-len(pending);nfail=0
    print(f'Configuration {cid}; existing {done}; pending {len(pending)}; workers {a.workers}',flush=True)
    log=ROOT/'work/production_events.jsonl'
    with ProcessPoolExecutor(max_workers=a.workers) as pool,log.open('a') as stream:
        for future in as_completed([pool.submit(process,(q,cid)) for q in pending]):
            result=future.result();stream.write(json.dumps(result)+'\n');stream.flush();done+=1;nfail+=bool(result.get('failed'))
            if done%100==0 or result.get('failed'):
                print(f'{done}/{len(ids)} structures; failures this run {nfail}; elapsed {(time.monotonic()-start)/60:.1f} min; {result["qmof_id"]}',flush=True)
    if a.limit:return
    files=[parts/(q+'.parquet') for q in ids];assert all(f.exists() for f in files),'Missing structures; inspect production_events.jsonl, rerun to resume'
    dest=ROOT/'data/structural_features.parquet';temp=dest.with_suffix('.parquet.partial');writer=None;atoms=0;valid=0;flags=Counter()
    target=pq.read_table(ROOT/'data/charge_targets.parquet',columns=['qmof_id','atom_index','element']).to_pandas().set_index(['qmof_id','atom_index'])
    for k,(q,f) in enumerate(zip(ids,files)):
        table=pq.read_table(f);df=table.to_pandas();assert len(df)==int(cohort.iloc[k].n_atoms)
        assert df.atom_index.tolist()==list(range(len(df))) and df.qmof_id.eq(q).all()
        assert df.element.tolist()==target.loc[q].element.tolist(),'Source-site alignment mismatch'
        if writer is None:writer=pq.ParquetWriter(temp,table.schema,compression='zstd')
        writer.write_table(table.cast(writer.schema));atoms+=len(df);valid+=int(df.graph_valid.sum())
        for flag in ('ambiguous','graph_failure','zero_neighbors','metal_c_contact','metal_h_contact','metal_metal_contact','extended_metal_contact','cutoff_disagreement','radius_sensitivity'):
            flags[flag]+=int(df[flag].sum())
    writer.close();assert atoms==1420532
    # Never silently replace an accepted table from a different configuration.
    status=ROOT/'data/descriptor_manifest.json'
    if dest.exists() and status.exists():assert json.loads(status.read_text())['configuration_id']==cid
    temp.replace(dest)
    report=dict(structures=len(ids),atoms=atoms,valid_graph_atoms=valid,flags=dict(flags),configuration_id=cid,
                configuration=config,source_cache='work/source_cache_manifest.json',elapsed_this_invocation_seconds=time.monotonic()-start,
                charge_alignment='PASS all keys/elements/source site order',failures_this_invocation=nfail,
                table_sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
    status.write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
