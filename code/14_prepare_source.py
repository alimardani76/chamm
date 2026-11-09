"""Download exactly the pinned archive, with resume; stream geometry once to cache."""
import argparse, hashlib, json, time, zipfile
from pathlib import Path
import requests
ROOT=Path(__file__).resolve().parents[1]

def download():
    meta=json.loads((ROOT/'inputs/figshare_article_v18.json').read_text())
    f=next(x for x in meta['files'] if x['id']==59573735)
    dest=ROOT/'work/qmof_database.zip'; part=dest.with_suffix('.zip.partial')
    dest.parent.mkdir(exist_ok=True)
    if not dest.exists():
        for attempt in range(8):
            try:
                n=part.stat().st_size if part.exists() else 0
                with requests.get(f['download_url'],headers={'Range':f'bytes={n}-'} if n else {},stream=True,timeout=(30,90)) as r:
                    r.raise_for_status()
                    if n and r.status_code!=206: n=0
                    last=time.monotonic()
                    with part.open('ab' if n else 'wb') as out:
                        for block in r.iter_content(1024*1024):
                            out.write(block); n+=len(block)
                            if time.monotonic()-last>10:
                                print(f'Download {n/1e6:.1f}/{f["size"]/1e6:.1f} MB',flush=True); last=time.monotonic()
                assert n==f['size'], (n,f['size'])
                part.replace(dest); break
            except Exception as e:
                print('Download retry',attempt,repr(e),flush=True); time.sleep(min(3*(attempt+1),20))
        else: raise RuntimeError('Download unsuccessful; partial file retained')
    h=hashlib.md5()
    with dest.open('rb') as src:
        for b in iter(lambda:src.read(8*1024*1024),b''):h.update(b)
    assert dest.stat().st_size==f['size'] and h.hexdigest()==f['computed_md5']
    print('Publisher size and MD5 PASS',flush=True)
    return dest

def cache(archive):
    import csv, ijson
    cohort={r['qmof_id']:int(r['n_atoms']) for r in csv.DictReader((ROOT/'data/chamm_common_cohort.csv').open())}
    out=ROOT/'work/structures'; out.mkdir(exist_ok=True)
    manifest=ROOT/'work/source_cache_manifest.json'
    if manifest.exists():
        m=json.loads(manifest.read_text())
        assert m['source_md5']=='0d89aaf66f2c306e86e47fd91cb1346e' and m['structures']==9171 and m['atoms']==1420532
        assert len(list(out.glob('qmof-*.json')))==9171
        print('Reuse completed pinned geometry cache',flush=True);return
    start=time.monotonic(); seen=set(); total=0; atoms=0
    with zipfile.ZipFile(archive) as z, z.open('qmof_database/qmof_structure_data.json') as src:
        assert z.getinfo('qmof_database/qmof_structure_data.json').file_size==3153444757
        for r in ijson.items(src,'item',use_float=True):
            total+=1
            if r['qmof_id'] in cohort:
                q=r['qmof_id']; sites=r['structure']['sites']; assert len(sites)==cohort[q]
                # Deliberately exclude electronic properties from the geometry cache.
                geo={'qmof_id':q,'name':r['name'],'structure':{'lattice':r['structure']['lattice'],'sites':[
                    {k:v for k,v in s.items() if k in ('species','abc','xyz','label')} for s in sites]}}
                p=out/(q+'.json'); tmp=p.with_suffix('.tmp'); tmp.write_text(json.dumps(geo,separators=(',',':'))); tmp.replace(p)
                seen.add(q); atoms+=len(sites)
            if total%2000==0: print(f'Streamed {total}; cached {len(seen)} structures/{atoms} atoms',flush=True)
    assert total==20372 and seen==cohort.keys() and atoms==1420532
    report=dict(source_md5='0d89aaf66f2c306e86e47fd91cb1346e',structures=len(seen),atoms=atoms,total_source=total,seconds=time.monotonic()-start,properties_removed=True)
    (ROOT/'work/source_cache_manifest.json').write_text(json.dumps(report,indent=2));print(report,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--download-only',action='store_true'); a=p.parse_args()
    archive=download()
    if not a.download_only:cache(archive)
