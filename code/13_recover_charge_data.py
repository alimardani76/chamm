"""Recover the charge-only data directly from the verified QMOF archive.
No geometry reconstruction, connectivity inference, or model fitting occurs.
Requires ijson; standard-library modules handle all other operations.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import sys
import time
import zipfile

KEYS = ('pbe_ddec_charge', 'pbe_cm5_charge', 'pbe_bader_charge')
EXPECTED = dict(structures=20372, all_atoms=2312743, common_structures=9171,
                common_atoms=1420532, ddec_complete=20372, cm5_complete=20372,
                bader_complete=9171)


def file_hash(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with path.open('rb') as source:
        for block in iter(lambda: source.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def recover(archive, metadata, benchmark, outdir):
    import ijson
    article = json.loads(metadata.read_text(encoding='utf-8'))
    published = next(f for f in article['files'] if f['name']=='qmof_database.zip')
    assert archive.stat().st_size == published['size'], 'Archive size does not match publisher metadata'
    assert file_hash(archive, 'md5') == published['computed_md5'], 'Publisher MD5 mismatch'
    refs = json.loads(benchmark.read_text(encoding='utf-8'))
    reference = {r['qmof_id']:r for r in refs}
    assert len(reference)==len(refs)==34, 'Expected the recovered 34-structure benchmark'
    seen = set(); matched = set()
    counts = {key:0 for key in EXPECTED}
    outdir.mkdir(parents=True, exist_ok=False)
    charge_part = outdir/'chamm_atom_charges_recovered.csv.gz.partial'
    manifest_part = outdir/'chamm_common_cohort.csv.partial'
    start = time.monotonic()
    member = 'qmof_database/qmof_structure_data.json'
    with zipfile.ZipFile(archive) as z, z.open(member) as stream, \
            gzip.open(charge_part,'wt',newline='',encoding='utf-8',compresslevel=4) as charges, \
            manifest_part.open('w',newline='',encoding='utf-8') as manifest:
        writer = csv.writer(charges)
        writer.writerow(['qmof_id','name','atom_index','element','ddec6','cm5','bader'])
        cohort = csv.writer(manifest)
        cohort.writerow(['qmof_id','name','n_atoms'])
        for record in ijson.items(stream,'item',use_float=True):
            q = record['qmof_id']
            assert q not in seen, f'Duplicate QMOF ID: {q}'
            seen.add(q)
            sites = record['structure']['sites']
            assert sites, f'Empty structure: {q}'
            counts['structures'] += 1
            counts['all_atoms'] += len(sites)
            complete = []
            for key in KEYS:
                present = all(s.get('properties',{}).get(key) is not None for s in sites)
                if present:
                    assert all(math.isfinite(float(s['properties'][key])) for s in sites), f'Non-finite charge: {q}, {key}'
                complete.append(present)
            for key, value in zip(('ddec_complete','cm5_complete','bader_complete'), complete):
                counts[key] += int(value)
            if q in reference:
                assert record==reference[q], f'Retained benchmark differs from recovered source: {q}'
                matched.add(q)
            if all(complete):
                counts['common_structures'] += 1
                counts['common_atoms'] += len(sites)
                cohort.writerow([q,record.get('name'),len(sites)])
                for i,site in enumerate(sites):
                    species = site['species']
                    assert len(species)==1 and species[0]['occu']==1, f'Disordered charge-table site: {q}/{i}'
                    p = site['properties']
                    writer.writerow([q,record.get('name'),i,species[0]['element'],
                        p[KEYS[0]],p[KEYS[1]],p[KEYS[2]]])
            if counts['structures'] % 2000 == 0:
                print(f"Read {counts['structures']:,} structures; recovered {counts['common_atoms']:,} common-cohort atoms",flush=True)
    checks = dict(counts_match_original=counts==EXPECTED, benchmark_matches=len(matched),
                  benchmark_missing=sorted(reference.keys()-matched))
    report = dict(status='RECOVERED' if counts==EXPECTED and matched==reference.keys() else 'REVIEW_REQUIRED',
                  observed=counts,expected=EXPECTED,checks=checks,
                  seconds=time.monotonic()-start,python=sys.version,ijson=importlib.metadata.version('ijson'),
                  source=dict(article_id=article['id'],article_version=article['version'],
                    file_name=published['name'],file_id=published['id'],size_bytes=published['size'],
                    md5=published['computed_md5'],download_url=published['download_url'],member=member),
                  benchmark_sha256=file_hash(benchmark),script_sha256=file_hash(Path(__file__)),
                  note='Minimal charge-only reconstruction. This is not an identical copy of the lost CSV schema. Derived differences are recomputable from the three charges. Source site ordering and common-cohort membership are retained. No connectivity is validated by this recovery.')
    if report['status'] != 'RECOVERED':
        (outdir/'recovery_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        raise RuntimeError('Recovery checks differ from the original checkpoint; partial outputs retained, not finalized')
    charge = outdir/'chamm_atom_charges_recovered.csv.gz'
    manifest = outdir/'chamm_common_cohort.csv'
    charge_part.rename(charge); manifest_part.rename(manifest)
    report['outputs'] = {p.name:dict(bytes=p.stat().st_size,sha256=file_hash(p)) for p in (charge,manifest)}
    (outdir/'recovery_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2),flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--benchmark',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    a = p.parse_args()
    recover(a.archive,a.metadata,a.benchmark,a.outdir)


if __name__=='__main__':
    main()
