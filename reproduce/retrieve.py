"""Fetch the published production deposit, verify SHA256, and safely unpack it."""
from pathlib import Path
import argparse,hashlib,json,urllib.request,zipfile
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--archive',type=Path);p.add_argument('--verify-only',action='store_true');a=p.parse_args()
 m=json.loads((ROOT/'reproduce/deposit.json').read_text());archive=a.archive or ROOT/'downloads'/m['file']
 archive.parent.mkdir(exist_ok=True)
 if not archive.exists():
  part=archive.with_suffix('.partial')
  with urllib.request.urlopen(m['url'],timeout=120) as src,part.open('wb') as dst:
   while b:=src.read(8*1024*1024):dst.write(b)
  part.replace(archive)
 assert archive.stat().st_size==m['bytes'] and sha(archive)==m['sha256'],'Deposit integrity failure'
 if not a.verify_only:
  with zipfile.ZipFile(archive) as z:
   for n in z.namelist():
    target=(ROOT/n).resolve()
    assert ROOT.resolve() in target.parents,'Unsafe archive member'
   z.extractall(ROOT)
  for f in json.loads((ROOT/'reproduce/production_files.json').read_text()):
   path=ROOT/f['path'];assert path.stat().st_size==f['bytes'] and sha(path)==f['sha256'],f['path']
 print('PASS: production deposit SHA256'+(' and all extracted file hashes' if not a.verify_only else ''))
if __name__=='__main__':main()
