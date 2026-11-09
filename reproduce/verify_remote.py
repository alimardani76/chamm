"""Check the published asset metadata, server digest, and live download prefix."""
from pathlib import Path
import argparse,json,urllib.request
R=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--archive',type=Path);a=p.parse_args()
 m=json.loads((R/'reproduce/deposit.json').read_text())
 req=urllib.request.Request('https://api.github.com/repos/alimardani76/chamm/releases/tags/v1.0.0',headers={'Accept':'application/vnd.github+json','User-Agent':'ChaMM-reproduction'})
 with urllib.request.urlopen(req,timeout=60) as response:release=json.load(response)
 asset=next(x for x in release['assets'] if x['name']==m['file'])
 assert asset['state']=='uploaded' and asset['size']==m['bytes']
 if asset.get('digest'):assert asset['digest']=='sha256:'+m['sha256']
 req=urllib.request.Request(m['url'],headers={'Range':'bytes=0-4095','User-Agent':'ChaMM-reproduction'})
 with urllib.request.urlopen(req,timeout=60) as response:
  prefix=response.read(4096);status=response.status;cr=response.headers.get('Content-Range')
 assert prefix[:4]==b'PK\x03\x04'
 if a.archive:
  with a.archive.open('rb') as f:assert prefix==f.read(len(prefix))
 result=dict(status='PASS',release=release['html_url'],asset=asset['browser_download_url'],bytes=asset['size'],server_digest=asset.get('digest'),http_status=status,content_range=cr,prefix_bytes=len(prefix),scope='Public release metadata and live download-prefix check; full local archive and extracted files verified separately. Not a second complete network download.')
 (R/'reproduce/remote_asset_check.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
if __name__=='__main__':main()
