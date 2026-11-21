#!/usr/bin/env python3
"""Download source originals from NCBI GEO; not a claim they were uploaded to GitHub.

Requires: requests. Run separately if needed (12 original RNA/DNA files ~100MB).
Checks HTTP, compressed-format signature and records SHA256. Uses original GEO sample URLs.
"""
from pathlib import Path
import csv,hashlib,requests,sys,time
HERE=Path(__file__).resolve().parent
OUT=HERE/'downloaded_original_geo'
OUT.mkdir(exist_ok=True)
manifest=[]
with (HERE.parent.parent/'08_manifests/GEO_SOURCE_DOWNLOADS.csv').open(newline='',encoding='utf-8') as f:
    rows=list(csv.DictReader(f))
for row in rows:
    p=OUT/row['local_download_basename']
    if not p.is_file():
        print('Fetch:',row['GSM'])
        with requests.get(row['original_geo_url'],stream=True,timeout=(25,120),headers={'User-Agent':'HAR-Research-Reproducible-Sources/1.0'}) as r:
            r.raise_for_status()
            with p.open('wb') as output:
                for chunk in r.iter_content(1024*1024):
                    if chunk:output.write(chunk)
    with p.open('rb') as f:
        assert f.read(2)==b'\x1f\x8b', 'Not gzip: '+str(p)
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    manifest.append([row['GSM'],p.name,p.stat().st_size,h,row['original_geo_url']])
with (OUT/'SOURCES_SHA256.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f);w.writerow(['GSM','filename','bytes','SHA256','source_url']);w.writerows(manifest)
print('Downloaded/verified:',len(manifest),'GEO raw files; never treat these culture samples as distinct donors')
