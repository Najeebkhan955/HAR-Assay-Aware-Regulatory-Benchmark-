import csv, hashlib, json, os, pathlib, shutil, zipfile, textwrap, sys
root=pathlib.Path('/mnt/data')
source=root/'HAR_Complete_Stage1-to25_Available_Code_Data_Results_GitHub_Upload.zip'
stage=root/'HAR_GitHub_Recovery_And_Upload_v2'
if stage.exists(): shutil.rmtree(stage)
stage.mkdir()
with zipfile.ZipFile(source) as z:
    z.extractall(stage)

# Construct expected hashes from ORIGINAL machine-generated source manifests.
expected={}
with zipfile.ZipFile(root/'HAR_Stage21D_FullSixRaw_Corrected_Results.zip') as z:
    rows=list(csv.DictReader(z.read('S1_Six_Raw_Sources_QC.csv').decode('utf8').splitlines()))
    for r in rows: expected[r['filename']]=r['sha256']
with zipfile.ZipFile(root/'HAR_Stage24_WTc_FinalValidation_Results.zip') as z:
    rows=list(csv.DictReader(z.read('S2_SixGEO_SourceIntegrity.csv').decode('utf8').splitlines()))
    for r in rows: expected[r['file']]=r['sha256']
expected['GSE110758_HAR_design.tsv.gz']='c0552a2fc675f043017252a8a341d6e06cd91cedefdc7b65a0186118901a0d3f'
with (stage/'08_manifests'/'GEO_SOURCE_DOWNLOADS.csv').open() as f:
    urlmap={r['local_download_basename']:r for r in csv.DictReader(f)}
urlmap['GSE110758_HAR_design.tsv.gz']={'original_geo_url':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110758/suppl/GSE110758_HAR_design.tsv.gz','GSM':'GSE110758','sample_label':'MPRA HAR design','local_download_basename':'GSE110758_HAR_design.tsv.gz'}
assert set(expected)==set(urlmap) and len(expected)==13
# Existing file byte SHA verify, strong evidence.
raw=stage/'02_data/raw_GEO_original_uploads'
for f in raw.glob('*.tsv.gz'):
    assert hashlib.sha256(f.read_bytes()).hexdigest()==expected[f.name],f'Incorrect original raw file hash {f}'
with (stage/'08_manifests'/'GEO_EXPECTED_SOURCE_SHA256.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.writer(f);w.writerow(['GEO_accession','sample_label','raw_filename','source_url','expected_sha256','initially_present_in_archive'])
    for fname in sorted(urlmap):
        r=urlmap[fname];w.writerow([r['GSM'],r['sample_label'],fname,r['original_geo_url'],expected[fname],(raw/fname).exists()])

script=r'''#!/usr/bin/env python3
"""Recover only verified public GSE110758 original RNA/DNA barcode files.

Run from the repository's root, or simply double-click the Windows batch file.
Copies nothing to GitHub itself: the .bat wrapper calls Git after validation.
Expected SHA256 are from Stage21D and Stage24 executed GEO source manifests.
"""
import csv, gzip, hashlib, pathlib, sys, time, urllib.request, urllib.error

BASE=pathlib.Path(__file__).resolve().parents[1]
TABLE=BASE/'08_manifests'/'GEO_EXPECTED_SOURCE_SHA256.csv'
DEST=BASE/'02_data'/'raw_GEO_original_uploads'
DEST.mkdir(parents=True,exist_ok=True)
TIMEOUT=120

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as fp:
        for block in iter(lambda:fp.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()

def fetch(row):
    name=row['raw_filename'];target=DEST/name
    expected=row['expected_sha256'].lower()
    if target.is_file():
        actual=sha256(target)
        if actual!=expected:
            raise RuntimeError(f'ALREADY PRESENT BUT WRONG CHECKSUM: {name} {actual} expected {expected}. Resolve manually; original will not be overwritten.')
        print(f'[verified existing] {name} {target.stat().st_size:,} bytes',flush=True)
        return 'existing_verified'
    for attempt in range(1,4):
        part=target.with_suffix(target.suffix+'.incomplete')
        try:
            print(f'[download attempt {attempt}/3] {name}',flush=True)
            req=urllib.request.Request(row['source_url'], headers={'User-Agent':'HAR-Audit-Archive/1.0 (public GEO data, source cited)'})
            with urllib.request.urlopen(req,timeout=TIMEOUT) as response,open(part,'wb') as output:
                while True:
                    chunk=response.read(1024*1024)
                    if not chunk:break
                    output.write(chunk)
                    if output.tell()>80*1024*1024:
                        raise RuntimeError('Unexpected >80MB original source; check GEO source URL')
            if part.stat().st_size < 500 or part.open('rb').read(2) != b'\x1f\x8b':
                raise RuntimeError('Not a valid gzip source')
            actual=sha256(part)
            if actual!=expected:
                raise RuntimeError(f'SHA-256 mismatch for {name}: {actual} != {expected}')
            part.replace(target)
            print(f'[verified downloaded] {name} {target.stat().st_size:,} bytes',flush=True)
            return 'downloaded_verified'
        except Exception as exc:
            if part.exists():part.unlink()
            print('[failed]',str(exc),flush=True)
            if attempt<3:time.sleep(2*attempt)
            else:raise

def main():
    rows=list(csv.DictReader(TABLE.open(encoding='utf8',newline='')))
    if len(rows)!=13:raise RuntimeError(f'Expected exactly 13 source files, got {len(rows)}')
    results=[];failures=[]
    for row in rows:
        try:status=fetch(row)
        except Exception as exc:
            status='FAILED'; failures.append((row['raw_filename'],str(exc)))
        p=DEST/row['raw_filename']
        results.append([row['GEO_accession'],row['raw_filename'],status,p.stat().st_size if p.exists() else '',sha256(p) if p.exists() else '',row['expected_sha256'],row['source_url']])
    audit=BASE/'08_manifests'/'GEO_RAW_DOWNLOAD_VERIFICATION_RESULTS.csv'
    with audit.open('w',newline='',encoding='utf8') as f:
        w=csv.writer(f);w.writerow(['accession','file','status','bytes','actual_sha256','expected_sha256','source_url']);w.writerows(results)
    print(f'\nSources verified: {len(rows)-len(failures)}/{len(rows)}; audit={audit}',flush=True)
    if failures:
        print('Remaining failures:',flush=True)
        for name,err in failures:print(f'- {name}: {err}',flush=True)
        return 1
    return 0

if __name__=='__main__':sys.exit(main())
'''
(stage/'upload_tools'/'recover_raw_geo_and_verify.py').write_text(script,encoding='utf8')
launcher=r'''@echo off
cd /d "%~dp0"
where py >nul 2>&1
if %errorlevel% equ 0 (
    py -3 upload_tools\recover_raw_geo_and_verify.py
) else (
    where python >nul 2>&1
    if %errorlevel% neq 0 (
        echo Python 3 is required: https://www.python.org/downloads/windows/
        pause
        exit /b 1
    )
    python upload_tools\recover_raw_geo_and_verify.py
)
if errorlevel 1 (
    echo.
    echo Raw source download/checksum step was incomplete. No push started.
    echo You may rerun this script. It will not redownload validated files.
    pause
    exit /b 1
)
echo.
echo All 13 GEO files verified. Initiating GitHub push with Git Credential Manager...
call UPLOAD_TO_GITHUB_Windows.bat
'''
(stage/'RECOVER_RAW_AND_UPLOAD_Windows.bat').write_text(launcher,encoding='utf8')
notes='''# HAR research archive — recovery and GitHub upload

**Current status:** GitHub connector can read `Najeebkhan955/HAR-Research-`, but direct writes return GitHub HTTP 403 `Resource not accessible by integration`. This is a GitHub App authorization blocker, **not** a data-analysis error. Nothing has been uploaded by ChatGPT to GitHub.

## One-click workflow with your own GitHub login
1. Extract this ZIP to a local folder on Windows.
2. Install **Python 3** and **Git for Windows** if not installed. Git for Windows includes a browser-based credential manager.
3. Double-click `RECOVER_RAW_AND_UPLOAD_Windows.bat` at the top level.
4. The script downloads **nine missing original GEO gzip files**, verifies their exact expected SHA-256, checks four already-present originals, and writes `08_manifests/GEO_RAW_DOWNLOAD_VERIFICATION_RESULTS.csv`.
5. After **all 13 original reporter source files** pass, the batch file starts the prior safe Git push script. Complete the GitHub browser sign-in yourself. Never send credentials or personal-access tokens in chat.
6. Verify visible files at https://github.com/Najeebkhan955/HAR-Research- . The batch **never force-pushes**.

### Scientifically essential limitations
- Original Stage 1–6 standalone notebooks, full GSE270828 87,896-cell single-cell raw data, and author-defined original HAR target BED **were not present** in the uploaded local archive; they are not substituted with fabricated or reconstructed original files.
- GSE110758 raw barcode counts comprise 12 source files across HS1-11/WTc at N2 (three paired RNA/DNA cultures each) plus original HAR design (13 files). The collection before running the script includes **four** of these; the script retrieves the remaining nine directly from NCBI GEO, using the original analyses' recorded checksums.
- The repository also includes derived source/count evidence, historical ZIPs, reproducible scripts, validation results and audit trails. Final 10–12 publication supplements are **not yet frozen**; their folder is reserved.
- Original GEO/other third-party sources retain attribution. Source data are available via NCBI; check external source reuse terms and licenses before publishing redistributed material. For raw single-cell matrices too large for GitHub, retain original GEO accession or use archival repositories with stable DOIs, rather than committing multi-gigabyte files.
- No paper claims experimental TF binding or causal HAR–gene effects from these reporter constructs.

If browser login fails or GitHub refuses push, use the GitHub app settings to review installation access, or supply the **error message only**, not any token.
'''
(stage/'START_HERE_RECOVER_RAW_AND_UPLOAD.md').write_text(notes,encoding='utf8')
# Extend existing README with accurate instructions.
readme=stage/'README.md';readme.write_text(readme.read_text(encoding='utf8')+'\n\n## Current upload procedure and raw GEO recovery\n\nSee [`START_HERE_RECOVER_RAW_AND_UPLOAD.md`](START_HERE_RECOVER_RAW_AND_UPLOAD.md). User-authenticated Windows recovery + push will fetch/check nine missing original GSE110758 raw files; GitHub connector writes presently receive HTTP 403. The existing GitHub repository is not yet populated by these files.\n',encoding='utf8')
# Build complete checksum manifest EXCLUDING its own content, no recursion/circular integrity.
manifest=stage/'08_manifests'/'FILE_INDEX_SHA256.csv'
entries=[]
for f in stage.rglob('*'):
    if not f.is_file() or f==manifest: continue
    rel=f.relative_to(stage).as_posix()
    # Ignore platform metadata and intermediate files
    entries.append((rel,f.stat().st_size,hashlib.sha256(f.read_bytes()).hexdigest()))
with manifest.open('w',newline='',encoding='utf8') as fh:
    w=csv.writer(fh);w.writerow(['repository_path','bytes','sha256']);w.writerows(sorted(entries))
# Modify original stage completeness metadata without claiming new files already downloaded.
p=stage/'08_manifests'/'STAGE_COVERAGE_AND_COMPLETENESS.json'
obj=json.loads(p.read_text(encoding='utf8'))
obj['geo_reporter_raw_expected_including_design']=13
obj['geo_reporter_raw_present_prior_to_recovery']=4
obj['geo_reporter_raw_missing_prior_to_recovery']=9
obj['github_connector_write_status']='403 Resource not accessible by integration; not uploaded'
p.write_text(json.dumps(obj,indent=2),encoding='utf8')
# re-manifest changed JSON
def refresh_manifest():
    entries=[]
    for f in stage.rglob('*'):
        if not f.is_file() or f==manifest:continue
        entries.append((f.relative_to(stage).as_posix(),f.stat().st_size,hashlib.sha256(f.read_bytes()).hexdigest()))
    with manifest.open('w',newline='',encoding='utf8') as fh:
        w=csv.writer(fh);w.writerow(['repository_path','bytes','sha256']);w.writerows(sorted(entries))
    return len(entries)
count=refresh_manifest()
# Validate all checksums and compilation.
import py_compile
py_compile.compile(str(stage/'upload_tools'/'recover_raw_geo_and_verify.py'),doraise=True)
# Clean bytecode cached file
for pyc in (stage/'upload_tools').rglob('*.pyc'):pyc.unlink()
for d in list((stage/'upload_tools').rglob('__pycache__')):d.rmdir()
count=refresh_manifest()
with manifest.open() as f:
 rows=list(csv.DictReader(f))
assert len(rows)==count
for row in rows:
    f=stage/row['repository_path'];assert hashlib.sha256(f.read_bytes()).hexdigest()==row['sha256']
# Build zip with no path prefix, executable needed to run locally.
output=root/'HAR_Complete_Research_GitHub_RawRecovery_AutoUpload_v2.zip'
with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
    for f in sorted(stage.rglob('*')):
        if f.is_file():z.write(f,f.relative_to(stage).as_posix())
with zipfile.ZipFile(output) as z:
 assert z.testzip() is None
print(json.dumps({'new_archive':str(output),'archive_MB':round(output.stat().st_size/2**20,2),'file_count':count+1,'verified_existing_raw':len(list(raw.glob('*.tsv.gz'))),'raw_expected':len(expected),'raw_to_recover':len(expected)-len(list(raw.glob('*.tsv.gz'))),'readme':str(stage/'START_HERE_RECOVER_RAW_AND_UPLOAD.md')},indent=2))
