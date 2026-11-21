import json,zipfile,io,csv,nbformat,hashlib,pathlib
from nbformat.v4 import new_notebook,new_code_cell,new_markdown_cell
root=pathlib.Path('/mnt/data')
# Use the 180 original HAR identifiers from Stage 19B source crosswalk, not speculative names.
with open(root/'HAR_Stage19B_Outputs/S1_180HAR_Stage19B_DesignRule_Crosswalk.csv',encoding='utf-8-sig',newline='') as f:
    rows=list(csv.DictReader(f))
print('Input crosswalk headers:',rows[0].keys());print('input rows',len(rows))
key=next(k for k in rows[0] if k.lower() in ('screen_har','har','har_name'))
names=sorted(set(r[key] for r in rows if r[key]))
assert len(names)==180,len(names)

preamble='''# Stage 21B — Independent HAR lentiMPRA design-overlap preflight (NO Google Drive)

**Purpose:** Cross-check the original 180 GSE270828 HAR **names** against the published GSE110758 lentiMPRA design table (714 HARs). This notebook checks identifier availability only. It does **NOT** reproduce RNA/DNA enhancer activity, establish gene targets, validate TF binding, or certify genomic coordinate / ortholog equivalence.

- Original selected HARs: GSE270828, 180 names, from the verified internal crosswalk.
- External assay: GSE110758_HAR_design.tsv.gz from GEO (human/chimp neural progenitor lentiMPRA, oligo design metadata).
- All possible matches are **name-level hypotheses pending a curated ID and sequence concordance**.
- GEO archive is downloaded only if necessary; no Drive mount; no large RNA/DNA count matrices.
- Auto-outputs: `HAR_Stage21B_MPRA_Design_Preflight.zip` or `HAR_Stage21B_MPRA_Diagnostics.zip`, using plain Python/Colab Files sidebar download (avoids subprocess `google.colab.files.download()` error).
'''
names_literal=json.dumps(names,ensure_ascii=False)
code=f'''import csv, io, gzip, hashlib, json, re, time, traceback, zipfile, requests
from pathlib import Path
import pandas as pd

ROOT=Path('/content') if Path('/content').exists() else Path.cwd()
OUT=ROOT/'HAR_Stage21B_Preflight'
OUT.mkdir(exist_ok=True)
NAMES={names_literal}
assert len(NAMES)==180 and len(set(NAMES))==180
URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110758/suppl/GSE110758_HAR_design.tsv.gz'
DEST=ROOT/'GSE110758_HAR_design.tsv.gz'

def write_csv(path,headers,items):
    with open(path,'w',encoding='utf-8',newline='') as f:
        w=csv.writer(f);w.writerow(headers);w.writerows(items)

def pack(name):
    target=ROOT/name
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as z:
        for p in OUT.iterdir():
            if p.is_file():z.write(p,p.name)
    print('Saved:',target,'size:',target.stat().st_size)
    print('Colab left sidebar > Files > Download', target.name)
    return target

try:
    if not DEST.exists():
        errors=[]
        for attempt in range(1,4):
            try:
                with requests.get(URL,headers={{'User-Agent':'HAR-Research-Reproducibility/1.0'}},stream=True,timeout=(25,120)) as resp:
                    resp.raise_for_status()
                    with open(DEST,'wb') as f:
                        for chunk in resp.iter_content(chunk_size=1024*1024):
                            if chunk:f.write(chunk)
                if DEST.stat().st_size>20_000_000:
                    raise ValueError('Unexpectedly large file; inspect the source manually')
                if DEST.stat().st_size<1000 or DEST.open('rb').read(2)!=b'\\x1f\\x8b':
                    raise ValueError('Downloaded file is not a plausible gzip TSV')
                break
            except Exception as ex:
                errors.append(str(ex));print('Attempt',attempt,'failed:',ex)
                if DEST.exists():DEST.unlink()
                time.sleep(2)
        else:
            raise RuntimeError('GEO design download failed after 3 tries: '+repr(errors))
    payload=DEST.read_bytes()
    if payload[:2]!=b'\\x1f\\x8b':raise ValueError('Invalid GZIP signature')
    raw=gzip.decompress(payload)
    df=pd.read_csv(io.BytesIO(raw),sep='\\t',dtype=str,keep_default_na=False,low_memory=False,on_bad_lines='skip')
    if df.shape[1]<2 or df.empty:
        raise ValueError('Unexpected GEO design schema; do not infer a match from this file')
    print('GSE110758 Design:',df.shape,'Columns:',list(df.columns))
    # Only identifier/name-like columns; do not interpret nucleotide sequences as identifiers.
    patcols=[c for c in df.columns if any(token in str(c).lower() for token in ('name','har','element','id','oligo','annot','target','region','label'))
             and not any(token in str(c).lower() for token in ('seq','barcode','read','count'))]
    if not patcols:
        patcols=[c for c in df.columns if not any(tok in str(c).lower() for tok in ('seq','barcode','read','count','rna','dna'))][:min(10,len(df.columns))]
    pattern=re.compile(r'(?<![a-zA-Z0-9.])(?:'+'|'.join(sorted((re.escape(s) for s in NAMES),key=len,reverse=True))+r')(?![a-zA-Z0-9.])',re.I)
    matched={{s:[] for s in NAMES}}
    scans=df[patcols].fillna('').astype(str).agg('|'.join,axis=1)
    for i,t in enumerate(scans):
        hits={{h.casefold():h for h in pattern.findall(t)}}
        for h in hits.values():
            real=next((n for n in NAMES if n.casefold()==h.casefold()),None)
            if real is not None:matched[real].append(i)
    nmatched=sum(bool(v) for v in matched.values())
    write_csv(OUT/'S1_180HAR_NameOverlap_Preflight.csv',
       ['HAR','name_matched_in_GSE110758_design','n_design_rows_with_name','example_design_row_indices_0based','status'],
       [[n, bool(matched[n]),len(matched[n]),';'.join(map(str,matched[n][:5])),
       'NAME_LEVEL_ONLY_NOT_SEQUENCE_OR_ACTIVITY_VALIDATION' if matched[n] else 'NO_NAME_MATCH_NOT_PROOF_OF_NO_BIOLOGICAL_OVERLAP']
        for n in NAMES])
    write_csv(OUT/'S2_Source_Columns.csv',['column','description'],[[c,'Original GEO design field; no biological interpretation'] for c in df.columns])
    matched_indices=sorted({{i for ids in matched.values() for i in ids}})[:400]
    if matched_indices:
        df.iloc[matched_indices].head(400).to_csv(OUT/'S3_Example_Matched_Source_Rows.csv',index=False)
    summary={{'stage':'21B','status':'complete','external_GEO':'GSE110758','design_url':URL,
             'downloaded_bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest(),
             'design_rows':int(df.shape[0]),'design_columns':int(df.shape[1]),
             'identifier_columns_inspected':patcols,'HAR_names_from_GSE270828':len(NAMES),
             'strict_NAME_match_count':nmatched,'exact_name_match_is_NOT_interval_validation':True,
             'experimental_activity_computed':False,'independent_target_validation':False,
             'next':'Inspect actual design labels and sequences; verify biological/coordinate correspondence. Then obtain previously published validated RNA/DNA reporter-activity effect summaries before statistical benchmarking.'}}
    (OUT/'S0_Source_Status.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    pack('HAR_Stage21B_MPRA_Design_Preflight.zip')
    print(json.dumps(summary,indent=2))
except Exception as ex:
    (OUT/'S0_Diagnostics.txt').write_text(traceback.format_exc(),encoding='utf-8')
    (OUT/'S0_Error.json').write_text(json.dumps({{'status':'failed','error':str(ex),'design_url':URL}},indent=2),encoding='utf-8')
    pack('HAR_Stage21B_MPRA_Diagnostics.zip')
    print('Preflight failed; uploaded diagnostics will show why. No biological results were inferred.')
'''
# For review: ensure strings are valid, notebook metadata sane.
compile(code,'stage21b.py','exec')
nb=new_notebook(cells=[new_markdown_cell(preamble),new_code_cell(code)],metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'colab':{'name':'HAR_Stage21B_MPRA_Design_Preflight_NoDrive.ipynb'}})
out=root/'HAR_Stage21B_MPRA_Design_Preflight_NoDrive.ipynb'
nbformat.write(nb,out)
nbformat.validate(nb)
print('CREATED',out,out.stat().st_size)
