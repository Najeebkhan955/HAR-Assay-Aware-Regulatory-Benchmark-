from pathlib import Path
import base64,zlib,nbformat,hashlib
R=Path('/mnt/data/HAR_Stage16_Comparative_TF');nb=nbformat.v4.new_notebook()
files=['stage16_full180_no_drive.py','S1_All180_HAR_Sequence_Mapping_Readiness.csv','S2b_Human_NSC_Bait_Fragments_Embedded_for_Colab.csv']
P={x:base64.b64encode(zlib.compress((R/x).read_bytes(),9)).decode() for x in files}
head="""# Stage 16B — Comparative HAR–TF motif scoring | no Google Drive

**Goal:** extend Stage 16 preliminary 10-HAR motif pilot to up to 180 screened HARs using the publicly shared human–chimp HAR sequence catalog and the named Pal Capture-HiC bait fragments.

**Publication caution:** a named Capture-HiC bait is a restriction fragment, not the exact HAR sequence. Coordinate overlap is a provisional mapping; catalog identifiers `HAR_...` are not Noble-screen identifiers. Motif scoring evaluates sequence preference only, not ChIP occupancy, enhancer activity, causal HAR–gene linkage or confirmed human-specific TF binding. Some loci will be flagged as unmapped/ambiguous; do not force them into the analysis.

This notebook contains the 180-screen label manifest, exact named human bait fragments, and analysis script directly. No personal Google Drive account or mounting required. It downloads a small *public* third-party dataset pinned to a reproducible GitHub commit (2026-06-24), and produces CSV data, QC logs and a ZIP. It never needs the multi-GB single-cell expression file.
"""
nb.cells.append(nbformat.v4.new_markdown_cell(head))
nb.cells.append(nbformat.v4.new_code_cell('!pip -q install "biopython>=1.81" requests matplotlib numpy'))
nb.cells.append(nbformat.v4.new_markdown_cell('## 1. Unpack the embedded reproducibility inputs'))
nb.cells.append(nbformat.v4.new_code_cell('''import os,base64,zlib, pathlib, json, hashlib, requests, subprocess, zipfile, traceback, time
from pathlib import Path
ROOT=Path('/content/HAR_Stage16B_NoDrive'); ROOT.mkdir(exist_ok=True)
PAYLOADS = '''+repr(P)+'''
for filename,b64 in PAYLOADS.items():
    data=zlib.decompress(base64.b64decode(b64))
    (ROOT/filename).write_bytes(data)
    print(filename, len(data), 'bytes, sha256:',hashlib.sha256(data).hexdigest()[:16])
'''))
nb.cells.append(nbformat.v4.new_markdown_cell('## 2. Fetch pinned public BED / FASTA / JASPAR2024 MEME files'))
nb.cells.append(nbformat.v4.new_code_cell('''SOURCE_DIR=ROOT/'public_source';SOURCE_DIR.mkdir(exist_ok=True)
RESULTS=ROOT/'results';RESULTS.mkdir(exist_ok=True)
COMMIT='a9251f39976bbc15df1d8f0dfedaf9786819c298'
BASE=f'https://raw.githubusercontent.com/athenamarou/HAR-TFBS-Project/{COMMIT}/data/'
SOURCES={'hars_hg38.bed':105973,'hars_panTro6.bed':105765,'human_hars.fasta':945122,'chimp_hars.fasta':944312,'meme_file.txt':552150}
source_status=[];errors=[]
for name,expected_min_bytes in SOURCES.items():
    dst=SOURCE_DIR/name
    try:
        if not dst.exists() or dst.stat().st_size < expected_min_bytes*0.9:
            result=requests.get(BASE+name,timeout=(25,120))
            result.raise_for_status()
            if len(result.content)<expected_min_bytes*0.9:
                raise ValueError('Unexpected small response; check public repository')
            dst.write_bytes(result.content)
        record={'name':name,'bytes':dst.stat().st_size,'sha256':hashlib.sha256(dst.read_bytes()).hexdigest(),'source':BASE+name,'status':'OK'}
        print(name,record['bytes'],'bytes downloaded/verified')
    except Exception as exc:
        record={'name':name,'status':'ERROR','detail':str(exc)};errors.append(name)
        print('FAILED:',name,exc)
    source_status.append(record)
(RESULTS/'SOURCE_CHECKSUMS.json').write_text(json.dumps(source_status,indent=2))
print('Source errors:',errors)
'''))
nb.cells.append(nbformat.v4.new_markdown_cell('## 3. Strict mapping, alignment and six-TF PWM scoring'))
nb.cells.append(nbformat.v4.new_code_cell('''if errors:
    print('One or more source files could not be obtained. No motif statistics will be claimed. The final diagnostic ZIP will still be generated.')
    run_status={'status':'PUBLIC_SOURCE_DOWNLOAD_FAILED','missing_sources':errors}
else:
    command=['python',str(ROOT/'stage16_full180_no_drive.py'),'--input',str(SOURCE_DIR),'--bait',str(ROOT/'S2b_Human_NSC_Bait_Fragments_Embedded_for_Colab.csv'),'--screen',str(ROOT/'S1_All180_HAR_Sequence_Mapping_Readiness.csv'),'--output',str(RESULTS)]
    print('Running:', ' '.join(command))
    result=subprocess.run(command,capture_output=True,text=True)
    (RESULTS/'RUN_stdout.log').write_text(result.stdout)
    (RESULTS/'RUN_stderr.log').write_text(result.stderr)
    print('Exit code:',result.returncode)
    print(result.stdout[-2000:])
    print(result.stderr[-2500:])
    run_status={'status':'SUCCESS' if result.returncode==0 else 'RUN_FAILED','exit_code':result.returncode}
(RESULTS/'RUN_STATUS.json').write_text(json.dumps(run_status,indent=2))
'''))
nb.cells.append(nbformat.v4.new_markdown_cell('## 4. Inspect QC; export output ZIP (small files only)'))
nb.cells.append(nbformat.v4.new_code_cell('''if (RESULTS/'STATUS.json').exists():
    stats=json.loads((RESULTS/'STATUS.json').read_text())
    print(json.dumps(stats,indent=2))
else:
    print('Main STATUS.json unavailable. Review diagnostic logs in ZIP.')
zip_path=ROOT/'HAR_Stage16B_Ortholog_PWM_Results.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for file in RESULTS.rglob('*'):
        if file.is_file():z.write(file,arcname=file.relative_to(ROOT))
    for file in [ROOT/'stage16_full180_no_drive.py',ROOT/'S1_All180_HAR_Sequence_Mapping_Readiness.csv',ROOT/'S2b_Human_NSC_Bait_Fragments_Embedded_for_Colab.csv']:
        z.write(file,arcname='inputs/'+file.name)
print('ZIP ready:',zip_path,'size:',zip_path.stat().st_size)
try:
    from google.colab import files as colab_files
    colab_files.download(str(zip_path))
except Exception as e:
    print('Not running in Colab? ZIP remains saved at:',zip_path)
'''))
nb.metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python'}}
out=R/'HAR_Stage16B_Full180_Ortholog_TF_NoDrive_Colab.ipynb';nbformat.write(nb,str(out));print(out,'bytes',out.stat().st_size)
import ast
for c in nb.cells:
 if c.cell_type=='code' and not c.source.startswith('!'):ast.parse(c.source)
nbformat.validate(nb)
print('All code cells syntax-validated; runtime tests performed separately with synthetic data.')
