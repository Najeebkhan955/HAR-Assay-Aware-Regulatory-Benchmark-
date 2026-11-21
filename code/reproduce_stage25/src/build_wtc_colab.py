"""Build a self-contained, no-Google-Drive WTc N2 external line notebook from verified local source archives."""
import gzip, zipfile, hashlib, io, base64, json, nbformat
import pandas as pd
from pathlib import Path
root=Path(__file__).resolve().parents[1]
with zipfile.ZipFile(root/'data/source_archives/HAR_Stage21C_MPRA_HS111_N2_Results.zip') as z:
    mapping=z.read('S1B_Assayed_Barcode_Map.csv')
assert hashlib.sha256(mapping).hexdigest()=='6f30451ac623ed4948c8649aa3e9113f0d244b5d70db1f37fbb35314cf92cdeb'
with zipfile.ZipFile(root/'data/source_archives/HAR_Stage21D_FullSixRaw_Corrected_Results.zip') as z:
    df=pd.read_csv(io.BytesIO(z.read('S2_All_Cut_FullRaw_Sensitivity.csv')))
ref=df[(df.min_DNA_count_per_barcode==1)&(df.aggregate_pseudocount==1)]
assert len(ref)==342
baseline=ref[['HAR','cut','culture','corrected_HminusC_log2']].sort_values(['HAR','cut','culture']).to_csv(index=False).encode('utf-8')
script=(root/'src/stage23_wtc_n2.py').read_text()
assert 'if __name__' in script
cell="""# Independent human iPSC line WTc, N2; only public GEO source files downloaded.
# This cell embeds pinned Stage21C barcode-to-construct map and corrected HS1-11 baseline
# and performs the WTc analysis. No Google Drive authentication or subprocess is needed.
import pathlib, base64, gzip, runpy
ROOT=pathlib.Path('/content') if pathlib.Path('/content').exists() else pathlib.Path.cwd()
I=ROOT/'HAR_Stage23_WTc_Embedded';I.mkdir(exist_ok=True)
"""
for name,payload in [('S1B_Assayed_Barcode_Map.csv',mapping),('S2_HS111_Baseline_342.csv',baseline)]:
    comp=base64.b64encode(gzip.compress(payload,compresslevel=9,mtime=0)).decode('ascii')
    cell+=f"\n(I/'{name}').write_bytes(gzip.decompress(base64.b64decode('{comp}')))\n"
cell+="\n"+"S=ROOT/'stage23_wtc_n2.py'\nS.write_text("+repr(script)+",encoding='utf-8')\nrunpy.run_path(str(S),run_name='__main__')\n"
download="""# User-requested automatic download; main notebook cell, not a subprocess.
from pathlib import Path
root=Path('/content') if Path('/content').exists() else Path.cwd()
candidates=[root/'HAR_Stage23_WTc_N2_IndependentLine_Results.zip',root/'HAR_Stage23_WTc_N2_Diagnostics.zip']
found=next((p for p in candidates if p.exists()),None)
if not found:
    print('No ZIP found. Inspect previous cell or Colab Files sidebar.')
else:
    print('Ready:',found.name,'size',found.stat().st_size,'bytes')
    try:
        from google.colab import files
        files.download(str(found))
        print('Automatic download requested. If your browser blocks it, use the Colab Files sidebar.')
    except ImportError:
        print('Outside Colab. File available at:',found)
    except Exception as e:
        print('Automatic download blocked:',str(e)[:300]);print('Manual backup:',found)
"""
nb=nbformat.v4.new_notebook(cells=[
    nbformat.v4.new_markdown_cell('# Essential external validation: WTc/N2 human iPSC line\n\nThis runs **only** the missing independent-human-line lentiMPRA comparison against already-verified HS1-11/N2. It does not rerun earlier motif, guide or chromatin analyses. Human and chimp refer to **insert alleles**, not chimp host cells. No primary HAR BED identity or target gene is confirmed.\n\nThe notebook includes a checksum-verified existing 25,500-barcode-to-insert mapping and the archived 342 HS1-11 corrected contrasts; it downloads the six **WTc N2 GEO raw count files**, hashes them, parses all headerless records, computes descriptive per-culture effects, and compares only matching HAR/cut constructs. Actual WTc values cannot be known until the live GEO download runs.'),
    nbformat.v4.new_code_cell(cell),nbformat.v4.new_code_cell(download),
    nbformat.v4.new_markdown_cell('After download, upload `HAR_Stage23_WTc_N2_IndependentLine_Results.zip` here. On failure upload `HAR_Stage23_WTc_N2_Diagnostics.zip`. This is the only recommended additional Colab run; do not upload any personal files.')],
  metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python'}})
nbformat.validate(nb)
out=root/'notebooks/HAR_Stage23_WTc_N2_IndependentLine_AutoDownload_NoDrive.ipynb'
nbformat.write(nb,out)
for c in nb.cells:
    if c.cell_type=='code':compile(c.source,'cell','exec')
print('Notebook validated',out,out.stat().st_size,'bytes')
