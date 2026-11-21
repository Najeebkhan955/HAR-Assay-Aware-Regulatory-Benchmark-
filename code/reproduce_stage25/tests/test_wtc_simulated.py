"""Synthetic only: WTc code path smoke test without making biological claims or internet calls."""
import tempfile, zipfile, io, gzip, csv, hashlib, importlib.util, json
from pathlib import Path
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'src/stage23_wtc_n2.py';spec=importlib.util.spec_from_file_location('wtc_module',p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
with tempfile.TemporaryDirectory() as td:
 d=Path(td);m.ROOT=d;m.IN=d/'HAR_Stage23_WTc_Embedded';m.OUT=d/'HAR_Stage23_WTc_N2_Outputs';m.CACHE=d/'cache'
 for x in [m.IN,m.OUT,m.CACHE]:x.mkdir(exist_ok=True)
 with zipfile.ZipFile(ROOT/'data/source_archives/HAR_Stage21C_MPRA_HS111_N2_Results.zip') as z:
  raw=z.read('S1B_Assayed_Barcode_Map.csv')
 (m.IN/'S1B_Assayed_Barcode_Map.csv').write_bytes(raw)
 with zipfile.ZipFile(ROOT/'data/source_archives/HAR_Stage21D_FullSixRaw_Corrected_Results.zip') as z:
  df=pd.read_csv(io.BytesIO(z.read('S2_All_Cut_FullRaw_Sensitivity.csv')))
 df[(df.min_DNA_count_per_barcode==1)&(df.aggregate_pseudocount==1)][['HAR','cut','culture','corrected_HminusC_log2']].to_csv(m.IN/'S2_HS111_Baseline_342.csv',index=False)
 mapping=list(csv.DictReader(raw.decode('utf8').splitlines()))
 for rep,rg,dg in m.SAMPLES:
  for gsm,assay in [(rg,'RNA'),(dg,'DNA')]:
   pp=m.CACHE/f'{gsm}_{assay}.tsv.gz'
   with gzip.open(pp,'wt') as f:
    for row in mapping:
     # Synthetic allele signal varies by HAR and across 3 replicas.
     v=(6+(sum(map(ord,row['HAR']))%5)) if assay=='DNA' else (1+(sum(map(ord,row['HAR']))%9)+rep+(row['allele']=='chimp')*2)
     f.write(f"{row['Barcode']}\t{v}\t{row['insert']}\n")
 # Cache files are present, so no network expected.
 out=m.main(); status=json.loads((m.OUT/'S0_WTc_vs_HS111_Source_Status.json').read_text())
 assert status['status']=='complete',status
 assert status['n_shared_cuts_with_HS111']>10
 assert out.exists() and out.stat().st_size>1000
 print('PASS: offline WTc full processing path using SYNTHETIC counts only:',status['n_shared_cuts_with_HS111'],'shared cuts')
