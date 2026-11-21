"""Stage24 full independent-line flow smoke test with SYNTHETIC source count files."""
import tempfile,gzip,csv,importlib.util,io,zipfile,json
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('stage24_wtc',ROOT/'src/stage24_wtc_final.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
with tempfile.TemporaryDirectory() as td:
 d=Path(td);m.ROOT=d;m.INPUT=d/'HAR_Stage24_Embedded';m.OUT=d/'HAR_Stage24_WTc_Results';m.CACHE=d/'HAR_Stage24_WTc_GEO_Cache'
 for x in (m.INPUT,m.OUT,m.CACHE):x.mkdir()
 with zipfile.ZipFile(ROOT/'data/source_archives/HAR_Stage21C_MPRA_HS111_N2_Results.zip') as z:
  raw=z.read('S1B_Assayed_Barcode_Map.csv')
 (m.INPUT/'S1B_Assayed_Barcode_Map.csv').write_bytes(raw)
 with zipfile.ZipFile(ROOT/'data/source_archives/HAR_Stage21D_FullSixRaw_Corrected_Results.zip') as z:
  table=pd.read_csv(io.BytesIO(z.read('S2_All_Cut_FullRaw_Sensitivity.csv')))
 b=table[(table.min_DNA_count_per_barcode==1)&(table.aggregate_pseudocount==1)][['HAR','cut','culture','corrected_HminusC_log2']]
 b.to_csv(m.INPUT/'S2_HS111_Baseline_342.csv',index=False)
 rows=list(csv.DictReader(raw.decode('utf8').splitlines()))
 for rep,rg,dg in m.SAMPLES:
  for gsm,assay in [(rg,'RNA'),(dg,'DNA')]:
   with gzip.open(m.CACHE/f'{gsm}_{assay}.tsv.gz','wt') as f:
    for i,row in enumerate(rows):
     num=6+(i%13) if assay=='DNA' else (0 if i%17==0 else 9+(i%7)+(row['allele']=='chimp')*2+rep)
     f.write(f"{row['Barcode']}\t{num}\t{row['insert']}:001\n")
 result=m.main()
 status=json.loads((m.OUT/'S0_IndependentLine_Status.json').read_text())
 assert status['status']=='complete' and status['n_shared_cuts']==114 and status['n_shared_HARs']==66
 assert status['selected_insert_annotation_mismatches']==0
 assert len(pd.read_csv(m.OUT/'S2_SixGEO_SourceIntegrity.csv'))==6
 assert result.exists()
 print('PASS: synthetic Stage24 independent line end-to-end; no real WTc data asserted')
