#!/usr/bin/env python3
"""Pinned non-synthetic Stage25 derived-result contract.
Run after python src/stage25_crossline_reanalysis.py.
"""
from pathlib import Path
import json, zipfile, io, hashlib
import pandas as pd
root=Path(__file__).resolve().parents[1]
audit=json.loads((root/'results/stage25/S0_Independent_Numerical_Audit.json').read_text())
assert audit['status']=='PASS_derived_tables_independently_reproduced'
assert (audit['n_cut_culture_WTc'],audit['n_unique_cuts'],audit['n_HARs'])==(342,114,66)
assert (audit['n_shared_sign'],audit['n_all_six_cultures_same_sign'])==(86,48)
assert audit['priority_shared_cuts']==10
assert (audit['priority_crossline_same_sign'],audit['priority_all_six_cultures_stable'])==(7,3)
assert .54<audit['Pearson_crossline']<.55
assert .62<audit['Spearman_crossline']<.63
assert .67<audit['HAR_bootstrap_Pearson_95CI'][1]<.68
assert audit['HAR_cluster_block_signflip_B']==200000
assert audit['original_Noble_180HAR_source_interval_crosswalk_complete'] is False
assert audit['raw_results_zip_SHA256']==hashlib.sha256((root/'data/source_archives/HAR_Stage24_WTc_FinalValidation_Results.zip').read_bytes()).hexdigest()
with zipfile.ZipFile(root/'data/source_archives/HAR_Stage24_WTc_FinalValidation_Results.zip') as z:
  df=pd.read_csv(io.BytesIO(z.read('S5_IndependentLine_ConstructConcordance.csv')))
assert len(df)==114 and df.HAR.nunique()==66
assert len(pd.read_csv(root/'results/stage25/S2_Five_Priority_HAR_CrossLine_Checked.csv'))==10
print('PASS: real WTc-derived Stage25 contract (separate from synthetic parsing tests)')
