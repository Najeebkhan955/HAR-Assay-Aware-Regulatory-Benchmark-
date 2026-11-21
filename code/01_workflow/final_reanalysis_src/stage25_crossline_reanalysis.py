#!/usr/bin/env python3
"""Independent numerical reconstruction of WTc/N2 derived-data outputs.

Requires only pinned Stage24 WTc results ZIP and pinned Stage24 HS1-11 baseline CSV.
Does NOT independently fetch or verify the six raw GEO gz file bytes.
"""
from pathlib import Path
import io, json, zipfile, hashlib
import numpy as np
import pandas as pd
from scipy import stats

BASE=Path(__file__).resolve().parents[1]
WT=BASE/'data/source_archives/HAR_Stage24_WTc_FinalValidation_Results.zip'
HS=BASE/'data/derived_inputs/S2_HS111_Baseline_342.csv'
MAP=BASE/'data/derived_inputs/S1B_Assayed_Barcode_Map.csv'
OUT=BASE/'results/stage25';OUT.mkdir(parents=True,exist_ok=True)

with zipfile.ZipFile(WT) as z:
    assert not z.testzip(), 'ZIP CRC failure'
    tbl=lambda n:pd.read_csv(io.BytesIO(z.read(n)))
    s0=json.loads(z.read('S0_IndependentLine_Status.json'))
    a=tbl('S1_WTc_CutByCulture.csv');sources=tbl('S2_SixGEO_SourceIntegrity.csv')
    alleles=tbl('S3_WTc_Allele_ActivityAndDepth.csv');sens=tbl('S4_WTc_FullDepth_PseudocountSensitivity.csv')
    joined=tbl('S5_IndependentLine_ConstructConcordance.csv');prior=tbl('S6_Five_Priority_HARs.csv');depth=tbl('S7_WTc_DNAdepth_SensitivitySummary.csv')
assert hashlib.sha256(MAP.read_bytes()).hexdigest()=='6f30451ac623ed4948c8649aa3e9113f0d244b5d70db1f37fbb35314cf92cdeb'
assert len(a)==342 and len(joined)==114 and a.HAR.nunique()==66
assert a.groupby('culture').size().to_dict()=={1:114,2:114,3:114}
assert not a.duplicated(['HAR','cut','culture']).any()
assert len(alleles)==765 and len(sources)==6 and len(sens)==4089 and len(prior)==10
assert (sources.selected_insert_annotation_mismatches==0).all()
assert sources['sha256'].map(lambda x:len(x)==64 and set(x.lower())<=set('0123456789abcdef')).all()
assert sources['GSM'].nunique()==6 and set(sources.assay)=={'RNA','DNA'}
assert sources.groupby(['culture','assay']).size().eq(1).all()
assert all(s0[k] for k in ['raw_GEO_source_insert_annotations_checked','six_GEO_downloads_sha256_recorded'])
# Independent reconstruction of WTc H-C effects from allele totals
pa=alleles.pivot(index=['HAR','cut','culture'],columns='allele',values='log2_activity')
raw_effect=(pa['human']-pa['chimp']).dropna().sort_index()
reported=a.set_index(['HAR','cut','culture']).WTc_HminusC_log2.sort_index()
assert raw_effect.index.equals(reported.index)
assert np.max(np.abs(raw_effect-reported))<1e-12
# Reconstruct mean/median/reproducibility from independently pinned baseline
hs=pd.read_csv(HS)
assert len(hs)==342 and not hs.duplicated(['HAR','cut','culture']).any()
gh=hs.groupby(['HAR','cut']).corrected_HminusC_log2.agg(['count','mean','min','max'])
gw=a.groupby(['HAR','cut']).WTc_HminusC_log2.agg(['count','mean','min','max'])
assert gh.index.equals(gw.index)
agg=gh.join(gw,lsuffix='_HS111',rsuffix='_WTc').reset_index()
agg['HS111_all3_same_direction']=(agg.count_HS111.eq(3)&((agg.min_HS111>0)|(agg.max_HS111<0)))
agg['WTc_all3_same_direction']=(agg.count_WTc.eq(3)&((agg.min_WTc>0)|(agg.max_WTc<0)))
agg['same_direction_2lines']=agg.mean_HS111*agg.mean_WTc>0
agg['both_lines_same_direction_all3cultures']=agg.HS111_all3_same_direction&agg.WTc_all3_same_direction&agg.same_direction_2lines
comp=agg.merge(joined,on=['HAR','cut'],validate='one_to_one')
for j,k in [('mean_HS111','HS111_mean'),('mean_WTc','WTc_mean'),('min_HS111','HS111_min'),('min_WTc','WTc_min'),('max_HS111','HS111_max'),('max_WTc','WTc_max')]:
    assert np.allclose(comp[j],comp[k],atol=1e-12)
for k in ['same_direction_2lines','both_lines_same_direction_all3cultures','WTc_all3_same_direction','HS111_all3_same_direction']:
    assert (comp[k+'_x'].astype(int)==comp[k+'_y'].astype(int)).all() if k in ['same_direction_2lines','both_lines_same_direction_all3cultures','WTc_all3_same_direction','HS111_all3_same_direction'] and (k+'_x') in comp else True
assert comp['same_direction_2lines_x'].sum()==86
assert comp['both_lines_same_direction_all3cultures_x'].sum()==48
assert np.allclose(joined['WTc_mean'],a.groupby(['HAR','cut']).WTc_HminusC_log2.mean().reindex(pd.MultiIndex.from_frame(joined[['HAR','cut']])).to_numpy())
# Recalculate 5000-cluster CI using same independent implementation and recorded seed
rng=np.random.default_rng(20261120)
byhar=[t[['HS111_mean','WTc_mean']].to_numpy(float) for _,t in joined.groupby('HAR',sort=True)]
boot=np.zeros((5000,2))
for i in range(len(boot)):
    b=np.concatenate([byhar[k] for k in rng.integers(0,len(byhar),len(byhar))])
    boot[i,0]=np.corrcoef(b[:,0],b[:,1])[0,1]
    boot[i,1]=np.mean(b[:,0]*b[:,1]>0)
meta=s0['independent_line_cluster_bootstrap']
assert np.allclose([np.corrcoef(joined.HS111_mean,joined.WTc_mean)[0,1],(joined.HS111_mean*joined.WTc_mean>0).mean()], [meta['Pearson_cut_mean'],meta['mean_sign_agreement']],atol=1e-12)
assert np.allclose(np.percentile(boot[:,0],[2.5,97.5]),meta['Pearson_95pct_bootstrap_HAR_CI'],atol=1e-12)
assert np.allclose(np.percentile(boot[:,1],[2.5,97.5]),meta['mean_sign_agreement_95pct_HAR_bootstrap_CI'],atol=1e-12)
# Verify sensitivity baseline per-culture cuts and depth flip rates
sb=sens[(sens.min_DNA_count_per_barcode==1)&(sens.aggregate_pseudocount==1)].set_index(['HAR','cut','culture']).WTc_HminusC_log2
assert sb.index.is_unique and np.allclose(sb.reindex(reported.index),reported,atol=1e-12)
rows=[]
for cutoff in [5,10,20]:
    cur=sens[(sens.min_DNA_count_per_barcode==cutoff)&(sens.aggregate_pseudocount==1)].set_index(['HAR','cut','culture']).WTc_HminusC_log2
    flip=int((cur*sb.reindex(cur.index)<0).sum())
    assert len(cur)==int(depth.loc[depth.min_DNA_count_per_barcode==cutoff,'n_comparable'].iloc[0])
    assert flip==int(depth.loc[depth.min_DNA_count_per_barcode==cutoff,'n_sign_flips'].iloc[0])
    rows.append({'DNA_threshold':cutoff,'n_comparable':len(cur),'WTc_sign_flips':flip,'fraction_flipped':flip/len(cur)})
# External-line sign agreement conditional on HAR-cluster inversion; exploratory permutation
by_har=(joined.assign(sign_product=np.sign(joined.HS111_mean*joined.WTc_mean))
  .groupby('HAR',sort=True).sign_product.sum().to_numpy())
observed=float(by_har.sum()) # 86 - 28 = 58
rng=np.random.default_rng(20261121)
B=200000
extreme=0; extreme_abs=0
for start in range(0,B,20000):
    n=min(20000,B-start)
    sg=rng.integers(0,2,size=(n,len(by_har)),dtype=np.int8)*2-1
    x=(sg*by_har).sum(axis=1)
    extreme+=int((x>=observed).sum())
    extreme_abs+=int((np.abs(x)>=abs(observed)).sum())
p_one=(1+extreme)/(B+1);p_two=(1+extreme_abs)/(B+1)
# HAR-balanced agreement, cluster bootstrap, and two-line candidate descriptors
perhar=joined.assign(same=joined.same_direction_2lines.astype(float)) .groupby('HAR').agg(n_cuts=('cut','size'),n_same=('same','sum'),fraction_same=('same','mean')).reset_index()
priority=joined.loc[joined.HAR.isin(['2xHAR.70','2xHAR.310','2xHAR.234','2xHAR.377','2xHAR.3'])].copy()
priority['abs_delta_between_lines']=abs(priority.HS111_mean-priority.WTc_mean)
priority=priority[['HAR','cut','HS111_mean','WTc_mean','same_direction_2lines','HS111_all3_same_direction','WTc_all3_same_direction','both_lines_same_direction_all3cultures','abs_delta_between_lines']].sort_values(['HAR','cut'])
priority.to_csv(OUT/'S2_Five_Priority_HAR_CrossLine_Checked.csv',index=False)
perhar.to_csv(OUT/'S3_66HAR_HARBalanced_SignAgreement.csv',index=False)
pd.DataFrame(rows).to_csv(OUT/'S4_WTc_DepthSensitivity_Verified.csv',index=False)
# Compute within-line culture correlations; not independent donors
wide_w=a.pivot(index=['HAR','cut'],columns='culture',values='WTc_HminusC_log2')
wide_h=hs.pivot(index=['HAR','cut'],columns='culture',values='corrected_HminusC_log2')
within=[]
for line,frame in [('WTc',wide_w),('HS1-11',wide_h)]:
    for p1,p2 in [(1,2),(1,3),(2,3)]:
        within.append(dict(line=line,culture1=p1,culture2=p2,Pearson_r=float(frame[p1].corr(frame[p2])),Spearman_r=float(frame[p1].corr(frame[p2],method='spearman'))))
pd.DataFrame(within).to_csv(OUT/'S5_WithinLine_Culture_Correlations.csv',index=False)
summary={
'status':'PASS_derived_tables_independently_reproduced',
'limits':'Six WTc raw GEO gzip bytes NOT supplied in uploaded archive; checksums/source annotation are reported by Colab workflow, not independently rehashed here',
'raw_results_zip_SHA256':hashlib.sha256(WT.read_bytes()).hexdigest(),
'original_HS111_baseline_sha256':hashlib.sha256(HS.read_bytes()).hexdigest(),
'barcode_design_map_sha256':hashlib.sha256(MAP.read_bytes()).hexdigest(),
'number_GEO_samples_reported':len(sources),
'n_cut_culture_WTc':len(a),'n_unique_cuts':len(joined),'n_HARs':joined.HAR.nunique(),
'n_shared_sign':int(comp.same_direction_2lines_x.sum()),
'n_all_six_cultures_same_sign':int(comp.both_lines_same_direction_all3cultures_x.sum()),
'Pearson_crossline':float(joined.HS111_mean.corr(joined.WTc_mean)),
'Spearman_crossline':float(joined.HS111_mean.corr(joined.WTc_mean,method='spearman')),
'HAR_bootstrap_Pearson_95CI':meta['Pearson_95pct_bootstrap_HAR_CI'],
'HAR_bootstrap_sign_agreement_95CI':meta['mean_sign_agreement_95pct_HAR_bootstrap_CI'],
'HAR_balanced_sign_agreement_mean':float(perhar.fraction_same.mean()),
'HAR_cluster_block_signflip_B':B,
'HAR_cluster_block_signflip_one_sided_p_exploratory':p_one,
'HAR_cluster_block_signflip_two_sided_p_exploratory':p_two,
'priority_shared_cuts':len(priority),'priority_crossline_same_sign':int(priority.same_direction_2lines.sum()),
'priority_all_six_cultures_stable':int(priority.both_lines_same_direction_all3cultures.sum()),
'reporter_activity_is_not_endogenous_regulation':True,
'original_Noble_180HAR_source_interval_crosswalk_complete':False,
'GEO_data_HS111_WTc_HUMAN_host_alleles_human_vs_chimp':True,
}
(OUT/'S0_Independent_Numerical_Audit.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
print(json.dumps(summary,indent=2))
