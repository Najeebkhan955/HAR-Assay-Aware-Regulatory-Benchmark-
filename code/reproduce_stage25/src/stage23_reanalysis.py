#!/usr/bin/env python3
"""HAR cross-assay compatibility benchmark; source-only descriptive/exploratory analysis.
No new primary experimental observations. All random procedures deterministic and locus-level.
Run from repo root: python src/stage23_reanalysis.py
"""
from __future__ import annotations
import csv, io, json, math, hashlib, zipfile
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr
from collections import defaultdict
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'data/source_archives';OUT=ROOT/'results';OUT.mkdir(parents=True,exist_ok=True)
SEED=20261120

def zipcsv(zipname,entry):
    with zipfile.ZipFile(SRC/zipname) as z:
        return pd.read_csv(io.BytesIO(z.read(entry)))

def require(cond,msg):
    if not cond: raise AssertionError(msg)

def bh(p):
    p=np.asarray(p,float); n=len(p); ix=np.argsort(p); adj=np.empty(n,float)
    adj[ix]=np.minimum.accumulate((p[ix]*n/np.arange(1,n+1))[::-1])[::-1].clip(0,1)
    return adj

m=zipcsv('HAR_Stage20B_ContextMatched_Results.zip','S1_152_HAR_ContextMatched_Statistics.csv')
mapq=zipcsv('HAR_Stage16B_Ortholog_PWM_Results.zip','results/S1_All180_Sequence_Mapping_Audit.csv')
guide=zipcsv('HAR_Stage16B_Ortholog_PWM_Results.zip','inputs/S1_All180_HAR_Sequence_Mapping_Readiness.csv')
r=zipcsv('HAR_Stage21C_MPRA_HS111_N2_Results.zip','S3_Paired_Fragment_Activity_By_Replicate.csv')
design=zipcsv('HAR_Stage21C_MPRA_HS111_N2_Results.zip','S1_75HAR_Design_Allele_Cuts.csv')
raw=zipcsv('HAR_Stage21D_FullSixRaw_Corrected_Results.zip','S2_All_Cut_FullRaw_Sensitivity.csv')
require(len(m)==152 and m.screen_HAR.is_unique,'Stage20B source family invalid')
require(len(mapq)==180 and mapq.screen_HAR.is_unique,'Sequence mapping source invalid')
require(len(guide)==180 and guide.screen_HAR.is_unique,'Guide coverage source invalid')
require(len(r)==342 and len(r.groupby(['HAR','cut']))==114,'MPRA 114x3 source invalid')
require(all(s['replicate'].nunique()==3 for _,s in r.groupby(['HAR','cut'])),'MPRA incomplete culture triplets')
require(len(raw)==4059,'Full raw-count sensitivity source invalid')
r=r.sort_values(['HAR','cut','replicate']).reset_index(drop=True)

# Reconcile source counts and conditional null q-values.
require((m.context_BH_q_152>=0.05).all(),'Unexpected significant motif result')
require(np.max(np.abs(bh(m.context_matched_p)-m.context_BH_q_152))<1e-12,'Stage20B BH mismatch')
require(len(set(design.HAR))==75,'Design HAR label count invalid')

agg=(r.groupby(['HAR','cut'],as_index=False)
      .agg(mean_HminusC=('human_minus_chimp_log2_activity','mean'),sd_culture=('human_minus_chimp_log2_activity','std'),
           culture_1=('human_minus_chimp_log2_activity',lambda x: x.iloc[0]),
           culture_2=('human_minus_chimp_log2_activity',lambda x: x.iloc[1]),
           culture_3=('human_minus_chimp_log2_activity',lambda x: x.iloc[2])))
agg['consistent_direction_3cultures']=((agg[['culture_1','culture_2','culture_3']]>0).all(axis=1) |
                                    (agg[['culture_1','culture_2','culture_3']]<0).all(axis=1))
agg['abs_mean_HminusC']=agg.mean_HminusC.abs()
agg['stable_positive_cut']=agg.consistent_direction_3cultures & agg.mean_HminusC.gt(0)
agg['stable_negative_cut']=agg.consistent_direction_3cultures & agg.mean_HminusC.lt(0)
# Avoid treating reporter cuts as independent HAR loci.
loc=(agg.groupby('HAR',as_index=False).agg(n_cuts=('cut','nunique'),
       median_abs_cut_effect=('abs_mean_HminusC','median'),
       mean_abs_cut_effect=('abs_mean_HminusC','mean'),
       max_abs_cut_effect=('abs_mean_HminusC','max'),
       n_consistent_cuts=('consistent_direction_3cultures','sum'),
       n_stable_positive_cuts=('stable_positive_cut','sum'),
       n_stable_negative_cuts=('stable_negative_cut','sum'),
       min_cut_mean=('mean_HminusC','min'),max_cut_mean=('mean_HminusC','max')))
loc['opposing_cut_directions']=((loc.n_cuts>=2)&(loc.min_cut_mean<0)&(loc.max_cut_mean>0))
loc['sign_consistency_cut_fraction']=loc.n_consistent_cuts/loc.n_cuts
loc['opposing_cuts_each_consistent_3cultures']= (loc.n_cuts>=2)&loc.n_stable_positive_cuts.gt(0)&loc.n_stable_negative_cuts.gt(0)

# Coverage and unavailable evidence are explicitly NA (never imputed as zero).
all180=mapq[['screen_HAR','status','catalog_id_provisional','global_alignment_identity']].rename(columns={'screen_HAR':'HAR','status':'mapping_status'})
all180=all180.merge(guide[['screen_HAR','min_guide_batch_cells','screen_n_guides','passes_min20_guide_batch','passes_min50_guide_batch','all_3_context_aliases_available']],left_on='HAR',right_on='screen_HAR',how='left',validate='one_to_one').drop(columns='screen_HAR')
all180=all180.merge(m[['screen_HAR','observed_max_abs_pwm_delta','context_matched_p','context_BH_q_152','n_aligned_substitutions']],left_on='HAR',right_on='screen_HAR',how='left',validate='one_to_one').drop(columns='screen_HAR')
all180=all180.merge(loc,on='HAR',how='left',validate='one_to_one')
all180['in_MPRA_design_name']=all180.HAR.isin(design.HAR.unique())
all180['MPRA_paired_3cultures']=all180.n_cuts.notna()
all180['mapped_for_PWM']=all180.context_matched_p.notna()
all180['nonzero_motif_score']=all180.observed_max_abs_pwm_delta.gt(0) & all180.mapped_for_PWM
all180['assay_compatibility_gate']=np.select(
    [~all180.mapped_for_PWM & ~all180.MPRA_paired_3cultures,
     ~all180.mapped_for_PWM & all180.MPRA_paired_3cultures,
      all180.mapped_for_PWM & ~all180.MPRA_paired_3cultures,
      all180.mapped_for_PWM & all180.MPRA_paired_3cultures],
    ['neither_provisional_motif_nor_paired_reporter','reporter_only_missing_motif_mapping',
     'motif_mapped_no_paired_reporter','both_assays_available'],default='UNREACHABLE')
all180['mechanistic_proof']='NOT_ESTABLISHED'
all180=all180.sort_values('HAR')
all180.to_csv(OUT/'S1_180HAR_EvidenceCompatibility.csv',index=False)

# Reporter culture correlation: event-level but HAR-clustered bootstrap.
r=r.sort_values(['HAR','cut','replicate']).reset_index(drop=True)
w=r.pivot(index=['HAR','cut'],columns='replicate',values='human_minus_chimp_log2_activity').reset_index()
w.columns=['HAR','cut','R1','R2','R3']
rs=[]; rg=np.random.default_rng(SEED)
byhar=[frame[['R1','R2','R3']].to_numpy(float) for _,frame in w.groupby('HAR',sort=True)]
for a,b in [('R1','R2'),('R1','R3'),('R2','R3')]:
    pr=pearsonr(w[a],w[b]);sr=spearmanr(w[a],w[b])
    pts=[]
    for i in range(2500):
        ids=rg.integers(0,len(byhar),size=len(byhar))
        subset=np.concatenate([byhar[j] for j in ids],axis=0)
        idx1={'R1':0,'R2':1,'R3':2}[a];idx2={'R1':0,'R2':1,'R3':2}[b]
        if np.std(subset[:,idx1]) and np.std(subset[:,idx2]):pts.append(np.corrcoef(subset[:,idx1],subset[:,idx2])[0,1])
    rs.append(dict(culture_a=a,culture_b=b,n_cut_fragments=len(w),n_HAR_clusters=len(byhar),pearson_r=pr.statistic,
                   spearman_r=sr.statistic,cluster_bootstrap_Pearson_95CI_low=float(np.percentile(pts,2.5)),
                   cluster_bootstrap_Pearson_95CI_high=float(np.percentile(pts,97.5)),
                   note='Descriptive same-line cultured replicates; shared library and assay biases possible'))
pd.DataFrame(rs).to_csv(OUT/'S2_MPRA_Replicate_Reliability.csv',index=False)

# Evidence transfer: motif nonzero vs matched motif-zero; only shared 58 HARs.
common=all180.loc[all180.mapped_for_PWM & all180.MPRA_paired_3cultures].copy()
require(len(common)==58,'Expected shared 58-HAR assay cohort')
require(common.nonzero_motif_score.sum()==5,'Expected 5 signal-bearing reporter HARs')
common['n_cut_stratum']=common.n_cuts.astype(int).clip(upper=3)
obs=float(common.loc[common.nonzero_motif_score,'median_abs_cut_effect'].mean()-common.loc[~common.nonzero_motif_score,'median_abs_cut_effect'].mean())
# Stratified label permutation preserves the number of signal loci in each 1/2/3+ cut stratum.
strata=[np.where(common.n_cut_stratum.to_numpy()==s)[0] for s in (1,2,3)]
orig=common.nonzero_motif_score.to_numpy(bool)
nselect=[int(orig[ix].sum()) for ix in strata]
y=common.median_abs_cut_effect.to_numpy(float)
perms=np.empty(30000)
for k in range(len(perms)):
    mask=np.zeros(len(common),dtype=bool)
    for ix,ns in zip(strata,nselect):
        if ns:mask[rg.choice(ix,size=ns,replace=False)]=True
    perms[k]=y[mask].mean()-y[~mask].mean()
p2=(1+(np.abs(perms)>=abs(obs)).sum())/(1+len(perms))
pupper=(1+(perms>=obs).sum())/(1+len(perms))
# Simple group bootstrap, not a confidence interval valid for causal inference.
selected=y[orig];others=y[~orig]
boots=np.empty(10000)
for k in range(len(boots)):
    boots[k]=rg.choice(selected,len(selected),replace=True).mean()-rg.choice(others,len(others),replace=True).mean()
sp=spearmanr(common.observed_max_abs_pwm_delta,common.median_abs_cut_effect)
transfer=dict(shared_HARs=len(common),nonzero_PWM_HARs=int(orig.sum()),PWM_zero_HARs=int((~orig).sum()),
  primary_outcome='mean across HARs of median abs(three-culture mean human-minus-chimp per cut)',
  motif_nonzero_mean=float(selected.mean()),motif_zero_mean=float(others.mean()),observed_difference=obs,
  stratified_label_permutation_two_sided_p=float(p2),stratified_label_permutation_upper_p=float(pupper),
  perm_n=len(perms),permutation_label_strata='number of cuts 1,2,3+',
  bootstrap_difference_95CI=[float(x) for x in np.percentile(boots,[2.5,97.5])],
  spearman_continuous_PWM_delta_vs_reporter_abs_r=float(sp.statistic),spearman_exploratory_p=float(sp.pvalue),
  inference='Post hoc cross-assay association; no gene target, TF binding, or predictive superiority is established')
(OUT/'S3_MotifReporter_Transfer_Benchmark.json').write_text(json.dumps(transfer,indent=2))
common[['HAR','mapping_status','n_cuts','n_aligned_substitutions','observed_max_abs_pwm_delta','nonzero_motif_score','context_matched_p','context_BH_q_152','median_abs_cut_effect','max_abs_cut_effect','sign_consistency_cut_fraction','opposing_cut_directions','min_guide_batch_cells','passes_min20_guide_batch','all_3_context_aliases_available']].to_csv(OUT/'S4_58HAR_SharedAssays.csv',index=False)

# Fragment-level cut context vs depth threshold sensitivity.
raw['threshold_key']=raw.min_DNA_count_per_barcode.astype(str)+'/'+raw.aggregate_pseudocount.astype(str)
baseline=raw[(raw.min_DNA_count_per_barcode==1)&(raw.aggregate_pseudocount==1)].copy()
require(len(baseline)==342,'Raw sensitivity did not include all baseline triplets')
base=baseline.set_index(['HAR','cut','culture'])['corrected_HminusC_log2']
comp=raw[(raw.min_DNA_count_per_barcode.isin([5,10,20]))&(raw.aggregate_pseudocount==1)].copy()
comp['base']=comp.set_index(['HAR','cut','culture']).index.map(base)
comp['sign_changed']=(comp.corrected_HminusC_log2*comp.base<0)
threshold_summary=(comp.groupby('min_DNA_count_per_barcode',as_index=False)
                   .agg(n_comparable=('HAR','size'),n_sign_changed=('sign_changed','sum')))
threshold_summary['proportion_changed']=threshold_summary.n_sign_changed/threshold_summary.n_comparable
threshold_summary.to_csv(OUT/'S5_DNA_DepthSensitivity.csv',index=False)

# Multi-fragment heterogeneity; do NOT equate signed null with no biology.
hetero=loc.loc[loc.n_cuts>=2,['HAR','n_cuts','min_cut_mean','max_cut_mean','opposing_cut_directions','opposing_cuts_each_consistent_3cultures','n_consistent_cuts','n_stable_positive_cuts','n_stable_negative_cuts']].copy()
hetero['cross_cut_range_log2']=hetero.max_cut_mean-hetero.min_cut_mean
hetero.sort_values('cross_cut_range_log2',ascending=False).to_csv(OUT/'S6_MultiCut_Heterogeneity.csv',index=False)
agg.to_csv(OUT/'S7_114_Cut_3Culture_Summary.csv',index=False)

counts=dict(n_original_HAR=180,n_provisionally_mapped=int(all180.mapped_for_PWM.sum()),
  n_original_unmapped=int((~all180.mapped_for_PWM).sum()),n_GEO_design_name_HAR=int(all180.in_MPRA_design_name.sum()),
  n_MPRA_paired_HAR=int(all180.MPRA_paired_3cultures.sum()),n_both_motif_and_MPRA=int((all180.mapped_for_PWM&all180.MPRA_paired_3cultures).sum()),
  n_shared_nonzero_motif_HAR=int(common.nonzero_motif_score.sum()),
  n_all152_motif_context_BH_0p05=int((m.context_BH_q_152<0.05).sum()),
  n_MPRA_fragments=len(agg),n_consistent_3cultures=int(agg.consistent_direction_3cultures.sum()),
  n_multicut_loci=len(hetero),n_multicut_direction_opposing=int(hetero.opposing_cut_directions.sum()),
  n_multicut_stable_3culture_opposing=int(hetero.opposing_cuts_each_consistent_3cultures.sum()),
  total_fullraw_sensitivity_rows=len(raw),
  fullraw_5_10_20_pseudo1_comp_rows=len(comp),fullraw_5_10_20_pseudo1_sign_changes=int(comp.sign_changed.sum()))
summary=dict(stage='23',version='2026-11-20',status='locally_rerun_complete',counts=counts,
              motif_reporter_transfer=transfer,replicate_correlations=rs,
              stage20B_BH_recomputed=True,source_archive_checksums={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(SRC.glob('*.zip'))},
              scope_limits=['one HS1-11 N2 cell line assayed; three cultures not donor replicates',
               '152 primary-author interval identities remain unverified',
               '5 motif-nonzero loci in shared assay cohort; post-hoc and low power',
               'no independent TF occupancy or enhancer-gene causality demonstrated'])
(OUT/'S0_Stage23_Summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({'counts':counts,'transfer_difference':obs,'stratified_p_two_sided':p2,
                  'pearson':[(x['culture_a'],x['culture_b'],round(x['pearson_r'],4)) for x in rs]},indent=2))
