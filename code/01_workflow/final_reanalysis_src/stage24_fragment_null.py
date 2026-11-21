#!/usr/bin/env python3
"""Descriptive cut-grouping null: are opposite-sign cut patterns surprising vs random grouping?
Not a biological randomization test; exchangeability assumptions are purely methodological.
"""
import pandas as pd,numpy as np,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
R=ROOT/'results'
T=ROOT/'results/stage24'
T.mkdir(exist_ok=True)
a=pd.read_csv(R/'S7_114_Cut_3Culture_Summary.csv')
assert len(a)==114 and a.HAR.nunique()==66
sizes=a.groupby('HAR').size()
multi=sizes[sizes>=2]
assert len(multi)==28
b=a[a.HAR.isin(multi.index)].sort_values(['HAR','cut']).reset_index(drop=True)
starts=np.r_[0,np.cumsum(multi.values)[:-1]]
assert len(b)==sum(multi.values)
allvec=a[['mean_HminusC','stable_positive_cut','stable_negative_cut']].to_numpy()
subvec=b[['mean_HminusC','stable_positive_cut','stable_negative_cut']].to_numpy()

def tally(x):
    mean=x[:,0];sp=x[:,1].astype(bool);sn=x[:,2].astype(bool)
    pos=np.add.reduceat(mean>0,starts)>0
    neg=np.add.reduceat(mean<0,starts)>0
    stpos=np.add.reduceat(sp,starts)>0
    stneg=np.add.reduceat(sn,starts)>0
    return int((pos&neg).sum()),int((stpos&stneg).sum())
obs=tally(subvec)
rng=np.random.default_rng(20261120)
N=50000
outputs=[]
for label,source in [('multi_cut_pool_only',subvec),('all_114_cuts',allvec)]:
 arr=np.empty((N,2),np.int16)
 for k in range(N):
  if label=='multi_cut_pool_only':group=source[rng.permutation(len(source))]
  else: group=source[rng.choice(len(source),len(subvec),replace=False)]
  arr[k]=tally(group)
 for j,metric in enumerate(['opposite_mean_sign_loci','opposite_stable_3culture_sign_loci']):
  vals=arr[:,j]
  lower=(1+np.sum(vals<=obs[j]))/(N+1)
  upper=(1+np.sum(vals>=obs[j]))/(N+1)
  two=min(1,2*min(lower,upper))
  outputs.append(dict(null=label,metric=metric,observed=obs[j],n_multi_loci=len(multi),
                      permutations=N,null_mean=float(vals.mean()),null_sd=float(vals.std()),
                      null_95pct_low=float(np.percentile(vals,2.5)),null_95pct_high=float(np.percentile(vals,97.5)),
                      p_lower_tail=float(lower),p_upper_tail=float(upper),p_two_sided_approx=float(two),
                      interpretation='Descriptive grouping null, not biological significance'))
pd.DataFrame(outputs).to_csv(T/'S1_MultiCut_Exchangeability_Null.csv',index=False)
summary={'stage':'24','scope':'Can reported opposite-sign reporter-cut patterns exceed random regrouping expectations?',
         'data_source':'Archived Stage23 114 reporter cut mean/three-culture sign fields',
         'n_cut_events':len(a),'n_HARs':a.HAR.nunique(),'n_multi_HARs':len(multi),
         'n_multi_cut_events':len(b),'observed_opposite_mean_sign':obs[0],
         'observed_opposite_3culture_stable':obs[1],
         'null_method':'Shuffle cut vectors while preserving per-HAR number of cuts, either within multi-cut pool or from all 114; source original cuts exchangeability is NOT guaranteed',
         'findings':outputs}
(T/'S0_CutGrouping_Null_Status.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
