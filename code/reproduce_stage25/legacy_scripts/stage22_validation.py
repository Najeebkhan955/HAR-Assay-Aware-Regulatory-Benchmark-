import csv,zipfile,io,math,statistics,random,json,hashlib
from collections import defaultdict,Counter
from pathlib import Path
from scipy.stats import t as student_t
BASE=Path('/mnt/data');Z21=zipfile.ZipFile(BASE/'HAR_Stage21D_FullSixRaw_Corrected_Results.zip'); Z20=zipfile.ZipFile(BASE/'HAR_Stage20B_ContextMatched_Results.zip')
def read_zip(z,name):return list(csv.DictReader(io.StringIO(z.read(name).decode('utf-8'))))
def read_csv(path):return list(csv.DictReader(open(path,encoding='utf-8')))
raw=read_zip(Z21,'S2_All_Cut_FullRaw_Sensitivity.csv')
reps=read_csv(BASE/'S3_Paired_Fragment_Activity_By_Replicate.csv')
stage20=read_zip(Z20,'S1_152_HAR_ContextMatched_Statistics.csv')
source=read_csv(BASE/'Stage21D_Work/S5_20_ConstructAllele_Sequence_Alignment.csv')
summary21=read_csv(BASE/'Stage21D_Work/S1_114_CutReplicate_Exploratory_Intervals.csv')
print('source S1 columns',list(summary21[0]));print('rows',len(summary21))
print('S2 row counts by threshold and pseudo',Counter((int(x['min_DNA_count_per_barcode']),float(x['aggregate_pseudocount'])) for x in raw))
by_cut=defaultdict(list)
for r in reps:by_cut[(r['HAR'],r['cut'])].append(float(r['human_minus_chimp_log2_activity']))
assert len(by_cut)==114
full=defaultdict(list)
for r in raw:full[(r['HAR'],r['cut'],int(r['culture']),int(r['min_DNA_count_per_barcode']),float(r['aggregate_pseudocount']))].append(float(r['corrected_HminusC_log2']))
assert len(full)==len(raw)
base={k:v[0] for k,v in full.items() if k[-2:]==(1,1.0)}
print('baseline count',len(base),'expected',114*3)
for key, vals in by_cut.items():
 assert len(vals)==3
# Source S3 is ordered by HAR/cut/replicate; assert explicitly.
rep_order=defaultdict(list)
for r in reps: rep_order[(r['HAR'],r['cut'])].append(int(r['replicate']))
assert all(x == [1,2,3] for x in rep_order.values())
print('baseline vs previous abs difference',max(abs(base[(h,cut,i,1,1.0)]-v) for (h,cut),vs in by_cut.items() for i,v in enumerate(vs,1)))
# Group source keys and metadata
m20={x['screen_HAR']:x for x in stage20}
src={(r['HAR'],r['cut'],r['allele']):r for r in source}
records=[]
for (h,cut),xs in sorted(by_cut.items()):
 mu=statistics.mean(xs);sd=statistics.stdev(xs);p=2*student_t.sf(abs(mu)/(sd/math.sqrt(3)),df=2) if sd else (0.0 if mu else 1.0)
 minsig=0.25 if all(x>0 for x in xs) or all(x<0 for x in xs) else 1.0 # two-sided exact sign test n3: if 2+1 sign result =1
 threshold={}
 for t in [1,5,10,20]:
  vals=[base[(h,cut,i,1,1.0)] if t==1 else full.get((h,cut,i,t,1.0),[None])[0] for i in (1,2,3)]
  threshold[t]=vals
  # check each present
  if t==1: assert all(isinstance(x,float) for x in vals)
 sign_consistent=all(x>0 for x in xs) or all(x<0 for x in xs)
 changes=[int(v*base[(h,cut,i,1,1.0)]<0) for t in [5,10,20] for i,v in enumerate(threshold[t],1) if v is not None and base.get((h,cut,i,1,1.0)) not in (0,None)]
 missing=sum(v is None for t in [5,10,20] for v in threshold[t])
 r=m20.get(h,{})
 records.append(dict(HAR=h,cut=cut,n_cultures=3,mean_HminusC_log2=mu,sd_across_cultures=sd,all_3_same_direction=sign_consistent,nominal_student_t_p=p,exact_sign_2s_p=minsig,Stage20B_3mer_p=float(r.get('context_matched_p','nan')),Stage20B_BH_q=float(r.get('context_BH_q_152','nan')),sign_flips_under_DNA5_10_20=sum(changes),missing_threshold_cells=missing,n_tested_threshold_cells=9-missing))
# BH adjust 114 parametric exploratory
p_sorted=sorted(range(len(records)), key=lambda j:records[j]['nominal_student_t_p']); running=1.
for ri in reversed(range(len(records))):
 j=p_sorted[ri]; running=min(running,records[j]['nominal_student_t_p']*len(records)/(ri+1));records[j]['exploratory_t_BH_q_114']=running
print('t q<0.05',[(r['HAR'],r['cut'],r['exploratory_t_BH_q_114']) for r in records if r['exploratory_t_BH_q_114']<0.05]);print('sign consistent',sum(r['all_3_same_direction'] for r in records))
print('all-sign sign-test p',Counter(r['exact_sign_2s_p'] for r in records))
print('threshold flips',sum(r['sign_flips_under_DNA5_10_20'] for r in records),'missing',sum(r['missing_threshold_cells'] for r in records))
print('five')
priority=['2xHAR.70','2xHAR.310','2xHAR.234','2xHAR.377','2xHAR.3']
for r in records:
 if r['HAR'] in priority:print({k:round(v,5) if isinstance(v,float) else v for k,v in r.items() if k not in ('n_cultures','sd_across_cultures')})
# motif overlap Stage18 actual aligned positions in humans, index 1-based, restrict full match of payload; motif lengths JASPAR SOX11 8; SOX2 7? report start/end with uncertainty
motif_sites=[('2xHAR.70','SOX11',303),('2xHAR.310','SOX11',581),('2xHAR.377','SOX2',126)]
overlap=[]
for h,tf,pos in motif_sites:
 for r in source:
  if r['HAR']!=h or r['allele']!='human':continue
  start=int(r['catalog_alignment_start_0based']);end=int(r['catalog_alignment_end_0based'])
  # pos is aligned 1-based; treat catalog approximate start; gaps may affect 1 bp
  overlap.append(dict(HAR=h,TF=tf,motif_window_start_aligned_1based=pos,MPRA_cut=r['cut'],catalog_payload_start_0based=start,catalog_payload_end_0based_exclusive=end,motif_start_inside_reporter_payload=(start<=pos-1<end),payload_exact_fraction=float(r['fraction_identical_of_payload']),notes=('Alignment coordinate provisional; exact motif-window overlap needs allele-specific alignment and motif width'))) 
print('motif coverage:')
for r in overlap:print(r)
# locus-level negative control: 5 selected vs 61 other HAR with N2 reporter, median ABS effect of cuts
per_har=defaultdict(list)
for r in records:per_har[r['HAR']].append(abs(r['mean_HminusC_log2']))
agg={h:statistics.median(x) for h,x in per_har.items()}; selected=[agg[h] for h in priority];others=[v for h,v in agg.items() if h not in priority]
diff=statistics.mean(selected)-statistics.mean(others)
seed=2212026;rg=random.Random(seed);fullpop=list(agg.values());n=20000
ex=0
for _ in range(n):
 samp_idx=set(rg.sample(range(len(fullpop)),5));s=[fullpop[k] for k in samp_idx];o=[v for k,v in enumerate(fullpop) if k not in samp_idx]; st=statistics.mean(s)-statistics.mean(o)
 if st>=diff:ex+=1
pperm=(ex+1)/(n+1)
print('MPRA 5 vs others',len(agg),statistics.mean(selected),statistics.mean(others),diff,pperm)
# Write output
out=BASE/'Stage22_Scientific_Gate';out.mkdir(exist_ok=True)
def wr(n,arr):
 with open(out/n,'w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(arr[0]));w.writeheader();w.writerows(arr)
wr('S1_All114_Reporter_StatisticalGate.csv',records)
wr('S2_Motif_Site_Cut_Coverage.csv',overlap)
wr('S3_Five_Candidate_Reporter_Gate.csv',[r for r in records if r['HAR'] in priority]);
(out/'S0_Summary.json').write_text(json.dumps(dict(n_HAR=180,n_mpra_HAR=len(agg),tested_cuts=114,prioritized_HAR=5,Stage20B_BH_significant=0,exploratory_t_BH_significant=sum(r['exploratory_t_BH_q_114']<0.05 for r in records),min_exact_sign_test_p=0.25,MPRA_priority_vs_other_mean_median_abs_effect=diff,MPRA_locus_label_permutation_exploratory_p=pperm,permutation_reps=n,permutation_seed=seed,statistical_caveat='Post hoc non-independent, unequal cut counts; exploratory not confirmatory'),indent=2,default=float))
