"""HAR Stage20A exact conditional-null statistic, faster aligned-window scans.

Optimized using Numba to re-score only motif windows affected by randomized substitutions
or fixed indels; all unaffected matched human-chimp windows have score difference zero.
Input ortholog correspondences remain provisional, independent author BED not found.
This is an exploratory reference-base-conditional null, NOT dinucleotide-preserving.
"""
from pathlib import Path
import os, sys, csv, json, time, zipfile, random, hashlib
import numpy as np
from numba import njit
import HAR_Stage20A_ConditionalNull_NoDrive as original
from Bio.Seq import Seq

ROOT = Path('/content') if Path('/content').exists() else Path.cwd()
OUT = ROOT/'HAR_Stage20A_FullScreen_Results'
OUT.mkdir(exist_ok=True,parents=True)
N_PERM = int(os.environ.get('HAR_N_PERM', '4999'))
if N_PERM < 99:
    raise ValueError('N_PERM below 99 is insufficient for this research audit')
SEED = 20261120
CHECKPOINT = int(os.environ.get('HAR_CHECKPOINT', '1'))

@njit(cache=True)
def window_delta_affected(h, changed, n_changed, ref_score, pwm, widths, lows, highs, marked, starts, threshold):
    """Compute identical six-TF max statistic; unaffected windows contribute zero."""
    n = len(h)
    best = 0.0
    for tf in range(len(widths)):
        w = widths[tf]
        if w > n:
            continue
        limit = n - w
        n_starts = 0
        for z in range(n_changed):
            p = changed[z]
            low = max(0, p-w+1)
            high = min(limit, p)
            for st in range(low,high+1):
                if marked[st] == 0:
                    marked[st] = 1
                    starts[n_starts] = st
                    n_starts += 1
        for z in range(n_starts):
            st = starts[z]
            cf = ref_score[tf,0,st]
            cr = ref_score[tf,1,st]
            # Both strands of the reference have identical validity conditions.
            if not np.isfinite(cf):
                continue
            hf = 0.0
            hr = 0.0
            valid = True
            for j in range(w):
                f = h[st+j]
                r = h[st+w-1-j]
                if f < 0 or r < 0:
                    valid = False
                    break
                hf += pwm[tf,j,f]
                hr += pwm[tf,j,3-r]
            if not valid:
                continue
            scale = highs[tf] - lows[tf]
            hf = (hf - lows[tf]) / scale
            hr = (hr - lows[tf]) / scale
            if hf >= threshold or cf >= threshold:
                d = abs(hf-cf)
                if d > best: best = d
            if hr >= threshold or cr >= threshold:
                d = abs(hr-cr)
                if d > best: best = d
        for z in range(n_starts):
            marked[starts[z]] = 0
    return best


def score_data(chimp, models):
    n=len(chimp)
    tfs=list(models)
    widths=np.array([models[t][0].shape[0] for t in tfs],dtype=np.int32)
    maxw=int(widths.max())
    pwms=np.zeros((len(tfs),maxw,4),dtype=np.float64)
    los=np.array([models[t][1] for t in tfs],dtype=np.float64)
    his=np.array([models[t][2] for t in tfs],dtype=np.float64)
    reference=np.full((len(tfs),2,n),-np.inf,dtype=np.float64)
    codes=original.encoded(chimp)
    for t,tf in enumerate(tfs):
        w=widths[t]
        pwms[t,:w,:] = models[tf][0]
        sc=original.scan_scores(codes,models[tf])
        if sc.size:
            reference[t,:,:sc.shape[1]]=sc
    return reference,pwms,widths,los,his


def locus_simulations(human, chimp, models, rng, n_perms, observed):
    """Reproduce reference-compatible substitution-placement sampling with fixed indels."""
    h0=original.encoded(human)
    c=original.encoded(chimp)
    n=len(c)
    base=c.copy()
    gapcols=np.where((h0<0)|(c<0))[0]
    base[gapcols]=h0[gapcols]
    fixed=np.flatnonzero(base!=c).astype(np.int32)
    subs={}
    possible={}
    for x in range(4):
        products = [int(y) for y,r in zip(h0,c) if r==x and y>=0 and y!=x]
        if products:
            subs[x]=products
            possible[x]=[i for i,(r,y) in enumerate(zip(c,h0)) if r==x and y>=0]
            assert len(possible[x]) >= len(products)
    reference,pwm,widths,lo,hi=score_data(chimp,models)
    visited=np.zeros(n,dtype=np.uint8)
    starts=np.empty(n,dtype=np.int32)
    changed_init=np.flatnonzero(h0!=c).astype(np.int32)
    obs_fast=window_delta_affected(h0,changed_init,len(changed_init),reference,pwm,widths,lo,hi,visited,starts,0.9)
    if abs(obs_fast-observed)>1e-9:
        raise AssertionError(f'Optimized test statistic differs from original: fast {obs_fast}, original {observed}')
    if not subs:
        return 1.0,n_perms,float(observed),0
    count=0
    total=0.
    new_h=base.copy()
    changed=np.empty(n,dtype=np.int32)
    n_fixed=len(fixed)
    changed[:n_fixed]=fixed
    for it in range(n_perms):
        new_h[:]=base
        k=n_fixed
        for x,products in subs.items():
            slots=rng.sample(possible[x],len(products))
            newproducts=products.copy()
            rng.shuffle(newproducts)
            for pos,product in zip(slots,newproducts):
                new_h[pos]=product
                changed[k]=pos
                k+=1
        simulated=window_delta_affected(new_h,changed,k,reference,pwm,widths,lo,hi,visited,starts,0.9)
        if simulated >= observed-1e-12: count+=1
        total+=simulated
    return (count+1)/(n_perms+1),count,total/n_perms,len(subs)


def load_resume(path, npers):
    records=[]
    if path.is_file():
        with path.open(newline='',encoding='utf8') as f:
            records=list(csv.DictReader(f))
        if not records or 'empirical_p_within_HAR_sixTF_max' not in records[0]:
            raise ValueError('Existing S1 does not have the expected Stage20A schema')
        seen=set()
        for r in records:
            if r['screen_HAR'] in seen:raise ValueError('Duplicate HAR in resume S1')
            seen.add(r['screen_HAR'])
            if int(r['null_permutations'])!=npers:
                raise ValueError('Cannot resume from a run with a different number of permutations')
        print(f'Loaded {len(records)} already completed HARs from {path}',flush=True)
    return records


def main():
    human=original.read_fasta(ROOT/'STAGE20A_152_PROVISIONAL_HUMAN.fa')
    chimp=original.read_fasta(ROOT/'STAGE20A_152_PROVISIONAL_CHIMP.fa')
    assert len(human)==len(chimp)==152 and set(human)==set(chimp)
    with (ROOT/'STAGE20A_152_ALIGNMENT_QC.csv').open(newline='',encoding='utf-8-sig') as f:
        meta={r['screen_HAR']:r for r in csv.DictReader(f)}
    assert set(meta)==set(human)
    pwms=original.load_pwms(original.get_memes().read_text())
    existing=load_resume(OUT/'S1_HAR_152_SixTF_FamilyWise_ConditionalNull.csv',N_PERM)
    resumed_count=len(existing)
    done={r['screen_HAR'] for r in existing}
    stats=existing[:]
    tfstats=[]
    s2path=OUT/'S2_912_Pairwise_Aligned_TF_Descriptive.csv'
    if s2path.exists():
        with s2path.open(newline='',encoding='utf8') as f:
            tfstats=[r for r in csv.DictReader(f) if r['screen_HAR'] in done]
    stored={(r['screen_HAR'],r['TF']) for r in tfstats}
    assert not (len(tfstats)>len(done)*6)
    start=time.time()
    ordered=sorted(human)
    for i,k in enumerate(ordered,1):
        h=human[k];c=chimp[k]
        assert set(h)<=set('ACGT') and set(c)<=set('ACGT')
        ori=meta[k]['alignment_orientation']
        if ori=='reverse_complement':c=str(Seq(c).reverse_complement())
        elif ori!='forward':raise ValueError('Unexpected alignment orientation '+str(ori))
        h,c=original.aligned_strings(h,c)
        cached_ref={tf:original.scan_scores(original.encoded(c),model) for tf,model in pwms.items()}
        observed,top_tf,per=original.candidate_stats(h,c,pwms,reference_cache=cached_ref)
        for tf,data in per.items():
            if (k,tf) not in stored:
                tfstats.append({'screen_HAR':k,'TF':tf,'motif_id':original.TF_MOTIFS[tf],**data,
                    'interpretation':'position-matched exploratory PWM score; no measured TF occupancy'})
        if k in done:
            continue
        rng=random.Random(SEED+1009*i) # stable PER-LOCUS random stream for resumability
        p,extreme,mean_null,_=locus_simulations(h,c,pwms,rng,N_PERM,observed)
        sites=sum(x in 'ACGT' and y in 'ACGT' and x!=y for x,y in zip(h,c))
        stats.append({'screen_HAR':k,'catalog_id_provisional':meta[k]['catalog_id_provisional'],
            'sequence_provenance':'PROVISIONAL_SECONDARY_CATALOG_NOT_ORIGINAL_BED',
            'alignment_orientation':ori,'aligned_length':len(h),'aligned_substitutions':sites,
            'source_guide_distance_flag':'BEYOND_300BP_PROVISIONAL_CORE' if k in {'HACNS336','HACNS643','HAR169'} else 'NONE',
            'observed_max_abs_pwm_delta':round(observed,6),'top_TF_descriptive':top_tf or '',
            'null_mean_max_abs_pwm_delta':round(mean_null,6),'null_extreme_count':extreme,
            'null_permutations':N_PERM,'empirical_p_within_HAR_sixTF_max':p})
        done.add(k)
        if (len(done)%CHECKPOINT==0 or len(done)==152):
            stats.sort(key=lambda row:row['screen_HAR'])
            tfstats.sort(key=lambda row:(row['screen_HAR'],row['TF']))
            original.csvout(OUT/'S1_HAR_152_SixTF_FamilyWise_ConditionalNull.csv',stats)
            original.csvout(s2path,tfstats)
            elapsed=time.time()-start
            print(f'CHECKPOINT: {len(done)}/152 HARs completed | elapsed {elapsed:.1f}s | results in {OUT}',flush=True)
    stats.sort(key=lambda row:row['screen_HAR'])
    tfstats.sort(key=lambda row:(row['screen_HAR'],row['TF']))
    assert len(stats)==152 and len(tfstats)==912 and len({x['screen_HAR'] for x in stats})==152
    pvals=[float(r['empirical_p_within_HAR_sixTF_max']) for r in stats]
    q=original.adjust_bh(pvals)
    for r,v in zip(stats,q):r['BH_q_152_HAR_families']=float(v)
    rest=[r for r in stats if r['source_guide_distance_flag']=='NONE']
    q149=original.adjust_bh([float(r['empirical_p_within_HAR_sixTF_max']) for r in rest])
    for r,v in zip(rest,q149):r['BH_q_149_nonoutlier_sensitivity']=float(v)
    for r in stats:
        if r['source_guide_distance_flag']!='NONE':r['BH_q_149_nonoutlier_sensitivity']='EXCLUDED_FROM_SENSITIVITY'
    original.csvout(OUT/'S1_HAR_152_SixTF_FamilyWise_ConditionalNull.csv',stats)
    original.csvout(s2path,tfstats)
    status={
      'stage':'20A','status':'complete','runner':'numba_affected_windows_v2',
      'original_author_BED_confirmed':False,'provisional_pair_count':len(stats),
      'resumed_existing_loci':resumed_count,
      'descriptive_TF_records':len(tfstats),'n_permutations_per_locus':N_PERM,
      'random_seed_strategy':'20261120 + 1009 * sorted_har_1based_index',
      'minimum_p':min(pvals),'minimum_BH_q_152':float(min(q)),
      'count_BH_q_below_0p05':int(sum(q<.05)),
      'n_nonoutlier_har_sensitivity':len(rest),'minimum_q_149_nonoutlier':float(min(q149)),
      'TFs':original.TF_MOTIFS,'PWM_SHA256':original.EXPECTED_MEME_SHA,
      'threshold':0.9,
      'null':'reference-base-matched mutation relocation, fixed aligned indel columns, within-HAR six-TF maximum',
      'multiplicity':'BH across 152 HAR maximum-statistic tests; no separate 912-test motif q-values',
      'comparison_note':'Same statistic and conditional null as original Stage20A, but independent reproducible per-HAR RNG stream; Monte Carlo p-values may differ slightly from slow runner.',
      'limitations':['ortholog identities provisional','dinucleotide frequencies not preserved','non-HAR genomic null not evaluated','no phylogenetic derived-state inference','0.90 PWM scores not FIMO p-values','no experimental TF occupancy or HAR-to-gene validation']
    }
    (OUT/'S0_Run_QC_and_Limitations.json').write_text(json.dumps(status,indent=2))
    zip_path=ROOT/'HAR_Stage20A_FullScreen_Results.zip'
    with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
        for pth in OUT.glob('*'):z.write(pth,arcname=pth.name)
    print(json.dumps(status,indent=2))
    print('RESULT:',zip_path)
    # Subprocesses do not have an active IPython notebook kernel; never fail the analysis
    # because a browser download is not available from this process.
    try:
        from google.colab import files
        from IPython import get_ipython
        ip = get_ipython()
        if ip is not None and getattr(ip, 'kernel', None) is not None:
            files.download(str(zip_path))
        else:
            print('AUTODOWNLOAD_SKIPPED: result ZIP is ready; download from a notebook cell.', flush=True)
    except Exception as download_error:
        print('AUTODOWNLOAD_SKIPPED:', type(download_error).__name__, str(download_error), flush=True)

if __name__=='__main__':main()
