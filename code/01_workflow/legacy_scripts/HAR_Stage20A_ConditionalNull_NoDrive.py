"""Stage 20A: full 152 provisional HAR ortholog pairs, aligned-window motif-change
screen and within-HAR substitution-placement conditional null. Exploratory.

NOT FIMO calibration, not dinucleotide-preserving, not functional TF binding evidence.
Never label provisional third-party HAR coordinates author-primary validated.
"""
import csv, json, math, hashlib, os, re, time, zipfile, random
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import requests

ROOT=Path('/content') if Path('/content').exists() else Path.cwd()
OUT=ROOT/'HAR_Stage20A_FullScreen_Results'
OUT.mkdir(exist_ok=True,parents=True)
N_PERM=int(os.environ.get('HAR_N_PERM','4999'))
RNG_SEED=20261120
THRESHOLD=0.90
COMMIT='a9251f39976bbc15df1d8f0dfedaf9786819c298'
MEME_URL=f'https://raw.githubusercontent.com/athenamarou/HAR-TFBS-Project/{COMMIT}/data/meme_file.txt'
EXPECTED_MEME_SHA='ae4f77a17256c7b8b7b73402810c1ecf4fc411e9ad7832cdc87d586dee9a4a15'
TF_MOTIFS={'SOX2':'MA0143.5','PAX6':'MA0069.1','REST':'MA0138.3','SOX11':'MA0869.3','CTCF':'MA0139.2','NFIA':'MA0670.2'}
BASEMAP={'A':0,'C':1,'G':2,'T':3,'-':-1}


def read_fasta(path):
    seqs={}; k=None; block=[]
    for line in Path(path).read_text().splitlines():
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if k is not None: seqs[k]=''.join(block).upper()
            k=line[1:].split('|')[0];block=[]
        else: block.append(line)
    if k is not None: seqs[k]=''.join(block).upper()
    return seqs


def get_memes():
    path=ROOT/'stage20a_pinned_JASPAR2024_meme.txt'
    if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest()==EXPECTED_MEME_SHA:
        return path
    r=requests.get(MEME_URL,timeout=100,headers={'User-Agent':'HAR-stage20a-source-validation/1.0'})
    r.raise_for_status()
    if hashlib.sha256(r.content).hexdigest()!=EXPECTED_MEME_SHA:
        raise ValueError('PINNED PWM SOURCE CHECKSUM MISMATCH; abort')
    path.write_bytes(r.content)
    return path


def load_pwms(text):
    lines=text.splitlines(); parsed={}; want=set(TF_MOTIFS.values()); i=0
    while i<len(lines):
        if lines[i].startswith('MOTIF '):
            mid=lines[i].split()[1]
            if mid in want:
                j=i+1
                while j<len(lines) and not lines[j].startswith(('letter-probability matrix:', 'MOTIF ')): j+=1
                if j<len(lines) and lines[j].startswith('letter-probability matrix:'):
                    m=re.search(r'w=\s*(\d+)',lines[j]); width=int(m.group(1))
                    raw=np.array([[float(z) for z in lines[j+k+1].split()] for k in range(width)],dtype=np.float64)
                    if raw.shape!=(width,4):raise ValueError(f'Bad PWM {mid}')
                    pwm=np.log2(((raw+0.001)/1.004)/0.25)
                    lo=float(pwm.min(axis=1).sum()); hi=float(pwm.max(axis=1).sum())
                    parsed[mid]=(pwm,lo,hi)
                    i=j+width
        i+=1
    if set(parsed)!=want:raise RuntimeError('Missing motifs '+str(want-set(parsed)))
    return {name:parsed[tid] for name,tid in TF_MOTIFS.items()}


def aligned_strings(human,chimp):
    a=PairwiseAligner()
    a.mode='global';a.match_score=2;a.mismatch_score=-1;a.open_gap_score=-4;a.extend_gap_score=-1
    obj=a.align(human,chimp)[0];coordinates=obj.coordinates
    hh=[];cc=[]
    for i in range(coordinates.shape[1]-1):
        hs,he=int(coordinates[0,i]),int(coordinates[0,i+1]); cs,ce=int(coordinates[1,i]),int(coordinates[1,i+1])
        dh=he-hs;dc=ce-cs
        assert dh>=0 and dc>=0
        if dh and dc:
            assert dh==dc;hh.append(human[hs:he]);cc.append(chimp[cs:ce])
        elif dh:hh.append(human[hs:he]);cc.append('-'*dh)
        elif dc:hh.append('-'*dc);cc.append(chimp[cs:ce])
    human_al=''.join(hh);chimp_al=''.join(cc)
    assert len(human_al)==len(chimp_al)
    assert human_al.replace('-','')==human and chimp_al.replace('-','')==chimp
    return human_al,chimp_al


def encoded(seq):return np.fromiter((BASEMAP.get(x,-1) for x in seq),dtype=np.int8,count=len(seq))


def scan_scores(arr,info):
    """Two-strand relative PWM scores at SAME alignment coordinates; gap-spanning invalid."""
    pwm,lo,hi=info; w=pwm.shape[0]; n=len(arr)
    if n<w: return np.empty((2,0),dtype=np.float64)
    v=sliding_window_view(arr,w)
    valid=np.all(v>=0,axis=1)
    vals=v.clip(0,3).astype(np.intp)
    pos=np.arange(w)
    forward=pwm[pos,vals].sum(axis=1)
    reverse=pwm[pos,3-vals[:,::-1]].sum(axis=1)
    scores=np.stack(((forward-lo)/(hi-lo),(reverse-lo)/(hi-lo)))
    scores[:,~valid]=-math.inf
    return scores


def candidate_stats(h,c,models,reference_cache=None):
    """Within each HAR: max |normalized score contrast| if either species
       has a window crossing exploratory .90. Full six-TF family tested jointly.
    """
    hg=encoded(h); cg=encoded(c); per={}; allmax=0.; tf_best=None
    for tf,info in models.items():
        a=scan_scores(hg,info); b=(reference_cache[tf] if reference_cache is not None else scan_scores(cg,info))
        if a.shape!=b.shape or a.size==0:
            per[tf]={'max_abs_delta':None,'n_crossing_windows':None,'max_score_h':None,'max_score_c':None}
            continue
        eligible=np.isfinite(a)&np.isfinite(b)
        hits=eligible&((a>=THRESHOLD)|(b>=THRESHOLD))
        d=np.where(hits,np.abs(np.nan_to_num(a-b,nan=0,posinf=0,neginf=0)),0.0)
        maxdelta=float(d.max()) if d.size else 0.
        cross=int(np.sum(eligible&((a>=THRESHOLD)!=(b>=THRESHOLD))))
        only_h=np.where(np.isfinite(a),a,np.nan);only_c=np.where(np.isfinite(b),b,np.nan)
        maxh=float(np.nanmax(only_h)) if np.isfinite(a).any() else None
        maxc=float(np.nanmax(only_c)) if np.isfinite(b).any() else None
        per[tf]={'max_abs_delta':maxdelta,'n_crossing_windows':cross,'max_score_h':maxh,'max_score_c':maxc}
        if maxdelta>allmax:
            allmax=maxdelta;tf_best=tf
    return allmax,tf_best,per


def mutate_positions(h,c,rng):
    """Null: relocate each observed aligned chimp->human substitution to a
    random reference-compatible site. Preserve substitution types/counts,
    base-specific reference positions and indel alignment columns. This DOES
    NOT preserve dinucleotide frequencies, phylogeny, or local mutation rates.
    """
    ha=list(h);ca=list(c)
    substitutions=defaultdict(list)
    for i,(x,y) in enumerate(zip(ca,ha)):
        if x in 'ACGT' and y in 'ACGT' and x!=y:substitutions[x].append(y)
    z=ca.copy()
    for i,(x,y) in enumerate(zip(ca,ha)):
        if x=='-' or y=='-': z[i]=y # preserve all indels exactly
    for x,products in substitutions.items():
        options=[i for i,(reference,original) in enumerate(zip(ca,ha)) if reference==x and original in 'ACGT']
        if len(options)<len(products):raise ValueError('Incompatible null positions')
        slots=rng.sample(options,len(products)); perm=products[:];rng.shuffle(perm)
        for i,product in zip(slots,perm):z[i]=product
    return ''.join(z),sum(map(len,substitutions.values()))


def adjust_bh(pvals):
    n=len(pvals);order=np.argsort(pvals,kind='mergesort'); ranks=np.arange(1,n+1)
    sorted_q=np.array(pvals)[order]*n/ranks
    sorted_q=np.minimum.accumulate(sorted_q[::-1])[::-1]
    ret=np.empty(n);ret[order]=np.clip(sorted_q,0,1)
    return ret


def csvout(path,rows):
    if not rows:return
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),extrasaction='ignore');w.writeheader();w.writerows(rows)


def main():
    human=read_fasta(ROOT/'STAGE20A_152_PROVISIONAL_HUMAN.fa')
    chimp=read_fasta(ROOT/'STAGE20A_152_PROVISIONAL_CHIMP.fa')
    assert len(human)==len(chimp)==152 and set(human)==set(chimp)
    with open(ROOT/'STAGE20A_152_ALIGNMENT_QC.csv',newline='',encoding='utf-8-sig') as f:
        alignmeta={r['screen_HAR']:r for r in csv.DictReader(f)}
    assert len(alignmeta)==152
    pwms=load_pwms(get_memes().read_text())
    rng=random.Random(RNG_SEED)
    stats=[];tfstats=[]; started=time.time()
    for idx,k in enumerate(sorted(human),1):
        h=human[k];c=chimp[k]
        assert set(h)<=set('ACGT') and set(c)<=set('ACGT')
        ori=alignmeta[k]['alignment_orientation']
        if ori=='reverse_complement':c=str(Seq(c).reverse_complement())
        elif ori!='forward':raise ValueError('Unexpected alignment orientation '+str(ori))
        h,c=aligned_strings(h,c)
        cached_ref={tf:scan_scores(encoded(c),info) for tf,info in pwms.items()}
        observed,top_tf,per=candidate_stats(h,c,pwms,reference_cache=cached_ref)
        sites=sum(x in 'ACGT' and y in 'ACGT' and x!=y for x,y in zip(h,c))
        if sites==0:
            p=1.;ge=N_PERM;mean_null=0.
        else:
            ge=0;total=0.
            for b in range(N_PERM):
                new_h,_=mutate_positions(h,c,rng)
                sim,_,_=candidate_stats(new_h,c,pwms,reference_cache=cached_ref)
                ge+=sim>=observed-1e-12;total+=sim
            p=(ge+1)/(N_PERM+1);mean_null=total/N_PERM
        stats.append({'screen_HAR':k,'catalog_id_provisional':alignmeta[k]['catalog_id_provisional'],
            'sequence_provenance':'PROVISIONAL_SECONDARY_CATALOG_NOT_ORIGINAL_BED',
            'alignment_orientation':ori,'aligned_length':len(h),'aligned_substitutions':sites,
            'source_guide_distance_flag': 'BEYOND_300BP_PROVISIONAL_CORE' if k in {'HACNS336','HACNS643','HAR169'} else 'NONE',
            'observed_max_abs_pwm_delta':round(observed,6),'top_TF_descriptive':top_tf or '',
            'null_mean_max_abs_pwm_delta':round(mean_null,6),'null_extreme_count':ge,
            'null_permutations':N_PERM,'empirical_p_within_HAR_sixTF_max':p})
        for tf,data in per.items():
            tfstats.append({'screen_HAR':k,'TF':tf,'motif_id':TF_MOTIFS[tf],
                **data,'interpretation':'position-matched exploratory PWM score; no measured TF occupancy'})
        if idx%10==0 or idx==152:
            print(f'PROGRESS: {idx}/152 provisional loci assessed',flush=True)
            csvout(OUT/'S1_HAR_152_SixTF_FamilyWise_ConditionalNull.csv',stats)
            csvout(OUT/'S2_912_Pairwise_Aligned_TF_Descriptive.csv',tfstats)
    q=adjust_bh([r['empirical_p_within_HAR_sixTF_max'] for r in stats]);
    for r,v in zip(stats,q):r['BH_q_152_HAR_families']=float(v)
    sensitivity=[r for r in stats if r['source_guide_distance_flag']=='NONE']
    q149=adjust_bh([r['empirical_p_within_HAR_sixTF_max'] for r in sensitivity])
    for r,v in zip(sensitivity,q149):r['BH_q_149_nonoutlier_sensitivity']=float(v)
    for r in stats:
        if r['source_guide_distance_flag']!='NONE':r['BH_q_149_nonoutlier_sensitivity']='EXCLUDED_FROM_SENSITIVITY'
    csvout(OUT/'S1_HAR_152_SixTF_FamilyWise_ConditionalNull.csv',stats)
    csvout(OUT/'S2_912_Pairwise_Aligned_TF_Descriptive.csv',tfstats)
    status={'stage':'20A','status':'complete','original_author_BED_confirmed':False,
       'provisional_pair_count':len(stats),'descriptive_TF_records':len(tfstats),
       'n_permutations_per_locus':N_PERM,'random_seed':RNG_SEED,
       'min_p':min(r['empirical_p_within_HAR_sixTF_max'] for r in stats),
       'min_q':float(min(q)),'n_BH_q_below_0p05':int(sum(q<0.05)),
       'n_nonoutlier_har_sensitivity':len(sensitivity),'min_BH_q_nonoutlier_sensitivity':float(min(q149)),
       'TFs':TF_MOTIFS,'PWM_SHA256':EXPECTED_MEME_SHA,'threshold':THRESHOLD,
       'null':'within-HAR reference-base-matched mutation repositioning; six-motif maximum statistic',
       'multiplicity':'BH across 152 HAR-level maximum-statistic tests; no 908 single-motif confirmatory q values',
       'limitations':['source ortholog pairs not independently author-BED matched',
                      'reference-base preserving is NOT dinucleotide preserving',
                      'no phylogenetic ancestry or non-HAR matched controls',
                      'relative PWM>=0.90 exploratory and not FIMO-calibrated',
                      'observed conditional-null association is not binding, enhancer function or causal HAR-gene regulation']}
    (OUT/'S0_Run_QC_and_Limitations.json').write_text(json.dumps(status,indent=2))
    with zipfile.ZipFile(ROOT/'HAR_Stage20A_FullScreen_Results.zip','w',zipfile.ZIP_DEFLATED) as z:
        for item in OUT.iterdir():z.write(item,arcname=item.name)
    print(json.dumps(status,indent=2));print('RESULT:',ROOT/'HAR_Stage20A_FullScreen_Results.zip')
    try:
        from google.colab import files
        files.download(str(ROOT/'HAR_Stage20A_FullScreen_Results.zip'))
    except ImportError:pass

if __name__=='__main__':main()
