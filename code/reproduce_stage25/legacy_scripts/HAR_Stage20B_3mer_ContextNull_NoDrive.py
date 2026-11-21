"""HAR Stage 20B | pre-specified 3-mer reference-context substitution-placement null.

All source ortholog mappings are PROVISIONAL. This is a sensitivity comparison
against Stage 20A, not a dinucleotide-count-preserving shuffle or a matched
non-HAR genomic experiment. No novel biological mechanism is inferred.
"""
from __future__ import annotations
import csv, json, math, hashlib, random, os, time, traceback, zipfile, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio.Seq import Seq

ROOT=Path('/content') if Path('/content').exists() else Path.cwd()
sys.path.insert(0,str(ROOT))
import HAR_Stage20A_ConditionalNull_NoDrive as orig
import HAR_Stage20A_Fast_Resumable as fast

OUT=ROOT/'HAR_Stage20B_3mer_ContextNull_Results'
OUT.mkdir(parents=True,exist_ok=True)
N_PERM=int(os.environ.get('HAR20B_PERMS','4999'))
if not (999<=N_PERM<=200000):raise ValueError('Set HAR20B_PERMS >=999, ideally 4999 or more')
SEED=20261120
RESULT=OUT/'S1_152_HAR_ContextMatched_Statistics.csv'
FIELDS=['screen_HAR','catalog_id_provisional','source_identity_status',
    'observed_max_abs_pwm_delta','Stage20A_baseMatched_p','Stage20A_baseMatched_BH_q',
    'n_aligned_substitutions','n_movable_3mer_substitutions','n_fixed_no_3mer_substitutions',
    'n_context_strata','n_context_slots_total','n_strata_with_no_position_freedom',
    'n_null_permutations','extreme_null_count','context_matched_p','mean_null_max_delta',
    'context_BH_q_152','context_BH_q_149_nonoutlier','flag_source_guide_distance','null_interpretation']


def filehash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def csv_read(path):
    with open(path,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))


def csv_write(path, rows, fields):
    temporary=Path(str(path)+'.tmp')
    with open(temporary,'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
    temporary.replace(path)


def prepare_trimer_null(h_al, c_al):
    """Context defined on chimp as a sequence *reference*, not inferred ancestor.

    Substitution types (chimp x -> human y), and the chimp reference 3-mer at
    the original substitution, are preserved for movable substitutions.
    Indel columns and substitutions without a valid gap-free 3-mer are fixed.
    Candidate sites have the same chimp 3-mer, intact human center, and may be
    selected only once within their reference-context stratum.
    """
    h=orig.encoded(h_al); c=orig.encoded(c_al); n=len(c)
    base=c.copy()
    fixed=[]; groups=defaultdict(list)
    n_sub=0
    for i in range(n):
        if h[i]<0 or c[i]<0:
            base[i]=h[i]
            if h[i]!=c[i]:fixed.append(i)
        elif h[i]!=c[i]:
            n_sub+=1
            if i==0 or i==n-1 or c[i-1]<0 or c[i+1]<0:
                base[i]=h[i]
                fixed.append(i)
            else:
                groups[(int(c[i-1]),int(c[i]),int(c[i+1]))].append(int(h[i]))
    eligible=defaultdict(list)
    for i in range(1,n-1):
        if h[i]>=0 and c[i-1]>=0 and c[i]>=0 and c[i+1]>=0:
            eligible[(int(c[i-1]),int(c[i]),int(c[i+1]))].append(i)
    records=[]
    for context in sorted(groups):
        product=groups[context]; slots=eligible[context]
        if len(slots)<len(product):
            raise AssertionError('Not enough context-compatible positions')
        records.append((context,tuple(product),tuple(slots)))
    movable=sum(len(x[1]) for x in records)
    nfixed=n_sub-movable
    assert movable+nfixed==n_sub
    zero_freedom=sum(len(slots)==len(prods) and len(set(prods))==1 for _,prods,slots in records)
    return h,c,base,np.array(fixed,dtype=np.int32),records,dict(
        n_aligned_substitutions=n_sub,
        n_movable_3mer_substitutions=movable,
        n_fixed_no_3mer_substitutions=nfixed,
        n_context_strata=len(records),
        n_context_slots_total=sum(len(x[2]) for x in records),
        n_strata_with_no_position_freedom=zero_freedom)


def simulate_3mer(h,c,base,fixed,records,models,rng,nperms,observed):
    reference,pwm,widths,lo,hi=fast.score_data(''.join('ACGT'[x] if x>=0 else '-' for x in c),models)
    n=len(c); marked=np.zeros(n,dtype=np.uint8);starts=np.empty(n,dtype=np.int32)
    changed=np.flatnonzero(h!=c).astype(np.int32)
    obsfast=fast.window_delta_affected(h,changed,len(changed),reference,pwm,widths,lo,hi,marked,starts,0.9)
    if abs(obsfast-observed)>1e-9:
        raise AssertionError('Fast observed statistic does not match independent full scan')
    # The test statistic is nonnegative, hence 0 leads to an exact upper-tail p=1.
    if observed<=1e-12:return 1.0,nperms,0.0,True
    trial=np.empty(n,dtype=np.int8)
    changed=np.empty(n,dtype=np.int32)
    fixed_num=len(fixed);changed[:fixed_num]=fixed
    moreextreme=0;nullsum=0.0
    for iteration in range(nperms):
        trial[:]=base;k=fixed_num
        for _,products,options in records:
            picks=rng.sample(options,len(products))
            ys=list(products);rng.shuffle(ys)
            for pos,derived in zip(picks,ys):
                trial[pos]=derived;changed[k]=pos;k+=1
        stat=fast.window_delta_affected(trial,changed,k,reference,pwm,widths,lo,hi,marked,starts,0.9)
        moreextreme+=(stat>=observed-1e-12)
        nullsum+=stat
    return (moreextreme+1)/(nperms+1),moreextreme,nullsum/nperms,False


def validate_inputs():
    required=['STAGE20A_152_PROVISIONAL_HUMAN.fa','STAGE20A_152_PROVISIONAL_CHIMP.fa',
      'STAGE20A_152_ALIGNMENT_QC.csv','STAGE20A_S1_HAR_152.csv','STAGE20A_S2_912_TF.csv']
    absent=[x for x in required if not (ROOT/x).is_file()]
    if absent:raise FileNotFoundError('Notebook embedded inputs missing: '+str(absent))
    hashes={f:filehash(ROOT/f) for f in required}
    return hashes


def main():
    started=time.monotonic();hashes=validate_inputs()
    human=orig.read_fasta(ROOT/'STAGE20A_152_PROVISIONAL_HUMAN.fa')
    chimp=orig.read_fasta(ROOT/'STAGE20A_152_PROVISIONAL_CHIMP.fa')
    meta={r['screen_HAR']:r for r in csv_read(ROOT/'STAGE20A_152_ALIGNMENT_QC.csv')}
    previous={r['screen_HAR']:r for r in csv_read(ROOT/'STAGE20A_S1_HAR_152.csv')}
    assert len(human)==len(chimp)==len(meta)==len(previous)==152
    assert set(human)==set(chimp)==set(meta)==set(previous)
    pfile=orig.get_memes()  # original pinned JASPAR2024 source, SHA-256 checked
    models=orig.load_pwms(pfile.read_text())
    assert len(models)==6
    pwmhash=filehash(pfile)
    contract_path=OUT/'S0_Resume_Input_Contract.json'
    contract={'input_file_sha256':hashes,'pwm_sha256':pwmhash,'n_permutations':N_PERM,'test_version':'3mer-chimp-reference-v1'}
    if contract_path.exists():
        if json.loads(contract_path.read_text())!=contract:
            raise RuntimeError('Existing checkpoint input contract differs; move the prior output directory before starting a different run.')
    else:
        contract_path.write_text(json.dumps(contract,indent=2))
    # Resume only from completed 20B output with the current permutation count and same verified inputs.
    if RESULT.is_file():
        old=csv_read(RESULT)
        for row in old:
            count=int(row['n_null_permutations'])
            if count not in (0,N_PERM) or (count==0 and float(row['observed_max_abs_pwm_delta'])>0):
                raise RuntimeError('Existing 20B results have inconsistent N_PERM. Move the prior output folder to start a new run.')
    else:old=[]
    rows={r['screen_HAR']:r for r in old}
    assert len(rows)==len(old)
    ordered=sorted(human)
    preflight=[]
    for i,key in enumerate(ordered,1):
        m=meta[key];prior=previous[key]
        h=human[key];c=chimp[key]
        if m['alignment_orientation']=='reverse_complement':c=str(Seq(c).reverse_complement())
        elif m['alignment_orientation']!='forward':raise ValueError('Unknown orientation '+m['alignment_orientation'])
        h,c=orig.aligned_strings(h,c)
        ha,ca,base,fixed,groups,details=prepare_trimer_null(h,c)
        if int(prior['aligned_substitutions'])!=details['n_aligned_substitutions']:
            raise AssertionError('Stage20A count mismatch at '+key)
        if key in rows:continue
        observed,TF,parts=orig.candidate_stats(h,c,models)
        if abs(observed-float(prior['observed_max_abs_pwm_delta']))>0.000001:
            raise AssertionError('Stage20A motif score mismatch at '+key)
        rng=random.Random(SEED+35023*i)
        if observed<=1e-12:
            p,extreme,nullmean,exact_zero=1.0,'NOT_SIMULATED_EXACT_P1',None,True  # mathematically exact
        else:
            p,extreme,nullmean,exact_zero=simulate_3mer(ha,ca,base,fixed,groups,models,rng,N_PERM,observed)
        rows[key]={
          'screen_HAR':key,'catalog_id_provisional':m['catalog_id_provisional'],
          'source_identity_status':'PROVISIONAL_SECONDARY_CATALOG_NOT_ORIGINAL_AUTHOR_BED',
          'observed_max_abs_pwm_delta':round(observed,6),
          'Stage20A_baseMatched_p':prior['empirical_p_within_HAR_sixTF_max'],
          'Stage20A_baseMatched_BH_q':prior['BH_q_152_HAR_families'],
          **details,'n_null_permutations':0 if exact_zero else N_PERM,'extreme_null_count':extreme,
          'context_matched_p':p,'mean_null_max_delta':round(nullmean,8) if nullmean is not None else '',
          'context_BH_q_152':'PENDING_FULL_RUN','context_BH_q_149_nonoutlier':'PENDING_FULL_RUN',
          'flag_source_guide_distance':prior['source_guide_distance_flag'],
          'null_interpretation':'exact_p_1_for_zero_statistic' if exact_zero else 'chimp_reference_3mer_conditional_mutation_relocation'}
        csv_write(RESULT,[rows[k] for k in sorted(rows)],FIELDS)
        print(f'CHECKPOINT {len(rows)}/152: {key}, trimer_p={p:.5f}, variable_subs={details["n_movable_3mer_substitutions"]}; {time.monotonic()-started:.1f}s',flush=True)
    data=[rows[key] for key in ordered]
    assert len(data)==152
    q=orig.adjust_bh([float(r['context_matched_p']) for r in data]);
    for r,qv in zip(data,q):r['context_BH_q_152']=float(qv)
    subset=[r for r in data if r['flag_source_guide_distance']=='NONE']
    q149=orig.adjust_bh([float(r['context_matched_p']) for r in subset]); assert len(subset)==149
    for r,qv in zip(subset,q149):r['context_BH_q_149_nonoutlier']=float(qv)
    for r in data:
        if r['flag_source_guide_distance']!='NONE':r['context_BH_q_149_nonoutlier']='EXCLUDED'
    csv_write(RESULT,data,FIELDS)
    csv_write(OUT/'S2_16_Nonzero_Locus_Comparisons.csv',
        [r for r in data if float(r['observed_max_abs_pwm_delta'])>1e-12],FIELDS)
    (OUT/'S3_912_TF_Descriptive_Stage20A_NoNewTests.csv').write_bytes((ROOT/'STAGE20A_S2_912_TF.csv').read_bytes())
    assert sum(float(r['observed_max_abs_pwm_delta'])>0 for r in data)==16
    num_active=sum(float(r['observed_max_abs_pwm_delta'])>0 for r in data)
    summary={
       'stage':'20B','status':'complete','scientific_model':'chimp-reference trinucleotide-conditioned substitution placement; six-TF maximum',
       'original_author_exact_BED_verified':False,'ortholog_pairs_provisional':152,
       'num_TF_motifs':6,'motif_source':'pinned JASPAR2024 secondary catalog',
       'source_commit':orig.COMMIT,'motif_sha256':pwmhash,'input_file_sha256':hashes,
       'n_simulations_per_nonzero_locus':N_PERM,
       'number_nonzero_loci_simulated':num_active,'zero_loci_exact_p_1':152-num_active,
       'null_context':'chimp reference 3mer for every movable substitution; reference indels/edge-substitutions remain fixed',
       'input_not_ancestral':'Chimp genome reference used as context, NOT proof of ancestral state',
       'statistic':'max absolute difference in relative PWM score over same gap-free aligned motif windows, only where score>=0.9 in either species',
       'hypothesis_family':'152 within-locus 6-TF maxima; BH across all 152, not 912 per-TF confirmatory q values',
       'minimum_nominal_p':float(min(float(x['context_matched_p']) for x in data)),
       'minimum_BH_q_152':float(min(q)),'n_BH_q_below_0p05':int(sum(q<0.05)),
       'minimum_BH_q_149_sensitivity':float(min(q149)),'n_BH_q_below_0p05_149_sensitivity':int(sum(q149<0.05)),
       'n_loci_with_any_nonmovable_substitution':sum(int(r['n_fixed_no_3mer_substitutions'])>0 for r in data),
       'n_loci_with_no_3mer_relocation_possible':sum(int(r['n_movable_3mer_substitutions'])==0 for r in data),
       'limitations':['NOT dinucleotide-count-preserving', 'NOT matched non-HAR genomic controls',
          'NOT a phylogenetically informed ancestral-state test',
          'null may be degenerate at highly constrained short loci',
          'exploratory PWM>=0.90 is not FIMO p calibration','original HAR identity based on third-party sequence catalog',
          'no measured TF binding or new HAR-gene causal validation']}
    (OUT/'S0_Stage20B_RunStatus_and_Limitations.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    archive=ROOT/'HAR_Stage20B_ContextMatched_Results.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for f in sorted(OUT.glob('*')):z.write(f,arcname=f.name)
    print('FINISHED Stage20B. ZIP:',archive,flush=True)
    print(json.dumps(summary,indent=2),flush=True)
    return summary


if __name__=='__main__':
    main()
