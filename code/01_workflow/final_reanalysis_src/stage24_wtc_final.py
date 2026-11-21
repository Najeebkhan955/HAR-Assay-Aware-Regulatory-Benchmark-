#!/usr/bin/env python3
"""Stage24: source-verified, full six-file WTc/N2 lentiMPRA comparison.

GEO GSE110758, two separate HUMAN iPSC genetic backgrounds HS1-11 and WTc.
Human/chimp in this analysis denotes inserted reporter DNA alleles, not host species.
Only same construct/cut events are compared; no causal regulatory mechanism claimed.
"""
import csv, gzip, hashlib, json, math, re, time, traceback, zipfile
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path('/content') if Path('/content').exists() else Path.cwd()
INPUT=ROOT/'HAR_Stage24_Embedded'
OUT=ROOT/'HAR_Stage24_WTc_Results'
CACHE=ROOT/'HAR_Stage24_WTc_GEO_Cache'
for p in (INPUT,OUT,CACHE):p.mkdir(parents=True,exist_ok=True)
SAMPLES=[(1,'GSM3016026','GSM3016044'),(2,'GSM3016027','GSM3016045'),(3,'GSM3016028','GSM3016046')]
PRIORITY=('2xHAR.70','2xHAR.310','2xHAR.234','2xHAR.377','2xHAR.3')
EXPECTED_MAP_SHA='6f30451ac623ed4948c8649aa3e9113f0d244b5d70db1f37fbb35314cf92cdeb'
SEED=20261120

def url_for(gsm,assay,rep):
    return f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:7]}nnn/{gsm}/suppl/{gsm}_human-WTc-N2-{rep}-{assay}.tsv.gz'

def fetch(u,name):
    p=CACHE/name
    if p.exists() and p.stat().st_size>1000:
        with p.open('rb') as f:tag=f.read(2)
        if tag==b'\x1f\x8b':
            print('Using cached',name);return p
    import requests
    errors=[]
    for attempt in range(1,4):
        try:
            print('Fetching',name,'attempt',attempt)
            n=0
            with requests.get(u,stream=True,timeout=(25,120),headers={'User-Agent':'HAR-Stage24-Research-Reproducibility/1.0'}) as response:
                response.raise_for_status()
                with p.open('wb') as f:
                    for block in response.iter_content(1<<20):
                        if block:
                            n+=len(block)
                            if n>80_000_000:raise ValueError('Unexpectedly large GEO sample file')
                            f.write(block)
            if n<1000:raise ValueError('File too small')
            with p.open('rb') as f:
                if f.read(2)!=b'\x1f\x8b':raise ValueError('Missing gzip magic')
            with gzip.open(p,'rb') as f:
                if not f.read(160):raise ValueError('Gzip has no data')
            return p
        except Exception as ex:
            errors.append(repr(ex));print('Retry after error:',str(ex)[:230])
            p.unlink(missing_ok=True);time.sleep(2)
    raise RuntimeError('Cannot download '+name+': '+repr(errors))

def load_design():
    p=INPUT/'S1B_Assayed_Barcode_Map.csv'
    raw=p.read_bytes();sha=hashlib.sha256(raw).hexdigest()
    if sha!=EXPECTED_MAP_SHA:raise ValueError('Published selected design barcode map SHA mismatch')
    result={};bad=0
    for row in csv.DictReader(raw.decode('utf-8').splitlines()):
        barcode=row['Barcode'].strip().upper()
        if not re.fullmatch(r'[ACGT]{15}',barcode):raise ValueError('Bad mapped barcode')
        key=(row['HAR'],row['allele'],row['cut'],row['insert'])
        if barcode in result and result[barcode]!=key:bad+=1
        result[barcode]=key
    if bad:raise ValueError('Ambiguous barcode assignments: '+str(bad))
    if len(result)!=25500:raise ValueError('Unexpected selected barcode count')
    return result,sha

def parse_counts(p,mapping):
    """Parse all HEADERLESS lines and verify original source barcode->insert annotations."""
    counts=Counter();annotations={}
    reads=0;nrows=0;nselected=0;nzero=0;dupes=0
    with gzip.open(p,'rt',encoding='utf8') as f:
        for i,row in enumerate(csv.reader(f,delimiter='\t'),1):
            if len(row)!=3:raise ValueError(f'{p.name} invalid 3-column structure at line {i}')
            bc=row[0].strip().upper()
            if not re.fullmatch(r'[ACGT]{15}',bc):raise ValueError(f'Invalid barcode at line {i} {bc[:30]}')
            num=int(row[1])
            if num<0:raise ValueError('Negative read count')
            # Last colon separates original oligo barcode number, not hg19 chr:coordinates.
            fragment=row[2].rsplit(':',1)[0] if re.search(r':\d+$',row[2]) else row[2]
            if bc in mapping:
                if fragment!=mapping[bc][3]:
                    raise ValueError(f'Source insert ID disagrees with pinned design: {p.name} line {i}')
                nselected+=1
                if bc in annotations and annotations[bc]!=fragment:
                    raise ValueError('Repeated barcode assigned to conflicting inserts')
                annotations[bc]=fragment
            if bc in counts:dupes+=1
            counts[bc]+=num;reads+=num;nrows+=1
            if num==0:nzero+=1
    if nrows<1000 or reads<1000:raise ValueError('Incomplete or invalid GEO file')
    status=dict(file=p.name,headerless=True,parsed_records=nrows,all_reads=reads,selected_design_rows=nselected,
                unique_barcodes=len(counts),duplicate_barcode_rows=dupes,zero_count_rows=nzero,
                selected_insert_annotation_mismatches=0)
    return counts,status

def per_fragment(rna,dna,mapping,rep,threshold=1,pseudo=1):
    grouped=defaultdict(lambda:[0,0,0])
    for bc,key in mapping.items():
        depth=dna.get(bc,0)
        if depth<threshold:continue
        grouped[key][0]+=rna.get(bc,0);grouped[key][1]+=depth;grouped[key][2]+=1
    effects={};frag=[]
    eligible=defaultdict(dict)
    for (har,allele,cut,ins),(rc,dc,nbar) in grouped.items():
        if nbar<10:continue
        val=math.log2((rc+pseudo)/(dc+pseudo))
        eligible[(har,cut)].setdefault(allele,[]).append((val,nbar))
        if threshold==1 and pseudo==1:
            frag.append(dict(HAR=har,allele=allele,cut=cut,culture=rep,RNA_sum=rc,DNA_sum=dc,
                             DNApositive_barcodes=nbar,log2_activity=val))
    for (har,cut),sp in eligible.items():
        if len(sp.get('human',[]))==1 and len(sp.get('chimp',[]))==1:
            human=sp['human'][0];chimp=sp['chimp'][0]
            effects[(har,cut)]=(human[0]-chimp[0],human[1],chimp[1])
    return effects,frag

def bootstrap_between_lines(j,nboot=5000):
    """Resample HAR clusters, not individual cuts; descriptive sampling uncertainty."""
    rng=np.random.default_rng(SEED)
    byhar=[frame[['HS111_mean','WTc_mean']].to_numpy(float)
           for _,frame in j.groupby('HAR',sort=True)]
    raw=np.concatenate(byhar)
    vals=[]
    for k in range(nboot):
        ids=rng.integers(0,len(byhar),len(byhar))
        arr=np.concatenate([byhar[i] for i in ids],axis=0)
        corr=float(np.corrcoef(arr[:,0],arr[:,1])[0,1]) if (arr[:,0].std()>0 and arr[:,1].std()>0) else np.nan
        agreement=float(np.mean(arr[:,0]*arr[:,1]>0))
        vals.append((corr,agreement))
    valid=np.array(vals,float)
    corrs=valid[np.isfinite(valid[:,0]),0]
    return dict(n_HAR_cluster_bootstraps=nboot,
                Pearson_cut_mean=float(np.corrcoef(raw[:,0],raw[:,1])[0,1]) if raw.shape[0]>1 else None,
                Pearson_95pct_bootstrap_HAR_CI=[float(x) for x in np.percentile(corrs,[2.5,97.5])] if len(corrs)>50 else None,
                mean_sign_agreement=float(np.mean(raw[:,0]*raw[:,1]>0)),
                mean_sign_agreement_95pct_HAR_bootstrap_CI=[float(x) for x in np.percentile(valid[:,1],[2.5,97.5])],
                cautions='Across just two human iPSC lines; CIs reflect HAR clusters not population of donor lines')

def process():
    mapping,msha=load_design()
    h=INPUT/'S2_HS111_Baseline_342.csv'
    baseline=pd.read_csv(h)
    expected={'HAR','cut','culture','corrected_HminusC_log2'}
    if len(baseline)!=342 or not expected.issubset(baseline.columns):
        raise ValueError('Invalid independent HS1-11 N2 baseline')
    if baseline.duplicated(['HAR','cut','culture']).any():raise ValueError('Duplicate HS1-11 baseline identities')
    rows=[];source=[];frag_rows=[];sensitivity=[]
    for rep,rg,dg in SAMPLES:
        cache=[]
        for gsm,assay in [(rg,'RNA'),(dg,'DNA')]:
            url=url_for(gsm,assay,rep)
            file=fetch(url,f'{gsm}_{assay}.tsv.gz')
            vals,stat=parse_counts(file,mapping)
            cache.append(vals)
            source.append(dict(culture=rep,GSM=gsm,assay=assay,url=url,
                               sha256=hashlib.sha256(file.read_bytes()).hexdigest(),bytes=file.stat().st_size,**stat))
        rna,dna=cache
        baseline_effect,frag=per_fragment(rna,dna,mapping,rep,threshold=1,pseudo=1)
        for (har,cut),(delta,hn,cn) in baseline_effect.items():
            rows.append(dict(HAR=har,cut=cut,culture=rep,WTc_HminusC_log2=delta,
                             WTc_human_DNApositive_barcodes=hn,WTc_chimp_DNApositive_barcodes=cn))
        frag_rows.extend(frag)
        for thr in [1,5,10,20]:
            for pc in [0.5,1,2]:
                test,_=per_fragment(rna,dna,mapping,rep,thr,pc)
                for (har,cut),(delta,hn,cn) in test.items():
                    sensitivity.append(dict(HAR=har,cut=cut,culture=rep,min_DNA_count_per_barcode=thr,
                                            aggregate_pseudocount=pc,WTc_HminusC_log2=delta,
                                            WTc_human_DNApositive_barcodes=hn,WTc_chimp_DNApositive_barcodes=cn))
        pd.DataFrame(rows).to_csv(OUT/'S1_WTc_CutByCulture_PARTIAL.csv',index=False)
        print(f'WTc N2 culture {rep}: {len(baseline_effect)} usable cut comparisons; original source barcode/insert IDs verified')
    w=pd.DataFrame(rows)
    if w.empty:raise ValueError('No WTc reporter comparisons')
    pd.DataFrame(rows).to_csv(OUT/'S1_WTc_CutByCulture.csv',index=False)
    pd.DataFrame(source).to_csv(OUT/'S2_SixGEO_SourceIntegrity.csv',index=False)
    pd.DataFrame(frag_rows).to_csv(OUT/'S3_WTc_Allele_ActivityAndDepth.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(OUT/'S4_WTc_FullDepth_PseudocountSensitivity.csv',index=False)
    ha=(baseline.groupby(['HAR','cut'],as_index=False)
           .agg(HS111_n_cultures=('culture','nunique'),HS111_mean=('corrected_HminusC_log2','mean'),
                HS111_min=('corrected_HminusC_log2','min'),HS111_max=('corrected_HminusC_log2','max')))
    wa=(w.groupby(['HAR','cut'],as_index=False)
          .agg(WTc_n_cultures=('culture','nunique'),WTc_mean=('WTc_HminusC_log2','mean'),
               WTc_min=('WTc_HminusC_log2','min'),WTc_max=('WTc_HminusC_log2','max')))
    ha['HS111_all3_same_direction']=(ha.HS111_n_cultures.eq(3)&((ha.HS111_min>0)|(ha.HS111_max<0)))
    wa['WTc_all3_same_direction']=(wa.WTc_n_cultures.eq(3)&((wa.WTc_min>0)|(wa.WTc_max<0)))
    combined=ha.merge(wa,on=['HAR','cut'],how='outer',validate='one_to_one',indicator=True)
    combined['both_lines']=combined._merge.eq('both')
    combined['same_direction_2lines']=np.where(combined.both_lines,combined.HS111_mean*combined.WTc_mean>0,np.nan)
    combined['both_lines_same_direction_all3cultures']=np.where(combined.both_lines,
        combined.same_direction_2lines.fillna(False).astype(bool)&combined.HS111_all3_same_direction.fillna(False).astype(bool)&combined.WTc_all3_same_direction.fillna(False).astype(bool),np.nan)
    combined.to_csv(OUT/'S5_IndependentLine_ConstructConcordance.csv',index=False)
    combined.loc[combined.HAR.isin(PRIORITY)].to_csv(OUT/'S6_Five_Priority_HARs.csv',index=False)
    common=combined.loc[combined.both_lines].copy()
    bs=bootstrap_between_lines(common) if len(common)>5 else None
    sens=pd.DataFrame(sensitivity)
    base=sens[(sens.min_DNA_count_per_barcode==1)&(sens.aggregate_pseudocount==1)].set_index(['HAR','cut','culture'])['WTc_HminusC_log2']
    strict=sens[(sens.min_DNA_count_per_barcode.isin([5,10,20]))&(sens.aggregate_pseudocount==1)].copy()
    strict['baseline']=strict.set_index(['HAR','cut','culture']).index.map(base)
    strict['WTc_direction_flip_vs_1read']=(strict.WTc_HminusC_log2*strict.baseline<0)
    depth=(strict.groupby('min_DNA_count_per_barcode',as_index=False)
             .agg(n_comparable=('HAR','size'),n_sign_flips=('WTc_direction_flip_vs_1read','sum')))
    depth['flip_fraction']=depth.n_sign_flips/depth.n_comparable
    depth.to_csv(OUT/'S7_WTc_DNAdepth_SensitivitySummary.csv',index=False)
    summary=dict(stage='24-final-independent-line',status='complete',
        GEO='GSE110758',human_host_iPSC_lines=['HS1-11','WTc'],NPC_stage='N2',
        human_vs_chimp_are_INSERT_alleles=True,n_independent_human_line_backgrounds=2,
        n_cultures_per_line=3,n_selected_design_barcodes=len(mapping),
        n_WTc_comparison_cuts=int(len(wa)),n_WTc_3culture_cuts=int(wa.WTc_n_cultures.eq(3).sum()),
        n_shared_cuts=int(len(common)),n_shared_HARs=int(common.HAR.nunique()),
        n_cross_line_same_sign_cuts=int(common.same_direction_2lines.sum()),
        n_cross_line_stable_3culture_same_sign_cuts=int(common.both_lines_same_direction_all3cultures.sum()),
        n_prioritized_shared_cuts=int(common.HAR.isin(PRIORITY).sum()),
        raw_GEO_source_insert_annotations_checked=True,selected_insert_annotation_mismatches=0,
        six_GEO_downloads_sha256_recorded=True,
        original_screen_HAR_interval_identity_verified=False,
        formal_donor_population_inference=False,
        independent_line_cluster_bootstrap=bs,
        cautions=['Different human iPSC cell lines; no chimp host used in this comparison',
                  'MPRA 171bp inserts do not prove TF occupancy or endogenous HAR-gene link',
                  'Only two human cell line genetic backgrounds; do not treat cultures as donors',
                  'Original Noble HAR BED and 28 unresolved mappings still unverified'])
    (OUT/'S0_IndependentLine_Status.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2));return summary

def pack(name):
    path=ROOT/name
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for f in sorted(OUT.iterdir()):
            if f.is_file():z.write(f,f.name)
    print('RESULT ZIP:',path,'size:',path.stat().st_size)
    return path

def main():
    try:
        process();return pack('HAR_Stage24_WTc_FinalValidation_Results.zip')
    except Exception as exc:
        (OUT/'S0_FailureTraceback.txt').write_text(traceback.format_exc())
        (OUT/'S0_FailureStatus.json').write_text(json.dumps({'status':'FAILED','error':str(exc),'next':'Upload diagnostics ZIP; do not infer biology'},indent=2))
        return pack('HAR_Stage24_WTc_FinalValidation_Diagnostics.zip')

if __name__=='__main__':main()
