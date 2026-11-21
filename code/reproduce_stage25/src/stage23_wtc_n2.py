#!/usr/bin/env python3
"""GSE110758 WTc N2: independent human iPSC line lentiMPRA descriptive replication.
Only human and chimp *inserts* (alleles) are compared within each human line.
This does not validate gene targets or TF binding, nor replicate Noble perturbation.

Run in Colab via bundled notebook, or locally from repo root (requires internet).
"""
import csv, gzip, hashlib, json, math, re, time, traceback, zipfile
from pathlib import Path
from collections import defaultdict, Counter
import requests
import pandas as pd
import numpy as np

ROOT=Path('/content') if Path('/content').exists() else Path(__file__).resolve().parents[1]
IN=ROOT/'HAR_Stage23_WTc_Embedded';OUT=ROOT/'HAR_Stage23_WTc_N2_Outputs'
CACHE=ROOT/'HAR_Stage23_WTc_N2_Cache'
OUT.mkdir(exist_ok=True);CACHE.mkdir(exist_ok=True)
SAMPLES=[(1,'GSM3016026','GSM3016044'),(2,'GSM3016027','GSM3016045'),(3,'GSM3016028','GSM3016046')]
PRIORITY={'2xHAR.70','2xHAR.310','2xHAR.234','2xHAR.377','2xHAR.3'}
EXPECTED_MAP_SHA='6f30451ac623ed4948c8649aa3e9113f0d244b5d70db1f37fbb35314cf92cdeb'

def url_for(gsm,assay,rep):
    return f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:7]}nnn/{gsm}/suppl/{gsm}_human-WTc-N2-{rep}-{assay}.tsv.gz'

def fetch(u,name):
    path=CACHE/name
    if path.exists() and path.stat().st_size>1000:
        with path.open('rb') as f:
            if f.read(2)==b'\x1f\x8b':print('Cached',name);return path
    err=''
    for attempt in range(1,4):
        try:
            print('GEO',name,'attempt',attempt)
            with requests.get(u,stream=True,timeout=(25,120),headers={'User-Agent':'HAR-assay-aware-benchmark/1.0'}) as r:
                r.raise_for_status();n=0
                with path.open('wb') as f:
                    for block in r.iter_content(1024*1024):
                        if block:
                            f.write(block);n+=len(block)
                            if n>80_000_000:raise ValueError('Unusually large sample file')
            if n<200 or path.read_bytes()[:2]!=b'\x1f\x8b':raise ValueError('Not a plausible GEO GZIP')
            with gzip.open(path,'rb') as f:
                if len(f.read(150))==0:raise ValueError('Empty GZIP')
            return path
        except Exception as e:
            err=f'{type(e).__name__}: {e}';print('Retry',err[:200]);path.unlink(missing_ok=True);time.sleep(2)
    raise RuntimeError(f'{name}: {err}')

def parse_counts(path):
    """Original GEO sample files have THREE HEADERLESS columns: barcode/count/insert ID."""
    values=Counter();nrows=0;reads=0
    with gzip.open(path,'rt',encoding='utf-8') as fh:
        for i,row in enumerate(csv.reader(fh,delimiter='\t'),1):
            if len(row)<3:raise ValueError(f'{path.name} malformed TSV row {i}')
            bc=row[0].strip().upper()
            if re.fullmatch(r'[ACGT]{15}',bc) is None:raise ValueError(f'Invalid barcode at row {i}: {bc!r}')
            n=int(row[1]);
            if n<0:raise ValueError('Negative count')
            values[bc]+=n;nrows+=1;reads+=n
    if nrows<1000 or reads<=0:raise ValueError('Empty or truncated count source '+str(path))
    return values,nrows,reads

def load_mapping():
    path=IN/'S1B_Assayed_Barcode_Map.csv';
    data=path.read_bytes();h=hashlib.sha256(data).hexdigest()
    if h!=EXPECTED_MAP_SHA: raise ValueError('Embedded Stage21C barcode map checksum mismatch')
    mapping={};groups=defaultdict(set)
    for x in csv.DictReader(data.decode('utf-8').splitlines()):
        bc=x['Barcode'];key=(x['HAR'],x['allele'],x['cut'],x['insert'])
        if bc in mapping and mapping[bc]!=key:raise ValueError('Ambiguous design barcode')
        mapping[bc]=key;groups[key[:3]].add(x['insert'])
    return mapping,groups,h

def process():
    mapping,designgroups,msha=load_mapping()
    hpath=IN/'S2_HS111_Baseline_342.csv'
    ref=pd.read_csv(hpath)
    assert len(ref)==342 and set(['HAR','cut','culture','corrected_HminusC_log2']).issubset(ref.columns)
    rows=[];qc=[];urls=[];insert_stats=[]
    for rep,rg,dg in SAMPLES:
        sources=[]
        for gsm,kind in [(rg,'RNA'),(dg,'DNA')]:
            u=url_for(gsm,kind,rep);f=fetch(u,f'{gsm}_{kind}.tsv.gz')
            vals,n,tot=parse_counts(f)
            sources.append((vals,n,tot))
            urls.append(dict(culture=rep,assay=kind,GSM=gsm,url=u,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),bytes=f.stat().st_size,rows=n,total_counts=tot))
        RNA,nr,tr=sources[0];DNA,nd,td=sources[1]
        counts=defaultdict(lambda:[0,0,0]);
        for bc,key in mapping.items():
            nc=DNA.get(bc,0)
            if nc<1:continue
            counts[key][0]+=RNA.get(bc,0);counts[key][1]+=nc;counts[key][2]+=1
        alleles=[]
        for (h,allele,cut,insert),(rna,dna,barcodes) in counts.items():
            # Aggregate pseudocount=1 on summed RNA and DNA; full raw includes zero-RNA.
            if barcodes<10:continue
            effect=math.log2((rna+1)/(dna+1))
            alleles.append((h,allele,cut,insert,rna,dna,barcodes,effect))
            insert_stats.append(dict(HAR=h,allele=allele,cut=cut,insert=insert,culture=rep,
                                     RNA_sum=rna,DNA_sum=dna,DNApositive_barcodes=barcodes,
                                     log2_RNA_DNA_pseudocount1=effect))
        eligible=defaultdict(dict)
        for h,allele,cut,insert,rna,dna,n,eff in alleles:
            k=(h,cut);eligible[k].setdefault(allele,[]).append((eff,insert,n))
        npaired=0
        for (h,cut),allele in eligible.items():
            if len(allele.get('human',[]))!=1 or len(allele.get('chimp',[]))!=1:continue
            a=allele['human'][0];b=allele['chimp'][0]
            rows.append(dict(HAR=h,cut=cut,culture=rep,WTc_HminusC_log2=a[0]-b[0],
                             WTc_human_DNApositive_barcodes=a[2],WTc_chimp_DNApositive_barcodes=b[2]))
            npaired+=1
        qc.append(dict(culture=rep,RNA_GSM=rg,DNA_GSM=dg,RNA_total_counts=tr,DNA_total_counts=td,
                       raw_RNA_rows=nr,raw_DNA_rows=nd,paired_HAR_cuts=npaired,
                       compared_insert_alleles='human and chimp constructs in human WTc background'))
        print('WTc N2 culture',rep,'paired cuts',npaired)
        pd.DataFrame(rows).to_csv(OUT/'S1_WTc_N2_CutEffects_PARTIAL.csv',index=False)
    w=pd.DataFrame(rows)
    if w.empty: raise ValueError('No usable same-cut human-chimp comparisons in WTc')
    w.to_csv(OUT/'S1_WTc_N2_CutEffects.csv',index=False)
    pd.DataFrame(qc).to_csv(OUT/'S2_WTc_Sample_QC.csv',index=False)
    pd.DataFrame(urls).to_csv(OUT/'S3_WTc_GEO_Source_Manifest.csv',index=False)
    pd.DataFrame(insert_stats).to_csv(OUT/'S4_WTc_Insert_ReadDepth.csv',index=False)
    a=w.groupby(['HAR','cut'],as_index=False).agg(WTc_n_cultures=('culture','nunique'),WTc_mean=('WTc_HminusC_log2','mean'),
                                                    WTc_min=('WTc_HminusC_log2','min'),WTc_max=('WTc_HminusC_log2','max'))
    a['WTc_direction_consistent_all3']=(a.WTc_n_cultures.eq(3)&((a.WTc_min>0)|(a.WTc_max<0)))
    h=ref.groupby(['HAR','cut'],as_index=False).agg(HS111_n_cultures=('culture','nunique'),HS111_mean=('corrected_HminusC_log2','mean'),
                                                   HS111_min=('corrected_HminusC_log2','min'),HS111_max=('corrected_HminusC_log2','max'))
    h['HS111_direction_consistent_all3']=(h.HS111_n_cultures.eq(3)&((h.HS111_min>0)|(h.HS111_max<0)))
    j=h.merge(a,on=['HAR','cut'],how='outer',indicator=True,validate='one_to_one')
    j['two_human_lines_same_sign_of_mean']=np.where(j._merge.eq('both'),j.HS111_mean*j.WTc_mean>0,np.nan)
    j['independent_line_status']=np.where(j._merge.eq('both'),'IN_BOTH_LINES','AVAILABLE_IN_ONE_LINE_ONLY')
    j.to_csv(OUT/'S5_WTc_vs_HS111_N2_Cut_Comparison.csv',index=False)
    j.loc[j.HAR.isin(PRIORITY)].to_csv(OUT/'S6_Five_Priority_HARs_IndependentLine.csv',index=False)
    common=j[j._merge.eq('both')].copy()
    summary=dict(stage='23-WTc-N2',status='complete',independent_human_iPSC_line='WTc',
                 assay='lentiMPRA GSE110758 same HAR constructs measured in human WTc N2 neural progenitors',
                 cell_line_background_not_chimp=True,number_alleles='human and chimp insert sequences',
                 cultures_per_WTc_line=3,n_WTc_fragment_cuts=int(a.shape[0]),n_WTc_3culture_cuts=int(a.WTc_n_cultures.eq(3).sum()),
                 n_shared_cuts_with_HS111=int(len(common)),n_shared_HARs=int(common.HAR.nunique()),
                 n_shared_same_sign_mean=int(common.two_human_lines_same_sign_of_mean.sum()),
                 shared_fraction_same_sign_mean=(float(common.two_human_lines_same_sign_of_mean.mean()) if len(common) else None),
                 n_five_priority_cuts_in_both=int(common.HAR.isin(PRIORITY).sum()),
                 mapping_SHA256=msha,primary_HAR_coordinate_validation=False,
                 formal_multi_donor_inference=False,
                 cautions=['two human cell lines are different individuals but not enough for population-level inference',
                 'original HAR screen remains in one genetic background','RNA/DNA pseudocount=1 and DNApositive>=10 per allele cut',
                 'no TF binding or endogenous HAR-gene causality determined'])
    if len(common)>5:
        from scipy.stats import spearmanr,pearsonr
        summary['interline_Pearson_cut_mean']=float(pearsonr(common.WTc_mean,common.HS111_mean).statistic)
        summary['interline_Spearman_cut_mean']=float(spearmanr(common.WTc_mean,common.HS111_mean).statistic)
    (OUT/'S0_WTc_vs_HS111_Source_Status.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2));return summary

def pack(name):
    p=ROOT/name
    with zipfile.ZipFile(p,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for f in sorted(OUT.iterdir()):
            if f.is_file():z.write(f,f.name)
    print('ZIP',str(p),p.stat().st_size)
    return p

def main():
    try:
        process();return pack('HAR_Stage23_WTc_N2_IndependentLine_Results.zip')
    except Exception as e:
        (OUT/'S0_Run_Error_Traceback.txt').write_text(traceback.format_exc())
        (OUT/'S0_Run_Error_Status.json').write_text(json.dumps({'status':'ERROR','error':str(e)},indent=2))
        return pack('HAR_Stage23_WTc_N2_Diagnostics.zip')
if __name__=='__main__':main()
