# Stage 21C — replicated HS1-11 N2 lentiMPRA activity reanalysis
# Strict data provenance: GEO GSE110758; measured barcode RNA/DNA counts, not TF occupancy.
import re, io, csv, json, time, hashlib, zipfile, traceback, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import requests

ROOT = Path('/content') if Path('/content').exists() else Path.cwd()
CACHE = ROOT / 'HAR_Stage21C_Source'
OUT = ROOT / 'HAR_Stage21C_Results'
CACHE.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)
HAR_NAMES = ['2xHAR.116', '2xHAR.128', '2xHAR.13', '2xHAR.132', '2xHAR.135', '2xHAR.165', '2xHAR.17', '2xHAR.177', '2xHAR.180', '2xHAR.183', '2xHAR.186', '2xHAR.20', '2xHAR.213', '2xHAR.215', '2xHAR.223', '2xHAR.234', '2xHAR.235', '2xHAR.238', '2xHAR.24', '2xHAR.240', '2xHAR.256', '2xHAR.258', '2xHAR.274', '2xHAR.287', '2xHAR.3', '2xHAR.30', '2xHAR.305', '2xHAR.310', '2xHAR.313', '2xHAR.33', '2xHAR.332', '2xHAR.333', '2xHAR.341', '2xHAR.346', '2xHAR.352', '2xHAR.377', '2xHAR.385', '2xHAR.39', '2xHAR.413', '2xHAR.419', '2xHAR.435', '2xHAR.437', '2xHAR.446', '2xHAR.471', '2xHAR.487', '2xHAR.49', '2xHAR.491', '2xHAR.492', '2xHAR.509', '2xHAR.51', '2xHAR.510', '2xHAR.531', '2xHAR.536', '2xHAR.62', '2xHAR.65', '2xHAR.67', '2xHAR.70', '2xHAR.82', '2xHAR.87', '2xHAR.99', 'HACNS10', 'HACNS105', 'HACNS125', 'HACNS129', 'HACNS149', 'HACNS15', 'HACNS155', 'HACNS174', 'HACNS188', 'HACNS196', 'HACNS198', 'HACNS205', 'HACNS238', 'HACNS239', 'HACNS255', 'HACNS259', 'HACNS278', 'HACNS295', 'HACNS325', 'HACNS33', 'HACNS332', 'HACNS336', 'HACNS368', 'HACNS377', 'HACNS38', 'HACNS383', 'HACNS420', 'HACNS423', 'HACNS426', 'HACNS439', 'HACNS442', 'HACNS454', 'HACNS457', 'HACNS464', 'HACNS480', 'HACNS497', 'HACNS50', 'HACNS52', 'HACNS522', 'HACNS544', 'HACNS549', 'HACNS551', 'HACNS563', 'HACNS567', 'HACNS579', 'HACNS583', 'HACNS584', 'HACNS587', 'HACNS588', 'HACNS607', 'HACNS635', 'HACNS636', 'HACNS643', 'HACNS656', 'HACNS659', 'HACNS67', 'HACNS671', 'HACNS674', 'HACNS682', 'HACNS684', 'HACNS686', 'HACNS694', 'HACNS698', 'HACNS703', 'HACNS71', 'HACNS712', 'HACNS719', 'HACNS722', 'HACNS734', 'HACNS735', 'HACNS754', 'HACNS76', 'HACNS761', 'HACNS762', 'HACNS763', 'HACNS772', 'HACNS780', 'HACNS786', 'HACNS798', 'HACNS805', 'HACNS839', 'HACNS851', 'HACNS858', 'HACNS860', 'HACNS865', 'HACNS870', 'HACNS880', 'HACNS901', 'HACNS909', 'HACNS915', 'HACNS924', 'HACNS925', 'HACNS927', 'HACNS938', 'HACNS96', 'HACNS962', 'HACNS964', 'HACNS968', 'HACNS975', 'HACNS976', 'HACNS99', 'HAR106', 'HAR108', 'HAR116', 'HAR149', 'HAR154', 'HAR156', 'HAR160', 'HAR165', 'HAR169', 'HAR17', 'HAR172', 'HAR181', 'HAR184', 'HAR190', 'HAR198', 'HAR202', 'HAR25', 'HAR67', 'HAR93']
MOTIF_CANDIDATES = {'2xHAR.234': {'Stage20B_3mer_p': 0.3652, 'Stage20B_BH_q': 1.0}, '2xHAR.3': {'Stage20B_3mer_p': 0.4434, 'Stage20B_BH_q': 1.0}, '2xHAR.310': {'Stage20B_3mer_p': 0.1526, 'Stage20B_BH_q': 1.0}, '2xHAR.377': {'Stage20B_3mer_p': 0.4416, 'Stage20B_BH_q': 1.0}, '2xHAR.70': {'Stage20B_3mer_p': 0.1188, 'Stage20B_BH_q': 1.0}}
EXPECTED_DESIGN_SHA256 = 'c0552a2fc675f043017252a8a341d6e06cd91cedefdc7b65a0186118901a0d3f'
DESIGN_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110758/suppl/GSE110758_HAR_design.tsv.gz'
# GEO sample IDs, same iPSC line and same N2 stage. These are three CULTURE replicates, not 3 donors.
SAMPLES = [
    (1, 'GSM3016032', 'GSM3016050'),
    (2, 'GSM3016033', 'GSM3016051'),
    (3, 'GSM3016034', 'GSM3016052'),
]
MANIFEST = []

def official_url(gsm, assay, rep):
    # e.g. GSM3016032_human-HS1-11-N2-1-RNA.tsv.gz
    fn = f'{gsm}_human-HS1-11-N2-{rep}-{assay}.tsv.gz'
    base = gsm[:7]+'nnn'  # GSM3016nnn
    return f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{base}/{gsm}/suppl/{fn}'

def download_checked(url, destination, max_bytes=100_000_000):
    destination=Path(destination)
    if destination.is_file() and destination.stat().st_size > 200:
        with destination.open('rb') as f:
            if f.read(2)==b'\x1f\x8b':
                print('Cached:',destination.name)
                return destination
    error=None
    for attempt in range(1,4):
        try:
            print('Fetching:',url, 'attempt',attempt)
            with requests.get(url,stream=True,timeout=(25,120),headers={'User-Agent':'HAR-Stage21C-Research/1.0'}) as r:
                r.raise_for_status()
                total=0
                with open(destination,'wb') as f:
                    for chunk in r.iter_content(1024*1024):
                        if not chunk: continue
                        f.write(chunk); total+=len(chunk)
                        if total>max_bytes: raise ValueError('File exceeded configured safe size')
            if total<200: raise ValueError('Downloaded content too small')
            with destination.open('rb') as f:
                if f.read(2)!=b'\x1f\x8b': raise ValueError('Not a GZIP file')
            with gzip.open(destination,'rb') as f:
                if len(f.read(100))==0: raise ValueError('GZIP content empty')
            return destination
        except Exception as e:
            error=e; print('Retryable download error:',type(e).__name__,str(e)[:250]);
            if destination.exists(): destination.unlink()
            time.sleep(2)
    raise RuntimeError(f'Download failed {url}: {error}')

def manifest_file(p,source):
    p=Path(p); hash=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): hash.update(chunk)
    rec={'filename':p.name,'download_url':source,'size_bytes':p.stat().st_size,'sha256':hash.hexdigest()}
    MANIFEST.append(rec)
    return rec

def barcode_counts(path,assay):
    # GEO processed TSV should contain barcode and molecule/read count; infer column labels conservatively.
    d=pd.read_csv(path,sep='\t',compression='gzip',dtype=str,keep_default_na=False,on_bad_lines='error')
    if d.empty: raise ValueError(f'Empty {assay} table')
    cols=list(d.columns)
    bc_candidates=[c for c in cols if 'barcode' in c.lower() or c.strip().lower() in {'bc','tag'}]
    if not bc_candidates:
        bc_candidates=[c for c in cols if d[c].astype(str).head(400).str.fullmatch('[ACGT]{15}',case=False).mean()>.90]
    if len(bc_candidates)!=1:
        raise ValueError(f'Could not uniquely identify barcode column from {cols}. First rows: {d.head(3).to_dict("records")}')
    bc=bc_candidates[0]
    numeric_cols=[]
    for c in cols:
        if c==bc: continue
        v=pd.to_numeric(d[c],errors='coerce')
        if v.notna().mean()>.98 and (v.dropna()>=0).all(): numeric_cols.append(c)
    if len(numeric_cols)!=1:
        raise ValueError(f'Expected one numeric count column, found {numeric_cols} from {cols}')
    val=numeric_cols[0]
    x=d[[bc,val]].copy(); x.columns=['Barcode','count']
    x['Barcode']=x['Barcode'].str.upper().str.strip()
    x['count']=pd.to_numeric(x['count'],errors='raise')
    x=x[x['Barcode'].str.fullmatch('[ACGT]{15}') & x['count'].gt(0)].copy()
    x=x.groupby('Barcode',as_index=False,sort=False)['count'].sum()
    if len(x)<1000 or x['count'].sum()<=0:
        raise ValueError(f'{assay}: unexpectedly few valid barcodes or zero total counts')
    print(assay, 'rows',len(d),'valid barcodes',len(x),'sum',int(x['count'].sum()),'columns',cols)
    return x

def process():
    manifest_file(download_checked(DESIGN_URL,CACHE/'GSE110758_HAR_design.tsv.gz',max_bytes=25_000_000),DESIGN_URL)
    if MANIFEST[-1]['sha256']!=EXPECTED_DESIGN_SHA256:
        raise ValueError('GEO design checksum changed relative to Stage21B. Review source version before computing.')
    df=pd.read_csv(CACHE/'GSE110758_HAR_design.tsv.gz',sep='\t',compression='gzip',dtype=str,keep_default_na=False,on_bad_lines='error')
    assert set(['#Name','Barcode','Sequence']).issubset(df.columns), list(df.columns)
    assert len(df)==244000, f'Design row count changed: {len(df)}'
    # Only the exact 180 specified input labels; no fuzzy name matches.
    names_set=set(HAR_NAMES)
    df['HAR']=df['#Name'].str.split('__species=',regex=False).str[0]
    df=df[df['HAR'].isin(names_set)].copy()
    df['allele']=df['#Name'].str.extract(r'__species=([^_]+)',expand=False).str.lower()
    df['cut']=df['#Name'].str.extract(r'__cut=(.*?)(?:__|$)',expand=False).str.lower()
    df['insert']=df['#Name'].str.replace(r':[0-9]+$','',regex=True)
    df['Barcode']=df['Barcode'].str.upper().str.strip()
    if df['allele'].isna().any() or df['cut'].isna().any():
        raise ValueError('Unexpected design identifier: some species/cut are missing; inspect GSE metadata')
    df=df[df['allele'].isin(['human','chimp']) & df['Barcode'].str.fullmatch('[ACGT]{15}')].copy()
    coll=df.groupby('Barcode')['insert'].nunique()
    bad=set(coll[coll>1].index); print('Ambiguous multi-insert barcodes excluded:',len(bad))
    df=df[~df['Barcode'].isin(bad)].drop_duplicates(subset=['Barcode','insert'])
    assert df['Barcode'].is_unique, 'Expected one distinct insert per retained barcode'
    print('Design selected:',len(df),'barcodes;',df['HAR'].nunique(),'HAR names;',df['insert'].nunique(),'inserts')
    design_summary=(df.groupby(['HAR','allele','cut'])['Barcode'].nunique().reset_index(name='design_n_barcodes'))
    design_summary.to_csv(OUT/'S1_75HAR_Design_Allele_Cuts.csv',index=False)
    df[['HAR','allele','cut','insert','Barcode']].to_csv(OUT/'S1B_Assayed_Barcode_Map.csv',index=False)
    rstats=[]; insert_all=[]; pair_all=[]
    for rep,rna_gsm,dna_gsm in SAMPLES:
        rna_url=official_url(rna_gsm,'RNA',rep)
        dna_url=official_url(dna_gsm,'DNA',rep)
        rna_path=download_checked(rna_url,CACHE/f'{rna_gsm}_RNA.tsv.gz')
        dna_path=download_checked(dna_url,CACHE/f'{dna_gsm}_DNA.tsv.gz')
        manifest_file(rna_path,rna_url);manifest_file(dna_path,dna_url)
        rn=barcode_counts(rna_path,'RNA'); dn=barcode_counts(dna_path,'DNA')
        total_rna=float(rn['count'].sum());total_dna=float(dn['count'].sum())
        paired=rn.merge(dn,on='Barcode',how='inner',suffixes=('_rna','_dna'),validate='one_to_one')
        paired=paired.merge(df[['HAR','allele','cut','insert','Barcode']],on='Barcode',how='inner',validate='one_to_one')
        paired['rna_cpm']=1e6*paired['count_rna']/total_rna
        paired['dna_cpm']=1e6*paired['count_dna']/total_dna
        g=paired.groupby(['HAR','allele','cut','insert'],as_index=False).agg(
            n_paired_barcodes=('Barcode','nunique'),rna_cpm=('rna_cpm','sum'),dna_cpm=('dna_cpm','sum'))
        g['activity_log2RNA_over_DNA']=np.log2(g['rna_cpm']/g['dna_cpm'])
        g['replicate']=rep;g['cell_line']='HS1-11';g['cell_stage']='N2'
        # Prevent single-barcode artifact from appearing as usable activity estimate.
        g['passes_min10_paired_barcodes']=g['n_paired_barcodes'].ge(10)
        insert_all.append(g)
        # Aggregate at fragment/cut pairing level. Multiple insert designs at same cut -> disallow averaging.
        ok=g[g['passes_min10_paired_barcodes']].copy()
        uniq=ok.groupby(['HAR','allele','cut']).size().rename('n_insert_designs').reset_index()
        ok=ok.merge(uniq,on=['HAR','allele','cut'])
        ok=ok[ok['n_insert_designs'].eq(1)].copy()
        hu=ok[ok.allele.eq('human')].drop(columns=['allele']).rename(columns={
            'activity_log2RNA_over_DNA':'human_log2_activity','n_paired_barcodes':'human_paired_barcodes'})
        ch=ok[ok.allele.eq('chimp')].drop(columns=['allele']).rename(columns={
            'activity_log2RNA_over_DNA':'chimp_log2_activity','n_paired_barcodes':'chimp_paired_barcodes'})
        x=hu[['HAR','cut','replicate','human_log2_activity','human_paired_barcodes']].merge(
           ch[['HAR','cut','replicate','chimp_log2_activity','chimp_paired_barcodes']],on=['HAR','cut','replicate'],how='inner',validate='one_to_one')
        x['human_minus_chimp_log2_activity']=x['human_log2_activity']-x['chimp_log2_activity']
        pair_all.append(x)
        rstats.append({'replicate':rep,'RNA_GSM':rna_gsm,'DNA_GSM':dna_gsm,'total_RNA_counts':total_rna,'total_DNA_counts':total_dna,
                       'RNA_barcodes':len(rn),'DNA_barcodes':len(dn),'RNA_DNA_shared_barcodes':len(rn.merge(dn,on='Barcode')),
                       'design_mapped_paired_barcodes':len(paired),'HARs_with_pairable_human_chimp_cut':x.HAR.nunique(),
                       'paired_fragment_cuts':len(x)})
        print('Replicate',rep,':',len(x),'paired fragment/cut effects for',x.HAR.nunique(),'HARs')
        pd.concat(insert_all,ignore_index=True).to_csv(OUT/'S2_Insert_Activity_By_Replicate_PARTIAL.csv',index=False)
        pd.concat(pair_all,ignore_index=True).to_csv(OUT/'S3_Paired_Fragment_Activity_By_Replicate_PARTIAL.csv',index=False)
    all_insert=pd.concat(insert_all,ignore_index=True)
    all_pair=pd.concat(pair_all,ignore_index=True)
    all_insert.to_csv(OUT/'S2_Insert_Activity_By_Replicate.csv',index=False)
    all_pair.to_csv(OUT/'S3_Paired_Fragment_Activity_By_Replicate.csv',index=False)
    pd.DataFrame(rstats).to_csv(OUT/'S4_Replicate_QC.csv',index=False)
    # Balanced 3-replicate paired observations and unique species cuts only.
    agg=(all_pair.groupby(['HAR','cut'],as_index=False).agg(
        n_replicates=('replicate','nunique'),mean_HminusC_log2=('human_minus_chimp_log2_activity','mean'),
        min_HminusC_log2=('human_minus_chimp_log2_activity','min'),max_HminusC_log2=('human_minus_chimp_log2_activity','max')))
    agg['all_three_replicates']=agg['n_replicates'].eq(3)
    agg['direction_same_all_three']=agg.apply(lambda row: bool((row['min_HminusC_log2']>0) or (row['max_HminusC_log2']<0)) if row['all_three_replicates'] else False,axis=1)
    agg['Stage20B_3mer_p']=agg.HAR.map({k:v['Stage20B_3mer_p'] for k,v in MOTIF_CANDIDATES.items()})
    agg['Stage20B_BH_q']=agg.HAR.map({k:v['Stage20B_BH_q'] for k,v in MOTIF_CANDIDATES.items()})
    agg['analysis_status']='EXPLORATORY_THREE_CULTURE_REPLICATES_NO_FORMAL_SIGNIFICANCE'
    agg.to_csv(OUT/'S5_75HAR_Replicate_Summary_By_Cut.csv',index=False)
    agg[agg.HAR.isin(MOTIF_CANDIDATES)].to_csv(OUT/'S6_Five_Motif_Candidate_MPRA_Summary.csv',index=False)
    pd.DataFrame(MANIFEST).to_csv(OUT/'S7_Source_Manifest.csv',index=False)
    status={'stage':'21C','status':'complete','source':'GSE110758','cell_line':'HS1-11','stage_NPC':'N2',
            'replicates_3_independent_cultures_in_same_line':True,'two_DNA_species_are_INSERT_alleles_not_cell_genetic_backgrounds':True,
            'design_SHA256':EXPECTED_DESIGN_SHA256,'original_Noble_HAR_coordinate_crosswalk_confirmed':False,
            'activity_has_formal_FDR':False,'TF_occupancy_confirmed':False,
            'n_HAR_exact_design_names':int(df.HAR.nunique()),'n_HAR_with_at_least_one_3rep_paired_cut':int(agg[agg.all_three_replicates].HAR.nunique()),
            'n_motif_candidates_with_3rep_paired_cut':int(agg[(agg.all_three_replicates)&(agg.HAR.isin(MOTIF_CANDIDATES))].HAR.nunique())}
    (OUT/'S0_Status_and_Limitations.json').write_text(json.dumps(status,indent=2))
    print(json.dumps(status,indent=2))

def zip_outputs(name):
    target=ROOT/name
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(OUT.glob('*')):
            if p.is_file():z.write(p,p.name)
    print('ZIP available:',target,'bytes',target.stat().st_size)
    print('Download manually in Colab left-side Files panel to avoid subprocess kernel crashes.')
    return target

try:
    process()
    zip_outputs('HAR_Stage21C_MPRA_HS111_N2_Results.zip')
except Exception as e:
    (OUT/'S0_Run_Error_Traceback.txt').write_text(traceback.format_exc())
    (OUT/'S0_Run_Error_Status.json').write_text(json.dumps({'status':'failed','stage':'21C','error':str(e)},indent=2))
    zip_outputs('HAR_Stage21C_Diagnostics_and_Partial.zip')
    print('Detailed error:',traceback.format_exc())
