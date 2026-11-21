# Stage 21D completion workflow: six original GEO RNA/DNA files, three cultures, no Google Drive.
# Uses only publicly available sources and the embedded, checksum-verified Stage 21C barcode-to-insert map.
from pathlib import Path
from collections import Counter, defaultdict
import csv, gzip, json, io, math, zipfile, hashlib, re, time, traceback
import requests

ROOT = Path('/content') if Path('/content').exists() else Path.cwd()
CACHE = ROOT/'HAR_Stage21D_GEO_Raw'
OUT = ROOT/'HAR_Stage21D_FullRaw_Output'
CACHE.mkdir(exist_ok=True);OUT.mkdir(exist_ok=True)
BASE = ROOT/'HAR_Stage21D_Baseline_Embedded.zip'
SOURCE_MANIFEST = [{'filename': 'GSE110758_HAR_design.tsv.gz', 'download_url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110758/suppl/GSE110758_HAR_design.tsv.gz', 'size_bytes': '3474773', 'sha256': 'c0552a2fc675f043017252a8a341d6e06cd91cedefdc7b65a0186118901a0d3f'}, {'filename': 'GSM3016032_RNA.tsv.gz', 'download_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3016nnn/GSM3016032/suppl/GSM3016032_human-HS1-11-N2-1-RNA.tsv.gz', 'size_bytes': '5142648', 'sha256': 'fd2f8a64151ab00817f89071e3992308079acde79d48300a3075032c8b08eb8b'}, {'filename': 'GSM3016050_DNA.tsv.gz', 'download_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3016nnn/GSM3016050/suppl/GSM3016050_human-HS1-11-N2-1-DNA.tsv.gz', 'size_bytes': '4977755', 'sha256': '8587175f2fd30e8512f2a2733a3531bb2c683fb7922058ca8dc0a92c10049c54'}, {'filename': 'GSM3016033_RNA.tsv.gz', 'download_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3016nnn/GSM3016033/suppl/GSM3016033_human-HS1-11-N2-2-RNA.tsv.gz', 'size_bytes': '5096019', 'sha256': '9e00e5883cbe4c20be95b4f142f3498df5e20fef832d9c47203031ba95c5c16e'}, {'filename': 'GSM3016051_DNA.tsv.gz', 'download_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3016nnn/GSM3016051/suppl/GSM3016051_human-HS1-11-N2-2-DNA.tsv.gz', 'size_bytes': '4972599', 'sha256': 'b6f943a49b5ff3d28f0bd288db3f4194a646482a4db6690bff90a2c990041802'}, {'filename': 'GSM3016034_RNA.tsv.gz', 'download_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3016nnn/GSM3016034/suppl/GSM3016034_human-HS1-11-N2-3-RNA.tsv.gz', 'size_bytes': '5036064', 'sha256': 'f1d066b7456161b230e6b72625af40d43e29042401ff2a655f12f6d35ca9ab9b'}, {'filename': 'GSM3016052_DNA.tsv.gz', 'download_url': 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3016nnn/GSM3016052/suppl/GSM3016052_human-HS1-11-N2-3-DNA.tsv.gz', 'size_bytes': '4908953', 'sha256': '49b6229b9b4ca5cd357cc094aedf46a95d6ed9e5f2f18184961afd679dbaed44'}]
PAIRS = [(1,'GSM3016032_RNA.tsv.gz','GSM3016050_DNA.tsv.gz'),
         (2,'GSM3016033_RNA.tsv.gz','GSM3016051_DNA.tsv.gz'),
         (3,'GSM3016034_RNA.tsv.gz','GSM3016052_DNA.tsv.gz')]
CAND = {'2xHAR.70','2xHAR.310','2xHAR.234','2xHAR.377','2xHAR.3'}

def sha256_file(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for buf in iter(lambda:f.read(1024*1024),b''):h.update(buf)
    return h.hexdigest()

def get_source(item):
    target=CACHE/item['filename']
    if target.exists() and sha256_file(target)==item['sha256']:
        print('Verified cached:',target.name);return target
    for k in range(3):
        try:
            print('Downloading',item['filename'],'attempt',k+1)
            with requests.get(item['download_url'],timeout=(25,120),stream=True) as r:
                r.raise_for_status()
                with target.open('wb') as f:
                    total=0
                    for buf in r.iter_content(1<<20):
                        if not buf:continue
                        f.write(buf);total+=len(buf)
                        if total>40_000_000:raise ValueError('Unexpected download size')
            if target.stat().st_size!=int(item['size_bytes']) or sha256_file(target)!=item['sha256']:
                raise ValueError('GEO source changed / checksum mismatch')
            return target
        except Exception as e:
            print('Retry error',str(e)[:180]);target.unlink(missing_ok=True);time.sleep(2)
    raise RuntimeError('Failed source '+item['filename'])

def parse_count_tsv(path):
    # CRITICAL: GEO files have no header. Capture row zero as data.
    counts=Counter();rows=0;sumcounts=0
    with gzip.open(path,'rt') as f:
        for row in csv.reader(f,delimiter='\t'):
            if len(row)<3:raise ValueError('Malformed count row '+repr(row))
            bc=row[0].strip().upper()
            if re.fullmatch('[ACGT]{15}',bc) is None:raise ValueError('Invalid barcode '+repr(bc))
            v=int(row[1]);assert v>=0
            rows+=1;sumcounts+=v;counts[bc]+=v
    return counts,rows,sumcounts

def writecsv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def workflow():
    assert BASE.exists(),'Embedded baseline missing: rerun notebook from first cell'
    with zipfile.ZipFile(BASE) as z:
        mp=list(csv.DictReader(io.StringIO(z.read('S1B_Assayed_Barcode_Map.csv').decode())))
        original=list(csv.DictReader(io.StringIO(z.read('S3_Paired_Fragment_Activity_By_Replicate.csv').decode())))
    assert len(mp)==25500 and len(original)==342
    bcs_bygrp=defaultdict(list)
    for row in mp:bcs_bygrp[(row['HAR'],row['cut'],row['allele'])].append(row['Barcode'])
    orig={(r['HAR'],r['cut'],int(r['replicate'])):float(r['human_minus_chimp_log2_activity']) for r in original}
    source_rows=[];output=[];max_error=0.
    for rep,rname,dname in PAIRS:
        ritem=next(x for x in SOURCE_MANIFEST if x['filename']==rname)
        ditem=next(x for x in SOURCE_MANIFEST if x['filename']==dname)
        rp=get_source(ritem);dp=get_source(ditem)
        r,Rrows,Rtot=parse_count_tsv(rp)
        d,Drows,Dtot=parse_count_tsv(dp)
        source_rows += [{'replicate':rep,'assay':'RNA','filename':rname,'rows':Rrows,'library_counts':Rtot,'sha256':sha256_file(rp)},
                        {'replicate':rep,'assay':'DNA','filename':dname,'rows':Drows,'library_counts':Dtot,'sha256':sha256_file(dp)}]
        # Recompute Stage21C positive-positive original estimate with proper headerless parsing.
        base={}
        for key,bcs in bcs_bygrp.items():
            good=[bc for bc in bcs if r.get(bc,0)>0 and d.get(bc,0)>0]
            if len(good)>=10:
                base[key]=math.log2(sum(r[bc] for bc in good)/sum(d[bc] for bc in good))
        for har,cut,rr in sorted({(h,c,rr) for h,c,rr in orig if rr==rep}):
            k1=(har,cut,'human');k2=(har,cut,'chimp')
            new=base[k1]-base[k2]
            max_error=max(max_error,abs(new-orig[(har,cut,rep)]))
        print('Culture',rep,'raw RNA/DNA rows',Rrows,Drows,'21C recompute max error so far',max_error)
        if max_error>1e-9:raise ValueError('Stage21C paired-contrast mismatch: review mapping before proceeding')
        for depth in (1,5,10,20):
            for pseudocount in (0.5,1,2):
                by_allele={}
                for (har,cut,alle),bcs in bcs_bygrp.items():
                    good=[bc for bc in bcs if d.get(bc,0)>=depth]
                    if len(good)<10:continue
                    sumr=sum(r.get(bc,0) for bc in good)
                    sumd=sum(d[bc] for bc in good)
                    zero=sum(r.get(bc,0)==0 for bc in good)
                    by_allele[(har,cut,alle)]=(math.log2((sumr+pseudocount)/(sumd+pseudocount)),len(good),zero)
                for (har,cut,rr),old in sorted(orig.items()):
                    if rr!=rep:continue
                    h=by_allele.get((har,cut,'human'));c=by_allele.get((har,cut,'chimp'))
                    if h is None or c is None:continue
                    effect=h[0]-c[0]
                    output.append({'HAR':har,'cut':cut,'culture':rep,'min_DNA_count_per_barcode':depth,
                                   'aggregate_pseudocount':pseudocount,'DNApositive_H_barcodes':h[1],
                                   'DNApositive_C_barcodes':c[1],'RNAzero_H_barcodes':h[2],
                                   'RNAzero_C_barcodes':c[2],
                                   'corrected_HminusC_log2':effect,'Stage21C_positive_only_log2':old,
                                   'delta_from_Stage21C':effect-old})
        writecsv(OUT/'S2_All_Cut_FullRaw_Sensitivity_PARTIAL.csv',output)
    assert len(source_rows)==6
    writecsv(OUT/'S1_Six_Raw_Sources_QC.csv',source_rows)
    writecsv(OUT/'S2_All_Cut_FullRaw_Sensitivity.csv',output)
    focus=[r for r in output if r['HAR'] in CAND and r['min_DNA_count_per_barcode']==1 and r['aggregate_pseudocount']==1]
    assert len(focus)==30
    writecsv(OUT/'S3_FiveCandidate_FullRaw_Sensitivity.csv',focus)
    out={'stage':'21D_FULL_RAW','run':'completed','n_GEO_count_files':6,'n_selected_GEO_HARs':75,
         'n_prior_pairable_cuts':114,'n_culture_replicates':3,
         'first_row_headerless_corrected':True,'n_output_sensitivity_rows':len(output),
         'max_Stage21C_2026_paired_contrast_recalc_abs_error':max_error,
         'tested_DNA_depths':[1,5,10,20],'tested_aggregate_pseudocounts':[0.5,1,2],
         'Noble_primary_HAR_BED_identity_confirmed':False,
         'warning':'This is one human iPSC line with 3 cultures, not three independent human donors. Comparisons are exploratory.'}
    (OUT/'S0_FullRaw_Status.json').write_text(json.dumps(out,indent=2))
    print('SUCCESS',json.dumps(out,indent=2))

def package(name):
    target=ROOT/name
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(OUT.glob('*')):
            if p.is_file():z.write(p,p.name)
    print('Results:',target,'bytes',target.stat().st_size)
    print('Download manually from the Colab Files sidebar (no files.download inside Python subprocess).')

try:
    workflow()
    package('HAR_Stage21D_FullSixRaw_Corrected_Results.zip')
except Exception as e:
    (OUT/'S0_Diagnostics.txt').write_text(traceback.format_exc())
    package('HAR_Stage21D_FullSixRaw_Diagnostics.zip')
    print('Error:',type(e).__name__,str(e))
