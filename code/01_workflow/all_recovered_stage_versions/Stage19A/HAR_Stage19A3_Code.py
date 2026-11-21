# HAR Stage 19A3 | No-Google-Drive genome-flank guide rescue
# Install dependency if needed: pip install requests
# Run: python HAR_Stage19A3_Code.py

import csv, json, time, zipfile, hashlib, requests
from pathlib import Path
from collections import defaultdict

RECORDS = json.loads('[{"HAR": "2xHAR.183", "chrom": "chr2", "end0": 11391826, "guides": [{"guide_id": "2xHAR_20", "sequence": "GGGTCCAGGCACACTCGGTC"}], "start0": 11391805}, {"HAR": "2xHAR.223", "chrom": "chr14", "end0": 33661912, "guides": [{"guide_id": "2xHAR_29", "sequence": "GCACATTTACATACTCTCGA"}, {"guide_id": "2xHAR_30", "sequence": "GCCACCATATCTAGACGGCA"}], "start0": 33661114}, {"HAR": "2xHAR.274", "chrom": "chr1", "end0": 243476360, "guides": [{"guide_id": "2xHAR_46", "sequence": "ACATGGAACTATCAGACGCA"}], "start0": 243476286}, {"HAR": "2xHAR.435", "chrom": "chr2", "end0": 204082416, "guides": [{"guide_id": "2xHAR_82", "sequence": "GTTAGCCATCATCTCCCGTG"}], "start0": 204082251}, {"HAR": "2xHAR.510", "chrom": "chr5", "end0": 171667741, "guides": [{"guide_id": "2xHAR_101", "sequence": "CACCCCGCATGAATGCGGCT"}], "start0": 171667678}, {"HAR": "HACNS188", "chrom": "chr13", "end0": 104144430, "guides": [{"guide_id": "HACNS_16", "sequence": "TGTGAAACTGGCGACTTAAG"}], "start0": 104144174}, {"HAR": "HACNS336", "chrom": "chr13", "end0": 43315916, "guides": [{"guide_id": "HACNS_41", "sequence": "CATGCTGCTAGAGCGCAATG"}], "start0": 43315612}, {"HAR": "HACNS50", "chrom": "chr7", "end0": 25861122, "guides": [{"guide_id": "HACNS_72", "sequence": "GCAAATGCACTGGACCGCTT"}], "start0": 25860952}, {"HAR": "HACNS544", "chrom": "chr14", "end0": 81681713, "guides": [{"guide_id": "HACNS_77", "sequence": "ATGTATGTGATGAGTCGCGT"}], "start0": 81681397}, {"HAR": "HACNS643", "chrom": "chr13", "end0": 94948541, "guides": [{"guide_id": "HACNS_103", "sequence": "AACGGATTATGCAACGCAAC"}], "start0": 94947886}, {"HAR": "HACNS659", "chrom": "chr6", "end0": 22118012, "guides": [{"guide_id": "HACNS_108", "sequence": "TAATCAGTAGAGCACGACAA"}], "start0": 22117752}, {"HAR": "HACNS712", "chrom": "chr4", "end0": 168743936, "guides": [{"guide_id": "HACNS_129", "sequence": "TGTGAGAATCGTGCGTGTCC"}], "start0": 168743672}, {"HAR": "HACNS719", "chrom": "chr2", "end0": 175569365, "guides": [{"guide_id": "HACNS_133", "sequence": "GGAGGTAGGGGGTACCCATA"}], "start0": 175569173}, {"HAR": "HACNS76", "chrom": "chr9", "end0": 1202104, "guides": [{"guide_id": "HACNS_144", "sequence": "ACATATTTAGGTTTGATCGT"}], "start0": 1201875}, {"HAR": "HACNS763", "chrom": "chr13", "end0": 106294777, "guides": [{"guide_id": "HACNS_149", "sequence": "AGATGTCAAACTTCCTCGAG"}], "start0": 106294649}, {"HAR": "HACNS772", "chrom": "chr17", "end0": 65439366, "guides": [{"guide_id": "HACNS_151", "sequence": "CCGCTCGATCACCCCCTGAT"}], "start0": 65438982}, {"HAR": "HACNS858", "chrom": "chr11", "end0": 119980586, "guides": [{"guide_id": "HACNS_165", "sequence": "TCCCCCGTCAGATTTCCGGC"}], "start0": 119980415}, {"HAR": "HACNS880", "chrom": "chr5", "end0": 126924809, "guides": [{"guide_id": "HACNS_173", "sequence": "TGTGAGAGGCTTAATACGTG"}], "start0": 126924602}, {"HAR": "HACNS962", "chrom": "chr20", "end0": 46225479, "guides": [{"guide_id": "HACNS_192", "sequence": "TGACTGTCAAGAAAGGGCGC"}], "start0": 46225170}, {"HAR": "HAR108", "chrom": "chr16", "end0": 73482925, "guides": [{"guide_id": "HAR_3", "sequence": "CTGAGGCCCGATCCCTCCGA"}, {"guide_id": "HAR_4", "sequence": "CCCGGAGGAAGAGAACGCGT"}], "start0": 73482771}, {"HAR": "HAR156", "chrom": "chr8", "end0": 115779602, "guides": [{"guide_id": "HAR_12", "sequence": "TCCTACTCCCTATTCACCCG"}], "start0": 115779442}, {"HAR": "HAR169", "chrom": "chr2", "end0": 176851576, "guides": [{"guide_id": "HAR_17", "sequence": "GGTGGTTTAGCACCTCAATA"}, {"guide_id": "HAR_18", "sequence": "TTATTCCAATGGCCATGTCG"}], "start0": 176851415}, {"HAR": "HAR181", "chrom": "chr11", "end0": 114928663, "guides": [{"guide_id": "HAR_23", "sequence": "TAGACTGAAATCATCTGCGT"}], "start0": 114928473}, {"HAR": "HAR202", "chrom": "chr14", "end0": 33576566, "guides": [{"guide_id": "HAR_33", "sequence": "AAGTGAAAGATACTCGGGAC"}], "start0": 33576363}, {"HAR": "HAR67", "chrom": "chr13", "end0": 71799565, "guides": [{"guide_id": "HAR_37", "sequence": "GCAGTAGGCATTATGACTAC"}], "start0": 71799450}]')
EXPECTED_SHA256 = '4cd882c8a981a1bfc56ae26701316f49ce6f33bda0aac73a86978988e1d212ae'
assert hashlib.sha256(json.dumps(RECORDS,ensure_ascii=True,sort_keys=True).encode()).hexdigest()==EXPECTED_SHA256
assert len(RECORDS)==25 and sum(len(r["guides"]) for r in RECORDS)==28
OUT=Path("/content/HAR_Stage19A3_results") if Path("/content").exists() else Path("HAR_Stage19A3_results")
OUT.mkdir(exist_ok=True,parents=True)
FLANK_BP=2000
print("Embedded source rows:",len(RECORDS),"HARs,",sum(len(r["guides"]) for r in RECORDS),"unmatched sgRNAs")


# UCSC reference genome retrieval; retries and explicit failure statuses.
S=requests.Session()
S.headers.update({'User-Agent':'HAR-Stage19A3-Research-Crosscheck/1.0'})
CACHE=OUT/'UCSC_25_Slices.jsonl'
cache={}
if CACHE.exists():
    for line in CACHE.read_text().splitlines():
        try:
            obj=json.loads(line)
            if obj.get('dna'): cache[obj['key']]=obj
        except Exception: pass

def fetch(r):
    left=max(0,r['start0']-FLANK_BP)
    right=r['end0']+FLANK_BP
    key=f"{r['chrom']}:{left}-{right}"
    if key in cache:return cache[key]
    error=''
    for nap in (0,2,6,12):
        if nap:time.sleep(nap)
        try:
            res=S.get('https://api.genome.ucsc.edu/getData/sequence',params={
                'genome':'hg38','chrom':r['chrom'],'start':left,'end':right}, timeout=50)
            res.raise_for_status(); j=res.json()
            dna=j['dna'].upper()
            if len(dna)!=right-left:raise ValueError('unexpected sequence length')
            obj={'key':key,'left':left,'right':right,'dna':dna,'url':res.url}
            cache[key]=obj
            with CACHE.open('a') as f:f.write(json.dumps(obj)+'\n')
            return obj
        except Exception as e:error=f'{type(e).__name__}: {e}'
    return {'key':key,'left':left,'right':right,'dna':None,'error':error}

for i,r in enumerate(RECORDS,1):
    r['reference']=fetch(r)
    if i%5==0 or i==len(RECORDS):print(f'Slices {i}/{len(RECORDS)}. Missing: {sum(x.get("reference",{}).get("dna") is None for x in RECORDS[:i])}')
    if i%6==0:time.sleep(.3)


COMP=str.maketrans('ACGT','TGCA')
def rc(s):return s.translate(COMP)[::-1]
def scan(seq,guide):
    if not seq or len(guide)!=20 or set(guide.upper())-set('ACGT'):return []
    hits=[]
    for strand,pattern in [('+',guide.upper()),('-',rc(guide.upper()))]:
        offset=0
        while True:
            at=seq.find(pattern,offset)
            if at<0:break
            pam=seq[at+20:at+23] if strand=='+' else rc(seq[at-3:at]) if at>=3 else ''
            hits.append((at,strand,pam,bool(len(pam)==3 and pam[1:]=='GG')))
            offset=at+1
    return hits

rows=[];status=[]
for r in RECORDS:
    ref=r['reference'];seq=ref.get('dna');left=ref['left'];start=r['start0'];end=r['end0']
    successes=0
    for item in r['guides']:
        hits=scan(seq,item['sequence'])
        # Require a complete 20-mer target to be inside the original 250-bp flank window.
        prior_hits=[at for at,_,_,_ in hits if at>=max(0,start-250)-left and at+20<=end+250-left]
        if prior_hits:raise AssertionError('Stage19A2 mismatch: guide appeared within original ±250bp window: '+item['guide_id'])
        if hits:successes+=1
        if not hits:
            rows.append({'HAR':r['HAR'],'guide_id':item['guide_id'],'sequence_20mer':item['sequence'],
                'reference_status':'FETCHED' if seq else 'FAILED','guide_hit_status':'NO_HIT_2KB' if seq else 'NO_REFERENCE',
                'n_hits_2kb':0,'start0':None,'end_exclusive':None,'strand':'','gap_to_core_bp':None,
                'relation':'','PAM_3mer':'','NGG_PAM':False,'UCSC_url':ref.get('url',''), 'error':ref.get('error','')})
        for at,strand,pam,ngg in hits:
            abs_start=left+at;abs_end=abs_start+20
            gap=max(start-abs_end,abs_start-end,0)
            relation='INSIDE' if abs_start>=start and abs_end<=end else 'UPSTREAM' if abs_end<=start else 'DOWNSTREAM' if abs_start>=end else 'CROSSES'
            rows.append({'HAR':r['HAR'],'guide_id':item['guide_id'],'sequence_20mer':item['sequence'],
                'reference_status':'FETCHED','guide_hit_status':'MULTI_HIT' if len(hits)>1 else 'EXACT_20MER',
                'n_hits_2kb':len(hits),'start0':abs_start,'end_exclusive':abs_end,'strand':strand,'gap_to_core_bp':gap,
                'relation':relation,'PAM_3mer':pam,'NGG_PAM':ngg,'UCSC_url':ref['url'],'error':''})
    status.append({'HAR':r['HAR'],'chrom':r['chrom'],'catalog_start0':start,'catalog_end_exclusive':end,
        'n_previously_unmatched_guides':len(r['guides']),'guides_recovered_2kb':successes,
        'reference_fetched':seq is not None,'exact_original_HAR_boundaries':'NOT_VERIFIED'})

def write_csv(path,data):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
write_csv(OUT/'S1_28_Guide_2kb_Search.csv',rows)
write_csv(OUT/'S2_25_HAR_2kb_Summary.csv',status)
manifest={'stage':'19A3','genome':'hg38','flank_bp_each_side':FLANK_BP,
    'n_guide_inputs':28,'n_har_inputs':25,
    'n_guide_recovered':sum(s['guides_recovered_2kb'] for s in status),
    'n_guide_remain_unresolved':28-sum(s['guides_recovered_2kb'] for s in status),
    'n_reference_failures':sum(not s['reference_fetched'] for s in status),
    'original_har_BED_confirmed':False,
    'embedded_input_sha256':EXPECTED_SHA256,
    'caution':'Local guide proximity is not exact author HAR identity, regulatory target evidence, or genome-wide off-target specificity.'}
(OUT/'STAGE19A3_STATUS.json').write_text(json.dumps(manifest,indent=2))
result=OUT.parent/'HAR_Stage19A3_Expanded_Flank_Results.zip'
with zipfile.ZipFile(result,'w',zipfile.ZIP_DEFLATED,compresslevel=8) as z:
    for n in ['S1_28_Guide_2kb_Search.csv','S2_25_HAR_2kb_Summary.csv','STAGE19A3_STATUS.json']:
        z.write(OUT/n,arcname=n)
print(json.dumps(manifest,indent=2))
try:
    from google.colab import files
    files.download(str(result))
except ImportError:
    print('Download the ZIP at:',result)

