import csv,json,hashlib,zipfile,io,statistics
from pathlib import Path
from collections import defaultdict
from artifact_tool import Workbook, SpreadsheetFile
ROOT=Path('/mnt/data'); OUT=ROOT/'Stage22_Scientific_Gate'; OUT.mkdir(exist_ok=True)
with zipfile.ZipFile(ROOT/'HAR_Stage21D_FullSixRaw_Corrected_Results.zip') as z:
 sourceqc=list(csv.DictReader(io.StringIO(z.read('S1_Six_Raw_Sources_QC.csv').decode())))
def csvrows(name):return list(csv.DictReader(open(OUT/name,encoding='utf-8')))
rows=csvrows('S1_All114_Reporter_StatisticalGate.csv'); sites=csvrows('S2_Motif_Site_Cut_Coverage.csv'); five=csvrows('S3_Five_Candidate_Reporter_Gate.csv')
rawreps=list(csv.DictReader(open(ROOT/'S3_Paired_Fragment_Activity_By_Replicate.csv',encoding='utf-8')))
repmap=defaultdict(dict)
for r in rawreps:repmap[(r['HAR'],r['cut'])][int(r['replicate'])]=float(r['human_minus_chimp_log2_activity'])
assert len(rows)==114 and len(five)==10 and len(sites)==6 and len(sourceqc)==6 and len(repmap)==114
srcs=[]
for r in sourceqc:
 p=ROOT/r['filename']; matched=(hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']) if p.exists() else None
 srcs.append([r['replicate'],r['assay'],r['filename'],int(r['rows']),int(r['library_counts']),r['sha256'], 'LOCALLY_SHA256_VERIFIED' if matched else 'COLAB_MANIFEST_ONLY_RAW_BYTES_NOT_AVAILABLE', 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110758'])
assert sum(s[6]=='LOCALLY_SHA256_VERIFIED' for s in srcs)==3
sourceGEO='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110758'
sourceHAR='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270828'
wb=Workbook.create()
NAVY='#132C45'; BLUE='#0C7A88'; PALE='#EAF3F6'; WHITE='#FFFFFF'
def sheet(name,title,note,headers,records,widths=None):
 sh=wb.worksheets.add(name); n=len(headers)
 sh.merge_cells(f'A1:{chr(64+n)}1');sh.get_range('A1').values=[[title]]
 sh.get_range(f'A1:{chr(64+n)}1').format={'fill':NAVY,'font':{'color':WHITE,'bold':True,'size':13},'row_height':33}
 sh.merge_cells(f'A2:{chr(64+n)}2');sh.get_range('A2').values=[[note]]
 sh.get_range(f'A2:{chr(64+n)}2').format={'fill':PALE,'font':{'color':'#264653','size':10},'row_height':36,'wrap_text':True}
 sh.get_range_by_indexes(3,0,1,n).values=[headers]
 sh.get_range_by_indexes(3,0,1,n).format={'fill':BLUE,'font':{'color':WHITE,'bold':True},'row_height':32,'wrap_text':True}
 for offset in range(0,len(records),150): sh.get_range_by_indexes(4+offset,0,len(records[offset:offset+150]),n).values=records[offset:offset+150]
 sh.freeze_panes.freeze_rows(4)
 sh.get_range_by_indexes(0,0,max(len(records)+4,7),n).format.column_width=18
 for j,w in enumerate(widths or []):sh.get_range_by_indexes(0,j,max(len(records)+4,7),1).format.column_width=w
 return sh
s=wb.worksheets.add('OVERVIEW')
s.merge_cells('A1:F1');s.get_range('A1').values=[['STAGE 22 | HAR SCIENTIFIC VALIDITY GATE']]
s.get_range('A1:F1').format={'fill':NAVY,'font':{'color':WHITE,'size':15,'bold':True},'row_height':37}
s.get_range('A3:C3').values=[['Evidence component','Verified result','Interpretation']]
s.get_range('A3:C3').format={'fill':BLUE,'font':{'color':WHITE,'bold':True}}
kpis=[
 ['HARs targeted in Noble screen',180,'Original primary study: GSE270828'],
 ['Provisional ortholog pairs',152,'Secondary catalog correspondence, primary BED unresolved'],
 ['Context-null FDR-significant HARs',0,'Stage 20B BH q<0.05; none'],
 ['HAR names in independent lentiMPRA design',75,'Name overlap only'],
 ['HARs with three-culture MPRA cuts',66,'Same host human HS1-11 line'],
 ['Paired reporter cuts with 3 cultures',114,'Three culture averages; not 114 independent biological repeats'],
 ['Cuts with same reporter sign in 3 cultures',68,'Descriptive reproducibility'],
 ['Min 2-sided exact sign p with 3 cultures',0.25,'Cannot claim sign-test significance'],
 ['Stage 21D exploratory t-FDR positives',1,'2xHAR.385 cut 2; fragile inferential interpretation'],
 ['Original Noble HAR BED confirmed',0,'Essential unresolved provenance gate'],
 ['Motif site directly overlapped by reported MPRA cuts',2,'Two of three Stage 18 motif cases tested by MPRA; 2xHAR.377 partial mapping'],
 ]
s.get_range_by_indexes(3,0,len(kpis),3).values=kpis
s.get_range('A:A').format.column_width=48;s.get_range('B:B').format.column_width=18;s.get_range('C:C').format.column_width=82
s.get_range('A18').values=[['DECISION']];s.get_range('A18').format={'font':{'bold':True,'color':NAVY,'size':12}}
s.merge_cells('A19:F20');s.get_range('A19').values=[['Suitable as a carefully scoped, multi-layer reproducibility and evidence-calibration study. Not yet a confirmed HAR-TF-gene mechanism. Required for stronger biological claims: primary HAR locus crosswalk + independent reporter/TF validation.']]
s.get_range('A19:F20').format={'fill':PALE,'wrap_text':True,'row_height':31,'font':{'color':'#253B4B'}}
s.get_range('A23:B25').values=[['Source','URL'],['Noble/GSE270828',sourceHAR],['MPRA/GSE110758',sourceGEO]]
s.get_range('B:B').format.column_width=83

reporter=[]
for r in rows:
 key=(r['HAR'],r['cut']);vals=repmap[key]
 reporter.append([r['HAR'],r['cut'],vals[1],vals[2],vals[3],None,
                  r['all_3_same_direction']=='True',float(r['exact_sign_2s_p']),float(r['nominal_student_t_p']),float(r['exploratory_t_BH_q_114']),(float(r['Stage20B_BH_q']) if r['Stage20B_BH_q'].lower() not in ('nan','') else None),int(r['sign_flips_under_DNA5_10_20']),int(r['missing_threshold_cells']),sourceGEO])
sh=sheet('S1_Reporter114','S1 | 114 REPORTER CUTS: THREE CULTURES',
'Frozen exploratory p/q diagnostics; means calculated in worksheet; exact sign-test p is non-parametric, 3 cultures per cut, single iPSC line.',
['HAR','cut','rep1 H-C','rep2 H-C','rep3 H-C','mean H-C (formula)','3/3 same sign','exact sign p','exploratory t p','exploratory t BH q','Stage20B motif BH q','sign flips under DNA thresholds','missing sensitivity cells','Source URL'],reporter,[19,14,16,16,16,20,16,15,18,20,20,27,25,55])
# Formula columns are editable and calculate mean from 3 experimentally sampled culture contrasts
sh.get_range('F5').formulas=[['=AVERAGE(C5:E5)']]; sh.get_range(f'F5:F{len(reporter)+4}').fill_down()
sh.get_range(f'C5:F{len(reporter)+4}').set_number_format('0.000000')
sh.get_range(f'H5:K{len(reporter)+4}').set_number_format('0.000000')

siteout=[]
for r in sites:
 st=int(r['catalog_payload_start_0based']);en=int(r['catalog_payload_end_0based_exclusive']);pos=int(r['motif_window_start_aligned_1based']);within=st<=pos-1<en
 siteout.append([r['HAR'],r['TF'],pos,r['MPRA_cut'],st,en, bool(within),float(r['payload_exact_fraction']),
                 'motif-start inside; full width/author boundaries need confirmation' if within else 'assayed cut does not include predicted motif site',sourceGEO])
sheet('S2_MotifVsMPRA','S2 | MOTIF-TO-REPORTER SITE CROSSWALK',
'Coordinates represent provisional secondary catalog local positions; MPRA library uses hg19 labels. Not a primary-author genomic crosswalk.',
['HAR','TF','aligned motif start 1-based','MPRA cut','payload start 0-based','payload end exclusive','motif start inside','exact payload fraction','Evidence classification','Source URL'],siteout,[19,14,28,14,25,25,22,23,68,50])

fivelook={(r['HAR'],r['cut']):r for r in five}
fivereg=[]
for k in [('2xHAR.70','1_of_3'),('2xHAR.70','2_of_3'),('2xHAR.70','3_of_3'),('2xHAR.310','1_of_2'),('2xHAR.310','2_of_2'),('2xHAR.377','1_of_1'),('2xHAR.234','1_of_2'),('2xHAR.234','2_of_2'),('2xHAR.3','1_of_2'),('2xHAR.3','2_of_2')]:
 r=fivelook[k];h,cut=k
 mflag='Stage18 predicted SOX11 site outside both MPRA cuts' if h=='2xHAR.310' else ('Stage18 motif site within cut1; cut3 is non-overlapping' if h=='2xHAR.70' else ('Stage18 SOX2 start inside cut; catalog payload partial' if h=='2xHAR.377' else 'Stage20B motif event location not linked to cut here'))
 fivereg.append([h,cut,float(r['mean_HminusC_log2']),r['all_3_same_direction']=='True',float(r['exploratory_t_BH_q_114']),float(r['Stage20B_3mer_p']),(float(r['Stage20B_BH_q']) if r['Stage20B_BH_q'].lower() not in ('nan','') else None),int(r['sign_flips_under_DNA5_10_20']),mflag,sourceGEO])
sheet('S3_FiveHAR','S3 | FIVE SHORTLISTED HARs',
'Candidate status exploratory. No five-candidate MPRA cut or Stage 20B motif comparison achieved FDR q<0.05.',
['HAR','MPRA cut','mean H-C log2','3/3 same sign','MPRA t BH q (exploratory)','3mer context p','motif BH q','depth-induced sign flips','TF-motif site assessment','Source URL'],fivereg,[19,15,19,20,26,18,16,23,67,50])

sheet('S4_RawSources','S4 | SIX SOURCE COUNT FILES AND CHECKSUMS',
'Source manifest from Stage 21D. Three local uploaded raw files independently checked; three others checksum recorded from completed Colab run, not re-fetched here.',
['culture','assay','GEO basename','rows (headerless)','library counts','SHA-256','verification level','GEO accession'],srcs,[11,12,31,19,19,71,50,55])

claims=[
['The original 180-HAR perturbation dataset was reprocessed','SUPPORTED','Stage12–15 raw-count reconciliation','Use one iPSC source line and batch/guide as distinct technical units',sourceHAR],
['152 human–chimp HAR orthologs are author-verified','NOT ESTABLISHED','Stage19 mapping from secondary catalog; author BED absent','Acquire primary original 180-HAR BED and assembly crosswalk',sourceHAR],
['HAR TF motif divergence passes genome-wide FDR','NOT ESTABLISHED','Stage20B all 152 BH q>=0.05','Avoid TF binding/derived lineage-specific claims',sourceHAR],
['The 5 priority HARs show confirmed differential MPRA','NOT ESTABLISHED','10 cuts in 3 cultures, none corrected t BH q<0.05','Original count-based model & independent lines/stages required',sourceGEO],
['One 2xHAR.385 cut has an exploratory t-test BH q<0.05','EXPLORATORY ONLY','3-culture parametric t q=0.0179; exact sign p=0.25','Do not present as independent donor/TF validation',sourceGEO],
['2xHAR.310 reporter validates predicted SOX11 motif','UNSUPPORTED - WRONG SITE','Predicted SOX11 site starts around 581, assayed cuts cover 11-208','Assay motif-overlapping interval or omit mechanism claim',sourceGEO],
['2xHAR.70 MPRA cuts support one directional HAR effect','NOT ESTABLISHED','Cut 1 -0.279 vs Cut 3 +0.301, only cut1 overlaps predicted SOX11 site','Treat fragment-cut results separately',sourceGEO],
['HAR–gene causal targets are experimentally established by names','UNSUPPORTED','Chromatin table name representation only, guide RNA contrasts descriptive','Gene-resolved contacts + guide concordance + independent perturbation',sourceHAR],
['Sequence/pseudobulk/MPRA reliability benchmark is reproducible','SUPPORTED WITH LIMITS','Independent checks, full raw MPRA sensitivity and saved hashes','Keep exact source files, scripts, checksums and selected sample universe',sourceGEO],
]
sheet('S5_Claims','S5 | CLAIMS THAT A JOURNAL REVIEWER MAY CHECK',
'Promote only claims that are directly supported by uploaded evidence; avoid inferring causality from associations or motif predictions.',
['Potential manuscript claim','Audit status','Evidence','Needed change or next step','Source URL'],claims,[72,35,86,91,55])

print('first row audit',wb.inspect({'kind':'table','range':'S3_FiveHAR!A4:G7','include':'values','table_max_rows':4,'table_max_cols':7}).ndjson[:900])
print('formula check',wb.inspect({'kind':'match','search_term':'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A','options':{'use_regex':True,'max_results':30}}).ndjson[:400])
path=ROOT/'HAR_Stage22_Scientific_Evidence_Gate.xlsx';SpreadsheetFile.export_xlsx(wb).save(str(path))
print('exported',path,path.stat().st_size)
