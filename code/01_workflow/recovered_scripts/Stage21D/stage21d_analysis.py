#!/usr/bin/env python3
"""Stage 21D reproducible audit from supplied Stage21C output, 3 uploaded raw count files and Stage19A FASTA.
No internet, no donor-level claims, no recalculation from unprovided RNA2/DNA2/DNA3.
Outputs: CSVs, QC status, markdown report, figures and notebook-derived reusable building blocks.
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv, gzip, io, math, json, re, zipfile, hashlib
from statistics import mean, stdev
import numpy as np
from scipy import stats
from Bio.Align import PairwiseAligner

ROOT=Path('/mnt/data')
OUT=ROOT/'Stage21D_Work'
ZIP=ROOT/'HAR_Stage21C_MPRA_HS111_N2_Results.zip'
FAS=ROOT/'HAR_Stage19A_Colab_Source_Verified_Results.zip'
DESIGN=ROOT/'GSE110758_HAR_design.tsv.gz'
CAND=['2xHAR.70','2xHAR.310','2xHAR.234','2xHAR.377','2xHAR.3']


def readzip_csv(name):
    with zipfile.ZipFile(ZIP) as z:
        return list(csv.DictReader(io.StringIO(z.read(name).decode('utf-8'))))

def savecsv(name,rows,cols=None):
    if cols is None:cols=list(rows[0]) if rows else []
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)

def num(s):return float(s)
def bh(ps):
    n=len(ps); order=sorted(range(n),key=lambda j:ps[j]); q=[None]*n;best=1.0
    for j in range(n-1,-1,-1):
        k=order[j];best=min(best,ps[k]*n/(j+1));q[k]=best
    return q

def rawcounts(path):
    count=Counter();tot=0;n=0;bad=0
    with gzip.open(path,'rt') as f:
        for row in csv.reader(f,delimiter='\t'):
            if len(row)<2:bad+=1;continue
            bc=row[0].strip().upper()
            if not re.fullmatch('[ACGT]{15}',bc):bad+=1;continue
            try:v=int(row[1])
            except:bad+=1;continue
            if v<0: raise ValueError('Negative count')
            count[bc]+=v; tot+=v;n+=1
    assert not bad, (path,bad)
    return count,tot,n

s3=readzip_csv('S3_Paired_Fragment_Activity_By_Replicate.csv')
s2=readzip_csv('S2_Insert_Activity_By_Replicate.csv')
s4=readzip_csv('S4_Replicate_QC.csv')
s5=readzip_csv('S5_75HAR_Replicate_Summary_By_Cut.csv')
barcodes=readzip_csv('S1B_Assayed_Barcode_Map.csv')
assert len(s3)==342 and len(s2)==765 and len(s5)==114 and len(barcodes)==25500

# Independent three-culture summaries, exploratory t intervals and appropriate multiplicity awareness.
groups=defaultdict(dict)
for row in s3:
    groups[(row['HAR'],row['cut'])][int(row['replicate'])]=num(row['human_minus_chimp_log2_activity'])
assert len(groups)==114 and all(set(x)=={1,2,3} for x in groups.values())
allcuts=[]
for (har,cut),dic in sorted(groups.items()):
    a=np.array([dic[i] for i in (1,2,3)],float)
    avg=float(a.mean()); sd=float(a.std(ddof=1)); sem=sd/(3**.5)
    crit=float(stats.t.ppf(.975,df=2))
    if sd==0:
        pval=1.0 if avg==0 else None
    else:
        pval=float(2*stats.t.sf(abs(avg/sem),2))
    allcuts.append({'HAR':har,'cut':cut,'n_cultures':3,'rep1':a[0],'rep2':a[1],'rep3':a[2],
                    'mean_HminusC_log2':avg,'sd_between_cultures':sd,'se_between_cultures':sem,
                    'exploratory_t95_lower':avg-crit*sem,'exploratory_t95_upper':avg+crit*sem,
                    'exploratory_t_p':pval,'same_direction_3_of_3':bool(np.all(a>0) or np.all(a<0)),
                    'exact_two_sided_sign_test_p':0.25 if np.all(a>0) or np.all(a<0) else 1.0,
                    'analysis_label':'EXPLORATORY_N3_ONE_I PSC_LINE'.replace(' ','')})
q=bh([x['exploratory_t_p'] if x['exploratory_t_p'] is not None else 1.0 for x in allcuts])
for x,v in zip(allcuts,q):x['exploratory_BH_q_114_cuts']=v
savecsv('S1_114_CutReplicate_Exploratory_Intervals.csv',allcuts)
six=[x for x in allcuts if x['HAR'] in CAND]
assert len(six)==10
savecsv('S2_Five_HAR_Ten_Cut_Exploratory_Results.csv',six)

# Pseudocount sensitivity of 21C positive-positive selected barcode aggregated counts.
qc={int(x['replicate']):x for x in s4}
sp={}
for x in s2:
    rep=int(x['replicate']);key=(x['HAR'],x['cut'],rep,x['allele'])
    tot_r=num(qc[rep]['total_RNA_counts']);tot_d=num(qc[rep]['total_DNA_counts'])
    sp[key]={'n':int(x['n_paired_barcodes']),'R':num(x['rna_cpm'])*tot_r/1e6,'D':num(x['dna_cpm'])*tot_d/1e6,
             'rep':rep}
old={(x['HAR'],x['cut'],int(x['replicate'])):num(x['human_minus_chimp_log2_activity']) for x in s3}
sens=[]; basemax=0
for (har,cut),reps in sorted(groups.items()):
    for rep in sorted(reps):
        h=sp[(har,cut,rep,'human')];c=sp[(har,cut,rep,'chimp')]
        for alpha in (0,0.5,1,2,5):
            diff=math.log2((h['R']+alpha)/(h['D']+alpha))-math.log2((c['R']+alpha)/(c['D']+alpha))
            if alpha==0:basemax=max(basemax,abs(diff-old[(har,cut,rep)]))
            sens.append({'HAR':har,'cut':cut,'replicate':rep,'aggregate_count_pseudocount':alpha,
                         'positive_barcodes_human':h['n'],'positive_barcodes_chimp':c['n'],
                         'HminusC_log2_positive_only':diff,
                         'delta_from_original_21C':diff-old[(har,cut,rep)],
                         'scope':'Positive-RNA-and-DNA barcodes only; does not add zero-RNA barcodes'})
assert basemax<1e-10,basemax
savecsv('S3_114Cut_3Rep_AggregatePseudocount_Sensitivity.csv',sens)

# Raw replicate1: strict headerless parse, zero-RNA sensitivity, depth filters.
R,Rtotal,Rn=rawcounts(ROOT/'GSM3016032_RNA.tsv.gz')
D,Dtotal,Dn=rawcounts(ROOT/'GSM3016050_DNA.tsv.gz')
assert (Rtotal,Dtotal,Rn,Dn)==(27936486,8685572,201098,197187)
assert len(set(x['Barcode'] for x in barcodes))==25500
bygrp=defaultdict(list)
for x in barcodes:bygrp[(x['HAR'],x['cut'],x['allele'])].append(x['Barcode'])
rawrows=[];bcore={}
for depth in (1,5,10,20):
    for alpha in (0.5,1,2):
        out={}
        for key,bcs in bygrp.items():
            eligible=[b for b in bcs if D.get(b,0)>=depth]
            if not eligible:continue
            rsum=sum(R.get(b,0) for b in eligible);dsum=sum(D.get(b,0) for b in eligible)
            nr=sum(R.get(b,0)>0 for b in eligible)
            out[key]={'R':rsum,'D':dsum,'n':len(eligible),'zero_RNA':len(eligible)-nr,'act':math.log2((rsum+alpha)/(dsum+alpha))}
        for har,cut in groups:
            hh=out.get((har,cut,'human'));ch=out.get((har,cut,'chimp'))
            if not (hh and ch):continue
            if hh['n']<10 or ch['n']<10:continue
            diff=hh['act']-ch['act']
            row={'HAR':har,'cut':cut,'replicate':1,'min_DNA_reads_per_barcode':depth,'aggregate_pseudocount':alpha,
                 'n_DNA_positive_human':hh['n'],'n_DNA_positive_chimp':ch['n'],
                 'n_missing_RNA_human':hh['zero_RNA'],'n_missing_RNA_chimp':ch['zero_RNA'],
                 'HminusC_including_RNA_zeros':diff,'original_RNA_and_DNA_positive_HminusC':old[(har,cut,1)],
                 'delta_vs_original':diff-old[(har,cut,1)],'scope':'RAW_REPLICATE_1_ONLY'}
            rawrows.append(row)
            if depth==1 and alpha==1:bcore[(har,cut)]=row
savecsv('S4_RawRep1_ZeroRNA_and_DNADepth_Sensitivity.csv',rawrows)
assert len(bcore)==114

# Conservatively map the actual published MPRA construct insert core to prior hg38/panTro6 sequence catalog.
# This is sequence alignment only; does not establish original author's target BED coordinates.
with zipfile.ZipFile(FAS) as z:
    seq={}
    for alle,fn in [('human','S7_152_Human_CatalogSequences_PROVISIONAL.fa'),('chimp','S8_152_Chimp_CatalogSequences_PROVISIONAL.fa')]:
        lines=z.read(fn).decode().splitlines()
        for i in range(0,len(lines),2):seq[(lines[i][1:].split('|')[0],alle)]=lines[i+1]
construct={}
with gzip.open(DESIGN,'rt') as f:
    for row in csv.DictReader(f,delimiter='\t'):
        nm=row['#Name'];har=nm.split('__',1)[0]
        if har not in CAND:continue
        alle=re.search('__species=([^_]+)',nm)
        cut=re.search('__cut=([^_]+)',nm)
        if not (alle and cut):continue
        key=(har,alle.group(1),cut.group(1))
        if key in construct:
            if construct[key]['Sequence'][15:186]!=row['Sequence'][15:186]:raise ValueError('Inconsistent biological 171bp payload by barcode '+str(key))
        else:construct[key]=row
assert len(construct)==20

aligner=PairwiseAligner(mode='local',match_score=2,mismatch_score=-2,open_gap_score=-6,extend_gap_score=-1)
trans=str.maketrans('ACGT','TGCA')
alignment_rows=[]
for key,row in sorted(construct.items()):
    har,alle,cut=key
    ref=seq[(har,alle)]
    # In this specific 171bp GSE110758 construct design, sequence has a 15bp constant 5' primer.
    payload=row['Sequence'][15:186]
    assert len(row['Sequence'])==230 and len(payload)==171
    best=None
    for strand,template in [('+',ref),('-',ref.translate(trans)[::-1])]:
        a=aligner.align(template,payload)[0]
        idxr,idxq=a.aligned
        n_ref=int(sum(e-s for s,e in idxr));n_query=int(sum(e-s for s,e in idxq))
        matched=0
        for (rs,rend),(qs,qe) in zip(idxr,idxq):
            matched+=sum(u==v for u,v in zip(template[rs:rend],payload[qs:qe]))
        res={'HAR':har,'allele':alle,'cut':cut,'assembly_catalog':'hg38' if alle=='human' else 'panTro6',
             'MPRA_coordinate_label_hg19':re.search(r'__coordinates_tested_in_hg19=\[([^\]]+)\]',row['#Name']).group(1),
             'MPRA_total_oligo_length':len(row['Sequence']),'MPRA_171bp_payload_length':len(payload),
             'catalog_length':len(template),'sequence_orientation':strand,
             'aligned_MPRA_bases':n_query,'identical_MPRA_bases':matched,
             'MPRA_payload_coverage':n_query/len(payload),
             'fraction_identical_of_payload':matched/len(payload),
             'fraction_identical_of_aligned':matched/n_query if n_query else 0,
             'catalog_alignment_start_0based':int(idxr[0][0]) if len(idxr) else None,
             'catalog_alignment_end_0based':int(idxr[-1][1]) if len(idxr) else None}
        if best is None or (matched,n_query)>(best['identical_MPRA_bases'],best['aligned_MPRA_bases']):best=res
    best['interpretation']='sequence_correspondence_only_not_author_coordinate_confirmation'
    alignment_rows.append(best)
savecsv('S5_20_ConstructAllele_Sequence_Alignment.csv',alignment_rows)

# Quantify cut heterogeneity as range of cut means within HAR (not independent samples).
het=[]
for h in CAND:
    cc=[x for x in six if x['HAR']==h]
    vals=np.array([x['mean_HminusC_log2'] for x in cc])
    het.append({'HAR':h,'n_distinct_cuts':len(cc),'means_by_cut':'; '.join(f"{x['cut']}={x['mean_HminusC_log2']:+.6f}" for x in cc),
                'range_across_cut_means':float(vals.max()-vals.min()),
                'opposing_cut_directions':bool(np.any(vals>0) and np.any(vals<0)),
                'stage20B_BH_q':1.0,
                'interpretation':'Cut-specific; do not pool as independent biological replicates'})
savecsv('S6_FiveCandidate_Cut_Heterogeneity.csv',het)

# Sensitivity outcomes across cuts, not candidate assertions.
from collections import Counter
nominal=sum(x['exploratory_t_p'] is not None and x['exploratory_t_p']<.05 for x in allcuts)
fdr=sum(x['exploratory_BH_q_114_cuts']<.05 for x in allcuts)
consistent=sum(x['same_direction_3_of_3'] for x in allcuts)
diff_r1=[abs(x['delta_vs_original']) for x in rawrows if x['min_DNA_reads_per_barcode']==1 and x['aggregate_pseudocount']==1]
qc_report={
 'Stage21C_count_file_parser':'Three uploaded count files headerless; read as data header=None',
 'source_design_SHA256':hashlib.sha256(DESIGN.read_bytes()).hexdigest(),
 'stage21C_zip_SHA256':hashlib.sha256(ZIP.read_bytes()).hexdigest(),
 'recomputed_R1_RNA_rows':Rn,'recomputed_R1_DNA_rows':Dn,
 'recomputed_R1_RNA_counts':Rtotal,'recomputed_R1_DNA_counts':Dtotal,
 'source_corrected_R3_RNA_counts':rawcounts(ROOT/'GSM3016034_RNA.tsv.gz')[1],
 'three_culture_114_cut_summary_rows':len(allcuts),
 'all_cuts_3_of_3_direction':consistent,
 'cut_level_nominal_t_p_lt_0.05':nominal,
 'cut_level_BH_q_lt_0.05':fdr,
 'warning_t_tests':'Exploratory n=3 culture means, not donor replication, t-assumptions cannot be tested',
 'raw_R1_DNApositive_zeroRNA_sensitivity_max_absolute_delta_from_21C':max(diff_r1),
 'raw_R1_DNApositive_zeroRNA_sensitivity_median_absolute_delta_from_21C':float(np.median(diff_r1)),
 'three_rep_aggregate_pseudocount_max_error_at_zero':basemax,
 'candidate_allele_construct_matches':len(alignment_rows),
 'remaining_raw_count_files_not_supplied':['GSM3016033_RNA','GSM3016051_DNA','GSM3016052_DNA'],
 'author_HAR_coordinate_crosswalk_confirmed':False,
 'full_zero_RNA_3rep_analysis_completed':False,
 'raw_R1_only_sensitivity_completed':True,
 'no_experimental_new_TF_occupancy_or_gene_target_validation':True
}
(OUT/'S0_QC_and_Scope.json').write_text(json.dumps(qc_report,indent=2))
print(json.dumps(qc_report,indent=2))
print('\nCANDIDATE SUMMARY')
for x in six:
 print(x['HAR'],x['cut'],f"mean={x['mean_HminusC_log2']:+.6f}",f"p={x['exploratory_t_p']:.5g}",f"q={x['exploratory_BH_q_114_cuts']:.5g}",f"R1zero_delta={bcore[(x['HAR'],x['cut'])]['delta_vs_original']:+.5f}")
print('\nALIGNMENT')
for x in alignment_rows:print(x['HAR'],x['allele'],x['cut'],x['identical_MPRA_bases'],x['aligned_MPRA_bases'],x['fraction_identical_of_payload'])
