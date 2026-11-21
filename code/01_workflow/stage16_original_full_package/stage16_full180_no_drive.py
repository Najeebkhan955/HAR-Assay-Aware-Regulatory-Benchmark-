"""Stage 16 no-Google-Drive HAR comparative motif preflight.

Public source is a third-party 2026 HAR catalog mirroring sequence/orthologs.
Uses unique hg38 interval overlap to link a screened HAR alias to a catalog row.
This does NOT assert equivalence of original HAR names or validate TF binding.
"""
import argparse, csv, json, math, collections, hashlib, re, pathlib
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TF_MOTIFS={'SOX2':'MA0143.5','PAX6':'MA0069.1','REST':'MA0138.3','SOX11':'MA0869.3','CTCF':'MA0139.2','NFIA':'MA0670.2'}
THRESHOLDS=(0.80,0.85,0.90)

def read_csv(path):
 with open(path,encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
def write_csv(path,rows,headers=None):
 if not rows: return
 with open(path,'w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=headers or list(rows[0]),extrasaction='ignore');w.writeheader();w.writerows(rows)
def read_bed(path):
 out={}
 for l in open(path,encoding='utf-8'):
  if not l.strip() or l.startswith('#'):continue
  p=l.rstrip().split('\t')
  if len(p)<4:continue
  out[p[3]]={'chrom':p[0],'start':int(p[1]),'end':int(p[2]),'id':p[3]}
 return out
def read_fasta(path):
 records={};k=None;buf=[]
 for line in open(path,encoding='utf-8'):
  line=line.strip()
  if not line:continue
  if line.startswith('>'):
   if k: records[k]=''.join(buf).upper()
   k=line[1:].split()[0];buf=[]
  else:buf.append(line)
 if k:records[k]=''.join(buf).upper()
 return records

def load_pwm(path):
 mats={}; want=set(TF_MOTIFS.values());lines=open(path,encoding='utf-8').read().splitlines();i=0
 while i<len(lines):
  if lines[i].startswith('MOTIF '):
   mid=lines[i].split()[1]
   if mid in want:
    j=i+1
    while j<len(lines) and not lines[j].startswith('letter-probability matrix:') and not lines[j].startswith('MOTIF '):j+=1
    if j<len(lines) and lines[j].startswith('letter-probability matrix:'):
     m=re.search(r'w=\s*(\d+)',lines[j]);w=int(m.group(1)) if m else 0
     arr=[list(map(float,lines[j+q+1].strip().split())) for q in range(w)]
     if len(arr)==w and all(len(x)==4 for x in arr): mats[mid]=arr
     i=j+w
  i+=1
 if set(mats)!=want:raise ValueError('Missing pinned JASPAR motif IDs: '+str(want-set(mats)))
 return mats

def pwm_scan(seq,rows):
 """Exploratory relative PWM score, not FIMO p-value or binding occupancy.

 p'= (p+0.001)/(1+0.004); equal background 0.25. Both strands.
 Each strand's windows scanned; N windows discarded. Palindromic sites may count twice.
 """
 pw=np.log2((np.asarray(rows,dtype=float)+0.001)/1.004/0.25)
 lo=pw.min(axis=1).sum();hi=pw.max(axis=1).sum();k=len(rows)
 best=-math.inf;cnt={t:0 for t in THRESHOLDS}
 for s in (seq,str(Seq(seq).reverse_complement())):
  for i in range(len(s)-k+1):
   x=s[i:i+k]
   if any(a not in 'ACGT' for a in x):continue
   score=sum(float(pw[j,'ACGT'.index(c)]) for j,c in enumerate(x))
   rel=(score-lo)/(hi-lo)
   best=max(best,rel)
   for t in THRESHOLDS:
    if rel>=t:cnt[t]+=1
 return (best if math.isfinite(best) else None),cnt

def alignment_quality(a,b):
 """Align both orientations, returning relative identity on aligned non-gap columns.
 No phylogenetic ancestral states are inferred."""
 aln=PairwiseAligner();aln.mode='global';aln.match_score=2;aln.mismatch_score=-1;aln.open_gap_score=-4;aln.extend_gap_score=-1
 fwd=aln.align(a,b)[0];rc=str(Seq(b).reverse_complement());rev=aln.align(a,rc)[0]
 choose=rev if rev.score>fwd.score else fwd
 try:
  co=choose.counts();denom=co.identities+co.mismatches+co.gaps
  return ('reverse_complement' if choose is rev else 'forward',round(co.identities/denom,5) if denom else None)
 except Exception:
  return ('reverse_complement' if choose is rev else 'forward',None)

def main():
 p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--bait',required=True);p.add_argument('--screen',required=True);p.add_argument('--output',required=True)
 args=p.parse_args();src=pathlib.Path(args.input);dst=pathlib.Path(args.output);dst.mkdir(parents=True,exist_ok=True)
 bed_h=read_bed(src/'hars_hg38.bed'); bed_c=read_bed(src/'hars_panTro6.bed');fa_h=read_fasta(src/'human_hars.fasta');fa_c=read_fasta(src/'chimp_hars.fasta');pwms=load_pwm(src/'meme_file.txt')
 assert len(bed_h)>=3200 and len(bed_c)>=3200 and len(fa_h)>=3200 and len(fa_c)>=3200,'Public source dataset incomplete'

 # FASTA headers encode exact 0-based half-open chromosome coordinates.
 def fasta_for(bed,fasta,vid):
  b=bed[vid];k=f"{b['chrom']}:{b['start']}-{b['end']}";s=fasta.get(k)
  return s if s and len(s)==b['end']-b['start'] else None
 by_har=collections.defaultdict(list)
 for r in read_csv(args.bait):
  if r['context']!='hNSCs' or r['assembly']!='hg38':continue
  by_har[r['screen_HAR']].append((r['chrom'],int(r['start_0based']),int(r['end_exclusive'])))
 # Deduplicate: fragment may be repeated after combining aliases.
 by_har={k:list(set(v)) for k,v in by_har.items()}
 screen=read_csv(args.screen); allmap=[]; paired=[]; motif=[]
 for row in screen:
  name=row['screen_HAR'];frags=by_har.get(name,[])
  candidates=[]
  for vid,r in bed_h.items():
   ov=sum(max(0,min(r['end'],e)-max(r['start'],s)) for c,s,e in frags if c==r['chrom'])
   if ov>0:
    frac=ov/(r['end']-r['start'])
    candidates.append((vid,ov,frac))
  eligible=[z for z in candidates if z[2]>=.90]
  base={'screen_HAR':name,'number_named_hNSC_bait_fragments':len(frags),'any_overlapping_catalog_records':len(candidates),'candidate_records_covering_GE90pct':len(eligible)}
  if not frags: base['status']='no_named_hNSC_bait';allmap.append(base);continue
  if len(eligible)!=1:base['status']='ambiguous' if len(eligible)>1 else 'no_catalog_element_covered_GE90pct';allmap.append(base);continue
  vid,ov,frac=eligible[0];base.update({'catalog_id_provisional':vid,'overlap_bp':ov,'overlap_fraction':round(frac,5)})
  if vid not in bed_c:base['status']='catalog_chimp_liftover_unmapped';allmap.append(base);continue
  hs=fasta_for(bed_h,fa_h,vid);cs=fasta_for(bed_c,fa_c,vid)
  if not hs or not cs:base['status']='sequence_length_or_header_failure';allmap.append(base);continue
  orientation,ident=alignment_quality(hs,cs)
  if ident is None or ident<0.7:base['status']='ortholog_alignment_identity_below_70pct_or_unavailable';allmap.append(base);continue
  base['status']='passed_provisional_catalog_match';base['orientation']=orientation;base['global_alignment_identity']=ident;allmap.append(base)
  cp={'screen_HAR':name,'catalog_id_provisional':vid,'h_chr':bed_h[vid]['chrom'],'h_start':bed_h[vid]['start'],'h_end':bed_h[vid]['end'],'c_chr':bed_c[vid]['chrom'],'c_start':bed_c[vid]['start'],'c_end':bed_c[vid]['end'],'human_length_bp':len(hs),'chimp_length_bp':len(cs),'human_GC_fraction':round(sum(x in 'GC' for x in hs)/len(hs),4),'chimp_GC_fraction':round(sum(x in 'GC' for x in cs)/len(cs),4),'alignment_orientation':orientation,'global_alignment_identity':ident,'interpretation':'provisional catalog ortholog pair; not experimental binding validation'};paired.append(cp)
  for tf,mid in TF_MOTIFS.items():
   sh,ch=pwm_scan(hs,pwms[mid]);sc,cc=pwm_scan(cs,pwms[mid])
   z={'screen_HAR':name,'catalog_id_provisional':vid,'TF':tf,'motif':mid,'human_best_relative_score':round(sh,5) if sh is not None else None,'chimp_best_relative_score':round(sc,5) if sc is not None else None,'delta_human_minus_chimp_best_relative_score':round(sh-sc,5) if sh is not None and sc is not None else None}
   for thr in THRESHOLDS:
    label=str(thr).replace('.','p');z[f'human_hits_{label}']=ch[thr];z[f'chimp_hits_{label}']=cc[thr];z[f'delta_hits_{label}']=ch[thr]-cc[thr]
   z['inference']='exploratory relative PWM; best site may not be homologous; occupancy unmeasured'
   motif.append(z)
 write_csv(dst/'S1_All180_Sequence_Mapping_Audit.csv',allmap)
 write_csv(dst/'S2_Aligned_Ortholog_Sequence_QC.csv',paired)
 write_csv(dst/'S3_SixTF_Motif_Score_Comparisons.csv',motif)
 status=collections.Counter(z['status'] for z in allmap)
 json.dump({'n_screened':len(screen),'mapping_status_counts':dict(status),'n_paired_orthologs':len(paired),'n_motif_comparisons':len(motif),'JASPAR_motif_version':'2024 mirrored at pinned third-party commit','limitations':['genomic restriction bait is broader than HAR','public catalog IDs are not original nomenclature','PWM threshold sensitivity','best-region motif score not homologous site comparison','no direct TF occupancy or causal target validation']},open(dst/'STATUS.json','w'),indent=2)
 if motif:
  hnames=[r['screen_HAR'] for r in paired];tfnames=list(TF_MOTIFS)
  M=np.full((len(hnames),len(tfnames)),np.nan)
  for r in motif:M[hnames.index(r['screen_HAR']),tfnames.index(r['TF'])]=r['delta_human_minus_chimp_best_relative_score']
  ix=np.argsort(np.nanmax(np.abs(M),axis=1))[::-1][:25];M=M[ix];labels=[hnames[i] for i in ix]
  plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','font.size':9})
  fig,ax=plt.subplots(figsize=(9,6.5));mx=max(.05,np.nanmax(abs(M)));im=ax.imshow(M,cmap='BrBG',vmin=-mx,vmax=mx,aspect='auto');ax.set_yticks(range(len(labels)),labels);ax.set_xticks(range(len(tfnames)),tfnames);ax.set_title('Exploratory human–chimp PWM score changes (top 25 loci)',loc='left',fontweight='bold');fig.colorbar(im,ax=ax,label='Δ best normalized PWM score');fig.tight_layout();fig.savefig(dst/'F1_Sequence_Only_Motif_Score_Sensitivity.png',dpi=180);plt.close(fig)
 print('STATUS:',json.dumps(dict(status)),'n_motif',len(motif))

if __name__=='__main__':main()
