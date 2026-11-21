from pathlib import Path
import csv, json, re, collections, statistics, math, hashlib, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT=Path('/mnt/data'); OUT=ROOT/'HAR_Stage16_Comparative_TF';OUT.mkdir(exist_ok=True)
S15=ROOT/'HAR_Stage15_Reference_Gated_Benchmark'/'S4_180HAR_Benchmark_Eligibility_and_Evidence.csv'
URL='https://github.com/athenamarou/HAR-TFBS-Project/tree/a9251f39976bbc15df1d8f0dfedaf9786819c298/data'

def rd(path,sep=','):
 with open(path,newline='',encoding='utf-8') as f:return list(csv.DictReader(f,delimiter=sep))
def wr(name,rows,fields=None):
 if not rows: return
 if not fields:fields=list(rows[0])
 with open(OUT/name,'w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)

rows=rd(S15); names={r['HAR'] for r in rows}
# Pal chromatin bait fragments, explicitly NOT the precise HAR regulatory intervals.
by_name=collections.defaultdict(dict); contexts={}
for ctx,fn in [('hNSCs','HAR_interactions_hNSCs.txt'),('human_neurons','HAR_interactions_human_neurons.txt'),('cNSCs','HAR_interactions_cNSCs.txt')]:
 dat=rd(ROOT/fn,sep='\t')
 for r in dat:
  toks={re.sub(r'_(sec|tert)$','',s.strip()) for s in r['Bait_name'].split(',')}
  for n in toks & names:
   # deduplicate; CHiCAGO scores are contact scores, not sequence or binding evidence.
   loc=(r['Bait_chr'],int(r['Bait_start']),int(r['Bait_end']),r['Bait_name'])
   by_name[(n,ctx)][loc]=True
print('Number of exact alias-matched HAR names per context:',{ctx:sum(bool(by_name.get((n,ctx))) for n in names) for ctx in ['hNSCs','human_neurons','cNSCs']})

bait=[]
for n in sorted(names):
 for ctx in ['hNSCs','human_neurons','cNSCs']:
  for chrom,start,end,rawname in sorted(by_name.get((n,ctx),{})):
   bait.append({'screen_HAR':n,'context':ctx,'assembly':'panTro6' if ctx=='cNSCs' else 'hg38','chrom':chrom,'start_0based':start,'end_exclusive':end,'bait_raw_alias':rawname,'coordinates_meaning':'capture_HiC_restriction_fragment_NOT_precise_HAR'})
wr('S2_All_Pal_Bait_Fragments_for_Screen.csv',bait)
human_baits=[r for r in bait if r['context']=='hNSCs']
wr('S2b_Human_NSC_Bait_Fragments_Embedded_for_Colab.csv',human_baits)

# Raw, preliminary independently computed pilot results from publicly accessible source records.
# Record IDs HAR_# are local sequential catalog IDs, not the original nomenclature (e.g., HAR181).
pilot=[
 ('HAR181','HAR_459','chr11',114928473,114928663,185),
 ('HAR202','HAR_788','chr14',33576363,33576566,203),
 ('2xHAR.238','HAR_1521','chr2',121075207,121075564,357),
 ('2xHAR.13','HAR_1742','chr20',37851093,37851271,178),
 ('HACNS549','HAR_2551','chr6',19976197,19977166,961),
 ('HACNS238','HAR_1581','chr2',157785403,157785570,159),
 ('HAR67','HAR_675','chr13',71799450,71799565,115),
 ('HAR25','HAR_2301','chr4',181350707,181351729,1018),
 ('2xHAR.437','HAR_1028','chr16',7395670,7395836,166),
 ('2xHAR.446','HAR_491','chr11',131652951,131652976,25)
]
maprows=[]
for n,pid,chr,st,en,ov in pilot:
 r=next(x for x in rows if x['HAR']==n)
 maprows.append({'screen_HAR':n,'candidate_catalog_ID_not_original_HAR_name':pid,'human_chr_hg38':chr,'human_start_0based':st,'human_end_exclusive':en,'overlap_bp_with_named_hNSC_bait_fragments':ov,'overlap_fraction_of_catalog_element':round(ov/(en-st),4),'guide_min_cells_per_batch':int(r['min_guide_batch_cells']),'guide_min20':r['eligible_min20'],'guide_min50':r['eligible_min50'],'catalog_mapping':'unique_overlap_in_pilot_catalog; provisional; not independent positional confirmation','sequence_source':URL,'experimental_TF_binding_validated':False})
wr('S3_TenHAR_Provisional_Catalog_Matches.csv',maprows)

# Six TF expression real raw NTC counts from the Stage 11 RNA pseudobulk data (Stage 12 same-group reconciliation).
tfs=[('SOX2','MA0143.5',829,1181,484,104.65,182.49,77.28),('PAX6','MA0069.1',104,177,154,13.13,27.35,24.59),('REST','MA0138.3',570,541,440,71.96,83.59,70.25),('SOX11','MA0869.3',903,1276,614,114.00,197.17,98.03),('CTCF','MA0139.2',637,622,525,80.42,96.11,83.82),('NFIA','MA0670.2',12,25,15,1.51,3.86,2.39)]
wr('S4_NonTargetingControl_NeuralTF_Expression.csv',[{'TF':n,'JASPAR2024_motif_ID':mid,'rep1_raw_counts':a,'rep2_raw_counts':b,'rep3_raw_counts':c,'rep1_CPM':v1,'rep2_CPM':v2,'rep3_CPM':v3,'mean_CPM':round((v1+v2+v3)/3,3),'biological_note':'TF transcription detectable; does not establish DNA binding'} for n,mid,a,b,c,v1,v2,v3 in tfs])

# Source-observed PWM best scores with pseudocount 0.001, p_bg=0.25, both strands.
# From exact source data (public mirrored BED/FASTA + MEME) pulled with Github connector 2026-11-19.
DATA='''HAR181 HAR_459 SOX2 MA0143.5 190 190 0.7935 0.7935 0 0 0 0
HAR181 HAR_459 PAX6 MA0069.1 190 190 0.8358 0.7612 0.0746 0 0 0
HAR181 HAR_459 REST MA0138.3 190 190 0.6606 0.6606 0 0 0 0
HAR181 HAR_459 SOX11 MA0869.3 190 190 0.8355 0.8355 0 0 0 0
HAR181 HAR_459 CTCF MA0139.2 190 190 0.8412 0.8412 0 0 0 0
HAR181 HAR_459 NFIA MA0670.2 190 190 0.8333 0.8333 0 0 0 0
HAR202 HAR_788 SOX2 MA0143.5 203 203 0.9057 0.9057 0 1 1 0
HAR202 HAR_788 PAX6 MA0069.1 203 203 0.7699 0.7771 -0.0072 0 0 0
HAR202 HAR_788 REST MA0138.3 203 203 0.6577 0.6577 0 0 0 0
HAR202 HAR_788 SOX11 MA0869.3 203 203 0.8145 0.8145 0 0 0 0
HAR202 HAR_788 CTCF MA0139.2 203 203 0.7858 0.7858 0 0 0 0
HAR202 HAR_788 NFIA MA0670.2 203 203 0.8333 0.8333 0 0 0 0
2xHAR.238 HAR_1521 SOX2 MA0143.5 357 356 0.8183 0.861 -0.0427 0 0 0
2xHAR.238 HAR_1521 PAX6 MA0069.1 357 356 0.8129 0.831 -0.0181 0 0 0
2xHAR.238 HAR_1521 REST MA0138.3 357 356 0.7396 0.7013 0.0383 0 0 0
2xHAR.238 HAR_1521 SOX11 MA0869.3 357 356 0.8219 0.8219 0 0 0 0
2xHAR.238 HAR_1521 CTCF MA0139.2 357 356 0.796 0.796 0 0 0 0
2xHAR.238 HAR_1521 NFIA MA0670.2 357 356 0.8333 0.8333 0 0 0 0
2xHAR.13 HAR_1742 SOX2 MA0143.5 178 178 0.8072 0.8072 0 0 0 0
2xHAR.13 HAR_1742 PAX6 MA0069.1 178 178 0.8219 0.8219 0 0 0 0
2xHAR.13 HAR_1742 REST MA0138.3 178 178 0.6594 0.6201 0.0393 0 0 0
2xHAR.13 HAR_1742 SOX11 MA0869.3 178 178 0.7991 0.7991 0 0 0 0
2xHAR.13 HAR_1742 CTCF MA0139.2 178 178 0.8258 0.8412 -0.0154 0 0 0
2xHAR.13 HAR_1742 NFIA MA0670.2 178 178 0.8333 0.8333 0 0 0 0
HACNS549 HAR_2551 SOX2 MA0143.5 969 968 0.9057 0.9057 0 1 1 0
HACNS549 HAR_2551 PAX6 MA0069.1 969 968 0.822 0.822 0 0 0 0
HACNS549 HAR_2551 REST MA0138.3 969 968 0.7097 0.7097 0 0 0 0
HACNS549 HAR_2551 SOX11 MA0869.3 969 968 0.9403 0.9403 0 1 1 0
HACNS549 HAR_2551 CTCF MA0139.2 969 968 0.8139 0.8365 -0.0226 0 0 0
HACNS549 HAR_2551 NFIA MA0670.2 969 968 0.8333 0.8333 0 0 0 0
HACNS238 HAR_1581 SOX2 MA0143.5 167 167 0.9051 0.9051 0 2 2 0
HACNS238 HAR_1581 PAX6 MA0069.1 167 167 0.8064 0.8064 0 0 0 0
HACNS238 HAR_1581 REST MA0138.3 167 167 0.7983 0.7377 0.0606 0 0 0
HACNS238 HAR_1581 SOX11 MA0869.3 167 167 0.9403 0.9403 0 1 1 0
HACNS238 HAR_1581 CTCF MA0139.2 167 167 0.7086 0.7086 0 0 0 0
HACNS238 HAR_1581 NFIA MA0670.2 167 167 0.8333 0.8333 0 0 0 0
HAR67 HAR_675 SOX2 MA0143.5 115 115 0.9321 0.9321 0 1 1 0
HAR67 HAR_675 PAX6 MA0069.1 115 115 0.7517 0.7517 0 0 0 0
HAR67 HAR_675 REST MA0138.3 115 115 0.6599 0.6599 0 0 0 0
HAR67 HAR_675 SOX11 MA0869.3 115 115 0.7699 0.7699 0 0 0 0
HAR67 HAR_675 CTCF MA0139.2 115 115 0.7767 0.7767 0 0 0 0
HAR67 HAR_675 NFIA MA0670.2 115 115 0.8333 0.8333 0 0 0 0
HAR25 HAR_2301 SOX2 MA0143.5 1022 1021 0.9557 0.9557 0 1 1 0
HAR25 HAR_2301 PAX6 MA0069.1 1022 1021 0.8322 0.8322 0 0 0 0
HAR25 HAR_2301 REST MA0138.3 1022 1021 0.6872 0.6872 0 0 0 0
HAR25 HAR_2301 SOX11 MA0869.3 1022 1021 1 1 0 1 1 0
HAR25 HAR_2301 CTCF MA0139.2 1022 1021 0.7983 0.7983 0 0 0 0
HAR25 HAR_2301 NFIA MA0670.2 1022 1021 0.8333 0.8333 0 0 0 0
2xHAR.437 HAR_1028 SOX2 MA0143.5 166 166 0.7944 0.8061 -0.0117 0 0 0
2xHAR.437 HAR_1028 PAX6 MA0069.1 166 166 0.7774 0.7606 0.0168 0 0 0
2xHAR.437 HAR_1028 REST MA0138.3 166 166 0.6929 0.6929 0 0 0 0
2xHAR.437 HAR_1028 SOX11 MA0869.3 166 166 0.8561 0.8561 0 0 0 0
2xHAR.437 HAR_1028 CTCF MA0139.2 166 166 0.7714 0.7714 0 0 0 0
2xHAR.437 HAR_1028 NFIA MA0670.2 166 166 0.8333 0.6667 0.1666 0 0 0
2xHAR.446 HAR_491 SOX2 MA0143.5 25 25 0.7007 0.7007 0 0 0 0
2xHAR.446 HAR_491 PAX6 MA0069.1 25 25 0.5785 0.5764 0.0021 0 0 0
2xHAR.446 HAR_491 REST MA0138.3 25 25 0.5001 0.4965 0.0036 0 0 0
2xHAR.446 HAR_491 SOX11 MA0869.3 25 25 0.6923 0.6923 0 0 0 0
2xHAR.446 HAR_491 CTCF MA0139.2 25 25 0.6331 0.6724 -0.0393 0 0 0
2xHAR.446 HAR_491 NFIA MA0670.2 25 25 0.6667 0.6667 0 0 0 0'''
allmot=[]
for line in DATA.splitlines():
 p=line.split();n,cat,tf,mid=p[:4];hl,cl=map(int,p[4:6]);x,y,d=map(float,p[6:9]);nh,nc,dn=map(int,p[9:12]);allmot.append({'screen_HAR':n,'catalog_ID':cat,'TF':tf,'JASPAR2024_MOTIF':mid,'human_length_bp':hl,'chimp_length_bp':cl,'human_best_relative_score':x,'chimp_best_relative_score':y,'delta_human_minus_chimp_best_score':d,'human_hits_relscore_GE_0p90':nh,'chimp_hits_relscore_GE_0p90':nc,'difference_0p90_hit_counts':dn,'site_homology_status':'not_established; best-of-region score only','binding_validation':'none; in_silico_PWM'} )
wr('S5_TenHAR_SixTF_PWM_Exploratory_Pilot.csv',allmot)

out=[]
for r in rows:
 n=r['HAR']; h=by_name.get((n,'hNSCs'),{}); hu=by_name.get((n,'human_neurons'),{}); c=by_name.get((n,'cNSCs'),{})
 out.append({'screen_HAR':n,'HAR_family':'HACNS' if n.startswith('HACNS') else '2xHAR' if n.startswith('2xHAR') else 'HAR','screen_n_guides':int(r['n_guides']),'min_guide_batch_cells':int(r['min_guide_batch_cells']),'two_guide_complete':r['eligible_2guides_3batches'],'passes_min20_guide_batch':r['eligible_min20'],'passes_min50_guide_batch':r['eligible_min50'],'human_NSC_named_bait_intervals':len(h),'human_neuron_named_bait_intervals':len(hu),'chimp_NSC_named_bait_intervals':len(c),'all_3_context_aliases_available':bool(h and hu and c),'pilot_unique_bed_overlap':n in {z[0] for z in pilot},'evidence_status':'guide_QC_and_naming_only; no validated TF binding / HAR–gene relationship'})
wr('S1_All180_HAR_Sequence_Mapping_Readiness.csv',out)

source=[
 {'source':'Stage 12/15 original processed guide & HAR pseudobulk','URL':'provided by user in conversation','assembly':'human hg38','relevance':'180 HAR guide coverage and RNA counts','limitations':'three batches from a single cell line; no independent donors'},
 {'source':'Pal et al 2025 Dryad chromatin contact BEDPE-derived table','URL':'https://doi.org/10.5061/dryad.c59zw3rkj','assembly':'human hg38; chimp panTro6','relevance':'bait-fragment intervals and chromatin-context presence','limitations':'bait is restriction fragment, not exact sequence boundary; contact not causal gene assignment'},
 {'source':'Marougka 2026 public sequence catalog (mirrors public genomes and Cui et al HAR intervals)','URL':URL,'assembly':'human hg38; chimp panTro6','relevance':'provisional overlap-derived human and chimp sequences; motif input','limitations':'third-party generated IDs HAR_1 etc are not original HAR accession names; must map coordinates and verify ortholog orientation'},
 {'source':'JASPAR2024 vertebrate CORE MEME mirrored in Marougka repository','URL':URL+'/meme_file.txt','assembly':'N/A','relevance':'6 selected motif matrices for proof of concept','limitations':'PWM scores cannot establish TF occupancy; version pinned to mirrored input'},
 {'source':'JASPAR current official downloads','URL':'https://jaspar2026.elixir.no/downloads/','assembly':'N/A','relevance':'independent motif source for later sensitivity analysis','limitations':'2026 motifs not necessarily equivalent to 2024 motifs'},
 {'source':'Noble et al open code','URL':'https://github.com/NoonanLab/Noble_et_al_2024','assembly':'human hg38','relevance':'Perturb-seq perturbation and statistical design','limitations':'published DEG effect table not retrieved; screen naming is not exact coordinate key'},
]
wr('S6_Independent_Source_Register.csv',source)

# Final visualizations: original QC and newly computed pilot motif scores.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
fig,ax=plt.subplots(figsize=(10.4,4.6))
categories=['2 guides in 3 batches','≥20 cells / guide / batch','≥50 cells / guide / batch','hNSC bait alias present','all 3 context bait aliases']
vals=[sum(x['two_guide_complete']=='True' for x in out),sum(x['passes_min20_guide_batch']=='True' for x in out),sum(x['passes_min50_guide_batch']=='True' for x in out),sum(x['human_NSC_named_bait_intervals']>0 for x in out),sum(x['all_3_context_aliases_available'] for x in out)]
y=np.arange(len(vals));ax.barh(y,vals,color=['#006F83','#087E87','#5996A5','#658F83','#9FB7BC'],height=.64)
for i,v in enumerate(vals):ax.text(v+2,i,f'{v}/180',va='center',fontweight='semibold',fontsize=11)
ax.set_yticks(y,labels=categories);ax.invert_yaxis();ax.set_xlim(0,198);ax.set_xlabel('HARs meeting descriptive preflight criterion')
ax.set_title('Screen-wide sequence-analysis preflight: what is ready?',fontsize=14,fontweight='bold',loc='left',pad=18)
ax.grid(axis='x',alpha=.16);ax.set_axisbelow(True)
fig.tight_layout();fig.savefig(OUT/'F1_180HAR_Sequence_Mapping_Preflight.png',dpi=230,bbox_inches='tight');plt.close(fig)

har_order=[x[0] for x in pilot];tf_order=[x[0] for x in tfs]
M=np.array([[next(x['delta_human_minus_chimp_best_score'] for x in allmot if x['screen_HAR']==h and x['TF']==t) for t in tf_order] for h in har_order])
fig,ax=plt.subplots(figsize=(10.8,5.7));mx=max(.08,np.max(abs(M)))
im=ax.imshow(M,cmap='BrBG',norm=TwoSlopeNorm(vmin=-mx,vcenter=0,vmax=mx),aspect='auto')
ax.set_yticks(range(len(har_order)),har_order);ax.set_xticks(range(len(tf_order)),tf_order)
for i in range(M.shape[0]):
 for j in range(M.shape[1]):
  v=M[i,j];ax.text(j,i, '0.000' if abs(v)<.00005 else f'{v:+.3f}',ha='center',va='center',fontsize=9,color='black' if abs(v)<mx*.55 else '#222222')
ax.set_title('Human − chimp motif-score shifts | 10 uniquely overlapping HAR loci',fontweight='bold',loc='left',fontsize=13,pad=14)
cb=fig.colorbar(im,ax=ax,shrink=.8,pad=.02);cb.set_label('Δ best normalized PWM score (descriptive)')
ax.set_xticks(np.arange(-.5,len(tf_order),1),minor=True);ax.set_yticks(np.arange(-.5,len(har_order),1),minor=True);ax.grid(which='minor',linewidth=.8,color='#dddddd');ax.tick_params(which='minor',bottom=False,left=False)
fig.text(.12,.03,'In-silico sequence preferences only • JASPAR2024 mirrored motifs • 0.90 relative-score site-count differences = 0 for this six-TF pilot',fontsize=9,color='#4e5963')
fig.tight_layout(rect=(0,.05,1,1));fig.savefig(OUT/'F2_10HAR_SixTF_Motif_Score_Heatmap.png',dpi=230,bbox_inches='tight');plt.close(fig)

fig,ax=plt.subplots(figsize=(10,4.9));ovr=[float(x['overlap_fraction_of_catalog_element']) for x in maprows]
ax.barh(range(10),[100*v for v in ovr],color='#368F97',height=.64)
ax.set_yticks(range(10),har_order);ax.invert_yaxis();ax.set_xlim(0,106)
for i,v in enumerate(ovr):ax.text(100*v+0.7,i,f'{100*v:.1f}%',va='center',fontsize=9)
ax.set_xlabel('Catalog sequence covered by named hNSC bait fragments (%)');ax.set_title('Provisional sequence correspondence: unique genomic overlap',loc='left',fontweight='bold',pad=14)
fig.text(.13,.01,'NOT direct original HAR-ID equivalence; capture-HiC bait may span outside the HAR element.',fontsize=9,color='#4e5963')
fig.tight_layout(rect=(0,.04,1,1));fig.savefig(OUT/'F3_Provisional_Catalog_Overlap_Fractions.png',dpi=230,bbox_inches='tight');plt.close(fig)

print('Preflight metrics',dict(zip(categories,vals)))
print('Pilot absolute nonzero best-score shifts',sum(abs(x['delta_human_minus_chimp_best_score'])>.00005 for x in allmot),'of',len(allmot))
json.dump({'HAR_total':180,'bait_context_name_coverage':{ctx:sum(bool(by_name.get((n,ctx))) for n in names) for ctx in ['hNSCs','human_neurons','cNSCs']},'pilot_matched_catalog_count':len(pilot),'pilot_TF_gene_count':len(allmot),'nonzero_best_PWM_deltas':sum(abs(x['delta_human_minus_chimp_best_score'])>.00005 for x in allmot),'observed_0p90_hit_count_differences':sum(x['difference_0p90_hit_counts']!=0 for x in allmot),'evidence_level':'preliminary computational motif comparison; no causal TF binding/functional confirmation'},open(OUT/'S0_Stage16_Audit.json','w'),indent=2)
