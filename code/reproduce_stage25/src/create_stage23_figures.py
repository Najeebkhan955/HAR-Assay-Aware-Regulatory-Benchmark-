"""Generate exact-data graphical mind map and figures; clean white publication style."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from graphviz import Digraph
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results';F=ROOT/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','savefig.dpi':300})
NAVY='#172B4D';BLUE='#2374AB';TEAL='#0C877A';ORANGE='#DB7D2F';GREY='#8694A5';RED='#BF5B5B';PALE='#F3F7FA'
S=json.load(open(R/'S0_Stage23_Summary.json'));C=S['counts']

# Research decision mind map; source .gv is editable/reproducible.
g=Digraph('HAR23',format='svg',graph_attr={'rankdir':'LR','bgcolor':'white','pad':'0.3','nodesep':'0.42','ranksep':'0.95','splines':'spline','fontname':'DejaVu Sans'},node_attr={'shape':'box','style':'rounded,filled','color':'#D4DFEA','penwidth':'1','fillcolor':'#F4F8FB','fontname':'DejaVu Sans','fontsize':'12','fontcolor':NAVY,'margin':'0.19,0.14'},edge_attr={'color':'#91A4B8','penwidth':'1.8','arrowsize':'0.75'})
g.node('core','HAR regulatory evidence\nASSAY-AWARE BENCHMARK',fillcolor=NAVY,fontcolor='white',fontsize='15',color=NAVY,penwidth='2')
g.node('pert','Perturb-seq\n180 HARs; one iPSC background\nTwo-guide agreement weak')
g.node('map','Sequence orthology\n152 provisional, 28 unresolved\nOriginal BED identity pending')
g.node('motif','Six TF motifs\n16 nonzero / 152\n0 context-null FDR discoveries')
g.node('report','Human vs chimp INSERT MPRA\n75 named, 66 paired; 114 cuts\n3 cultures, HS1-11 N2')
g.node('fragment','Within-HAR fragment discordance\n20/28 multi-cut HARs opposite signs\n11 retain opposing 3-culture signs',fillcolor='#EAF5F5',color=TEAL)
g.node('transfer','Cross-assay transfer test\n58 shared HARs; 5 motif nonzero\nNo reporter effect advantage; p=0.62',fillcolor='#EFF4FC',color=BLUE)
g.node('genome','Hard validation gates\nAuthor BED + cut-specific motif overlap\nNo confirmed HAR-TF-gene causality',fillcolor='#FFF2EF',color=ORANGE)
g.node('wtc','NEXT: WTc N2 independent human line\nThree paired RNA/DNA cultures\nNeeds one live GEO/Colab run',fillcolor='#FFF7EA',color=ORANGE)
g.node('repo','OPEN RESEARCH RESOURCE\nAudit tables + tests + CSVs\nGitHub-ready, source hashes',fillcolor='#F0F7EC',color=TEAL)
for a,b in [('core','pert'),('core','map'),('core','motif'),('core','report'),('report','fragment'),('motif','transfer'),('report','transfer'),('map','genome'),('transfer','genome'),('transfer','wtc'),('wtc','repo'),('genome','repo')]:g.edge(a,b)
g.save(str(F/'F1_Research_Novelty_Mindmap.gv'))
g.render(filename='F1_Research_Novelty_Mindmap',directory=str(F),cleanup=True,format='svg')
g.render(filename='F1_Research_Novelty_Mindmap',directory=str(F),cleanup=True,format='png')

# Figure 2: honest coverage intersection (no misleading waterfall for nonnested subsets).
fig,ax=plt.subplots(figsize=(9.5,4.4));fig.subplots_adjust(left=.3,right=.95,top=.79,bottom=.16)
parts=[('Mapped + paired MPRA',58,TEAL),('Mapped; reporter absent',94,BLUE),('Reporter; mapping absent',8,ORANGE),('Neither provisional map nor reporter',20,GREY)]
y=np.arange(len(parts))[::-1]
ax.barh(y,[t[1] for t in parts],height=.63,color=[t[2] for t in parts],edgecolor='none')
ax.set_yticks(y,[t[0] for t in parts],fontsize=10);ax.set_xlim(0,112);ax.set_xlabel('HAR loci (denominator: original 180 screened HARs)');ax.set_title('Assay-compatible coverage, not a sequential funnel',loc='left',pad=23,fontweight='bold',fontsize=14,color=NAVY)
for yi,(_,n,_) in zip(y,parts):ax.text(n+1.2,yi,str(n),va='center',fontweight='bold',fontsize=12)
fig.text(.3,.89,'152/180 provisional orthologs   |   75/180 design names   |   66/180 paired MPRA loci',fontsize=10,color='#567185')
fig.text(.3,.04,'Unmapped and unassayed observations are NA, not negative biological results.',fontsize=9,color=GREY)
for ext in ('png','svg'):fig.savefig(F/f'F2_180HAR_Assay_Compatibility.{ext}',bbox_inches='tight')
plt.close(fig)

# Figure 3 reproducible MPRA contrast vs culture, per HAR cut.
q=pd.read_csv(R/'S7_114_Cut_3Culture_Summary.csv');a=pd.read_csv(R/'S1_180HAR_EvidenceCompatibility.csv');pos=set(a[a.nonzero_motif_score].HAR)
fig,axs=plt.subplots(1,3,figsize=(13,4),sharex=True,sharey=True)
for ax, (x,y,pr) in zip(axs,[('culture_1','culture_2',S['replicate_correlations'][0]['pearson_r']),('culture_1','culture_3',S['replicate_correlations'][1]['pearson_r']),('culture_2','culture_3',S['replicate_correlations'][2]['pearson_r'])]):
    rest=q[~q.HAR.isin(pos)];sel=q[q.HAR.isin(pos)]
    ax.scatter(rest[x],rest[y],s=22,c=BLUE,alpha=.52,linewidths=0)
    ax.scatter(sel[x],sel[y],s=42,c=ORANGE,alpha=.85,linewidths=.5,edgecolors='white')
    ax.axhline(0,color=GREY,lw=.6);ax.axvline(0,color=GREY,lw=.6)
    ax.plot([-.65,.65],[-.65,.65],ls='--',lw=.75,color=GREY)
    ax.set_xlim(-.62,.62);ax.set_ylim(-.62,.62)
    ax.set_xlabel(x.replace('culture_','Culture ')+' H−C log2 activity');ax.set_ylabel(y.replace('culture_','Culture ')+' H−C log2 activity')
    ax.text(.04,.94,f'Pearson r = {pr:.3f}',transform=ax.transAxes,ha='left',va='top',color=NAVY,fontsize=10,fontweight='bold')
fig.suptitle('Three-culture reproducibility of allele-contrast measurements',y=1.02,color=NAVY,fontsize=14,fontweight='bold')
fig.text(.03,-.04,'Each point = one HAR × reporter cut; 114 cuts across 66 loci. Orange = loci with nonzero six-TF PWM statistic.\nCulture replicates share the same HS1-11 iPSC background; correlations do not confer independent-donor validation.',fontsize=9,color='#637282')
for ext in ('png','svg'):fig.savefig(F/f'F3_MPRA_ThreeCulture_ReplicateAgreement.{ext}',bbox_inches='tight')
plt.close(fig)

# Figure 4 intra-HAR cut-specific effects. Loci with both signs, not effect-size CIs.
multi=pd.read_csv(R/'S6_MultiCut_Heterogeneity.csv');d=multi[multi.opposing_cut_directions].sort_values('cross_cut_range_log2',ascending=True)
fig,ax=plt.subplots(figsize=(10,8.2));y=np.arange(len(d))
ax.hlines(y,d.min_cut_mean,d.max_cut_mean,colors=[TEAL if stable else GREY for stable in d.opposing_cuts_each_consistent_3cultures],lw=3,alpha=.9)
ax.scatter(d.min_cut_mean,y,s=35,marker='<',c=BLUE,zorder=3);ax.scatter(d.max_cut_mean,y,s=35,marker='>',c=ORANGE,zorder=3)
ax.axvline(0,color='#9BA8B3',ls=':',lw=1)
ax.set_yticks(y,d.HAR);ax.set_xlabel('Three-culture mean human − chimp log2 activity, across different 171-bp cuts')
ax.set_title('One HAR label can contain opposing fragment-level allele effects',loc='left',fontweight='bold',color=NAVY,fontsize=13,pad=17)
fig.text(.13,.025,'20/28 multi-cut loci span zero; 11/28 contain cuts that each keep opposite directions in all three cultures.\nEndpoints are minimum/maximum CUT EFFECTS (not a confidence interval).',fontsize=9,color='#5D7082')
fig.subplots_adjust(left=.24,bottom=.13,right=.97,top=.89)
for ext in ('png','svg'):fig.savefig(F/f'F4_MultiCut_AlleleEffect_Heterogeneity.{ext}',bbox_inches='tight')
plt.close(fig)

# Figure 5 five priority loci with all 10 distinct cut effects across 3 cultures.
focus=q[q.HAR.isin({'2xHAR.70','2xHAR.310','2xHAR.234','2xHAR.377','2xHAR.3'})].copy().sort_values(['HAR','cut']);y=np.arange(len(focus))[::-1]
fig,ax=plt.subplots(figsize=(9,5.2));
for i,(_,r) in enumerate(focus.iterrows()):
    yi=y[i];v=[r.culture_1,r.culture_2,r.culture_3]
    ax.plot([min(v),max(v)],[yi,yi],lw=2,color='#BCCBD6',zorder=1)
    ax.scatter(v,[yi]*3,c=[BLUE,TEAL,ORANGE],s=19,zorder=3)
    ax.scatter(r.mean_HminusC,yi,c=NAVY,s=65,marker='D',edgecolor='white',linewidth=.4,zorder=4)
ax.axvline(0,ls=':',lw=1,color=GREY)
ax.set_yticks(y,[f'{rr.HAR}  {rr.cut}' for _,rr in focus.iterrows()]);ax.set_xlabel('Human − chimp log2 reporter activity')
ax.set_title('Five nominated HARs: keep all 10 assayed fragment identities',loc='left',color=NAVY,fontweight='bold',fontsize=13,pad=16)
fig.text(.2,.02,'Small colored dots = three same-line cultures; dark diamonds = means. No formal FDR-positive candidate among these five.',fontsize=9,color='#5D7082')
fig.subplots_adjust(left=.26,bottom=.17,top=.88)
for ext in ('png','svg'):fig.savefig(F/f'F5_FiveHAR_10Cut_EffectDirections.{ext}',bbox_inches='tight')
plt.close(fig)
print('Figure package written to',F)
