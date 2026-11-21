from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
out=Path('/mnt/data/HAR_Stage19_Figures'); out.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':14,'axes.labelsize':11,'figure.facecolor':'white','axes.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white','axes.titleweight':'bold','axes.edgecolor':'#CBD5E1','xtick.color':'#455667','ytick.color':'#455667','text.color':'#203347','axes.labelcolor':'#203347'})
blue='#2779A7'; teal='#0E9488'; coral='#D96E59'; slate='#8793A1'; light='#E7F1F4'; purple='#7263A3'
def lab(ax, txt): ax.text(-0.11, 1.09,txt,transform=ax.transAxes,weight='bold',fontsize=15,color='#162F49',va='top')
def save(fig,stem):
 fig.savefig(out/(stem+'.png'),dpi=330,bbox_inches='tight',pad_inches=0.23); plt.close(fig)

# Figure 1 - Data universe and filter retention
fig,ax=plt.subplots(1,2,figsize=(13.8,5.2),gridspec_kw={'width_ratios':[1.12,1]})
fig.suptitle('Dataset provenance and guide-coverage constraints',fontsize=17,weight='bold',y=1.03)
A=ax[0]; A.barh(['PRISM-processed cells','Single-guide selected cells','Non-targeting selected cells'],[87896,63335,1027],color=[slate,blue,teal],height=.53)
A.invert_yaxis(); A.set_xlim(0,100000); A.set_xlabel('Cells (count)'); A.grid(axis='x',alpha=.13)
for i,v in enumerate([87896,63335,1027]): A.text(v+1500,i,f'{v:,}',va='center',fontsize=11,weight='bold')
A.text(0.01,-.24,'19,698 genes  |  1,110 guide-by-batch pseudobulks  |  3 batches / 1 iPSC background',transform=A.transAxes,fontsize=9.1,color='#5A6878')
lab(A,'A')
B=ax[1]; labels=['Two-guide HARs','Minimum 20 cells','Minimum 50 cells']; val=[179,145,45];
B.barh(labels,val,color=[slate,blue,teal],height=.52); B.invert_yaxis(); B.set_xlim(0,198);B.set_xlabel('HARs satisfying coverage rule');B.grid(axis='x',alpha=.13)
for i,v in enumerate(val): B.text(v+4,i,f'{v}/179',va='center',weight='bold')
B.text(0,-.25,'145 and 45 counts use the 179 two-guide HAR denominator; single-guide HACNS255 excluded.',transform=B.transAxes,fontsize=9,color='#5A6878')
lab(B,'B'); fig.subplots_adjust(left=.15,right=.98,bottom=.20,wspace=.47,top=.84);save(fig,'Figure1_Source_and_Coverage')

# Figure 2 - calibration
fig,ax=plt.subplots(1,2,figsize=(14,5.9),gridspec_kw={'width_ratios':[1,1.05]})
fig.suptitle('Guide-reproducibility calibration and positive-control sensitivity',fontsize=16,weight='bold',y=1.02)
A=ax[0]; x=np.arange(2); ww=.34
A.bar(x-ww/2,[.1298,-.0005],ww,color=blue,label='Same-HAR guides')
A.bar(x+ww/2,[.1363,.0017],ww,color=slate,label='Coverage-matched different HARs')
A.axhline(0,color='#789',lw=.9);A.set_xticks(x,['Uncentered','Gene-response centered']);A.set_ylabel('Median Spearman correlation');A.set_ylim(-.025,.17);A.legend(loc='upper right',fontsize=9,frameon=False);A.grid(axis='y',alpha=.15)
for xi,y in zip([-.17,.17,.83,1.17],[.1298,.1363,-.0005,.0017]): A.text(xi,y+(.003 if y>=0 else -.009),f'{y:+.4f}',fontsize=9,ha='center')
lab(A,'A')
B=ax[1]; gens=['DLGAP1','TFRC','UBQLN2','RAB1A','GRN']; vals=[-.55683,-.27683,-.23495,-.06612,.08065]
B.barh(gens,vals,color=[coral if z<0 else blue for z in vals],height=.52);B.invert_yaxis();B.axvline(0,color='#6B7784',lw=1)
B.set_xlim(-.70,.24); B.set_xlabel('Mean log₂(CPM+1) difference vs batch-matched NTC')
for i,v in enumerate(vals): B.text(v+(-.02 if v<0 else .014),i,f'{v:+.3f}',va='center',ha='right' if v<0 else 'left',weight='bold',fontsize=10)
B.grid(axis='x',alpha=.15);lab(B,'B')
fig.text(.5,.01,'Exploratory CPM contrasts; no donor-level replication, validated target suppression or formal gene-wise DE significance.',ha='center',fontsize=10,color='#566879')
fig.subplots_adjust(left=.10,right=.985,bottom=.17,wspace=.42,top=.83);save(fig,'Figure2_Guide_Control_Calibration')

# Figure 3 - HAR181
fig,ax=plt.subplots(1,2,figsize=(14,5.5),gridspec_kw={'width_ratios':[1.18,1]})
fig.suptitle('Predefined HAR181 reference panel: descriptive transcript readouts',fontsize=16.3,weight='bold',y=1.02)
A=ax[0]; genes=['ABI2','CADM1','PARD3','PCM1','SETD2']; means=[-.196797,-.060492,-.032859,-.076489,.039545]
A.barh(genes,means,color=[coral if z<0 else blue for z in means],height=.55);A.invert_yaxis();A.axvline(0,color='#6B7784',lw=1);A.set_xlim(-.27,.10);A.set_xlabel('Guide/batch mean log₂(CPM+1) difference')
for i,v in enumerate(means): A.text(v-.008 if v<0 else v+.006,i,f'{v:+.3f}',ha='right' if v<0 else 'left',va='center',fontsize=10,weight='bold')
A.grid(axis='x',alpha=.13);lab(A,'A')
B=ax[1]; ranks=[6,25,91,94,158]; colors=[purple]*len(ranks)
B.hlines(range(5),1,ranks,color='#DCE6E9',lw=6);B.scatter(ranks,range(5),color=purple,s=95,zorder=4)
B.set_yticks(range(5),genes);B.invert_yaxis();B.set_xlim(0,185);B.set_xlabel('Descriptive rank: most-negative first (out of 179 HARs)');B.xaxis.set_major_locator(MaxNLocator(5));B.grid(axis='x',alpha=.13)
for i,r in enumerate(ranks): B.text(r+4,i,f'{r}/179',va='center',fontsize=10,weight='bold')
lab(B,'B');fig.text(.5,.012,'The original paper’s PARD3 apical-localization phenotype was not retested by these RNA summaries.',ha='center',fontsize=9.5,color='#566879')
fig.subplots_adjust(left=.11,right=.985,bottom=.20,top=.83,wspace=.35);save(fig,'Figure3_HAR181_Reference_Benchmark')

# Figure 4 - motif findings & selected candidate null
fig,ax=plt.subplots(1,3,figsize=(16.1,5.5),gridspec_kw={'width_ratios':[.8,1,1.25]})
fig.suptitle('Provisional ortholog mapping, aligned motif events, and exploratory calibration',fontsize=16.2,weight='bold',y=1.03)
A=ax[0]; A.bar(['Provisional\nmatched','Unresolved'],[152,28],color=[teal,slate],width=.65)
for i,v in enumerate([152,28]):A.text(i,v+4,str(v),ha='center',weight='bold')
A.set_ylim(0,180);A.set_ylabel('Original screen HAR IDs');A.grid(axis='y',alpha=.12); lab(A,'A')
B=ax[1];tf=['SOX2','SOX11','NFIA','PAX6','REST','CTCF'];ev=[9,4,5,0,0,0]
B.barh(tf,ev,color=[blue,teal,purple,slate,slate,slate]);B.invert_yaxis();B.set_xlim(0,11);B.set_xlabel('Collapsed aligned-site events');B.xaxis.set_major_locator(MaxNLocator(integer=True));B.grid(axis='x',alpha=.12)
for i,v in enumerate(ev): B.text(v+.23,i,str(v),va='center',weight='bold',fontsize=10)
lab(B,'B')
C=ax[2]; cands=['2xHAR.70 / SOX11','HACNS149 / SOX2','2xHAR.310 / SOX11','HACNS38 / SOX11','HACNS772 / SOX2','2xHAR.377 / SOX2']
pvals=[.027888,.043825,.123506,.231076,.239044,.350598];q=[.131474,.131474,.247012,.286853,.286853,.350598]
y=np.arange(6);C.scatter(pvals,y,s=80,color=blue,label='Nominal permutation p');C.scatter(q,y,s=80,color=coral,marker='D',label='Selected-six BH q')
C.axvline(.05,ls='--',lw=1.2,color='#8F98A2');C.set_xlim(0,.41);C.set_yticks(y,cands);C.invert_yaxis();C.set_xlabel('Exploratory event-count probability');C.grid(axis='x',alpha=.12);C.legend(frameon=False,loc='upper right',fontsize=8.5)
lab(C,'C');fig.text(.5,.008,'No selected candidate reached q < 0.05; the six tests were nominated from earlier motif screening.',ha='center',fontsize=9.7,color='#536475')
fig.subplots_adjust(left=.06,right=.99,bottom=.17,top=.82,wspace=.5);save(fig,'Figure4_Motif_Ortholog_Null')
for p in sorted(out.glob('*.png')):
 from PIL import Image
 im=Image.open(p); print(p.name,im.size,p.stat().st_size)
