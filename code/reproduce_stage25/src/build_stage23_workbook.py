"""Build editable Stage23 scientific evidence workbook via artifact_tool (not openpyxl)."""
import csv,json
from pathlib import Path
from artifact_tool import Workbook, SpreadsheetFile
ROOT=Path(__file__).resolve().parents[1];RES=ROOT/'results';SRC=ROOT/'data/source_archives'
wb=Workbook.create()
HEADER={'fill':'#183550','font':{'color':'#FFFFFF','bold':True},'row_height':26,'vertical_alignment':'center'}

def a1col(n):
    s=''
    while n:
        n,rem=divmod(n-1,26);s=chr(65+rem)+s
    return s

def rows_csv(path):
    with open(path,newline='',encoding='utf-8') as fh:
        d=list(csv.reader(fh))
    result=[]
    for row in d:
        out=[]
        for v in row:
            if v=='':out.append(None)
            elif v=='True':out.append(True)
            elif v=='False':out.append(False)
            else:
                try:out.append(float(v) if '.' in v or 'e' in v.lower() else int(v))
                except ValueError:out.append(v)
        result.append(out)
    return result

def add(name,path,cols=None):
    data=rows_csv(path)
    if cols:data=[[v for v in r[:cols]] for r in data]
    sh=wb.worksheets.add(name);n=len(data);c=len(data[0]);last=a1col(c)
    sh.get_range(f'A1:{last}{n}').values=data
    sh.get_range(f'A1:{last}1').format=HEADER
    sh.get_range(f'A1:{last}1').format.row_height=30
    sh.get_range(f'A1:{last}{n}').format.row_height=19
    sh.get_range(f'A1:{last}1').format.row_height=30
    sh.freeze_panes.freeze_rows(1)
    sh.get_range(f'A1:{last}{n}').format.vertical_alignment='center'
    for i,h in enumerate(data[0],1):
        s=str(h).lower();width=29 if any(x in s for x in ('status','interpretation','note','mapping','gate','context')) else (24 if any(x in s for x in ('har','cut','motif','score')) else 18)
        sh.get_range(f'{a1col(i)}:{a1col(i)}').format.column_width=width
    sh.get_range('A:A').format.column_width=17
    print('Sheet',name,n-1,'rows',c,'columns')
    return sh

sumr=json.loads((RES/'S0_Stage23_Summary.json').read_text());c=sumr['counts'];t=sumr['motif_reporter_transfer']
readme=wb.worksheets.add('00_README_AND_KEY')
items=[['STAGE 23 | ASSAY-AWARE HAR EVIDENCE BENCHMARK','VALUE','INTERPRETATION / LIMIT'],
       ['Reference date','2026-11-20','All data already provided; this is independent secondary analysis'],
       ['Original Noble HARs',c['n_original_HAR'],'Original screen, one iPSC background'],
       ['Provisional mapped orthologs',c['n_provisionally_mapped'],'NOT author-verified genomic boundaries'],
       ['Original names found in MPRA design',c['n_GEO_design_name_HAR'],'Name match only'],
       ['HARs with paired HS1-11/N2 MPRA',c['n_MPRA_paired_HAR'],'3 cultures, same human line'],
       ['HARs with both motif and paired MPRA',c['n_both_motif_and_MPRA'],'Matched assay subcohort'],
       ['HARs with nonzero motif statistic in matched subcohort',c['n_shared_nonzero_motif_HAR'],'None survived 152-HAR BH correction'],
       ['Fragment comparisons',c['n_MPRA_fragments'],'Cuts are nested within HARs'],
       ['Same direction in all 3 cultures',c['n_consistent_3cultures'],'Descriptive, not independent donor replication'],
       ['HARs with more than one reporter fragment',c['n_multicut_loci'],'Comparisons must remain cut-specific'],
       ['HARs with opposite mean cut effects',c['n_multicut_direction_opposing'],'Not necessarily all strong differences'],
       ['HARs with 3-culture-stable opposite-sign cuts',c['n_multicut_stable_3culture_opposing'],'Two opposing cuts each keep their sign across all 3 cultures'],
       ['Motif-positive vs motif-zero effect difference',round(t['observed_difference'],6),'Mean within motif-positive minus motif-zero, HAR level'],
       ['Stratified-label permutation p (two-sided)',round(t['stratified_label_permutation_two_sided_p'],5),'Post-hoc; preserves number-of-cuts stratum'],
       ['Stage 20B BH significant loci',c['n_all152_motif_context_BH_0p05'],'Family = 152 provisionally mapped loci'],
       ['Stricter DNA threshold sign flips',f"{c['fullraw_5_10_20_pseudo1_sign_changes']}/{c['fullraw_5_10_20_pseudo1_comp_rows']}",'Correlated sensitivity cells, not experiments'],
       ['Upcoming external-line test','WTc N2','One essential Colab notebook, not yet run in this workspace'],
       ['Claims that are NOT supported','New HAR–TF–gene causal mechanism','No primary author BED or TF occupancy/target validation'],
       ['Output scope','Secondary computational benchmarking','Do not call this prospective experimental validation'],
       ['Sources','GSE270828; GSE110758; Pal 2025; pinned catalog','Provenance details in repository docs and S7'],
       ['README','Missing source is NA','No effects or p-values should be imputed as zero'],]
n=len(items);readme.get_range(f'A1:C{n}').values=items
readme.get_range('A1:C1').format=HEADER;readme.get_range('A:A').format.column_width=53;readme.get_range('B:B').format.column_width=23;readme.get_range('C:C').format.column_width=76
readme.get_range(f'A2:C{n}').format.row_height=23
readme.freeze_panes.freeze_rows(1)
add('S1_180_AssayGate',RES/'S1_180HAR_EvidenceCompatibility.csv')
add('S2_58_SharedHARs',RES/'S4_58HAR_SharedAssays.csv')
add('S3_114_MPRAcuts',RES/'S7_114_Cut_3Culture_Summary.csv')
add('S4_MultiCut',RES/'S6_MultiCut_Heterogeneity.csv')
add('S5_DNA_Sensitivity',RES/'S5_DNA_DepthSensitivity.csv')
add('S6_ReplicateReliab',RES/'S2_MPRA_Replicate_Reliability.csv')
readme2=wb.worksheets.add('S7_TransferBenchmark')
entries=[['Metric','Value','Definition / scope'],
         ['Shared HARs',t['shared_HARs'],'58 loci with both 3-mer motif and paired reporter'],
         ['Motif-positive',t['nonzero_PWM_HARs'],'Nonzero six-TF observed max |ΔPWM|, not significant'],
         ['Motif-zero',t['PWM_zero_HARs'],'Context-conditioned statistic exactly zero under fixed threshold'],
         ['Motif-positive mean absolute effect',t['motif_nonzero_mean'],'Per-HAR median of absolute fragment mean, then group mean'],
         ['Motif-zero mean absolute effect',t['motif_zero_mean'],'Same definition'],
         ['Difference',t['observed_difference'],'Motif-positive minus motif-zero'],
         ['Permutation p two-sided',t['stratified_label_permutation_two_sided_p'],'30,000 shuffles; preserve number-of-cut strata'],
         ['Permutation p upper',t['stratified_label_permutation_upper_p'],'Preliminary directional check, not an adjusted p'],
         ['Bootstrap lower (95%)',t['bootstrap_difference_95CI'][0],'Descriptive uncertainty, resample within groups'],
         ['Bootstrap upper (95%)',t['bootstrap_difference_95CI'][1],'Descriptive uncertainty, resample within groups'],
         ['Continuous motif/reporter Spearman rho',t['spearman_continuous_PWM_delta_vs_reporter_abs_r'],'Post hoc 58-HAR exploratory association'],
         ['Continuous motif/reporter nominal p',t['spearman_exploratory_p'],'Unadjusted exploratory statistic'],
         ['Scientific interpretation','NO PREDICTIVE ADVANTAGE ESTABLISHED','No mechanistic HAR-TF-gene validation established']]
readme2.get_range(f'A1:C{len(entries)}').values=entries
readme2.get_range('A1:C1').format=HEADER
readme2.get_range('A:A').format.column_width=45;readme2.get_range('B:B').format.column_width=21;readme2.get_range('C:C').format.column_width=83
readme2.get_range(f'A2:C{len(entries)}').format.row_height=22;readme2.freeze_panes.freeze_rows(1)

# File checksums and public URLs in plain-text cells.
src=wb.worksheets.add('S8_SourceRegister')
srows=[['Source','Accession or pinned revision','Public link','Usage and caveat'],
 ['Original HAR Perturb-seq','GSE270828','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270828','180 HAR; one genetic background and 3 experimental batches'],
 ['Independent lentiMPRA','GSE110758 / GSE110760','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110758','WTc and HS1-11 human lines; N2/N3. This report currently HS1-11 N2 only'],
 ['3D chromatin','Pal et al. 2025','https://pmc.ncbi.nlm.nih.gov/articles/PMC11928272/','Bait-name representation is not validated physical target loop'],
 ['Secondary genomic sequence catalog','commit a9251f39976bbc15df1d8f0dfedaf9786819c298','https://github.com/athenamarou/HAR-TFBS-Project/tree/a9251f39976bbc15df1d8f0dfedaf9786819c298','Provisional mapped ortholog sequences; author BED unresolved'],
 ['Source workbooks','Stages 12, 13, 15, 17, 18','Bundled in repository data/prior_workbooks','User-provided prior derived research'],
 ['Prior validation archive','Stage20B','Bundled in repository data/source_archives','152-locus context-conditioned null; all q=1'],
 ['Reporter outputs','Stage21C/Stage21D','Bundled in repository data/source_archives','Assay-matched reporter effects and full count sensitivity'],
 ['Next independent-line validation','WTc N2','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3016026','Notebook available but source count reanalysis not completed here']]
src.get_range(f'A1:D{len(srows)}').values=srows;src.get_range('A1:D1').format=HEADER
for col,width in [('A',36),('B',55),('C',85),('D',87)]:src.get_range(f'{col}:{col}').format.column_width=width
src.get_range(f'A2:D{len(srows)}').format.row_height=25;src.freeze_panes.freeze_rows(1)
# Focused verification (no silent truncated layouts)
print('Workbook sheets',len(wb.worksheets.items))
print(wb.inspect({'kind':'table','range':'00_README_AND_KEY!A1:C8','include':'values','table_max_rows':8,'table_max_cols':3}).ndjson[:1600])
errors=wb.inspect({'kind':'match','search_term':'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A','options':{'use_regex':True,'max_results':200},'summary':'scientific workbook error scan'})
print('error scan',errors.ndjson[:800])
out=ROOT/'HAR_Stage23_Novelty_and_Scientific_Benchmark_S1-S8.xlsx'
SpreadsheetFile.export_xlsx(wb).save(str(out));print('EXPORTED',out,out.stat().st_size)
