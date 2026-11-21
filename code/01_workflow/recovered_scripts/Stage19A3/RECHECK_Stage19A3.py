"""Reproduce the independent local Stage19A3 ZIP consistency tests.
Place alongside these three source ZIP archives and run `python RECHECK_Stage19A3.py`.
This deliberately does not call UCSC API or equate proximal guides with original HAR identity.
"""
import csv, io, zipfile, json, statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def table(archive,name):
    with zipfile.ZipFile(ROOT/archive) as z:
        return list(csv.DictReader(io.StringIO(z.read(name).decode('utf-8-sig'))))
A2='HAR_Stage19A2_Flank_Guide_Results.zip'
A3='HAR_Stage19A3_Expanded_Flank_Results.zip'
A19='HAR_Stage19A_Colab_Source_Verified_Results.zip'
a2=table(A2,'S2_Original_GuideFlank_Alignment_Details.csv')
a2sum=table(A2,'S1_152_Assembly_and_GuideFlank_Crosscheck.csv')
a3=table(A3,'S1_28_Guide_2kb_Search.csv')
a3sum=table(A3,'S2_25_HAR_2kb_Summary.csv')
source=table(A19,'GSE270828_feature_README.csv')
org={r['id']:r for r in source};cores={r['HAR']:r for r in a2sum}
miss=[r for r in a2 if int(r['n_matches_in_plusminus250_window'])==0]
assert len(a2)==304 and len(miss)==28 and len(a3)==28 and len(cores)==152
assert len(set(r['guide_id'] for r in a2))==304
assert set((r['HAR'],r['guide_id']) for r in miss)==set((r['HAR'],r['guide_id']) for r in a3)
assert len(a3sum)==25 and sum(int(r['guides_recovered_2kb']) for r in a3sum)==28
for r in a3:
    c=cores[r['HAR']]; a,b=int(r['start0']),int(r['end_exclusive'])
    assert b-a==20 and r['sequence_20mer']==org[r['guide_id']]['sequence']
    assert int(r['n_hits_2kb'])==1 and r['guide_hit_status']=='EXACT_20MER'
    assert r['reference_status']=='FETCHED' and r['NGG_PAM']=='True'
    assert len(r['PAM_3mer'])==3 and r['PAM_3mer'].endswith('GG')
    s,e=int(c['core_start0']),int(c['core_end_exclusive'])
    assert int(r['gap_to_core_bp'])==max(s-b,a-e,0)
    assert r['relation']==('UPSTREAM' if b<=s else 'DOWNSTREAM' if a>=e else 'CROSSES')
print(json.dumps({
 'stage':'19A3','locally_verified':True,
 'provisional_catalog_mapped_hars':len(cores),
 'guides_in_initial_250bp':len(a2)-len(miss),
 'guides_recovered_additional_2kb':len(a3),
 'total_guides_locally_matched_2kb':len(a2),
 'newly_rescued_hars':len(a3sum),
 'ngg_all_new_recovered':True,
 'median_new_gap_bp':statistics.median(int(r['gap_to_core_bp']) for r in a3),
 'original_author_exact_HAR_BED_verified':False,
 'independent_UCSC_API_redownload':False
}, indent=2))
