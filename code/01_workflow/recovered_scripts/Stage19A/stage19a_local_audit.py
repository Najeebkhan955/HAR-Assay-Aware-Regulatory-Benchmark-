#!/usr/bin/env python3
"""Stage 19A: Offline internal consistency audit of provided Stage 16 and 16B archives.
Never claims original author exact HAR element or guide-coordinate confirmation.
"""
import csv, json, hashlib, collections, pathlib
ROOT=pathlib.Path(__file__).resolve().parent
OUT=ROOT/'output'
OUT.mkdir(exist_ok=True)
def load(p):
    with open(p, newline='', encoding='utf-8') as f:return list(csv.DictReader(f))
def save(name, records, columns=None):
    with (OUT/name).open('w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns or list(records[0]) if records else columns)
        w.writeheader(); w.writerows(records)
def overlap_union(fragments,chrom,start,end):
    intervals=sorted({(max(start,int(f['start_0based'])),min(end,int(f['end_exclusive']))) for f in fragments if f['chrom']==chrom and min(end,int(f['end_exclusive']))>max(start,int(f['start_0based']))})
    merged=[]
    for a,b in intervals:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    return sum(b-a for a,b in merged)

screen=load(ROOT/'Stage16/S1_All180_HAR_Sequence_Mapping_Readiness.csv')
fragments=load(ROOT/'Stage16/S2b_Human_NSC_Bait_Fragments_Embedded_for_Colab.csv')
audit=load(ROOT/'results/S1_All180_Sequence_Mapping_Audit.csv')
paired=load(ROOT/'results/S2_Aligned_Ortholog_Sequence_QC.csv')
motif=load(ROOT/'results/S3_SixTF_Motif_Score_Comparisons.csv')
manifest=json.loads((ROOT/'results/SOURCE_CHECKSUMS.json').read_text())
assert len(screen)==len(audit)==180 and len(paired)==152 and len(motif)==912
by_bait=collections.defaultdict(list)
for z in fragments:
    assert z['context']=='hNSCs' and z['assembly']=='hg38'
    by_bait[z['screen_HAR']].append(z)
by_screen={z['screen_HAR']:z for z in screen}
by_audit={z['screen_HAR']:z for z in audit}
by_pair={z['screen_HAR']:z for z in paired}
assert len(by_screen)==len(by_audit)==180 and len(by_pair)==152
assert set(by_screen)==set(by_audit) and set(by_pair)<=set(by_audit)
assert len({z['catalog_id_provisional'] for z in paired})==152
assert set(z['screen_HAR'] for z in motif)==set(by_pair)
assert collections.Counter(x['status'] for x in audit)=={'passed_provisional_catalog_match':152,'no_catalog_element_covered_GE90pct':22,'no_named_hNSC_bait':6}
assert collections.Counter(x['TF'] for x in motif)=={k:152 for k in ('SOX2','PAX6','REST','SOX11','CTCF','NFIA')}
checks=[]; rows=[]; seqrows=[]; unresolved=[]
for name in sorted(by_audit):
    a=by_audit[name]; s=by_screen[name];p=by_pair.get(name)
    if a['status']=='passed_provisional_catalog_match':
        assert p
        start,end=int(p['h_start']),int(p['h_end'])
        recomputed=overlap_union(by_bait[name],p['h_chr'],start,end)
        recorded=int(float(a['overlap_bp']))
        match=recomputed==recorded
        frac=round(recomputed/(end-start),5)
        frac_ok=abs(frac-float(a['overlap_fraction']))<0.000012
        checks.append(match and frac_ok)
        assert match and frac_ok, f"overlap mismatch: {name}"
        baits=sorted({(x['chrom'],x['start_0based'],x['end_exclusive']) for x in by_bait[name]})
        classification='PROVISIONAL_CATALOG_MATCH_NOT_PRIMARY_CONFIRMED'
        reason='Unique qualifying secondary-catalog candidate; primary author element coordinates absent'
    else:
        recomputed=''; frac=''; baits=sorted({(x['chrom'],x['start_0based'],x['end_exclusive']) for x in by_bait.get(name,[])})
        classification='UNRESOLVED_NO_HNSC_BAIT' if a['status']=='no_named_hNSC_bait' else 'UNRESOLVED_NO_90PCT_CATALOG_COVERAGE'
        reason=('No exact-name human NSC bait record; do not assign zero motif effect' if a['status']=='no_named_hNSC_bait' else 'Secondary catalog has no qualifying >=90% candidate for named bait; cannot impute ortholog')
    row={
    'screen_HAR':name,'stage16B_status':a['status'],'stage19A_status':classification,
    'primary_author_exact_HAR_interval_verified':'NO','primary_author_interval_source':'NOT_PROVIDED_IN_ARCHIVES',
    'independent_bait_to_catalog_overlap_rechecked': 'YES' if p else 'NOT_APPLICABLE',
    'mapping_confidence':'PROVISIONAL' if p else 'UNRESOLVED',
    'catalog_ID_not_original_name':p['catalog_id_provisional'] if p else '',
    'human_hg38_chrom':p['h_chr'] if p else '', 'human_start_0based':int(p['h_start']) if p else '', 'human_end_exclusive':int(p['h_end']) if p else '',
    'chimp_panTro6_chrom':p['c_chr'] if p else '', 'chimp_start_0based':int(p['c_start']) if p else '', 'chimp_end_exclusive':int(p['c_end']) if p else '',
    'human_length_bp':int(p['human_length_bp']) if p else '', 'chimp_length_bp':int(p['chimp_length_bp']) if p else '',
    'alignment_orientation':p['alignment_orientation'] if p else '',
    'pairwise_alignment_identity':float(p['global_alignment_identity']) if p else '',
    'fragment_union_overlap_bp':recomputed,'fragment_union_overlap_fraction':frac,
    'number_distinct_hNSC_bait_fragments':len(baits),
    'n_catalog_candidates_GE90pct':int(a['candidate_records_covering_GE90pct']),
    'n_original_guides':int(s['screen_n_guides']),
    'min_guide_batch_cells':int(s['min_guide_batch_cells']),
    'pass_20_cells_guide_batch':s['passes_min20_guide_batch'],
    'pass_50_cells_guide_batch':s['passes_min50_guide_batch'],
    'present_name_in_3_chromatin_contexts':s['all_3_context_aliases_available'],
    'original_HAR_to_catalog_equivalence':'NOT_ESTABLISHED','TF_occupancy_validation':'NOT_PERFORMED',
    'notes':reason,
    'source_perturbation':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270828',
    'source_chromatin':'https://doi.org/10.5061/dryad.c59zw3rkj',
    'source_catalog':'https://github.com/athenamarou/HAR-TFBS-Project/tree/a9251f39976bbc15df1d8f0dfedaf9786819c298'
    }
    rows.append(row)
    if p:
        seqrows.append({k:row[k] for k in ('screen_HAR','catalog_ID_not_original_name','human_hg38_chrom','human_start_0based','human_end_exclusive','chimp_panTro6_chrom','chimp_start_0based','chimp_end_exclusive','human_length_bp','chimp_length_bp','alignment_orientation','pairwise_alignment_identity','fragment_union_overlap_fraction','primary_author_exact_HAR_interval_verified','original_HAR_to_catalog_equivalence','source_catalog')})
    else:unresolved.append(row)

save('S1_All180_Stage19A_Crosswalk.csv',rows)
save('S2_152_Provisional_Ortholog_Coordinates.csv',seqrows)
save('S3_28_Unresolved_HARs.csv',unresolved)
# Reusable 0-based half-open BEDs; names are screened alias + provisional catalog ID, NOT canonical identity
with (OUT/'S4_152_Provisional_hg38.bed').open('w') as f:
    for r in sorted(seqrows,key=lambda z:(z['human_hg38_chrom'],int(z['human_start_0based']))):
        f.write(f"{r['human_hg38_chrom']}\t{r['human_start_0based']}\t{r['human_end_exclusive']}\t{r['screen_HAR']}|{r['catalog_ID_not_original_name']}|PROVISIONAL\n")
with (OUT/'S5_152_Provisional_panTro6.bed').open('w') as f:
    for r in sorted(seqrows,key=lambda z:(z['chimp_panTro6_chrom'],int(z['chimp_start_0based']))):
        f.write(f"{r['chimp_panTro6_chrom']}\t{r['chimp_start_0based']}\t{r['chimp_end_exclusive']}\t{r['screen_HAR']}|{r['catalog_ID_not_original_name']}|PROVISIONAL\n")
# The direct source catalog sequences are NOT in the provided ZIP; BED manifest only.
selected=['2xHAR.70','HACNS772','2xHAR.377','HACNS149','HACNS38','2xHAR.310','HAR181']
flags=[]
for name in selected:
    if name not in by_audit:continue
    r=next(z for z in rows if z['screen_HAR']==name)
    flags.append({k:r[k] for k in ('screen_HAR','stage19A_status','catalog_ID_not_original_name','human_hg38_chrom','human_start_0based','human_end_exclusive','chimp_panTro6_chrom','chimp_start_0based','chimp_end_exclusive','min_guide_batch_cells','pass_20_cells_guide_batch','pass_50_cells_guide_batch','primary_author_exact_HAR_interval_verified','original_HAR_to_catalog_equivalence')})
save('S6_Six_Motif_Candidates_and_HAR181_Flags.csv',flags)
summary={
 'original_HAR_screened':180,'independent_original_author_exact_coordinate_confirmed':0,
 'provisional_secondary_catalog_mapped':152,'unresolved':28,
 'unresolved_no_named_hNSC_bait':6,'unresolved_no_90pct_catalog_coverage':22,
 'provisional_catalog_candidate_ids_unique':len({z['catalog_id_provisional'] for z in paired}),
 'provisional_overlap_independently_rechecked_count':len(checks),
 'provisional_overlap_disagreements':sum(not z for z in checks),
 'source_FASTA_absent_from_uploaded_archives':True,
 'original_author_exact_HAR_BED_absent_from_uploaded_archives':True,
 'primary_source_guide_design_code_only_no_original_BED_retrieved':True,
 'inference':'The arithmetic of the 152 provisional overlaps is independently reproduced within provided records; this is not independent verification of experimental HAR identity, orthology ancestry or causal biological function.',
 'next':'Run companion no-Drive Colab to fetch/check pinned BED+FASTA and GEO feature README. Obtain exact original author HAR BED for independent element identity comparison.'
}
(OUT/'S0_Stage19A_Validation_Summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
