#!/usr/bin/env python3
"""Human Accelerated Region evidence compatibility query (secondary benchmark).

Data input: Stage23 results/S1_180HAR_EvidenceCompatibility.csv;
Stage23 results/S7_114_Cut_3Culture_Summary.csv;
Stage24 S2_255_Source_Construct_AnnotationQC.csv.
Does NOT establish endogenous genomic interval identity, gene targets or TF occupancy.
"""
import argparse,csv,json,math
from pathlib import Path

def as_type(s):
    if s is None or s=='':return None
    if isinstance(s,(float,int,bool)):return s
    if s.lower() in ('true','false'):return s.lower()=='true'
    try:
        x=float(s)
        if math.isfinite(x):return x
    except (ValueError,TypeError):pass
    return s

def records(p):
    if not p.exists():raise FileNotFoundError('Missing required published analysis file: '+str(p))
    with p.open(newline='') as f:return list(csv.DictReader(f))

def make_report(har,stage23_dir,construct_csv):
    root=Path(stage23_dir)
    all180=records(root/'S1_180HAR_EvidenceCompatibility.csv')
    matches=[x for x in all180 if x['HAR'].lower()==har.lower()]
    if len(matches)!=1:
        raise ValueError(f'HAR name {har!r} not found in 180-screen target list')
    e=matches[0];h=e['HAR']
    cuts=[x for x in records(root/'S7_114_Cut_3Culture_Summary.csv') if x['HAR']==h]
    inserts=[x for x in records(Path(construct_csv)) if x['HAR']==h]
    status={k:as_type(e.get(k)) for k in
        ['mapping_status','catalog_id_provisional','min_guide_batch_cells','screen_n_guides',
         'passes_min20_guide_batch','passes_min50_guide_batch','all_3_context_aliases_available',
         'observed_max_abs_pwm_delta','context_matched_p','context_BH_q_152',
         'in_MPRA_design_name','MPRA_paired_3cultures','n_cuts','opposing_cut_directions',
         'assay_compatibility_gate']}
    flags=[]
    if status['mapping_status']!='passed_provisional_catalog_match':
        flags.append('PRIMARY_IDENTITY_UNVERIFIED: no provisional ortholog motif comparison')
    else:flags.append('PRIMARY_IDENTITY_UNVERIFIED: provisional catalog mapping is not original-author BED confirmation')
    if status['context_BH_q_152'] is None:
        flags.append('MOTIF_TEST_NOT_AVAILABLE, not q=1 or a zero effect')
    elif status['context_BH_q_152']>=0.05:
        flags.append('NO_FDR_SIGNIFICANT_6TF_MOTIF_RESULT')
    if not status['MPRA_paired_3cultures']:
        flags.append('NO_PAIRED_3CULTURE_REPORTER_CONTRAST')
    if status['opposing_cut_directions'] is True:
        flags.append('OPPOSITE_SIGN_REPORTER_CUTS: avoid whole-HAR signed average')
    if any(x['possible_indel_or_coordinate_disagreement']=='True' for x in inserts):
        flags.append('SOURCE_ANNOTATION_SPAN_LENGTH_MISMATCH: do not directly project sequence offsets')
    if not status['all_3_context_aliases_available']:
        flags.append('NOT_ALL_3_CHROMATIN_CONTEXT_NAME_COVERAGES')
    return {'HAR':h,'source_basis':'180-Screen secondary evidence; source GEO GSE270828 & GSE110758',
            'evidence_status':status,'reporter_cuts':[
                {k:as_type(x[k]) for k in ['cut','mean_HminusC','consistent_direction_3cultures','culture_1','culture_2','culture_3']} for x in cuts],
            'GSE110758_design_construct_annotations':len(inserts),
            'number_of_span_length_discordant_constructs':sum(x['possible_indel_or_coordinate_disagreement']=='True' for x in inserts),
            'warnings':flags,
            'causal_HAR_to_TF_to_gene_mechanism':'NOT_ESTABLISHED',
            'not_a_ranking':'Missingness preserved. Reports evidence availability and restrictions, not a candidate prioritization score.'}

def main():
    ap=argparse.ArgumentParser(description='Query evidence compatibility for one screened HAR (descriptive, not causal).')
    ap.add_argument('--har',required=True,help='e.g. 2xHAR.70')
    ap.add_argument('--stage23-results',default='results',help='directory containing Stage23 S1 and S7 CSVs')
    ap.add_argument('--construct-qc',default='results/stage24/S2_255_Source_Construct_AnnotationQC.csv',help='Stage24 construct annotation QC CSV')
    ap.add_argument('--out',help='output JSON path (default stdout)')
    args=ap.parse_args()
    obj=make_report(args.har,args.stage23_results,args.construct_qc)
    content=json.dumps(obj,indent=2,allow_nan=False)
    if args.out:Path(args.out).write_text(content,encoding='utf8');print(args.out)
    else:print(content)
if __name__=='__main__':main()
