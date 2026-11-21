"""Input/output validation smoke tests; no claims about biological significance."""
from pathlib import Path
import json,zipfile,io
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];O=ROOT/'results'
S=json.loads((O/'S0_Stage23_Summary.json').read_text());C=S['counts']
assert (C['n_original_HAR'],C['n_provisionally_mapped'],C['n_original_unmapped'])==(180,152,28)
assert (C['n_MPRA_paired_HAR'],C['n_both_motif_and_MPRA'],C['n_shared_nonzero_motif_HAR'])==(66,58,5)
assert (C['n_MPRA_fragments'],C['n_consistent_3cultures'])==(114,68)
assert (C['n_multicut_loci'],C['n_multicut_direction_opposing'],C['n_multicut_stable_3culture_opposing'])==(28,20,11)
assert (C['fullraw_5_10_20_pseudo1_comp_rows'],C['fullraw_5_10_20_pseudo1_sign_changes'])==(1011,81)
assert len(pd.read_csv(O/'S1_180HAR_EvidenceCompatibility.csv'))==180
assert len(pd.read_csv(O/'S4_58HAR_SharedAssays.csv'))==58
assert len(pd.read_csv(O/'S6_MultiCut_Heterogeneity.csv'))==28
assert len(pd.read_csv(O/'S7_114_Cut_3Culture_Summary.csv'))==114
assert abs(S['motif_reporter_transfer']['stratified_label_permutation_two_sided_p']-0.6222459251358288)<1e-12
assert C['n_all152_motif_context_BH_0p05']==0
print('PASS: benchmark data invariants; no missingness imputation detected in 180HAR mapping table')
