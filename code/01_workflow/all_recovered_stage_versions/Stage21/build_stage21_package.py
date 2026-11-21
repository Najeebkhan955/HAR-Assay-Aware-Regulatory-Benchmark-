from pathlib import Path
import zipfile,hashlib,csv,io,json,statistics
root=Path('/mnt/data')
zip20=root/'HAR_Stage20B_ContextMatched_Results.zip'
diag=root/'HAR_Stage20B_Diagnostics_and_Partial.zip'
with zipfile.ZipFile(zip20) as z:
 s1=list(csv.DictReader(io.StringIO(z.read('S1_152_HAR_ContextMatched_Statistics.csv').decode())))
 s2=list(csv.DictReader(io.StringIO(z.read('S2_16_Nonzero_Locus_Comparisons.csv').decode())))
 s3=list(csv.DictReader(io.StringIO(z.read('S3_912_TF_Descriptive_Stage20A_NoNewTests.csv').decode())))
 status=json.loads(z.read('S0_Stage20B_RunStatus_and_Limitations.json'))
 ziphashes={n:hashlib.sha256(z.read(n)).hexdigest() for n in z.namelist()}
with zipfile.ZipFile(diag) as z:assert all(hashlib.sha256(z.read(n)).hexdigest()==ziphashes[n] for n in z.namelist())
assert len(s1)==152 and len(s2)==16 and len(s3)==912
assert len({r['screen_HAR'] for r in s1})==152
assert sum(int(r['n_null_permutations'])==4999 for r in s1)==16
assert sum(int(r['n_null_permutations'])==0 for r in s1)==136
assert min(float(r['context_matched_p']) for r in s1)==0.1188
assert all(float(r['context_BH_q_152'])==1 for r in s1)
assert all(float(r['context_BH_q_149_nonoutlier'])==1 for r in s1 if r['context_BH_q_149_nonoutlier'] not in ('EXCLUDED',''))
assert all(abs((int(r['extreme_null_count'])+1)/5000-float(r['context_matched_p']))<1e-10 for r in s2)

S1_by={r['screen_HAR']:r for r in s1}
focus=['2xHAR.377','HACNS149','HACNS880','2xHAR.70','2xHAR.310','HACNS772']
def f4(v):return f'{float(v):.4f}'
table='\n'.join([f"| {h} | {f4(S1_by[h]['Stage20A_baseMatched_p'])} | {f4(S1_by[h]['context_matched_p'])} | {f4(S1_by[h]['context_BH_q_152'])} |" for h in focus])
report=f'''# HAR Stage 20B — Independent verification and Stage 21A evidence integration

**Date:** 20 November 2026  
**Scope:** Uploaded `HAR_Stage20B_ContextMatched_Results.zip` and `HAR_Stage20B_Diagnostics_and_Partial.zip`; cross-reference with prior verified Stage 15 and Stage 17 project workbooks.

## 1. Run integrity and numerical audit

- Both uploaded archives contain **exactly the same five output files at matching SHA-256 content hashes**. The diagnostics ZIP is not an incomplete alternative; the Stage 20B recorded status is `complete`.
- There are **152 unique provisionally paired human–chimp HARs**, with **912 descriptive HAR × TF records** (six motifs per HAR). No new TF-specific significance model was fitted.
- Exactly **16 HAR loci have a nonzero observed threshold-gated PWM statistic** and were simulated with **4,999 conditional permutations each**. The other **136** have a nonnegative test statistic of zero, giving an exact upper-tail Monte Carlo *p*=1 under the defined test; treating those as untested biological negatives would be incorrect.
- The 16 reported empirical p-values exactly satisfy **(extreme_null_count + 1)/5,000**. An independent Benjamini–Hochberg calculation across the **152 HAR-level six-TF maxima** matches the results: minimum nominal **p=0.1188**, minimum adjusted **q=1.000**, and **0/152 HARs q<0.05**.
- The sensitivity analysis retaining **149** HARs after excluding three source/guide-distance outlier loci also has **no BH-FDR discoveries** (all recorded adjusted q values 1.000).
- The tested outputs leave **28 of the 180 original HARs unmapped**. These 28 should be represented as **missing sequence evidence, not zero motif effects or q=1**.

## 2. Changes in nominal evidence under chimp-reference 3-mer conditioning

| Previously examined HAR | Stage 20A base-matched p | Stage 20B 3-mer p | Stage 20B BH q (152) |
|---|---:|---:|---:|
{table}

For example, **2xHAR.377** has nominal *p*=0.0228 under the Stage 20A placement model but *p*=0.4416 when substitutions are relocated among positions sharing the chimpanzee reference trinucleotide context; **HACNS149** changes from 0.0498 to 1.0000. The disparity shows **sensitivity to null-model assumptions**, not proof that any DNA site does or does not bind a TF.

## 3. Cross-layer evidence triage — existing datasets only

The Stage 21A crosswalk combines the 180 original perturbation identifiers with Stage 17 provisional sequence associations, Stage 15 guide/batch QC and name-level chromatin representation, and Stage 20B statistical summaries. All 180 identifiers match between the prior Stage 15/17 tables; all 152 Stage 20B loci align to these identifiers.

- **16** loci have a nonzero observed computational motif statistic.
- **12/16** meet the descriptive minimum of **20 cells per guide per experimental batch**; **5/16** meet the stricter minimum of 50.
- **14/16** have the HAR name represented in all three chromatin-contact context tables (human NSCs, human neurons and chimp NSCs).
- **11/16** satisfy both the ≥20 cell criterion and three-context name representation. These are **follow-up data-availability flags only**, not independent support for a causal gene target.
- The three experimental batches arise from one iPSC genetic background; prior centered guide-to-guide correlations were not globally higher than coverage-matched different-HAR guide pairs. A guide-design gene annotation is **not** an experimentally validated HAR target, and the source contact tables' name overlap alone does **not** establish HAR–gene chromatin loops.

## 4. Model and provenance limitations

1. The 152 ortholog pairs rely on a **third-party sequence catalog**, and the original author's exact HAR BED intervals remain unverified.
2. The 3-mer null preserves each movable substitution's **chimpanzee-reference three-base neighborhood and substitution identity**, but does **not** preserve exact dinucleotide counts, local mutation rate spectra or phylogenetic ancestral state. Indels and substitutions without a valid local three-base context are held fixed.
3. Relative PWM score ≥0.90 is an **exploratory motif-score threshold** rather than an empirically calibrated FIMO *p*-value or observed protein occupancy.
4. The familywise maximum of six TF scores per HAR and BH across 152 HAR families should not be described as 912 independently significant TF-specific tests.
5. No independent donor/line perturbation replication, TF occupancy measurement or new HAR–TF–gene functional mechanism was established in Stages 20A–21A.
6. A null-model failure to reject the hypothesis cannot be interpreted as evidence that a locus is nonfunctional.

## 5. What to do next: independent functional evidence, not another permutation loop

**Stage 21B:** Use the public human/chimp neural-progenitor lentiMPRA study **GSE110758** first as an **identity-overlap feasibility check** (study design includes 714 HARs). The companion notebook downloads only `GSE110758_HAR_design.tsv.gz` and checks literal names; this does **not** analyze enhancer activity or prove overlapping functional loci. Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110758 .

If genuine overlapping HARs are found after **original identifier, sequence and genome-build harmonization**, obtain published reporter-activity contrasts or source-appropriate RNA/DNA readouts. Prespecify species, neural developmental stage and assay endpoints; compare recorded human–chimp enhancer-activity direction/magnitude to previously observed motif scores and separately to guide-level expression/contacts. Evaluate concordance only with sufficient matched labels, calibrated negative regions and appropriate multiple-testing procedures. Where no overlap exists, report non-comparability—not failure of validation.

Alternative source: **GSE292632** provides comparative CRISPRi of 106 cis-regulated open chromatin regions in human/chimp neural progenitors; those regions must **not** be assumed to be the same 180 HARs. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292632 . Additional independent published HAR experiments: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222113 .

**Decision:** Retain all Stage 20B results, avoid claiming FDR-significant motif rewiring, and advance to well-defined external evidence comparison only after a verified locus crosswalk. Do not pursue protein–ligand docking for this noncoding regulatory question without a new justified molecular mechanism.

## 6. Deliverables

- `HAR_Stage21A_CrossLayer_Evidence_S1-S6.xlsx`: 180-locus crosswalk, 16 nonzero motif loci, six previously shortlisted cases, Stage 20A/B sensitivity, source inventory and full 912 descriptive TF records.
- `HAR_Stage21B_MPRA_Design_Preflight_NoDrive.ipynb`: external GEO design-file name-overlap feasibility; generates independent results ZIP or diagnostics and makes **no** biological activity claim.
- Original Stage 20B user-produced results archive retained in the package.

**Numerical assurance:** Both uploaded ZIPs identical internally; 152 unique HARs; 16 × 4,999 permutation statistics; full-family BH recalculated independently; 0 discoveries. The source workbook merger is a separate **data-provenance reconciliation**, not novel experimental validation.
'''
report_path=root/'HAR_Stage20B_Independent_Audit_and_Stage21_Plan.md'
report_path.write_text(report,encoding='utf-8')
manifest={'stage':'21A','source_files':{},'numerical_checks':{
'har_original':180,'har_provisional':152,'har_missing':28,'active_har':16,'perms_per_active_har':4999,'zero_stat_p_exact1':136,'min_p20B':0.1188,'min_q20B':1.0,'n_q_below_005':0,'data_available_flags_for_active':11},
'constraints':["Loci mapped only to secondary catalog, not primary author's exact HAR BED",'No causal validated HAR-gene or TF occupancy','Stage21B design feasibility is NOT enhancer activity analysis']}
for p in [zip20,diag,root/'HAR_Stage21A_CrossLayer_Evidence_S1-S6.xlsx',report_path,root/'HAR_Stage21B_MPRA_Design_Preflight_NoDrive.ipynb',root/'files/HAR_Stage15_Editable_Supplementary.xlsx',root/'files/HAR_Stage17_Editable_Research_Supplementary.xlsx']:
 manifest['source_files'][p.name]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'size_bytes':p.stat().st_size}
manifest_path=root/'HAR_Stage21_Integrity_Manifest.json'
manifest_path.write_text(json.dumps(manifest,indent=2),encoding='utf-8')
package=root/'HAR_Stage20B_Reviewed_Stage21A_and_Next_Stage21B.zip'
with zipfile.ZipFile(package,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in [zip20,root/'HAR_Stage21A_CrossLayer_Evidence_S1-S6.xlsx',report_path,root/'HAR_Stage21B_MPRA_Design_Preflight_NoDrive.ipynb',manifest_path]:z.write(p,p.name)
print('PACKAGE',package,package.stat().st_size)
for p in [report_path,root/'HAR_Stage21A_CrossLayer_Evidence_S1-S6.xlsx',root/'HAR_Stage21B_MPRA_Design_Preflight_NoDrive.ipynb',manifest_path]:
 print('verified',p.exists(),p.name,p.stat().st_size)
