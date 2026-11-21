import nbformat as nbf
from pathlib import Path
import hashlib, zipfile, io, base64, json, csv
ROOT=Path('/mnt/data')
inputs={
 'STAGE20A_152_PROVISIONAL_HUMAN.fa':ROOT/'stage20afast/STAGE20A_152_PROVISIONAL_HUMAN.fa',
 'STAGE20A_152_PROVISIONAL_CHIMP.fa':ROOT/'stage20afast/STAGE20A_152_PROVISIONAL_CHIMP.fa',
 'STAGE20A_152_ALIGNMENT_QC.csv':ROOT/'stage20afast/STAGE20A_152_ALIGNMENT_QC.csv',
 'STAGE20A_S1_HAR_152.csv':ROOT/'Stage20A_Recovered/S1_HAR_152_SixTF_FamilyWise_ConditionalNull.csv',
 'STAGE20A_S2_912_TF.csv':ROOT/'Stage20A_Recovered/S2_912_Pairwise_Aligned_TF_Descriptive.csv',
 'HAR_Stage20A_ConditionalNull_NoDrive.py':ROOT/'stage20afast/HAR_Stage20A_ConditionalNull_NoDrive.py',
 'HAR_Stage20A_Fast_Resumable.py':ROOT/'stage20afast/HAR_Stage20A_Fast_Resumable.py',
 'HAR_Stage20B_3mer_ContextNull_NoDrive.py':ROOT/'HAR_Stage20B_3mer_ContextNull_NoDrive.py',
}
for a,p in inputs.items(): assert p.is_file(),p
names=list(inputs)
inp_hash={k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in inputs.items()}
with io.BytesIO() as buf:
 with zipfile.ZipFile(buf,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
  for name,p in inputs.items():z.writestr(name,p.read_bytes())
 encoded=base64.b64encode(buf.getvalue()).decode()
assert len(encoded)<190000,len(encoded)

nb=nbf.v4.new_notebook(metadata={"colab":{"name":"HAR_Stage20B_3mer_Context_Null_NoDrive.ipynb","provenance":[]},"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"}})
nb.cells.append(nbf.v4.new_markdown_cell('''# Stage 20B — Reference-3-mer-conditioned HAR motif sensitivity analysis

**Goal:** independent context-conditioned reanalysis of the completed Stage 20A results: 152 **provisional** human–chimp HAR ortholog pairs × six JASPAR2024 TF motifs. Uses exactly the same Stage 20A statistic and BH correction across **152 HAR families**, but replaces the reference-base-only substitution-placement null with a **chimp-reference trinucleotide (3-mer) context–matched** substitution-placement null.

**Status:** This is a *computational sensitivity test*, **not** an empirically validated human-lineage gain of TF binding. The chimp reference is a comparator sequence, not proven ancestral state. It does **not** preserve the final sequence's dinucleotide frequencies, provide matched non-HAR genomic negative controls, establish TF occupancy, or resolve original-author HAR interval provenance. TF motif relative-score cutoff ≥0.90 is exploratory, not a FIMO-calibrated p-value.

**Why faster:** Stage 20A showed 136 of 152 HARs have a zero-valued nonnegative maximum PWM-difference statistic. These receive mathematically exact p=1; only the 16 nonzero-statistic HARs need permutation simulation. The familywise BH comparison still includes **all 152 HARs**.

**How to run:** `Runtime → Run all` on Google Colab. There is no Google Drive mount or file upload. All previous input sequences, 152-HAR alignment metadata and Stage 20A result CSVs are embedded and SHA256-checked. The six motif matrices are read from the same checksum-pinned public JASPAR2024 MEME source as Stage 20A (internet required). Default **4,999 permutations per nonzero HAR**. Results checkpoint after *each HAR*; if Colab interrupts, run the notebook again in the same runtime to resume. The last cell triggers browser download **from the notebook process**, avoiding the previous subprocess/kernel error.'''))
nb.cells.append(nbf.v4.new_code_cell('''# Step 1 — Restore and verify the embedded original research files
from pathlib import Path
import hashlib,base64,io,zipfile,json
ROOT=Path('/content') if Path('/content').exists() else Path.cwd()
PAYLOAD="""'''+encoded+'''"""
EXPECTED_HASHES = '''+repr(inp_hash)+'''
with zipfile.ZipFile(io.BytesIO(base64.b64decode(PAYLOAD))) as z:
    assert set(z.namelist())==set(EXPECTED_HASHES), 'Unexpected embedded file set'
    for name in z.namelist():
        raw=z.read(name)
        digest=hashlib.sha256(raw).hexdigest()
        assert digest==EXPECTED_HASHES[name], 'Corrupt embedded payload '+name
        target=ROOT/name
        if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest()==digest:
            print('Existing input OK:',name)
        else:
            target.write_bytes(raw)
            print('Restored:',name)
print('CHECKED:',len(EXPECTED_HASHES),'complete pinned research inputs')'''))
nb.cells.append(nbf.v4.new_code_cell('''# Step 2 — Install dependencies only if missing
import importlib.util,subprocess,sys
for mod,req in [('numpy','numpy>=1.24'),('Bio','biopython>=1.84'),('numba','numba>=0.59'),('requests','requests>=2.31')]:
    if importlib.util.find_spec(mod) is None:
        print('Installing',req)
        subprocess.check_call([sys.executable,'-m','pip','install','-q',req])
import numba,numpy,Bio,requests
print('Dependencies ready:', 'numpy',numpy.__version__,'numba',numba.__version__)'''))
nb.cells.append(nbf.v4.new_code_cell('''# Step 3 — Run Stage 20B in the notebook kernel (NOT a subprocess)
# 4,999 simulations per nonzero HAR; override via os.environ['HAR20B_PERMS'] if needed.
import sys,os,traceback,zipfile,importlib
from pathlib import Path
sys.path.insert(0,str(ROOT))
os.environ.setdefault('HAR20B_PERMS','4999')
try:
    import HAR_Stage20B_3mer_ContextNull_NoDrive as research
    research=importlib.reload(research)
    result = research.main()
    print('SUCCESS: Stage20B research ZIP created:',ROOT/'HAR_Stage20B_ContextMatched_Results.zip')
except Exception:
    traceback.print_exc()
    folder=ROOT/'HAR_Stage20B_3mer_ContextNull_Results'
    diagnostic_zip=ROOT/'HAR_Stage20B_Diagnostics_and_Partial.zip'
    with zipfile.ZipFile(diagnostic_zip,'w',zipfile.ZIP_DEFLATED) as z:
        if folder.exists():
            for p in folder.iterdir():
                if p.is_file():z.write(p,arcname=p.name)
    print('ERROR: Partial results preserved at',diagnostic_zip)
    raise'''))
nb.cells.append(nbf.v4.new_code_cell('''# Step 4 — Download from actual Colab notebook kernel, never inside a subprocess
from pathlib import Path
path=ROOT/'HAR_Stage20B_ContextMatched_Results.zip'
if path.is_file():
    print('Complete output:',path,'size:',path.stat().st_size,'bytes')
    try:
        from google.colab import files
        files.download(str(path))
    except ImportError:
        print('Non-Colab environment: retrieve ZIP from printed path')
else:
    print('No complete ZIP. Use the partial checkpoint from the next cell.')'''))
nb.cells.append(nbf.v4.new_code_cell('''# Optional recovery — use only if the run failed or the browser download was blocked.
from pathlib import Path
import zipfile
folder=ROOT/'HAR_Stage20B_3mer_ContextNull_Results'
backup=ROOT/'HAR_Stage20B_Diagnostics_and_Partial.zip'
with zipfile.ZipFile(backup,'w',zipfile.ZIP_DEFLATED) as z:
    if folder.exists():
        for p in folder.iterdir():
            if p.is_file():z.write(p,arcname=p.name)
print('Partial checkpoint:',backup,backup.stat().st_size,'bytes')
try:
    from google.colab import files
    files.download(str(backup))
except ImportError:pass'''))
notebook=ROOT/'HAR_Stage20B_3mer_Context_Null_NoDrive.ipynb'
nbf.write(nb,notebook)
print('Notebook',notebook,len(encoded),'base64 length',notebook.stat().st_size,'bytes')

qc=[]
for r in csv.DictReader(open(inputs['STAGE20A_S1_HAR_152.csv'])):
    qc.append({'HAR':r['screen_HAR'],'Stage20A_observed_max_pwm_delta':r['observed_max_abs_pwm_delta'],
               'Stage20A_p':r['empirical_p_within_HAR_sixTF_max'],
               'Stage20A_BH_q':r['BH_q_152_HAR_families'],
               'needs_stage20b_simulations':str(float(r['observed_max_abs_pwm_delta'])>0),
               'primary_coordinate_status':'PROVISIONAL'})
with (ROOT/'HAR_Stage20B_Stage20A_Preflight_152.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(qc[0]));w.writeheader();w.writerows(qc)

report=ROOT/'HAR_Stage20B_Design_and_QC_Report.md'
report.write_text('''# Stage 20B — Study design and independent code QC

## Evidence already available (Stage20A)

- 152 provisional ortholog pairs; 912 motif–locus descriptive records; 4,999 reference-base-matched mutations per HAR.
- Stage20A yielded 16 HARs with a nonzero six-TF maximum statistic and 136 with an exactly zero statistic; FDR significance was not reached (152-locus BH q minimum = 1.0).
- For any nonnegative test statistic T, P(T_null ≥ 0) = 1, so a zero observed statistic needs no simulation.

## New Stage20B method (not yet executed in Colab)

For each aligned substitution, compute its gap-free reference 3-mer using the chimpanzee sequence as an operational reference, without declaring chimp the ancestral allele. Keep the chimp→human substitution products within each exact chimpanzee 3-mer stratum; randomly choose eligible positions from the same 3-mer without replacement, preserving the mutation-type multiset within that stratum. Original insertion/deletion alignment columns and mutations lacking a complete reference 3-mer are kept fixed. Compare the same six-TF maximum normalized PWM-difference statistic used by Stage20A. Use a +1 corrected upper-tail Monte Carlo p-value after 4,999 permutations for each nonzero HAR, and BH correction across all 152 HAR families. The 149-HAR source-outlier sensitivity is also reported.

## Verification performed locally

- All 152 genuine sequence-pair alignments were recomputed with the original Biopython alignment procedure, reproducing the Stage20A substitution counts.
- 10 randomized reference-context placements per locus preserved substitution counts, valid alignment bases and gap handling (1,520 local structural test trials).
- The Numba accelerated statistic was cross-checked against the full-window PWM scorer using a synthetic six-motif fixture, including 499 nonzero-statistic permutations. Synthetic fixture scores are **test data only**, not HAR scientific findings.
- Source FASTA/metadata/Stage20A summaries embedded into the notebook with exact SHA-256 verification; six JASPAR matrices are downloaded from the original pinned source and independently checksum-checked at runtime.

## Major limitations

- This context-matched substitution-placement procedure does **not** preserve the complete dinucleotide composition of the resulting human-like sequence.
- Trinucleotide matching is made against a chimp reference, not a phylogenetic ancestral reconstruction. Exact author HAR intervals remain unverified.
- Substitutions at alignment ends or adjacent to reference gaps are held fixed; context strata with no alternative positions offer limited or no randomization. The output audits these restrictions.
- Relative PWM threshold 0.90 is not a calibrated TF binding probability or FIMO p-value.
- No matched non-HAR genomic controls can be validly constructed from the present unverified original HAR boundaries. Independent experimental binding and regulatory targets have not been confirmed.

**Next scientific decision:** Analyze completed Stage20B p, q, and context-relocation degeneracy; only then prioritize independent negative genomic controls or cross-source functional evidence. Do not claim binding or enhancer function based solely on predictions.
''')

zipout=ROOT/'HAR_Stage20B_NoDrive_Research_Code_Package.zip'
with zipfile.ZipFile(zipout,'w',zipfile.ZIP_DEFLATED) as z:
 for filename in [notebook,ROOT/'HAR_Stage20B_3mer_ContextNull_NoDrive.py', ROOT/'HAR_Stage20B_Stage20A_Preflight_152.csv',report]:
  z.write(filename,arcname=filename.name)
print('Zipped package',zipout,zipout.stat().st_size)
