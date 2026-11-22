# HAR Assay-Aware Regulatory Benchmark

**Repository overview.** Reproducible, assay-aware benchmarking of human accelerated region (HAR) regulatory evidence across primate neural development, integrating public CRISPRi perturbation data, comparative sequence and transcription-factor motif analyses, chromatin contacts, and human–chimp MPRA reporter assays.

**Research objective.** Evaluate which regulatory conclusions are supported across assays, distinguish construct-level reporter effects from endogenous HAR–gene mechanisms, and document limitations from genomic coordinate matching, fragment selection, replication, and statistical calibration.

**Status (22 November 2026).** Stage-wise recovered code archive is **online and checksum-verified**: 108 source-code files (42 Jupyter notebooks, 65 Python scripts, and 1 shell script). Stage 31 contains **17 prepared Excel supplements (S01–S17), 9 main and 10 Extended Data figure pairs**, and a checksum-audited local release ZIP. Stage 31 [release documentation](releases/stage31/RELEASE_NOTES.md) and a [token-free GitHub issue importer](.github/workflows/stage31-import-issue-attachment.yml) are public; **binary files remain pending until the ZIP import succeeds**. Attribution and redistribution requirements still apply.

**Key public accessions.** [GSE270828](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270828) and [GSE110758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110758).

**Scope limitation.** This project is an evidence benchmarking resource, not a claim of a newly experimentally confirmed causal HAR–TF–gene mechanism.

## Code archive — verified uploaded

- [Open the full source-code archive](code/)
- [Code guide, stage map and scientific limitations](code/README.md)
- [108-file SHA-256 and Git blob integrity manifest](code/CODE_HASH_MANIFEST.csv)
- [Additional recovered script provenance](code/CODE_PROVENANCE_ADDENDUM.csv)

Original interactive notebooks and subsequently produced analysis scripts are preserved as separate archival versions. Automated completeness is **limited to the recovered and accessible files**. Original code for some earlier stages is unavailable, and no missing script is passed off as an authentic original. External GEO raw datasets and final publication supplements are **not uploaded yet**.

### Repository attribution and reproducibility

The account `Najeebkhan955` maintains this repository. File origin is documented from recovered notebook uploads and historical analysis packages. Research conclusions must cite the underlying primary studies and reflect the actual data and analyses; GitHub source code by itself is not evidence of experimental validation, nor does a commit prove sole human authorship. Manuscript disclosure of computational assistance should follow the target journal's requirements.

## Data archive: Batch 01 uploaded and verified (21 November 2026)

The first user-provided input/output batch contains **19 original files**, representing **31,476,940 bytes** in total. All 19 have been uploaded either as intact files or numerically ordered binary parts, and their 53 Git blobs were checked against original source bytes (**53/53 exact matches**). This is a staged archive, not the final complete GSE270828 source-matrix deposit.

- [Batch 01 file provenance, SHA-256 and uploaded path ledger](provenance/batches/B01_2026-11-21/BATCH01_FILE_LEDGER.csv)
- [Batch 01 source details, attribution and reconstruction instructions](provenance/batches/B01_2026-11-21/README.md)
- [Restore all 19 originals and check SHA-256 locally](code/tools/verify_batch01_reconstruction.py)

Third-party GEO/CHiCAGO source files retain original attribution. `GSE270828_Full_Guide_Reference_FOR_REVIEW.zip` contains guide/source metadata rather than the complete single-cell expression matrix. Further Colab input/results batches and final S01–S12 binary uploads remain pending; the local Stage 29 release package has been prepared and its source-verification scope is documented below.

## Stage 29 — raw-source verification and final supplement audit (22 November 2026)

The **HS1-11** original HAR design + 6 culture source files passed **7/7 recorded SHA-256/size checks**, and **342/342** reconstructed reporter effects agreed with corrected reference values (maximum difference 9.71e-17). The **WTc** six original culture files passed **6/6 recorded SHA-256/size checks**; **342/342** effects agreed exactly (maximum difference zero). The resulting **114 matched cuts across 66 HARs** show Pearson r = **0.543083**, Spearman ρ = **0.627249**, **86/114** matching mean directions and **48/114** stable six-culture directions. These are two **human iPSC backgrounds**, with human/chimp denoting **reporter DNA alleles**.

The Stage 27 coordinate recheck reproduced **304/304** archived eligibility flags relative to *secondary-catalog* ±300 bp windows (300 inside, 4 outside). Of 180 screened HAR labels, 152 have **provisional** sequence crosswalks and 28 remain unresolved. Original Noble author **hg38/panTro6 target BED authentication remains incomplete**; four secondary-window deviations must not be called experimentally confirmed guide errors.

[Detailed Stage 29 source and deposition status](provenance/verification/STAGE29_FINAL_AUDIT_2026-11-22.md). The 12 curated Excel workbooks, 13 final PNG figure panels, vector PDFs and validation ZIP reports are assembled in an offline author release, but **their binary uploads have not been independently confirmed in this repository**. Do not cite this repository as a complete final-data deposit until the actual files are uploaded and hash-checked. Original source GEO data remain accessible through the accession links above.


## Stage 31 — final manuscript-aligned statistical release

The statistical manuscript was reconciled with supplementary S01–S17. Data-derived Figure 9 and updated Figures 6–7 use the final Stage 30 clustered estimates and calibrated nulls. The unpublished Word manuscript is **excluded** from this repository's public deposit. [Stage 31 release notes and authoritative SHA-256](releases/stage31/RELEASE_NOTES.md) explain provenance, manuscript caveats, and upload verification. The repository owner may use [issue #1](issues/1) to attach the verified 9.2 MB public-only ZIP (no unpublished manuscript); the repository [Actions importer](.github/workflows/stage31-import-issue-attachment.yml) checks the owner, exact SHA-256 and every source file before committing release materials. Do not claim this is a complete public release until the binary files can be opened under `supplementary/` and `figures/`.
