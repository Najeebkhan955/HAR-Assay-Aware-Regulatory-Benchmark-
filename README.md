# HAR Assay-Aware Regulatory Benchmark

**Repository overview.** Reproducible, assay-aware benchmarking of human accelerated region (HAR) regulatory evidence across primate neural development, integrating public CRISPRi perturbation data, comparative sequence and transcription-factor motif analyses, chromatin contacts, and human–chimp MPRA reporter assays.

**Research objective.** Evaluate which regulatory conclusions are supported across assays, distinguish construct-level reporter effects from endogenous HAR–gene mechanisms, and document limitations from genomic coordinate matching, fragment selection, replication, and statistical calibration.

**Status (21 November 2026).** Stage-wise recovered code archive is **online and checksum-verified**: 108 source-code files (42 Jupyter notebooks, 65 Python scripts, and 1 shell script). Public GEO raw data, derived results, images and 10–12 curated final supplements will be uploaded and audited in separate stages, after checking attribution and redistribution restrictions.

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

Third-party GEO/CHiCAGO source files retain original attribution. `GSE270828_Full_Guide_Reference_FOR_REVIEW.zip` contains guide/source metadata rather than the complete single-cell expression matrix. Further Colab input/results batches and the final curated 10–12 supplementary files remain pending.
