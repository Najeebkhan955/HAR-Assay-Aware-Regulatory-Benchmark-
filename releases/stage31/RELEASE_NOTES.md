# HAR supplementary statistical evidence and figure release

**Prepared:** 22 November 2026  
**Repository:** [HAR Assay-Aware Regulatory Benchmark](https://github.com/Najeebkhan955/HAR-Assay-Aware-Regulatory-Benchmark-)  
**Release type:** Retrospective computational benchmark; not an original wet-lab/causal validation  
**Deposit status:** Release documentation and importer are public; **binary supplements/figures are pending until the authenticated owner ZIP is committed and verified.**

## Prepared immutable package

`HAR_Supplementary_Data_and_Figures_Public_Release.zip`

- ZIP bytes: **9,185,548**
- ZIP SHA-256: `eaa92a8520171625f90fa8b45438c168a91772c7cf2e1038b9a4d878e0f1588c`
- ZIP integrity: **PASS**; exactly 60 public-safe members.
- Contents: **no unpublished manuscript**; **S01–S17** Excel workbooks; nine main figures each in PNG+PDF; ten Extended Data figures each in PNG+PDF; legends, QC/requirements files, and complete per-file SHA-256 manifest.
- Expected GitHub deposit: **60 public files**: 59 checked entries plus one SHA-256 manifest.

## Owner-controlled upload without a personal access token

**[Upload issue #1](https://github.com/Najeebkhan955/HAR-Assay-Aware-Regulatory-Benchmark-/issues/1)** is configured with an owner-only GitHub Actions import workflow. The repository owner attaches the exact prepared 9.2 MB public-only ZIP in a new issue comment and adds `PUBLISH SUPPLEMENTS` on its own line. The workflow rejects a different archive SHA-256, checks internal per-file checksums, excludes the unpublished DOCX, writes only whitelisted release paths, and commits to `main` using the repository's GitHub Actions token.

**[Importer workflow](../../.github/workflows/stage31-import-issue-attachment.yml)**

**Do not attach the older complete submission ZIP:** it contains an unpublished Word manuscript. Attach only the safe public-only ZIP above. The above is a prepared workflow, not proof of a successful upload. Check the `Actions` run and GitHub file tree to establish publication. Do not cite a complete data deposition until files exist.

## Scientific reproducibility summary

- Perturbation screen: **180 HAR labels**; 152 provisional ortholog mappings; 28 unmapped with available secondary sequence; original Noble target BED files unverified.
- Secondary-guide coordinate QC: 300/304 within ±300 bp of provisional targets; 4 outside, **not** confirmed design errors.
- Original WTc N2 reporter source checks: six compressed GEO files passed user-executed Colab SHA-256 integrity checks; 342/342 cut-by-culture effects reproduced.
- HS1-11 N2 source checks: seven original reporter files (design + six RNA/DNA) passed user-executed Colab checksum checks; 342/342 cut-by-culture effects reproduced.
- Matched reporter evidence: 114 cuts across 66 HARs; Pearson r = 0.543083, Spearman rho = 0.627249, Lin's CCC = 0.5242; 86/114 matching mean directions, 48/114 strictly consistent over six cultures. Technical cultures are not independent donor backgrounds.
- Six-TF locus tests: **0/152 BH-FDR positive**; 20/28 opposite-signed multi-cut HARs not enriched under permutation grouping null (Stage30 p approx 0.492). Findings should not be interpreted as validated TF occupancy or HAR→gene causal mechanisms.

## Primary public experiments

- [Noble perturbation GSE270828](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270828)
- [Human/chimp sequence lentiMPRA GSE110758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110758)

## Release and disclosure boundaries

Full original GEO single-cell matrices, the original Noble `hars_hg38.bed`/`hars_panTro6.bed`, and a repository DOI are **not** represented as uploaded or verified. Reuse of source assets requires attribution to primary studies. Author names, funders and journal-specific submissions must be author-confirmed. The journal-style public ZIP is separately available in the ongoing manuscript-preparation conversation prior to GitHub import.

## Journal-style file names

The public archive uses supplementary identifiers S01–S17 with descriptive discipline-specific filenames (e.g., `S06_Integrated_Source_Quality_Control.xlsx`, `S10_Independent_Human_Cell_Line_Replication.xlsx`). Historical processing-stage labels are not used in archive filenames. `FILE_NAME_CROSSWALK.csv` records every old and new filename and source SHA-256; each Excel workbook remains byte-identical to the reviewed input. `FILE_SHA256_MANIFEST.csv` covers all other public release files.
