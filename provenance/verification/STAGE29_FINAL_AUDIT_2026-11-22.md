# Stage 29 — Original-data and supplementary release audit (2026-11-22)

## Verified project scope

The secondary HAR study is an **assay-aware evidence benchmark**; it does not establish a new causal HAR–transcription-factor–gene mechanism.

- The original screen includes **180 HAR labels**. A secondary sequence/orthology catalog yields **152 provisional mappings** and **28 unresolved**. The original author BED files `hars_hg38.bed` and `hars_panTro6.bed` were **not recovered from checked public repositories**; original-target coordinate authentication remains **NOT VERIFIED**.
- Stage 27 independently recalculated **304/304** archived guide-to-secondary-window QC flags: **300** within ±300 bp of provisional cores, **4** outside. Those four flags are **not confirmed guide-design errors**.
- Original GSE110758 WTc N2 reporter sources: executed Colab validation reports **6/6 compressed RNA/DNA files** passing recorded SHA-256 and size checks; **342/342** reporter cut-by-culture effects matched the archived corrected Stage 24 reference exactly (maximum absolute difference **0**).
- Original GSE110758 HS1-11 N2 reporter sources: executed Colab validation reports **7/7 files** (shared HAR design plus 3 RNA and 3 DNA) passing recorded SHA-256 and size checks; **342/342** effects matched corrected reference within floating point error (maximum **9.71445146547012e−17**).
- Cross-background comparisons involve **two HUMAN iPSC backgrounds**, not human versus chimp host cells. Across **114 matched cuts** from **66 HARs**, Pearson `r=0.5430829699495423`, Spearman `rho=0.6272493469148053`; **86/114** matched two-line mean directions, **48/114** matched all six culture-level directions.
- These verify **existing reporter calculations**, not an independent biological experiment. The downloaded raw gzip bytes are not preserved inside the small shared diagnostic ZIPs; the Colab reports document source-file checks.

## Local release assembled (not yet deposited in this repository)

The author prepared a compact **Stage 29 release package** containing final Excel **S01–S12**, 8 main original-data PNG figures, 5 Extended Data PNG figures, data-derived vector PDFs, source diagnostics, Stage 27 coordinate review, a SHA-256 file manifest, and a readme. Archived Excel S09, S10, and S12 were copied and revised to correct *historical* source-verification/public-repository status, with a `FINAL_VERIFICATION_2026` tab; other workbooks were retained. The original large S12 source archive and full original GEO raw files are **not asserted as public deposits here**.

**Deposition gate:** The GitHub README and this verification record are deposited now. The binary S01–S12 workbooks, figures and downloadable verification ZIPs must be uploaded to the corresponding repository folders, followed by SHA-256 checks, before describing the GitHub repository as a complete release. This repository currently hosts earlier source code and a staged raw-data batch; it is **not yet the complete final supplemental archive**. An unpublished Word manuscript should not be claimed to be publicly deposited.

## Original sources and attribution

- Noble study/GSE270828: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270828
- lentiMPRA GSE110758: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110758
- Original Noble authors' code: https://github.com/NoonanLab/Noble_et_al_2024
- Author code archive: https://doi.org/10.5281/zenodo.12600310

## Recommended repository release paths

- `supplementary/S01_…xlsx` through `supplementary/S12_…xlsx`
- `figures/Figure_01_FINAL_ExactData.png` through `figures/Figure_08_FINAL_ExactData.png`
- `figures/Figure_ED01_…png` through `figures/Figure_ED05_…png`
- `verification/HAR_HS111_Diagnostic_Report_ALWAYS.zip`
- `verification/HAR_WTc_Diagnostic_Report_ALWAYS.zip`
- `verification/HAR_Stage27_Diagnostic_ALWAYS.zip`
- `STAGE29_SHA256_FILE_MANIFEST.csv`

**Limitations:** no authenticated Noble source BED, no newly established endogenous TF occupancy, no FDR-significant 152-HAR motif signal, and no automatic right to redistribute third-party raw GEO files.
