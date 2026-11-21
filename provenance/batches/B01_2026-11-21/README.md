# HAR input and results — Batch 01 archive ledger

This folder is the **provenance and upload status ledger** for the first 19 files provided by the study researcher on 21 November 2026.

- **19/19** local originals passed SHA-256 hashing and format checks (ZIP CRC or gzip integrity sample reads as applicable).
- **8/19** have already been directly committed as original binary files to the GitHub repository.
- **4/19** larger derived analysis/Colab ZIPs are awaiting an appropriately sized upload method and will not be falsely labeled uploaded.
- **7/19** external GEO/raw or comparative chromatin source files are registered with exact SHA-256 checksums and source links, **not rehosted** pending source-specific redistribution and citation review.

The complete details, category and current status for each file are in [BATCH01_FILE_LEDGER.csv](BATCH01_FILE_LEDGER.csv).

## Data identity and scientific cautions

The GSE270828 *Full Guide Reference* ZIP contains a feature README and guide metadata JSON; **it is not a full perturbation single-cell expression matrix**. Stage 7 and 8 archives contain diagnostics and rescue outputs rather than comprehensive final single-cell data.

Comparative CHiCAGO interaction tables are external published chromatin observations, not independently demonstrated endogenous target links for the user's selected HARs. Human–chimp MPRA barcode files reflect original GSE110758 reporter measurements; different fragment cuts and genomic coordinate systems need careful matching.

## Roadmap

1. Finish staged upload of four large derived ZIPs.
2. Record author/publisher reuse rights and citations for GEO/chromatin source tables. Wherever rehosting is restricted or unnecessary, preserve a checksum plus original source URL and automated downloader.
3. Continue batches in chronological order; place verified user-run outputs in `results/colab_runs/`, assistant-generated independent secondary audits in `results/independent_audits/`, and third-party metadata under `data/`.
4. Curate 10–12 publication supplementary files *after* full input/result reconciliation, not by simply renaming all historical ZIPs.

**No uploaded file is claimed to be generated from scratch by a single party.** Provenance is intended to make code and data reproducible and properly attributable.
