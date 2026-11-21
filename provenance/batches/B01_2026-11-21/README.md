# Batch 01 — Upload completion and source provenance (19 files)

**All 19 uploaded source files are represented on GitHub and independently integrity-verified.**

- **19 original user-provided files** classified, inspected and SHA-256-hashed.
- **53 Git blobs**, including split parts of large files, verified against locally computed original-file blob hashes: **53/53 exact**.
- All original ZIPs passed CRC checks. For gzipped GEO files, gzip format and first records were inspected (complete compressed bytes were verified through the Git blob hashes).
- Source catalogs and stage results are kept as distinctly attributed original Colab-run outputs, independent derived analyses, and external published experimental data.

[See 19-file authoritative status ledger](BATCH01_FILE_LEDGER.csv).

## Canonical archival layout

- `data/01_public_GEO/GSE110758/`: original publicly deposited MPRA experimental design and three uploaded RNA/DNA samples (split files).
- `data/02_source_reference_summaries/`: early GSE270828 feature and guide metadata export. **Not a raw scRNA expression matrix**.
- `data/03_pal2025_chromatin/`: three original publicly obtained CHiCAGO interaction tables (split files).
- `results/01_colab_original_archives/`: Stages 2, 6, 7, 8 and 9 user-run Colab results.
- `results/02_independent_audits/`: Stage 2, 4, 5 and 9 secondary audited analyses.
- `results/03_stage_workbooks/`: Stage 9 editable independent QC workbook.

Some earlier uploads exist under alias folders `results/colab_runs/`, `results/independent_audits/`, `data/01_external_sources/`; **do not double-count these as different experiments**. The paths above are the canonical current data locations and should be used for downstream reporting.

## How to reassemble the split original files

For each directory ending with `.parts`, concatenate `part-00001`, `part-00002`, ... in **numeric filename order**, without separators or line breaks. E.g., in Linux/macOS/Git Bash:

```bash
cat data/01_public_GEO/GSE110758/GSE110758_HAR_design.tsv.gz.parts/part-* > GSE110758_HAR_design.tsv.gz
sha256sum GSE110758_HAR_design.tsv.gz
```

Match the final original-file SHA-256 recorded in the [ledger](BATCH01_FILE_LEDGER.csv). Other files are reconstructed similarly. The part files are binary data, not tabular rows and must not be individually parsed as TSV. The exact path and SHA-256 values are part of the provenance contract.

**Important:** NCBI GEO and comparative 3D-chromatin tables are third-party source data. We preserve their original identity and source citation. Redistribution/derivative usage must comply with the primary source's policies and journal requirements; rehosting does not establish authorship or experimental ownership. GEO source record and Pal et al. attribution are recorded in the ledger.

## Scientific interpretation

- The guide-reference ZIP is source metadata, not a complete perturbation count matrix.
- CHiCAGO associations are chromatin interaction evidence, not validated HAR-to-gene causal assignments.
- GSE110758 MPRA data are reporter assays of constructed alleles; human/chimp allele direction is not endogenous human–chimp tissue expression.
- Stage 9 QC outputs include the cell metadata coverage and are not individual donor-level biological replicates.

## Next batch protocol

For each incoming batch: preserve original name in provenance ledger, compute SHA-256 and ZIP/GZIP integrity, record exact canonical GitHub path/part count, verify remote Git blob hashes, and distinguish user-run Colab outputs from subsequent independent analyses. Final curated supplementary S01–S12 will be assembled only after all stages are reconciled.
