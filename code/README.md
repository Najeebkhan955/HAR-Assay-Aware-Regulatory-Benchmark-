# HAR research — code archive

This folder contains **95 verified code files** recovered from the user's available Stage 1–25 research archives: **39 Jupyter notebooks**, **55 Python scripts** and **1 shell script**. Some files intentionally duplicate source bytes in archival and reproducible layouts.

## Find the appropriate code

| Path | Contents |
| --- | --- |
| `01_workflow/original_stage_notebooks/` | Recovered original notebooks: pilot, Stage 2, 7, 8, 9, 10/10C/10D and 12 |
| `01_workflow/stage16_original_full_package/` | Stage 16 comparative ortholog/TF pipeline code |
| `01_workflow/all_recovered_stage_versions/` | Original, fixed and intermediate scripts/notebooks, primarily Stages 19–23 |
| `01_workflow/archived_notebooks/` | Preserved archived notebooks |
| `01_workflow/final_reanalysis_src/` | Stage 23–25 benchmarking, QC, figure production, independent-line reanalysis and manifest scripts |
| `01_workflow/tests/` | Invariant and synthetic tests |
| `reproduce_stage25/` | Preserved executable Stage 25 folder layout, including notebooks, src, legacy scripts and tests |
| `02_data/raw_GEO_original_uploads/` | Downloader script for original GEO files; this code upload **does not include** original raw data |
| `upload_tools/` | Legacy Git helper script |

See [the SHA-256 and Git blob hash manifest](CODE_HASH_MANIFEST.csv) for exact paths, byte sizes and integrity digests.

### What this repository does and does not establish

This is an archival code upload, **not** a claim that all experiments were rerun or all original data were archived. Complete original Stage 1–6 source notebooks, some Stage 11–15 original code, the original target BED, and full single-cell raw matrices were not recovered from the available archive. Results and source data will be added in separately verified batches.

Reproducible assay-aware secondary analyses assess provisional ortholog mappings, TF-motif null calibration, guide reproducibility and human/chimp MPRA reporter-fragment agreement. **No new causal HAR–TF–gene mechanism is experimentally established.** Do not report missing/unedited analyses as completed. Original publication/provenance and third-party redistribution permissions must be respected.

### Integrity verification

The 95 files were compared byte-for-byte at the Git blob hash level against the original archived code before this README was uploaded: **95/95 exact matches**, including notebooks. Local JSON notebook validation and Python AST parsing passed; this does **not** mean every notebook has been executed against its live external inputs.
