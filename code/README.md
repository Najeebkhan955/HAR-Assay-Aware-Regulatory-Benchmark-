# HAR research — code archive

This folder contains **98 verified code files** recovered from the user's available Stage 1–25 research archives: **42 Jupyter notebooks**, **55 Python scripts** and **1 shell script**. Some files intentionally duplicate source bytes in archival and reproducible layouts.

## Find the appropriate code

| Path | Contents |
| --- | --- |
| `01_workflow/original_stage_notebooks/` | Recovered original notebooks: pilot, Stage 2, **Stage 3, Stage 6, Stage 7 (including Rescue)**, 8, 9, 10/10C/10D and 12 |
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

This is an archival code upload, **not** a claim that all experiments were rerun or all original data were archived. Complete original Stage 1, 4, 5 and other unrecovered early-stage source notebooks, some Stage 11–15 original code, the original target BED, and full single-cell raw matrices were not recovered from the available archive. Results and source data will be added in separately verified batches.

Reproducible assay-aware secondary analyses assess provisional ortholog mappings, TF-motif null calibration, guide reproducibility and human/chimp MPRA reporter-fragment agreement. **No new causal HAR–TF–gene mechanism is experimentally established.** Do not report missing/unedited analyses as completed. Original publication/provenance and third-party redistribution permissions must be respected.

### Integrity verification

The 98 files were compared byte-for-byte at the Git blob hash level against the original uploaded/archive source files: **98/98 exact matches**, including notebooks. Local JSON notebook validation and Python AST parsing passed; this does **not** mean every notebook has been executed against its live external inputs.

### Recovered code added 21 November 2026

The user supplied the original GSE270828 full-guide-reference notebook (Stage 3), public perturbation-effect source locator (Stage 6), and the Stage 7 RDS/GZIP rescue notebook. All three were added as exact original JSON notebooks, without normalizing their contents. Fifteen other supplied notebook copies were already byte-for-byte identical to the repository. The Stage 16B ZIP contained a Python script already present with an identical Git blob hash; only its results remain for the later **results/data upload** phase.

A Colab notebook that contains shell magics must be run in a compatible environment; generic Python AST parsing is not a valid runtime test for shell-magic cells.
