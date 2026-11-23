# Evidence-graded molecular benchmarks for colorectal cancer-on-a-chip

Code, frozen gene panels, protocols and metadata for:

> Mohammed AE, Shaikh S, Cetin Y. *Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics.* Bioinformatics Advances (submitted, 2026).

Release **v1.2** (November 2026) is the version described in the submitted manuscript, including the analyses added in response to pre-submission review. Archived at Zenodo: **[DOI]** · GitHub: **https://github.com/AbdurehmanEM/CRC_molecular_benchmarks**

## What this repository contains

| Folder | Content |
|---|---|
| `scripts/` | All R (and two Python) scripts, numbered in run order. `00_run_all.R` runs the main pipeline. |
| `scripts/revision/` | Post hoc analyses added in response to pre-submission review (scripts 45 and 50–67). |
| `frozen_panels/` | The two frozen panels exactly as used for external validation (`.rds`), their SHA-256 hashes, and a readable CSV of every gene, coefficient, intercept and scaling value. |
| `protocol/` | Dated protocol versions (v0.1–v2.2) and the amendment log. Supplementary Table S14 lists every change made after a related result was known. |
| `metadata/` | Input manifest with SHA-256 checksums of all original data files (Supplementary Table S16), dataset eligibility audit, and the literature-search logs (PubMed, Scopus). |
| `supplementary/` | Supplementary Tables (xlsx) and Supplementary Methods as submitted. |
| `environment/` | R version, package versions and `sessionInfo()` of the machine that produced the results. |
| `RELEASE_MANIFEST.csv` | SHA-256 and size of every file in this release. |
| `SCRIPTS.md` | One line per script, generated from each script's header. |

No patient-level data are redistributed. All data are public and are downloaded by the scripts (see below).

## Study in one paragraph

Ten transcriptional programs across five compartments (tumor-associated epithelium, fibroblasts, monocytes/macrophages, blood vascular endothelium, CD8 T cells) were defined before analysis; eight met gene-coverage criteria. For each program, the expressed genes were split into evaluation genes (which define the program target) and predictor genes (eligible for the panel), with no overlap. Sparse panels (relaxed elastic net, 1–4 genes; 30 genes in total) were chosen by donor-grouped cross-validation in the discovery cohort (GSE178341, 60 patients), frozen with a recorded SHA-256 hash, and applied without refitting to three independent single-cell cohorts (48 patients) and, after adjustment for spot composition, to 13 patients from three spatial cohorts. Each program received an evidence grade. The panels are **molecular benchmarks**: they have not been measured on any chip.

## Data sources (all public)

| Role | Dataset | Source |
|---|---|---|
| Discovery (scRNA-seq) | GSE178341 (Pelka et al., 2021) | GEO |
| Validation (scRNA-seq) | GSE132465 (Lee et al., 2020) | GEO |
| Validation (scRNA-seq) | GSE166555 (Uhlitz et al., 2021) | GEO |
| Validation (scRNA-seq) | GSE200997 (Khaliq et al., 2022) | GEO |
| Descriptive (scRNA-seq, <10 donors) | GSE144735 (Lee et al., 2020) | GEO |
| Spatial (Visium) | Valdeolivas et al., 2024 | Zenodo record 7739700 |
| Spatial (Visium HD) | GSE280318 (Oliveira et al., 2025) | GEO |
| Spatial (FFPE Visium) | GSE226997 (Park et al., 2024) | GEO (manual download step, see below) |

`metadata/input_manifest.csv` gives the SHA-256 of every original file used, so you can confirm that your downloads match ours.

## Requirements

- R 4.6.1 (R ≥ 4.3 should work); package versions in `environment/package_versions.csv`.
- Key packages: Seurat, edgeR, glmnet, Matrix, data.table, spacexr (Bioconductor version; provides `createRctd()`), metafor, ggplot2, patchwork, openxlsx, digest, rhdf5 or hdf5r.
- Python 3 with openpyxl (only for `65_merge_supplementary_workbook.py` and `67_scopus_screening_decisions.py`).
- About 64 GB RAM recommended for the full single-cell steps; the full pipeline takes roughly 12 hours.

## How to reproduce

1. Clone or download this repository and set the project root, either in `scripts/revision/50_revision_config.R` or with the environment variable `CRC_PROJECT_ROOT` (some older scripts also contain a `setwd()` line with the authors' path; change it to yours).
2. Run the download scripts (stages 1–5). **GSE226997** must be downloaded manually from GEO into `03_original_data/GSE226997/` (instructions in the script header).
3. Run `00_run_all.R`. Stages are skipped when their output already exists; each stage is logged with its duration.
4. Run the revision analyses in `scripts/revision/` in numerical order (see `scripts/revision/README_revision.md`).
5. Compare your outputs with ours using `26_compare_outputs.R`.

### Important: keep the supplied frozen panels

Re-running the panel-freezing stage regenerates `frozen_panel_v2.rds` with a new timestamp, which changes its SHA-256. Every validation script checks the hash and **stops by design** if it differs. To reproduce the published numbers, copy the files from `frozen_panels/` into `04_processed_data/frozen_panel/` and start the pipeline at the first validation stage.

| File | Programs | SHA-256 |
|---|---|---|
| `frozen_panel_v2.rds` | 7 programs, 26 genes | `a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb` |
| `frozen_panel_LYM_v1.rds` | CD8 IL2/STAT5, 4 genes | `7dc3256d97477a796fcacd1df71021740888f7978a01ea36caf70787da7f513f` |
| `frozen_panel_v1.rds` | superseded penalized version (protocol v0.8); not used for validation, kept for transparency | `f5be7339aa1a4e09b479335d6f4f4d50d3f680ac52648d1495053abbff72a0d4` |

## Main results (frozen discovery scaling, primary endpoint)

| Program | Evidence category |
|---|---|
| Endothelial TNF/NF-κB activation | Independently supported with spatial evidence (single external cohort, n = 11) |
| Epithelial hypoxia | Independently supported in single-cell data (not distinguishable from expression-matched random panels) |
| Epithelial apoptosis, myeloid IFN-γ response, myeloid inflammatory response, stromal TGF-β response, CD8 IL2/STAT5 | Exploratory |
| Endothelial angiogenesis | Unsupported |

Full numbers are in the Supplementary Tables.

## Known limitations of the code

- Scripts 18 and 18b are a mandatory pair (18 saves unlabelled RCTD weights; 18b labels them).
- The published RCTD reference selected cells in file order; script 42 rebuilds it with random selection as a sensitivity analysis.
- Some scripts were patched after producing their outputs during the study; the versions here are the corrected versions, and `26_compare_outputs.R` checks that they reproduce the published outputs.
- A complete re-run of the pipeline in a clean computing environment was not performed. Script 63 was run for its first part only, which records SHA-256 checksums of all input data and frozen panels and the R and package versions; its clean-room part is untested.

## Use of AI tools

We used Claude Opus 5.5 to assist with writing and debugging R code for data analysis and figure preparation. All the AI-generated code was reviewed, executed locally, and checked against intermediate outputs. Claude was also used to check the manuscript for consistency and numerical errors against analysis outputs and a source publication. Grammarly was used for language editing. We reviewed and verified all the text and take full responsibility for the content.

## Licence

All contents of this repository (code, frozen panels, tables and documentation) are released under the MIT licence (see `LICENSE`). The input datasets are not redistributed and remain under the terms of their original sources.

## Citation

Please cite the article and this archive (see `CITATION.cff`).

## Contact

Abdurehman Eshete Mohammed, Biotechnology Research Group, Climate and Life Sciences, Marmara Research Center (MAM), The Scientific and Technological Research Council of Türkiye (TÜBİTAK), Kocaeli, Türkiye. ORCID 0000-0002-9459-8378.
