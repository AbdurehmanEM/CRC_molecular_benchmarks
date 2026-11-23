# Work log — CRC Biomarker Project

Write 3 lines at the end of each day: **Ran** (what you ran) · **Decided** (decisions and why) · **Next** (tomorrow's first task).
Record problems and how you fixed them — they become Methods details and reviewer answers later.

---

## Thu 8 Oct 2026 — Day 1 (Phase 1: Feasibility)

**Ran**
- Created project folder `F:/CRC_Biomarker_Project 2/CRC_Biomarker_Project/` (folders 01–09), manuscript template `09_manuscript/manuscript_v0.1.docx`, protocol `01_protocol/protocol_v0.1.md`, audit spreadsheet `02_metadata/dataset_sample_eligibility_audit.xlsx`.
- Installed packages and recorded versions in `06_results/software_versions.csv`: R 4.6.1, Matrix 1.7.5, data.table 1.18.6.1, ggplot2 4.0.3, patchwork 1.3.2, Seurat 5.5.1, glmnet 5.1, digest 0.6.39, openxlsx 4.2.9, R.utils 2.13.0, edgeR 4.10.5, GEOquery 2.80.0, scDblFinder 1.26.7, SingleCellExperiment 1.34.0, msigdbr 26.1.1.
- Started download of GSE132465 (script 01).

**Problems and fixes**
- Working directory was one level too high (`.Rproj` in outer folder) → use `setwd("F:/CRC_Biomarker_Project 2/CRC_Biomarker_Project")` at the start of each session, or recreate the RStudio project in the inner folder.
- renv step in old setup script activated an empty project library → removed renv (`renv::deactivate(clean = TRUE)`), replaced script without renv. Lock file to be added before submission.

**Decided**
- R is the primary environment (Seurat, edgeR, glmnet).
- Dataset roles remain provisional until the audit is complete.
- GSE226997 has only 4 patients → identify a second spatial cohort before seeing any results.
- Working schedule: 7 h/day, 6 days/week; target submission in ~9 weeks.


**Pilot result (afternoon)**
- Download complete; 4 files, checksums in `02_metadata/GSE132465_checksums.csv` (count matrix 123.0 MiB, matches GEO).
- Pilot load v1 failed: temporary unzipped matrix (>4 GB) exceeded the FAT32 limit on F:, file truncated at line 33,593. Script 02 changed to stream directly from .gz; results discarded and re-run.
- Pilot load v2 passed: 4 donors (SMC01–04, tumor + adjacent mucosa), 15,245 cells, 33,694 genes, raw integer UMI counts (max 24,350), 0 duplicated genes, 100% cell–annotation match, donors recoverable, 315 MB in memory.
- Script 03: fixed integer/double type error in qc_summary; reproducibility confirmed (hash 073d78b16e7b68cd74dd4d9fec6eb993 on two runs).
- Compartment mapping fixed from the published subtype table: Endothelial = Tip-like, Stalk-like, Proliferative ECs; Lymphatic ECs and enteric glia excluded from core; Stromal = Stromal 1–3 + Myofibroblasts; Pericytes/smooth muscle = mural (context); epithelial malignant = CMS1–4 labels.
- Tumor cells per donor after QC flags (Endo / Epi / Myeloid / Stromal): SMC01 37/959/163/90; SMC02 32/244/610/62; SMC03 31/665/295/175; SMC04 69/173/807/36. All 4 donors ≥30 cells in all 4 compartments.
- All tumor epithelial cells carry CMS labels (published malignant annotation); no independent malignancy check yet.

**Decision: pilot PASSES → proceed to full 23-donor run.** Watch: endothelial counts are borderline (31–37 in 3 of 4 donors); expect some donors below 30 in the full cohort.

**Full GSE132465 run (late afternoon)**
- Full load: 63,689 cells (matches published total), 33,694 genes, 100% matched, raw integer counts, 1.45 GB. Checks saved as `06_results/pilot/full_load_checks_GSE132465.tsv`.
- Full QC reproducible (hash 562dd2e3517d9b721152b25208608142, 4 identical runs).
- Eligible tumour donors at ≥30 cells: epithelial 21, myeloid 20, stromal 11, endothelial 9 (all four: 9). At ≥10: 22 / 23 / 17 / 16.
- Authors pre-filtered cells (max mito ≈20% in all samples). Tumours SMC16/18/21/22: 4–6× higher UMIs/cell, epithelium-dominated → possible batch.

**Decision (protocol v0.2): Option A.** GSE178341 (Pelka, 62 patients) = discovery; GSE132465 = primary validation; GSE144735 (KUL, 6 patients) = secondary validation. Firewall: no program scores/models on validation cohorts until panel is frozen.

**Next**
- Run `04_download_GSE178341_GSE144735.R` (overnight if large); check file sizes vs 4 GB FAT32 limit.
- Read Lee 2020 Methods: filtering thresholds, chemistry of SMC16/18/21/22, treatment, MSI.
- Read Pelka 2021 Methods: CD45 sorting/fractions, chemistry, metadata columns.

**Cohort expansion (evening)**
- Protocol v0.3 → v0.4 (before opening any new data): validation = GSE132465, GSE144735, GSE200997, GSE166555 (~57 patients); spatial = Valdeolivas Visium (14 sections), GSE280318 (5 patients, Visium/Visium HD/Xenium + matched single-cell), GSE226997 (4, if storage allows). Joanito 2022 and GSE146771 considered and not included. Cohort list closed.
- Scripts ready: 04 (GSE178341 + GSE144735), 05 (Valdeolivas + GEO cohorts, skips files > 3.9 GB).

**Next (Day 2 start)**
1. `setwd("F:/CRC_Biomarker_Project 2/CRC_Biomarker_Project")`
2. Run script 04, then script 05; paste file lists/sizes and any "too large" lines.
3. While downloading: read Pelka 2021 and Lee 2020 Methods → Paper_notes; start literature search (group 2: CRC-on-chip).

**Next (was)**
- Confirm download finished (checksums printed; files in `03_original_data/GSE132465/`).
- Run `02_pilot_load_GSE132465.R`, paste output for review.
- Read Lee et al. 2020 (Nat Genet) → fill `Paper_notes` sheet.

---

## Fri 9 Oct 2026 — Day 2 (Phase 1: discovery audit)

**Ran**
- Script 04: GSE178341 (h5 1,203.6 MB, metatables, cluster table) and GSE144735 (raw UMI 56.6 MB, annotation) downloaded; all < 4 GB.
- Script 06: GSE178341 audit (metadata + h5 structure only).

**GSE178341 audit findings**
- 370,115 cells; 62 patients; 100 specimens = 64 tumour (C130, C171 have TA + TB) + 36 adjacent normal.
- h5 = 43,113 features × 370,115 cells, raw integer UMIs (slice min 1, max 3,094); barcodes match metadata 100%; 1,918 duplicated gene names → use Ensembl IDs, sum duplicates for symbol-level output.
- MMRd 34 / MMRp 28; MGH + DFCI; 10x v2 (266,604) + v3 (103,511 cells).
- Processing: unsorted 231,737; CD45+ MACS 122,766; mixed 10,670; Live MACS 4,942 → sorting rule needed.
- Treatment status not in metadata → verify in Pelka 2021.
- 3 patients M1 at surgery; histology includes mucinous/medullary.
- Full sparse matrix ≈ 9 GB in RAM (764 M non-zeros) → pseudobulk will be computed by reading the h5 in cell chunks.

**Decisions pending (must be fixed before expression analysis)**
- Sorting rule (all processing pooled vs unsorted only) — from script 07 counts.
- Specimen rule for C130/C171 (TA + TB pooled per donor proposed).
- Epithelial malignancy rule (Pelka clusters are lineage states, not CMS labels).
- Myeloid definition harmonised as monocytes + macrophages (THP-1 analogue); DCs/granulocytes as context.

**Script 07 (metadata counts) → protocol v0.5**
- Rule B (unsorted + LiveMACS) chosen as primary; rule A (all fractions) gains little (≥20 cells: stromal 29 vs 27, endothelial 47 vs 42). CD45+ MACS fractions contain 23,167 epithelial and 1,264 endothelial cells (impure). 2 donors have no unsorted fraction.
- Epithelial stem/TA-like clusters cE01–cE03 are 69–89% tumour-derived; mature lineages 6–28%.
- v0.5: pooled TA/TB; harmonised compartments; min 20 cells; chemistry-balanced CV.

**Script 08 (pseudobulk)**
- Problem: `Dim` slot must be integer → fixed with `as.integer()`.
- Eligible donors at ≥20 cells after QC: epithelial 60, stromal 27, myeloid 60, endothelial 40. Hash 831d778f3233a96f864f2f8f08ac10af.

**Protocol v0.6:** 8 programs (Hallmark/Reactome), coverage rule, random half-split partition, target definition.

**Script 09 (programs + partition)**
- First run matched 0 genes (ID format) → script made robust (version-stripped Ensembl / symbol fallback). Paste truncation twice → switched to download + file.copy method for long scripts.
- Coverage: 7/8 programs retained. STR_ECM 77.4% < 80% → not retained (pre-specified rule; not changed). END_ANGIO small (15 evaluation / 8 predictor genes).
- 129 predictor genes removed (evaluation genes of other programs). Partition hash 90280a13b2629844e6aa913559286ca6.

**Protocol v0.7:** sample-wise CPM normalisation (frozen-transferable); modelling, baselines, freezing, validation-scaling rules.

**Next**
- Run script 10 (panel development, discovery only); review CV performance vs baselines and selection stability.
