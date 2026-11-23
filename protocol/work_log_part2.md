
**GSE226997 + final 13-patient spatial synthesis (14 Oct 2026)**
- All 4 tar.gz downloaded (12.86/9.72/9.82/11.78 GB; parallel curl windows; tar -tzf listing clean, 21 entries each). Extracted only filtered_feature_bc_matrix.h5, tissue_positions_list.csv, scalefactors, lowres image.
- GEO metadata: 4 patients, Homo sapiens, primary CRC, FFPE Visium, Ajou University Hospital (AJOUIRB-OBS-2016-218). P1 F T3N0M0; P2 F T3N0M0; P3 M T3N1M1 (primary tissue); P4 M T3N1M0. **Treatment status not reported** → recorded as such; cohort-drop sensitivity available.
- Script 22 patched: per-section h5/pos_file columns in SECTIONS; `header = FALSE` for *_tissue_positions_list.csv; Read10X_h5 list handling. Spots after QC: P1 4,142 · P2 3,588 · P3 4,508 · P4 3,264.
- RCTD marker check: fibroblast 0.70–0.92, pericyte 0.58–0.74 OK; **endothelial only 0.47–0.51** (vs 0.74–0.88 in GSE280318); P3 is 89% malignant epithelium (myeloid 0.12, T/NK 0.14 weak).
- Probe-panel coverage: **MSX1 absent in P1/P2/P4** (angio evaluable in 10 of 13 patients); **HLA-B absent in all GSE226997 + GSE280318** (IFN-γ evaluable in 7).
- Script 23 fix: `rbindlist(fill = TRUE, lapply(...))` — MYE_IFNG has only 1 patient in the "additional cohorts only" slice → short row.

**FINAL spatial result — 13 patients (Valdeolivas 6 + GSE226997 4 + GSE280318 3)**
| Program | Pooled A (95% CI) | A | B | I² | Label |
|---|---|---|---|---|---|
| MYE_IFNG (n=7) | 0.241 (0.172–0.310) | 7/7 | 7/7 | 90.8 | Spatial support |
| END_TNFA_NFKB | 0.203 (0.133–0.273) | 13/13 | 13/13 | 97.0 | Spatial support |
| EPI_APOPTOSIS | 0.156 (0.105–0.207) | 13/13 | 9/13 | 94.6 | Spatial support |
| MYE_INFLAM | 0.147 (0.102–0.192) | 13/13 | 13/13 | 93.6 | Spatial support |
| STR_TGFB | 0.110 (0.072–0.149) | 13/13 | 11/13 | 92.1 | Spatial support |
| EPI_HYPOXIA | 0.055 (0.030–0.080) | 13/13 | 6/13 | 80.4 | Partial |
| END_ANGIO (n=10) | 0.015 (0.004–0.026) | 7/10 | 10/10 | 42.0 | Spatial support |
- Valdeolivas-only effects ≈ half the added-cohort effects (e.g. TNF/NF-κB 0.122 vs 0.271) → high I² is a magnitude/platform effect, direction replicates.
- Angiogenesis reaches both criteria only with the added cohorts, at exactly the 7/10 threshold, with B from endothelial-proportion quartiles and A not LOPO-robust → **evidence category remains Unsupported** (failed single-cell validation 2.11, 1.37–3.35). Script 20 confirmed all categories unchanged.

**Scripts/figures updated (14 Oct)**
- Script 23 writes `06_results/spatial_support_combined/`; derived `spatial_support_summary_main.csv` (main-analysis slice) now feeds script 20 and Figure 5.
- Script 21: Fig 5B/5C switched to the 13-patient combined files; panel C shape = cohort (criterion-B rule column no longer exists); `cpm` → `edgeR::cpm`. All 6 figures rebuilt as `*_updated` (PDF/TIFF 600 dpi/PNG). Figure 5 visually checked.

**Manuscript v0.4 (14 Oct 2026)**
- Methods: Spatial assessment rewritten for 3 cohorts (HD 7×7 → 56 µm pseudo-spots; spot QC; composition-quartile criterion B in unannotated cohorts; pooling of all 13 patients; coverage rule). Statistics subsection now describes the 12 sensitivity scenarios as **performed** (placeholder removed).
- Results: spatial subsection rewritten (13 patients, coverage limits, I² caveat, angiogenesis caveat); **new subsection "Sensitivity of the discovery model"** (82/85 fits MAE ratio < 1; gene identity not unique); evidence paragraph updated.
- Figure legends 1–6 written in full. **Table 2 inserted as a real Word table** (26 genes × gene/compartment/program/sign/evidence + interpretation limits).
- Data availability now lists GSE280318 and GSE226997.
- Remaining placeholders: hdf5r/rhdf5 versions, repository URL, Supplementary numbering, author/affiliation block.

**Next**
- Protocol v1.6 (document GSE226997 addition, treatment status not reported, coverage-driven exclusions).
- Literature comparison table from the key papers → Introduction ¶4 gap statement (Dai 2026 and Hachey 2021/2023 are the closest precedents and must be contrasted explicitly).
- Then: Introduction, Discussion (hMSC vs patient fibroblasts, THP-1 vs primary myeloid, HUVEC vs tumour endothelium, gene interchangeability), Conclusion, Abstract.

**Numerical audit of the manuscript (14 Oct 2026) — script 25_audit_manuscript_numbers.R**
Every numerical claim traced to its source file. Verified: frozen hash; 26 genes (4/4/4/2/4/4/4); evaluation-gene range 15–97; ECM coverage 0.774 (dropped); 129 removed predictors (1,186 expressed − 595 evaluation − 462 predictor); discovery donors 60/27/60/40; CV baselines (panel 0.41–0.84, equal-weight 0.43–0.84, conventional 0.54–1.02, larger 0.26–0.83, shrunk 0.84–0.96); compression cost IFN-γ 0.26 vs 0.41, inflam 0.27 vs 0.42, angio 0.75 vs 0.78, TGF-β 0.83 vs 0.82; selection stability 15/26 ≥75%, five genes at 1.00, ARID4B 0.26, FOSB 0.34; all per-cohort validation values and CIs; pooled SMC+KUL 0.46 (0.40–0.52) and 2.22 (1.49–3.35); all meta-analysis values, I² and pooled ρ (frozen); 15,904 spots / 12 sections / 6 patients; all 13-patient spatial values incl. Valdeolivas-only vs added-only; evidence table; LOPO (angio FALSE, rest TRUE); chip table 26 rows.
- **Calibration slopes verified** from `cv_performance_discovery*.csv`, column name is `cal_slope` (not `calib_slope`): v1 panel 1.45–5.61; v2 panel_refit 0.52–1.01. Both match the manuscript.
- **ERROR FIXED 1:** sensitivity fits were reported as 85 with 82 below 1; actual **78 fits, 75 below 1** (all-epithelial-states applies to 2 programs only; TGF-β not evaluable at ≥50 cells). Corrected in Results and Discussion. Exceptions unchanged (v2-only angio 1.115, TGF-β 1.204; partition seed 3 TGF-β 1.078).
- **ERROR FIXED 2:** HLA-B reported absent from "all GSE280318 and GSE226997 sections"; it is present in GSE226997 P3 (hence IFN-γ evaluable in 7 patients). Corrected to "all GSE280318 sections and three of four GSE226997 sections".
- **IMPRECISION FIXED 3:** offset diagnostic described as "−0.2 to −0.7"; actual values are negative in 24 of 28 cohort × program combinations (range −0.009 to −1.079, four slightly positive). Reworded.
- Noted, not changed: I² "84–98%" is 84.4–97.5 rounded; pooled ρ quoted in Results are the frozen-scaling values.

**Journal decision (14 Oct 2026)**
- Briefings in Bioinformatics **rejected as target**: review journal, "does not publish articles that are pure original research"; research-type categories capped at 5,000 words (main text was 6,854).
- Target is now **Bioinformatics Advances** (OUP, fully open access). Original Article: ≤8 pages, abstract ≤200 words with no citations or abbreviations, ≤5 keywords, ≤50 references, ≤6 figures/6 tables, author–date references, figures as separate files (vector PDF preferred; no .docx/.jpg), required "Software and Data availability" statement, supplementary files ≤15 MB. APC $2,726 (20% ISCB discount).
- **TÜBİTAK ULAKBİM covers OUP** (2026–2028, hybrid + fully OA, quota 400 for 2026, first-come). Springer R&P covers hybrid only; Wiley 2026 allowance exhausted. Library must confirm MAM eligibility and remaining quota.
- Files: `manuscript_BioinfAdv_v1.docx` (main text 5,093 words: Methods condensed to 834) and `supplementary_BioinfAdv_v1.docx` (full Methods + S1–S10, Figs S1–S2). Main text still ~500–1,000 words over 8 pages → plan: move Figures 2 and 3 to Supplementary.

**Open before submission**
- Reproducibility: pipeline never run end-to-end from clean state; scripts 21/22/23 patched in place after producing outputs; script 12 superseded by 12b; no renv lock; GSE226997 download not scripted. Plan: review scripts 14/15/17/18/18b/20/23, then `00_run_all.R` + README + renv.lock + clean-room re-run + GitHub/Zenodo archive.
- Literature comparison table → Introduction ¶4 and Discussion comparison paragraph (both bracketed).
- References author–date from Mendeley; author block/ORCID/CRediT/funding/AI disclosure; Figure 4 axis widen to 0.10; supplementary numbering.

**Code review of scripts 14, 15, 17, 18, 18b, 20, 23 (14 Oct 2026)**
- No computational errors in 14, 15, 17, 18, 18b, 23. Verified in detail: pseudobulk indicator-matrix sums; MAE-ratio baseline (mean|y| is correct because the discovery-mean prediction is 0 on the frozen scale); bootstrap quantile indexing; 0-based CSC row indices in the hand-built dgCMatrix (script 17); DerSimonian-Laird tau-squared; two-thirds thresholds via ceiling(); role=="main" exclusion of the same-patient section. Script 23's "Valdeolivas only" slice reproduces script 19 exactly (TNF/NF-kB 0.122, 0.070-0.174) — independent confirmation.
- **BUG FOUND AND FIXED (script 20):** `SP` had been repointed to the 13-patient combined summary but `AP` still read `spatial_support/criterionA_by_patient.csv` (6 Valdeolivas patients), so the leave-one-patient-out column described a 6-patient test beside a 13-patient pooled effect. Repointed to `criterionA_patients_combined.csv` and re-ran. **Result: A_robust_all_LOPO is now TRUE for all seven programs, including END_ANGIO** (previously FALSE on 6 patients). Evidence categories unchanged (LOPO is not part of the v1.3 rule). Manuscript corrected: the claim "criterion A for angiogenesis was also not robust to leaving out single patients" was false and has been replaced; the remaining angiogenesis caveats (smallest effect 0.015; exactly 7/10; criterion B from endothelial-proportion quartiles) stand. Protocol v1.6 amended.
- **Disclosure items for Supplementary Methods (not errors):**
  1. Script 14's cluster assignment compares marker sets on absolute mean log-expression, so the 1.10 "margin" is not comparable between sets of differing size and baseline expression; GSE200997 already contributes supporting evidence only. Optional fix: z-score the score matrix by column before max.col, then re-run 14 and 15.
  2. Script 17 assigns Malignant_epithelial only to cE01-cE03 in tumour specimens; other tumour epithelium falls to NA and is absent from the RCTD reference, which likely explains the weak malignant/normal epithelial separation (normal-epithelium rho 0.06-0.64).
  3. The RCTD "Myeloid" reference type merges monocytes, macrophages, DCs and granulocytes (protocol v1.2), broader than the monocyte/macrophage compartment used to train the myeloid panels; criterion B myeloid quartiles use the broader definition.
  4. Script 14 applies no doublet filtering although the protocol mentions doublet information — apply scDblFinder or drop the claim.
- **Latent fragilities to fix before archiving:** scripts 18 and 18b are a mandatory pair (18 saves unlabelled weights; the `labels_verified` guard in 19 prevents misuse); script 23 uses `sqrt(mean(se^2))` for multi-section patients (RMS, conservative; identity here since every added patient has one section); `eb` lacks the `is.finite(diff)` filter that `vb` has; `condition == "Tumor"` in 15 vs `grepl("tum", ignore.case=TRUE)` in 14.

**Reproducibility package (14 Oct 2026)**
- Added `00_run_all.R` (27 stages, sentinel-based skipping, resolves scripts by prefix and ignores `_before_*` patch backups, logs duration per stage), `26_compare_outputs.R` (field-by-field CSV comparison with numeric tolerance + hash comparison for other files), `README.md` (requirements, ~12 h runtime, 8 data sources, manual GSE226997 step, 18/18b pairing, frozen-panel hash trap, known limitations), `renv_setup.R`, `LICENSE` (MIT, name to be filled).
- **`05_download_spatial_and_GSE200997.R` was missing from the project** and was restored from the session copy. The dry run surfaced it as SCRIPT NOT FOUND; without the manifest this gap would have been invisible until someone tried to reproduce the spatial analysis.
- Dry run now reports all stages as "skip (output exists)" except stage 26 (read-only audit, marked `__always_run__`). Manifest validated against the real project.
- Note for the archive: re-running stage 10 regenerates `frozen_panel_v2.rds` and changes its SHA-256 (the file stores a timestamp), after which every downstream script halts by design. A clean re-run intended to reproduce published numbers must keep the supplied frozen panel and start at stage 11.

**Submission assembly (14 Oct 2026)**
- **Figures renumbered to 4 main + 4 supplementary** to fit the 8-page Bioinformatics Advances limit: Fig 1 cohorts, Fig 2 validation forest (was 4), Fig 3 spatial 13 patients (was 5), Fig 4 evidence summary (was 6). Old Figs 2 and 3 became S1 (program variation) and S2 (panel baselines); new S3 (RCTD marker check, all sections) and S4 (spatial support by cohort group) drawn by script 27. All in-text references and the four main legends rewritten; the `Supplementary Fig. S[x]` placeholder resolved to S3.
- `27_build_submission_files.R` added: widens the Figure 2 axis to 0.10, rebuilds figures, draws S3/S4, assembles `10_submission/Supplementary_Tables.xlsx` (17 sheets, S1a–S10d) and copies figures to `10_submission/figures/` under final names (PDF vector for submission; TIFF 600 dpi; PNG).
- Fix: S1 initially read the audit workbook's README sheet (3 rows); corrected to Study_level (10), Sample_level (29) and Paper_notes (8).
- **renv abandoned.** `renv::init(bare = TRUE)` created an empty project library and broke the session (same failure as day 1). `renv::deactivate()` restored it; no lockfile was written; `renv_setup.R` deleted. Version evidence is `06_results/software_versions.csv` plus the per-stage `sessionInfo_*.txt` files, which record what actually ran. README rewritten accordingly.
- **Fresh-session dependency defect found:** `21_figures.R` used sparse-matrix subsetting but never attached `Matrix`; it only worked because an earlier script had attached it. Fixed in the script and documented in the README. Other scripts may share this; a per-script fresh-session test would confirm.
- Outstanding and only the author can supply: references (none in the manuscript yet), author block/affiliations/ORCID/CRediT/funding/AI disclosure/ethics, and the two bracketed passages (Introduction gap, Discussion comparison) which depend on the **`Literature_comparison` sheet of the audit workbook — currently 0 rows**.
- Unresolved minor: Figure 2 still emits "log-10 transformation introduced infinite values" (3x) after the axis widening; no plotted interval appears to be lost, source of the non-positive value not yet identified.
