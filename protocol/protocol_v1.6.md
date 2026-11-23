# Study protocol — v1.6 (14 Oct 2026)

**Working title:** Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics

> **Amendment to v1.5.** All sections of v1.5 remain in force unchanged except §10d (spatial support), which is extended below, and §13 (new, sensitivity analyses). No change affects the frozen panel (`frozen_panel_v2.rds`, SHA-256 `a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb`), the programs, the gene partition or the validation criteria.

## 10d (extension) — GSE226997 as executed

**Eligibility, verified from GEO sample metadata before analysis.** Four samples, GSM7089855–58 (patients P1–P4), *Homo sapiens*, "Primary colorectal cancer", formalin-fixed paraffin-embedded Visium, collected at Ajou University Hospital under review-board protocol AJOUIRB-OBS-2016-218. P1 female T3 N0 M0; P2 female T3 N0 M0; P3 male T3 N1 **M1**; P4 male T3 N1 M0. All four are included:
- P3 had metastatic disease at diagnosis, but the profiled section is **primary tumour tissue**; the eligibility rule refers to the tissue site, not the stage. Recorded in the eligibility table.
- **Treatment status is not reported in GEO** for this cohort. It is recorded as "not reported" rather than assumed to be treatment-naive. Because §5 rule 2 specifies untreated primary CRC for the main analysis, the pre-specified handling is: this cohort is retained in the main pooled spatial analysis, and the Valdeolivas-only and added-cohorts-only analyses (already specified in v1.5) allow the reader to see the result with and without it. Any further claim restricted to treatment-naive tissue must use the Valdeolivas-only result.
- No pathologist annotations are available, so criterion B uses composition quartiles throughout, as specified in v1.5.

**Gene-coverage exclusions (consequence of the pre-specified evaluability rule, §10d of v1.1).** Targeted probe panels do not cover the whole transcriptome. Applying the rule unchanged:
- **MSX1** is absent from GSE226997 P1, P2 and P4 → END_ANGIO not evaluable in those three patients (evaluable in 10 of 13 patients overall).
- **HLA-B** is absent from all GSE280318 and all GSE226997 sections → MYE_IFNG not evaluable in either added cohort (evaluable in 7 of 13 patients overall, all but one from Valdeolivas).
No imputation was performed and no panel gene was substituted.

**Deconvolution quality, reported as a limitation.** RCTD marker-correlation checks pass in all four GSE226997 sections for fibroblast (ρ 0.70–0.92) and pericyte/smooth muscle (0.58–0.74), but **endothelial separation is weak (ρ 0.47–0.51)** compared with GSE280318 (0.74–0.88). Section P3 is estimated at 89% malignant epithelium, leaving little signal for the myeloid (0.12) and T/NK (0.14) components. Endothelial-program results in this cohort, and all P3 results for the non-epithelial compartments, are therefore interpreted with caution.

**Final spatial unit counts.** 13 patients: Valdeolivas 6 (12 sections, pathology-annotated), GSE226997 4 (4 sections, FFPE Visium), GSE280318 3 (3 Visium HD sections; the standard-Visium P2 section remains a within-patient sensitivity analysis only). "≥ 2/3 of patients" = ≥ 9 of 13, ≥ 7 of 10 for END_ANGIO and ≥ 5 of 7 for MYE_IFNG.

**Interpretation rule added for angiogenesis (stated before writing the Results, after the pooled numbers were seen; disclosed as post hoc).** END_ANGIO meets both spatial criteria only in the pooled 13-patient analysis, at exactly the 7-of-10 threshold, with criterion B resting on endothelial-proportion quartiles. Its criterion A effect is robust to leave-one-patient-out in the 13-patient analysis but is the smallest in the panel (0.015, 95% CI 0.004-0.026). Under the v1.3 evidence rule it remains **Unsupported**, because spatial evidence cannot promote a program that failed independent single-cell validation (frozen MAE ratio 2.11, 95% CI 1.37–3.35). This is stated explicitly in the Results rather than left implicit.

## 13 (new) — Pre-specified discovery sensitivity analyses, as executed

The sensitivity analyses specified in §6 (minimum cells), §8 (alternative gene partitions) and §8 note (TMM), plus chemistry and sorting variants, were run as **read-only re-runs of the discovery cross-validation**. The frozen panel was never modified; the script verifies its hash and reports, for each scenario, the MAE ratio against the discovery-mean baseline and how many frozen panel genes are re-selected.

Scenarios: minimum cells 10 / 30 / 50 (primary 20); five alternative random gene partitions (seeds 1–5); TMM normalisation; 10x v2 donors only; all sorting fractions pooled; all tumour-specimen epithelial states.

**Result.** The MAE ratio stays below 1 in 82 of 85 program × scenario fits. Exceptions: END_ANGIO 1.12 and STR_TGFB 1.20 restricted to 10x v2 donors, and STR_TGFB 1.08 under one alternative partition — all in programs already classified Unsupported or Exploratory. STR_TGFB is not evaluable at ≥50 cells (<10 donors).

**Interpretation, carried into the Discussion.** Program-level predictability by a ≤4-gene panel is robust to these choices; the **identity** of the selected genes is not (0–4 of 4 frozen genes re-selected, depending on scenario). The panels are therefore one of several near-equivalent representations of each program, which is a practical advantage for chip assays, where some transcripts are easier to measure than others.

*Known limitation:* alternative partitions re-split the retained gene pool only; the 129 genes removed because they are evaluation genes of another program are not restored, so the alternative splits are not fully independent of the original partition.

## Change log (continued from v1.5)

| Date | Change | Reason |
|---|---|---|
| 2026-10-14 | v1.6 (rev. same day: script 20 read leave-one-patient-out from the 6-patient Valdeolivas file while reporting the 13-patient pooled effect; repointed to the combined file, after which criterion A is LOPO-robust for all seven programs and the manuscript claim was corrected): GSE226997 executed as specified (treatment status not reported; P3 primary tissue with M1 stage included; coverage-driven non-evaluability of MSX1/HLA-B; weak endothelial deconvolution recorded); final 13-patient pooling; angiogenesis interpretation stated; discovery sensitivity analyses recorded as performed | Documents execution of v1.5 rules; the angiogenesis interpretation note is post hoc and disclosed |
