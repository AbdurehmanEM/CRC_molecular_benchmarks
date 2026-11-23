# Study protocol — v2.1 (20 Oct 2026)

**Title (revised):** Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics

**Previous titles:** v2.0 "Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics"; v0.1–v1.9 "Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics"

> **Amendment to v2.0.** It makes two changes: (A) an editorial title change and (B) a **post hoc sensitivity analysis** of the spatial deconvolution reference, added after all results were known. The primary analysis is unchanged: `frozen_panel_v2.rds` (`a14ccb45…e0eb`) and `frozen_panel_LYM_v1.rds` (`7dc3256d…513f`) are unmodified, the published spatial reference remains primary, and no evidence category changes.

## A. Title

"Molecular benchmarks" describes the output for its intended use: patient-derived, evidence-graded reference measurements against which a chip compartment can be compared. "Derived from single-cell and spatial transcriptomics" states the data source without listing compartments. To prevent overclaiming, the manuscript defines "benchmark" in the Introduction and Conclusion as a reference to compare against, **not** evidence that any existing chip meets it. The running head becomes "Evidence-graded benchmarks for CRC-on-a-chip".

## B. Sensitivity analysis of the RCTD reference (post hoc)

**Trigger.** Code review of 19 Oct 2026: script 17 samples 1,200 cells per type at random, re-sorts them by h5 position and keeps the first 800 that pass QC. The diagnostic `17b_check_reference_donor_balance.R` confirmed that donors stored late in the h5 file were under-represented: for every type, a donor's selection rate fell with its h5 position (ρ −0.47 to −0.73), and 35–55 of 58–62 eligible donors were represented (normal epithelium 34 of 36).

**Analysis.** `42_spatial_sensitivity_random_reference.R` rebuilt the reference with random selection (52–62 donors per type) and re-ran scripts 18, 18b, 19, 22, 23 and 36–38 into separate outputs (`*_randomsel`, `06_results/sens_randomref/`). Script 18b verified the labels. No published output was overwritten.

**Result (all spatial cohorts, main).** Pooled criterion A changed by at most 0.006, with near-identical CIs. Published versus corrected: END_ANGIO 0.015 vs 0.016; END_TNFA_NFKB 0.203 vs 0.203; EPI_APOPTOSIS 0.156 vs 0.157; EPI_HYPOXIA 0.055 vs 0.055; LYM_IL2_STAT5 0.146 vs 0.148; MYE_IFNG 0.241 vs 0.243; MYE_INFLAM 0.147 vs 0.141; STR_TGFB 0.110 vs 0.111. **One spatial label changed:** EPI_APOPTOSIS moved from Spatial support to Partial, because criterion B fell from 9/13 to 8/13 (threshold 9) in an added-cohort patient scored by composition quartiles; criterion A was unchanged (13/13). **No evidence category changed:** EPI_APOPTOSIS is Exploratory on its single-cell result.

**Disposition.** The published reference remains primary, as pre-specified. The sensitivity analysis is reported in Methods, Results and Limitations, including the apoptosis label change.

## Other disclosures added to the manuscript (no analysis change)

- No doublet filtering was applied.
- GSE166555 myeloid = "Monocyte 1/2" only; the cohort has no macrophage label.
- Genes were matched across cohorts by symbol.
- Fixed QC thresholds were used for the RCTD reference and for spots.
- RCTD marker check: 7 section × type flags, involving rare types and weak myeloid deconvolution in GSE226997 P4 (ρ 0.28).

## Change log (continued from v2.0)

| Date | Change | Reason |
|---|---|---|
| 2026-10-20 | v2.1: title revised to "Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics"; running head "Evidence-graded benchmarks for CRC-on-a-chip" | Names the output by its intended use; "benchmark" is defined in the text as a reference to compare against, not a demonstration of chip fidelity. Editorial only |
| 2026-10-20 | v2.1: post hoc spatial sensitivity analysis with a randomly selected RCTD reference | Code review found h5-order bias in reference cell selection (17b diagnostic). Result: effects within ±0.006, no category change; EPI_APOPTOSIS spatial label Spatial → Partial (criterion B 9/13 → 8/13). Published reference kept as primary |
