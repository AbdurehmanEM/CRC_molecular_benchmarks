# Study protocol — v2.0 (16 Oct 2026)

**Title (revised):** Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics

**Previous working title (v0.1–v1.9):** Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics

> **Amendment to v1.9.** Editorial only. The title of the manuscript was revised so that it describes the study as actually executed. **No analysis, dataset, program definition, gene partition, threshold, evidence rule or frozen panel is affected.** `frozen_panel_v2.rds` (`a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb`) and `frozen_panel_LYM_v1.rds` (`7dc3256d97477a796fcacd1df71021740888f7978a01ea36caf70787da7f513f`) are unmodified. Every reported number is unchanged.

## Reason for the revision

The original working title was fixed at project registration, before any data were analysed. Three features of the completed study are not described by it:

1. The study does not deliver a prioritised biomarker *panel* as its primary claim. It delivers an **evidence-graded assessment** of eight transcriptional programs, of which only one reached the highest evidence category, four are Exploratory, and three are Unsupported. A title promising "biomarker prioritization" overstates the result; "Evidence-Graded" states what the study actually produces.
2. The study acquired a **fifth (lymphoid) compartment** in v1.8–v1.9, so "tumour–stromal–myeloid–vascular" no longer enumerates the compartments analysed. The revised title names no compartment list and so does not go stale.
3. "Interpretable" described the model class (sparse linear panels) but not the study's main methodological asset, which is **pre-specification**: programs, eligibility, thresholds, evidence rules and the frozen panels were fixed in dated protocol versions before the corresponding data were opened. "Pre-Specified" names that asset. Interpretability is retained as a property of the panels and is stated in the Abstract and Methods rather than in the title.

## Scope of the change

Editorial, applied to the title string wherever it appears. The dated protocol files v0.1–v1.9 **retain the original working title** and are not retroactively edited: they are the pre-specification record, and altering dated documents after the fact would destroy the property that makes that record worth having. A reviewer who notices the discrepancy between the registration title and the submitted title will find this amendment, which is the intended behaviour.

| Artefact | Change |
|---|---|
| `09_manuscript/manuscript_BioinfAdv_v1.docx` | Title, running head |
| `09_manuscript/supplementary_BioinfAdv_v1.docx` | Title |
| `README.md` | Title line |
| Protocol v0.1–v1.9 | **Unchanged by design** |

## Change log (continued from v1.9)

| Date | Change | Reason |
|---|---|---|
| 2026-10-16 | v2.0: manuscript title revised to "Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics", replacing the registration-stage working title | The original title predates the analysis: it promises biomarker prioritisation, which the evidence grades do not support; it enumerates four compartments where five were analysed; and it names interpretability rather than pre-specification as the framework's defining property. Editorial only; no result, rule or frozen artefact changed. Earlier protocol versions keep the original title, as the pre-specification record must not be rewritten |
