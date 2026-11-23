# Study protocol — v1.9 (15 Oct 2026)

**Working title:** Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics

> **Amendment to v1.8.** Extends the lymphoid validation to the full set of single-cell validation cohorts (§14.10), and fixes the spatial assessment parameters for the lymphoid compartment (§14.11). Nothing else changes. `frozen_panel_v2.rds` (`a14cc…e0eb`) and `frozen_panel_LYM_v1.rds` (`7dc32…f513f`) are both unmodified, and **no gene, coefficient, intercept or threshold is altered by this amendment.**

## §14.10 — Extension of lymphoid validation to additional cohorts

### Timing disclosure — read this first

**This amendment was written after the GSE132465 lymphoid result was known** (MAE ratio 0.646, 95% CI 0.47–0.868, Supported; v1.8 §14.9). The decision to test further cohorts was therefore taken in knowledge of a favourable first result, and that fact is disclosed here and must be stated in the Methods.

What this does and does not compromise:

- **It does not compromise the panel.** `frozen_panel_LYM_v1.rds` was frozen with a recorded hash before GSE132465 was touched and is not modified. Each additional cohort is a genuine out-of-sample application of an unchanged model.
- **It does affect the choice of how many tests to run.** Having seen one success, the investigators chose to run more. The guard against this becoming selective reporting is §14.10.3.

### 14.10.1 Cohorts

The lymphoid panel is applied to every single-cell validation cohort already used in the four-compartment analysis, in addition to GSE132465: **GSE144735 (KUL), GSE166555, GSE200997**. A cohort is included if it carries a CD8 T-cell annotation separable from CD4, Treg, γδ and NK cells, and excluded otherwise, with the reason recorded.

### 14.10.2 Procedure, inherited unchanged

For each cohort: CD8 cells defined from that cohort's own annotation to match §14.2 as closely as its labels allow, with the mapping recorded before the model is run; donor-level pseudobulk at ≥20 cells; log₂-CPM (prior count 2); frozen scaling as primary and within-cohort scaling as secondary; MAE ratio against the zero baseline with 2,000-resample bootstrap CI; `MIN_DONORS` 10 and `MIN_EVAL_COVER` 0.80; missing panel genes render a cohort not evaluable, with no imputation. Pooling uses the same DerSimonian–Laird procedure as script 16.

### 14.10.3 Reporting commitment

**Every cohort attempted is reported**, including cohorts that fail, cohorts that are not evaluable, and cohorts whose MAE ratio exceeds 1. The pooled estimate is computed over all evaluable cohorts; no cohort may be dropped from the pool on the basis of its result. If the pooled estimate is worse than the GSE132465 result, the pooled estimate is the headline figure and the single-cohort result is not substituted for it.

### 14.10.4 Context for interpretation

Single-cohort validation is not in itself a deficiency in this study: both endothelial programs, including `END_TNFA_NFKB` — the only program classified as independently supported with spatial evidence — rest on GSE132465 alone (`k_cohorts` = 1). The extension is undertaken because the lymphoid compartment is a late addition and will attract the most scrutiny, not because one cohort would have been unacceptable.

### 14.10.5 Execution record (added after running)

CD8 definitions, fixed before each model run:

| Cohort | Definition | Outcome |
|---|---|---|
| GSE132465 | `Cell_type == "T cells" & Cell_subtype == "CD8+ T cells"` | 20 donors, **MAE ratio 0.646 (0.470–0.868), Supported** |
| GSE166555 | `minor_cell_type` matching `^CD8\+ T` (Tconv, Treg and others excluded) | 12 donors, **MAE ratio 0.214 (0.122–0.298), Supported** |
| GSE144735 (KUL) | `Cell_type == "T cells" & Cell_subtype == "CD8+ T cells"` | 369 cells over 6 patients; only **2 donors** reach ≥20 cells → **Not evaluable** |
| GSE200997 | — | **Excluded**: the cohort's own annotation resolves T and NK jointly (`ctx_T_NK`, 15,810 tumour cells); CD8 cannot be separated without constructing a new sub-annotation, which would be a methodological addition applied to this compartment alone |

**Pooled result (DerSimonian–Laird, frozen scaling, cohorts with ≥10 donors):**

| | k | Donors | Cohorts supported | Pooled MAE ratio | 95% CI | τ² | I² | Status |
|---|---|---|---|---|---|---|---|---|
| Frozen (primary) | 2 | 32 | 2 of 2 | **0.376** | 0.128–1.111 | 0.572 | **93.7%** | **Partially supported** |
| Within-cohort | 2 | 32 | 0 of 2 | 0.945 | 0.620–1.441 | 0.000 | 0.0% | Partially supported |

**Consequence, accepted under §14.10.3.** The pooled interval crosses 1, so `meta_status` is not "Supported", and under the v1.3 evidence rule `LYM_IL2_STAT5` is classified **Exploratory** — a weaker classification than the GSE132465 result alone would have produced. The single-cohort result is **not** substituted for the pooled figure.

**Interpretation.** Both evaluable cohorts are individually Supported and both point the same direction; the heterogeneity (I² 93.7%) concerns the magnitude of the benefit, not its sign. With k = 2, τ² is poorly estimated and the random-effects interval is correspondingly unstable. This is the same disposition already applied to `STR_TGFB` (k = 2, pooled 0.748, Exploratory).

**Structural asymmetry, for the limitations.** Programs with a single evaluable cohort — both endothelial programs, including `END_TNFA_NFKB` — are judged on that cohort and never incur a heterogeneity penalty, whereas programs with two or more are. This property of the evidence rule predates the lymphoid extension and applies to all eight programs; it is disclosed rather than corrected, since changing the rule after seeing results would be worse.

Offset diagnostic (frozen scale): mean predicted − mean observed = −0.168 (GSE132465) and −0.124 (GSE166555), consistent with the negative offsets seen across the other cohort × program combinations.

## §14.11 (new) — Spatial assessment parameters for the lymphoid compartment

Fixed before the spatial analysis was run. Criteria A and B otherwise follow v1.5 §10d unchanged, including the 8 × 8 array-block bootstrap (1,000 resamples), DerSimonian–Laird pooling across patients, and the ≥2/3-of-patients requirement.

**14.11.1 Composition covariate — declared proxy.** The RCTD reference built in script 17 has nine cell types, of which the lymphoid one is `T_NK`: CD8, CD4, Treg, γδ, NK and ILC pooled. No CD8-specific reference type exists. The compartment proportion used for the lymphoid criterion A and B is therefore the **T/NK proportion, which is broader than the CD8 population from which the panel was derived**. Rebuilding the reference with a CD8-specific type would alter the deconvolution weights for every existing program and invalidate the four-compartment spatial results, so it is not done. This proxy is a limitation of the spatial assessment for this compartment and must be stated in the Results and Discussion.

**14.11.2 Criterion B rule.** Lymphoid criterion B uses **T/NK proportion quartiles** (upper vs lower quartile within each section), mirroring the rule already applied to the endothelial compartment. The myeloid rule is deliberately **not** reused: it is based on "IC high/low" pathology labels and explicitly excludes spots annotated as lymphoid aggregates, which would be the most informative spots for a lymphoid program. Defining a new aggregate-based rule would require choosing a label vocabulary after inspecting the data, so it is not done. Quartiles are in any case the rule applied in the two cohorts that carry no pathology annotations.

**14.11.3 Classification is not affected by the outcome.** Under the v1.3 evidence rule, spatial evidence cannot promote a program that failed independent single-cell validation. `LYM_IL2_STAT5` is Exploratory on the pooled single-cell result (§14.10.5), so the spatial analysis **cannot change its category** and is reported descriptively. It is run nonetheless, for consistency with `END_ANGIO`, which has spatial support and remains Unsupported for the same reason.

## Change log (continued from v1.8)

| Date | Change | Reason |
|---|---|---|
| 2026-10-15 | v1.9 §14.10: lymphoid validation extended to GSE144735, GSE166555 and GSE200997, with an explicit timing disclosure and a commitment to report every cohort attempted | Strengthens the evidence for the one compartment added late in the study. Decided after the favourable GSE132465 result, which is disclosed; the frozen panel is unchanged, so each cohort remains a true out-of-sample test |
| 2026-10-15 | §14.10.5 execution record: GSE166555 Supported (0.214); KUL not evaluable (2 donors); GSE200997 excluded (T/NK not separable); **pooled 0.376 (0.128–1.111), I² 93.7% → Partially supported → program classified Exploratory** | The extension produced a weaker classification than the first cohort alone. Recorded and accepted under §14.10.3; the single-cohort figure is not substituted |
| 2026-10-15 | §14.11: lymphoid spatial parameters fixed — T/NK proportion as a declared proxy for CD8 composition, criterion B by T/NK quartiles, and acknowledgement that the result cannot alter the evidence category | Pre-specifies the two choices the spatial run forces, before it is run |
