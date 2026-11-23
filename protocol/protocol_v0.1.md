# Study protocol — v0.1 (draft, 8 Oct 2026)

**Working title:** Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics

> Version-control this file. Every change after the panel is frozen (Section 9) must be dated and justified in the change log at the bottom.

## 1. Research question

**Which reproducible tumor–stromal–myeloid–vascular programs in human colorectal cancer can be captured by a compact biomarker panel to guide cancer-on-chip assessment?**

## 2. Compartments and intended chip correspondence

| Human tissue compartment | Intended relevance to the chip |
|---|---|
| Malignant epithelial cells | Tumor compartment (HCT116, SW48) |
| Fibroblasts and related stromal states | Stromal biology represented partly by hMSC |
| Myeloid cells | Immune biology represented partly by THP-1 |
| Endothelial cells | Vascular compartment represented by HUVEC |

These are biological correspondences only. hMSC, THP-1 and HUVEC are not assumed to reproduce the full range of patient-derived fibroblast, myeloid or tumor-endothelial states.

## 3. Scope boundaries

- **In scope:** whether transcriptional programs in human CRC are reproducible across patients and cohorts, whether a compact gene panel predicts them, whether they have spatial support, and which future chip measurements are proposed.
- **Out of scope:** chip fidelity, drug-response prediction, protein-level validation, causal cell–cell communication. No wet-lab data are generated.
- Program scores are *transcriptional signatures* (e.g. "hypoxia-associated transcription"), not measurements of oxygen, apoptosis or vessel formation.

## 4. Candidate datasets (roles provisional until audit is complete)

| Accession | Provisional role | Status |
|---|---|---|
| GSE132465 | Pilot and discovery (single-cell) | Audit in progress |
| GSE178341 | Independent single-cell validation | Audit in progress; enrichment status to confirm |
| GSE226997 | Spatial support (Visium) | Only 4 patients — consider adding a second spatial cohort |

## 5. Eligibility rules (pre-specified)

1. Organism *Homo sapiens* for every sample.
2. Main analysis: untreated primary CRC tumor specimens.
3. Adjacent tissue is labelled "adjacent mucosa", never "healthy". It is excluded from the main panel analysis and may be used only in clearly labelled secondary analyses.
4. Raw UMI counts must be available for count-based analyses. Log-normalized/TPM matrices are not used as counts.
5. Donor identity must be recoverable for every cell.
6. Sorting/enrichment (e.g. CD45+) is recorded per sample. Cell-proportion comparisons are only made between samples with comparable enrichment.

## 6. Pseudobulk unit

Donor × specimen × compartment. Provisional minimum: **30 cells** per profile (sensitivity: 10, 20, 50). One primary tumor specimen per donor in the main analysis; if a donor has multiple tumor specimens, the pre-specified rule is: [choose before analysis — e.g. pool all tumor specimens per donor / take the specimen with most cells in the compartment].

## 7. Programs (candidates — at most two per compartment)

| Compartment | Candidate programs | Proposed source |
|---|---|---|
| Malignant epithelial | Hypoxia; apoptosis-related transcription | MSigDB Hallmark (record version) |
| Stromal | TGF-β response; ECM organization | MSigDB Hallmark / Reactome (record version) |
| Myeloid | Inflammatory response; interferon response | MSigDB Hallmark |
| Endothelial | Angiogenesis; inflammatory endothelial activation | MSigDB Hallmark / curated (justify) |

Retention rule: ≥ [e.g. 80%] of program genes detected in both discovery and validation data; biological justification recorded.

## 8. Leakage prevention

Two non-overlapping gene pools per program: **evaluation genes** (define the target score) and **candidate predictor genes** (eligible for the panel). No predictor gene may appear in any primary evaluation set. Saved to `04_processed_data/gene_partition.tsv` before modelling.

## 9. Modelling and freezing

- Elastic net (glmnet), donor-grouped cross-validation; all filtering, scaling and tuning inside training folds.
- Provisional panel size 2–4 genes per program.
- Baselines: discovery mean; equal-weight marker score; conventional markers; larger (unconstrained) model.
- Primary metric: MAE (prediction error). Supporting: Spearman ρ, calibration slope/intercept.
- **Freeze** (genes, ID rules, normalization, scaling values, coefficients, intercepts, target definitions, metrics, missing-gene rules) to `06_results/frozen_panel/` with a checksum *before* validation data are scored.

## 10. Pre-specified success criteria (fill in before modelling)

- Independent validation: panel MAE lower than the discovery-mean baseline, with [bootstrap 95% CI excluding no improvement / ≥ X% reduction].
- Spatial support: association in ≥ [3 of 4] patients after composition adjustment, using spatially aware inference.

## 11. Evidence categories for final panel

Independently supported with spatial evidence · Independently supported in single-cell data · Exploratory · Unsupported.

## 12. Go/no-go after the pilot

Proceed with GSE132465 as discovery only if: cell IDs match annotations; donor IDs recoverable; raw counts present; ≥ [N] donors reach 30 cells in each of the four compartments; scripts reproduce identical outputs on re-run.

## Change log

| Date | Change | Reason |
|---|---|---|
| 2026-10-08 | v0.1 created | — |
