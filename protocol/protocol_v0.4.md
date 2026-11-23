# Study protocol — v0.4 (8 Oct 2026)

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

## 4. Datasets and roles (v0.4 — fixed before any data from these cohorts were opened)

| Accession | Role | Status |
|---|---|---|
| GSE178341 (Pelka 2021; 62 patients, US) | **Discovery** (programs, gene partition, panel training, donor-grouped CV) | To download and audit; confirm CD45 enrichment per specimen |
| GSE132465 (Lee 2020; 23 patients, SMC Korea) | **Primary independent validation** | Loaded and QC'd; used only for feasibility counts, no program scores or models computed |
| GSE144735 (Lee 2020; 6 patients, KUL Belgium) | **Secondary validation** (reported separately and pooled with SMC as a sensitivity analysis) | To download and audit |
| GSE200997 (Khaliq 2022; 16 treatment-naive CRC + 7 adjacent normal, US) | **Validation 3** (included only if audit confirms raw counts, donor IDs, compartment annotations, no overlap) | To download and audit |
| GSE166555 (Uhlitz 2021; tumor + normal, ~12 patients, Germany) | **Validation 4** (same inclusion conditions) | To download and audit |
| Valdeolivas et al. 2024, Zenodo 7739700 (Visium; 14 sections, Belgian + Korean patients; pathologist spot annotations) | **Main spatial cohort** | To download and audit; check patient overlap with GSE132465/GSE144735 |
| GSE280318 (10x Genomics CRC study; subseries GSE280311/13/15; 5 patients; Visium, Visium HD, Xenium + matched single-cell from the same patients) | **Second spatial cohort** (matched single-cell used as the deconvolution reference) | To download and audit; files > 4 GB need a non-FAT32 drive |
| GSE226997 (Park et al.; 4 patients, Visium) | **Third spatial cohort** (if storage allows) | 41 GB RAW.tar exceeds FAT32 limit |

**Considered and not included (decided 8 Oct 2026, before any data were opened):** Joanito et al. 2022 (controlled Synapse access; likely includes SMC/KUL patients already used for validation); GSE146771 (CD45-based sorting, partly Smart-seq2; may be used only as a labelled myeloid sensitivity analysis); GSE161277, GSE221575, GSE231559 (few tumours each); metastasis cohorts and post-treatment cohorts (violate eligibility); T-cell-only datasets (missing compartments). No further cohorts will be added after model development begins.

**Reason for change (v0.1 → v0.2):** in GSE132465 only 11/23 donors (stromal) and 9/23 donors (endothelial) reach ≥30 tumour cells per compartment, too few for stable donor-grouped elastic-net training. GSE178341 is larger (62 patients). Roles were swapped before any expression-based analysis of either cohort, so GSE132465 remains untouched by model development.

**Overlap rule (v0.3):** discovery (GSE178341, US) must share no patients with any validation or spatial cohort. Spatial cohorts may share patients with validation cohorts (Valdeolivas Korean/Belgian samples vs SMC/KUL) because spatial analysis tests a different question; any overlap is checked, reported, and analysed as a sensitivity analysis excluding shared patients.

**Firewall rule:** no program scores, gene selection, model fitting or performance metrics may be computed on GSE132465 or GSE144735 until the panel is frozen (Section 9); this applies to all validation and spatial cohorts. Cell counts, QC metrics and annotation mapping already computed for feasibility are permitted.

**GSE132465 notes for validation preprocessing:** data were pre-filtered by the authors (max mitochondrial fraction ≈20% in every sample); tumours SMC16, SMC18, SMC21, SMC22 show 4–6× higher UMIs per cell and epithelium-dominated composition (possible chemistry/preparation batch — to verify in Lee 2020).

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
| 2026-10-08 | v0.4: added GSE166555 (validation 4); GSE280318 promoted to second spatial cohort; GSE226997 → third; exclusions listed; cohort list closed | Systematic re-check of public CRC datasets; changes made before any of these data were opened |
| 2026-10-08 | v0.3: added GSE200997 (tertiary validation), Valdeolivas Visium (main spatial), GSE226997 (second spatial), GSE280311 (optional); overlap rule | Increase independent validation and spatial patient numbers; added before any data from these cohorts were opened |
| 2026-10-08 | v0.2: GSE178341 → discovery; GSE132465 (+GSE144735) → validation; firewall rule added | GSE132465 has 9 (endothelial) and 11 (stromal) eligible donors at ≥30 cells; roles changed before any GSE178341 data or GSE132465 expression analyses |
