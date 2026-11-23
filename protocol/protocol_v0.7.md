# Study protocol — v0.7 (9 Oct 2026)

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

## 6. Pseudobulk unit and discovery-cohort rules (v0.5 — fixed from metadata counts only, before any expression analysis)

**Unit:** donor × compartment, tumour specimens only. Multiple tumour specimens from one donor (GSE178341: C130 TA/TB, C171 TA/TB) are **pooled** into one profile. Adjacent-normal specimens are not used in the primary analysis.

**Sorting rule (GSE178341):** primary analysis uses cells from *unsorted* and *LiveMACS* fractions only, matching the unsorted validation cohorts. CD45+ MACS and mixed fractions are excluded (they are not pure: they contain 23,167 epithelial and 1,264 endothelial tumour cells). Two donors with no unsorted fraction are therefore excluded from the primary analysis. Sensitivity analysis (all compartments): all fractions pooled (eligible at ≥20 cells: epithelial 62, stromal 29, myeloid 62, endothelial 47).

**Compartment definitions (harmonised across cohorts):**

| Compartment | GSE178341 (published clusters) | Validation cohorts (e.g. GSE132465) |
|---|---|---|
| Tumour-associated epithelium | Epithelial cells in tumour specimens, stem/TA-like states cE01, cE02, cE03 | Authors' malignant (CMS-labelled) epithelium in tumour samples |
| Stromal | clMidwayPr = Fibro (fibroblasts, CAFs, myofibroblasts) | Stromal 1–3 + Myofibroblasts |
| Myeloid | Monocytes + macrophages (cM01, cM02) | Myeloid monocyte/macrophage states (SPP1+, pro-inflammatory, proliferating); cDC and unknown excluded |
| Endothelial | Blood-vascular ECs (lymphatic excluded) | Tip-like, Stalk-like, Proliferative ECs |
| Context (kept, not modelled) | Lymphatic EC, pericytes, smooth muscle, Schwann, DCs, granulocytes, lymphocytes, plasma, mast | Same |

Sensitivity: all tumour-specimen epithelium. A copy-number-based malignancy check on a subset is planned as additional support; epithelial labels are reported as "tumour-associated", not proven malignant.

**Minimum cells:** 20 cells per donor × compartment profile (all compartments); sensitivity analyses at 10, 30 and 50. Eligible GSE178341 donors at ≥20 (primary rule): epithelial 60, stromal 27, myeloid 60, endothelial 42 (of 62); all four compartments 23.

**Cell QC:** authors' pre-filtering inherited; additional per-channel flags (3 MAD: low log library size, low log genes, high mitochondrial fraction) as in the validation cohorts. Thresholds and numbers removed reported per channel.

**Gene identifiers:** Ensembl IDs internally (GSE178341 has 1,918 duplicated gene symbols); symbol-level outputs sum duplicated symbols. Cross-cohort matching rule fixed before validation.

**Batch:** 10x chemistry (v2 / v3; one donor, C162, has both) and hospital (MGH / DFCI) recorded per donor. Donor-grouped CV folds are balanced by majority chemistry; a sensitivity analysis restricts to v2 donors.

## 7. Programs (v0.6 — fixed before any program score was computed)

| Program | Compartment | Gene set (MSigDB via msigdbr 26.1.1) | Described as |
|---|---|---|---|
| EPI_HYPOXIA | Tumour-associated epithelium | HALLMARK_HYPOXIA | Hypoxia-associated transcription |
| EPI_APOPTOSIS | Tumour-associated epithelium | HALLMARK_APOPTOSIS | Apoptosis-related transcription |
| STR_TGFB | Stromal | HALLMARK_TGF_BETA_SIGNALING | TGF-β response transcription |
| STR_ECM | Stromal | REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION | ECM-organisation transcription |
| MYE_INFLAM | Myeloid | HALLMARK_INFLAMMATORY_RESPONSE | Inflammatory-response transcription |
| MYE_IFNG | Myeloid | HALLMARK_INTERFERON_GAMMA_RESPONSE | IFN-γ response transcription |
| END_ANGIO | Endothelial | HALLMARK_ANGIOGENESIS | Angiogenesis-associated transcription |
| END_TNFA_NFKB | Endothelial | HALLMARK_TNFA_SIGNALING_VIA_NFKB | Inflammatory (TNF/NF-κB) activation transcription |

Genes matched by Ensembl ID. **Expressed** = CPM ≥ 1 in ≥ 25% of eligible discovery donors (≥ 20 cells) of that compartment. **Retention rule:** coverage ≥ 80% of the set's genes expressed, ≥ 10 evaluation genes and ≥ 5 predictor genes. Gene-set version recorded in `04_processed_data/programs/geneset_version.csv`; cached copy in `02_metadata/genesets/`.

## 8. Gene partition and leakage prevention (v0.6)

- Each program's expressed genes are split at random (seed 20261008) into two halves: **evaluation genes** (target) and **predictor genes** (panel-eligible).
- A predictor gene may not be an evaluation gene of **any** program; such genes are removed from the predictor pool.
- The panel for a program is selected only from that program's predictor genes.
- Saved before modelling: `04_processed_data/programs/gene_partition.tsv` + hash (`gene_partition.sha`).
- Sensitivity: five alternative random splits (seeds 1–5) to show results do not depend on one split.
- Limitation stated in the paper: predictor and evaluation genes come from the same assay and pathway, so biological correlation remains.

**Normalisation (v0.7 change, before modelling):** sample-wise log2-CPM (edgeR::cpm, prior count 2, library size only) instead of TMM, because TMM factors depend on the other samples in a cohort and cannot be frozen for application to new samples. TMM is a sensitivity analysis.

**Target score (per donor × program):** log2-CPM of evaluation genes; each gene z-scored using **discovery** mean and SD (stored and frozen); target = mean z-score. Inside cross-validation, z-scoring parameters are estimated in training folds only.

## 9. Modelling and freezing (v0.7)

**Programs modelled:** 7 of 8. STR_ECM was not retained (coverage 77.4% < 80% pre-specified threshold; decided by script 09 before any score was computed). It may be shown only as a labelled exploratory supplementary analysis.

- **Panel model:** elastic net (glmnet, α = 0.5, standardize = FALSE on z-scored predictors), restricted to the program's predictor genes; λ chosen by inner 5-fold CV among λ values giving 1–4 non-zero genes.
- **Outer validation:** 5-fold CV × 10 repeats; one profile per donor, so folds are donor-grouped; folds balanced by majority 10x chemistry; target and predictor z-score parameters, gene selection and λ estimated in training folds only.
- **Baselines:** (1) training-mean prediction; (2) equal-weight score of the 4 predictor genes most correlated with the target in training, linearly calibrated in training; (3) conventional markers (pre-specified list, restricted to the predictor pool), calibrated in training; (4) larger model = elastic net at λ.min without the gene cap.
- **Metrics:** primary MAE; secondary RMSE, Spearman ρ, calibration slope and intercept; MAE ratio vs training-mean baseline. Variation reported across the 10 CV repeats.
- **Selection stability:** frequency each predictor gene is selected across the 50 training fits.
- **Freeze:** final model refit on all eligible discovery donors; saved as `04_processed_data/frozen_panel/frozen_panel_v1.rds` with SHA-256, containing genes (Ensembl + symbol), coefficients, intercept, predictor and evaluation-gene means/SDs, normalisation, target definition, gene-ID rule, missing-gene rule (program not evaluable, no imputation) and metrics. **No change after validation data are scored**; any revision requires a new version and an untouched cohort.
- **Validation scaling:** primary analysis applies frozen discovery means/SDs; a pre-specified secondary analysis standardises each validation cohort internally (tests ranking ability independent of cohort shift).

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
| 2026-10-09 | v0.7: STR_ECM not retained (coverage 77%); sample-wise CPM normalisation; modelling, baselines, freezing and validation-scaling rules specified | After coverage check; before any model was fitted |
| 2026-10-09 | v0.6: programs, coverage rule, gene partition and target definition fixed | Before any program score or model was computed |
| 2026-10-09 | v0.5: discovery rules fixed (unsorted fractions; pooled TA/TB; harmonised compartments; epithelial = stem/TA-like in tumours; min 20 cells; batch handling) | Based on GSE178341 metadata counts only; no expression data analysed |
| 2026-10-08 | v0.4: added GSE166555 (validation 4); GSE280318 promoted to second spatial cohort; GSE226997 → third; exclusions listed; cohort list closed | Systematic re-check of public CRC datasets; changes made before any of these data were opened |
| 2026-10-08 | v0.3: added GSE200997 (tertiary validation), Valdeolivas Visium (main spatial), GSE226997 (second spatial), GSE280311 (optional); overlap rule | Increase independent validation and spatial patient numbers; added before any data from these cohorts were opened |
| 2026-10-08 | v0.2: GSE178341 → discovery; GSE132465 (+GSE144735) → validation; firewall rule added | GSE132465 has 9 (endothelial) and 11 (stromal) eligible donors at ≥30 cells; roles changed before any GSE178341 data or GSE132465 expression analyses |
