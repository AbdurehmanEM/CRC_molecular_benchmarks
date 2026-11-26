# Protocol version index

Twenty-three dated protocol versions, written during the study. Each version states what was fixed and when, relative to the data. **No version was edited after its date.** Versions v0.1–v1.9 therefore carry the original working title; the title was revised in v2.0 and v2.1, and the v2.1 title was restored in v2.2 after interim wording used during revision; the earlier files were deliberately left alone.

Frozen artefacts referenced throughout:
- `frozen_panel_v2.rds` — SHA-256 `a14ccb4583a901ec3fc1b0c4ba2afcf9e21de5b72150cb3598ce8a6f9da3e0eb` (7 programs, 26 genes, four compartments)
- `frozen_panel_LYM_v1.rds` — SHA-256 `7dc3256d97477a796fcacd1df71021740888f7978a01ea36caf70787da7f513f` (1 program, 4 genes, lymphoid compartment)
- `frozen_comparators_v2.3.rds` — SHA-256 `f9e6bb28187dfb275652f6c3e78cfb7441ee9a75592dca3b08a2e9d765b23b7f` (comparators C1/C2, 8 programs; v2.3)

| Version | Date | What it fixed | Timing relative to data |
|---|---|---|---|
| v0.1 | 8 Oct 2026 | Research question, compartments, scope boundaries, eligibility rules, candidate datasets | Before any data were opened |
| v0.2 | 8 Oct 2026 | GSE178341 → discovery, GSE132465 + GSE144735 → validation; firewall rule | Before any GSE178341 data or GSE132465 expression analysis |
| v0.3 | 8 Oct 2026 | Added GSE200997, Valdeolivas Visium, GSE226997; cohort-overlap rule | Before any data from these cohorts were opened |
| v0.4 | 8 Oct 2026 | Added GSE166555; spatial cohort order; exclusions listed; **cohort list closed** | Before any of these data were opened |
| v0.5 | 9 Oct 2026 | Pseudobulk unit, sorting rule, harmonised compartments, ≥20-cell minimum, batch handling | From metadata counts only; no expression analysed |
| v0.6 | 9 Oct 2026 | Eight programs, coverage rule, gene partition, target definition | Before any program score was computed |
| v0.7 | 9 Oct 2026 | STR_ECM dropped (coverage 77.4%); CPM normalisation; modelling, baselines, freezing, validation scaling | After coverage check; before any model was fitted |
| v0.8 | 9 Oct 2026 | Post-selection OLS refit of the elastic-net panel (`frozen_panel_v2`); v1 superseded | Discovery CV only; no validation data scored |
| v0.9 | 9 Oct 2026 | Validation criteria, baseline, evaluability rules; frozen hash recorded | Before scoring any validation cohort |
| v1.0 | 9 Oct 2026 | GSE200997 annotation method; cross-cohort meta-analysis method (timing disclosed) | Before any GSE200997 expression scoring |
| v1.1 | 9 Oct 2026 | Spatial criteria A and B, deconvolution, block bootstrap, spatial evidence labels | Before any spatial data were downloaded |
| v1.2 | 12 Oct 2026 | Criterion B region mapping from pathologist label names; RCTD reference definition | After label-name audit; before any spatial expression opened |
| v1.3 | 13 Oct 2026 | Evidence-category assignment rule; chip assay rule | After single-cell and spatial results — **disclosed as post hoc** |
| v1.4 | 13 Oct 2026 | Valdeolivas origin corrected; patient A798015 excluded as non-neoplastic | Clinical table read after the spatial run; eligibility rule unchanged |
| v1.5 | 13 Oct 2026 | GSE280318 and GSE226997 added to the spatial analysis; Visium HD aggregation; pooled spatial result | After file listing; before any expression from these cohorts opened |
| v1.6 | 14 Oct 2026 | GSE226997 as executed; coverage-driven non-evaluability; 13-patient pooling; discovery sensitivity analyses | Execution record; angiogenesis interpretation note disclosed as post hoc |
| v1.7 | 15 Oct 2026 | Epithelial chip correspondence amended to name one MMRd and one MMRp line | Proposal-level only; no analysis affected |
| v1.8 | 15 Oct 2026 | Lymphoid (CD8) compartment added; two pre-specified programs; separate frozen panel; partial-blinding disclosure | §§14.1–14.8 before any CD8 target computed; §14.9 is the execution record |
| v1.9 | 15 Oct 2026 | Lymphoid validation extended to all cohorts; lymphoid spatial parameters | Timing disclosed: written after the first lymphoid result was known |
| v2.0 | 16 Oct 2026 | Manuscript title revised | Editorial only; no result, rule or frozen artefact changed |
| v2.1 | 20 Oct 2026 | Title revised to "Evidence-Graded Molecular Benchmarks…"; spatial analysis repeated with a randomly drawn RCTD reference (h5-order bias found in code review) | **Post hoc, disclosed as such.** Frozen panels and primary spatial reference unchanged; no evidence category changed; apoptosis spatial label Spatial → Partial |
| v2.2 | 13 Nov 2026 | Title restored to "Evidence-Graded Molecular Benchmarks…" after the interim wording "candidate molecular references"; definition of "molecular benchmark" required in Introduction and Conclusion | Editorial only; no analysis, rule, number or frozen artefact changed |
| v2.3 | 23 Nov 2026 | External comparison of the frozen panels with two comparators (C1 compartment DE markers; C2 top-4 correlated genes), selected on discovery data and frozen (scripts 73a, 73b); grading-rule sensitivity under three alternative rules (script 74); execution records included | **Post hoc analysis, disclosed as such** (prompted by pre-submission review, after validation results were known); comparators frozen before any validation cohort was scored; no panel or evidence category changed |

## Companion records

| File | Contents |
|---|---|
| `work_log.md` | Daily log, 8–9 Oct 2026 (setup, pilot, discovery audit, programs, partition) |
| `work_log_part2.md` | Daily log, 14 Oct 2026 onward (GSE226997, 13-patient spatial synthesis, numerical audit, journal decision, code review, reproducibility package, submission assembly) |
| `roadmap_to_publication.md` | Original ten-phase plan, written 8 Oct 2026 before the work began |
| `code_review_2026-10-19.md` | Final code review of all analysis scripts (19 Oct 2026) and follow-ups; source of the v2.1 sensitivity analysis and the corrected scripts |

## How to read this record

Three kinds of entry appear, and they are not equivalent:

1. **Pre-specification** — a rule fixed before the relevant data were opened (v0.1–v1.2, v1.5, and §§14.1–14.8 of v1.8). These are the entries that make the validation claims meaningful.
2. **Execution record** — what happened when a pre-specified rule was applied (v1.6, §14.9 of v1.8, §14.10.5 of v1.9). These report outcomes, including unfavourable ones, without changing any rule.
3. **Post hoc, disclosed as such** — v1.3's evidence-category assignment rule, the angiogenesis interpretation note in v1.6 and the v2.1 reference sensitivity analysis were written after results were known. They are labelled in place rather than presented as pre-specified.

Where a pre-specified rule produced an unwelcome result, the result stands and the decision is recorded: STR_ECM dropped at 77.4% coverage (v0.7); LYM_CYTOTOX dropped at 78% coverage with an available substitute gene set explicitly refused (v1.8 §14.9); the lymphoid pooled estimate of 0.376 reported in place of the more favourable single-cohort 0.646 (v1.9 §14.10.5); angiogenesis left Unsupported despite meeting both spatial criteria (v1.6).
