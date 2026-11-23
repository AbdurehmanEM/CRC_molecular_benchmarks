# Study protocol — v2.2 (13 Nov 2026)

**Title:** Evidence-Graded Molecular Benchmarks for Colorectal Cancer-on-a-Chip Models Derived from Single-Cell and Spatial Transcriptomics

**Running head (44 characters):** Evidence-graded benchmarks for CRC-on-a-chip

**Previous titles:** interim revision wording (manuscript V9–V13, Nov 2026) "candidate molecular references"; v2.1 = same title as above; v2.0 "Evidence-Graded Biomarker Panels for Colorectal Cancer-on-a-Chip: A Pre-Specified Framework Integrating Single-Cell and Spatial Transcriptomics"; v0.1–v1.9 "Patient-Informed Biomarker Prioritization for Colorectal Cancer-on-a-Chip: An Interpretable Framework Integrating Single-Cell and Spatial Transcriptomics"

> **Amendment to v2.1.** Editorial only. The author restored the v2.1 title and returned to "molecular benchmarks" throughout. This replaces the interim wording "candidate (molecular) references" used during the pre-submission revision. No analysis, dataset, threshold, evidence rule or frozen panel changed: `frozen_panel_v2.rds` (`a14ccb45…e0eb`) and `frozen_panel_LYM_v1.rds` (`7dc3256d…513f`) are unmodified, and every number is unchanged. Protocol files v0.1–v2.1 are not edited.

## Safeguard against overclaiming (required in the text)
The reviewer's concern was that "benchmark" could imply wet-lab-validated chip performance. It is addressed by definition, not by avoiding the word.
- **Introduction, first use:** a molecular benchmark is defined as a patient-derived, evidence-graded reference measurement against which a chip compartment can be compared; it does not show that any existing chip meets it, and transferring it to qPCR or protein assays requires experimental validation.
- **Conclusion:** the limit is repeated in one sentence.

## Implementation (script `70_apply_benchmark_title.R`, an editing script for the submission documents; not included in this release)
| Artefact | Change |
|---|---|
| Manuscript, Supplementary Methods, cover letter | Title, running head and "candidate (molecular) reference(s)" → "molecular benchmark(s)" |
| Supplementary workbook | Cell text; formatting unchanged |
| Release `CRC_candidate_references_v1.1` → `CRC_molecular_benchmarks_v1.2` | README, CITATION.cff (version 1.2), .zenodo.json, txt/csv. v1.1 not released |
| Figure/table scripts | Scanned only; no figure label used the old wording |
| Protocol v0.1–v2.1 | Unchanged by design |

## Change log (continued from v2.1)
| Date | Change | Reason |
|---|---|---|
| 2026-11-13 | v2.2: title restored to the v2.1 "Evidence-Graded Molecular Benchmarks…" wording; "candidate molecular references" replaced by "molecular benchmarks" in all submission documents; definition required in Introduction and Conclusion | Author decision: "benchmark" names the intended use; the overclaiming concern is handled by explicit definition. Editorial only |
