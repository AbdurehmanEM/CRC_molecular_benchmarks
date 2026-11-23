# Roadmap to publication (v1, 8 Oct 2026)

Status: setup done (R 4.6.1, Seurat 5.5.1, edgeR 4.10.5, GEOquery 2.80.0). Next: download GSE132465 (script 01).

| Phase | Weeks | Scripts | Key output | Gate |
|---|---|---|---|---|
| 1 Feasibility | 1-2 | 01-03 | pilot QC, donor x compartment counts, audit | proceed only if donor-level analysis supported |
| 2 Full discovery QC + annotation | 3-4 | 03 full, 04 | QC plots, verified annotations, retained cells/donor | Methods: preprocessing, annotation |
| 3 Pseudobulk + programs + gene partition | 5-6 | 05-06 | per-compartment count matrices; program table; gene_partition.tsv | protocol v1.0 locked before modelling |
| 4 Panel development (discovery only) | 7-9 | 07 | genes, coefficients, CV MAE, baselines, stability | freeze panel (checksum) |
| 5 Independent validation (GSE178341) | 10-11 | 08 | external validation table + figure | report failures honestly |
| 6 Spatial support (GSE226997 + 2nd cohort) | 12-14 | 09 | maps, composition-adjusted, per-patient | spatial support only, no causal claims |
| 7 Final panel + chip assay table | 15 | 10 | evidence categories; proposed assays | |
| 8 Figures/tables | 16-17 | 11 | Fig 1-6, Supp tables | every number traceable |
| 9 Writing | 18-21 | - | Methods > Results > Intro > Discussion > Conclusion > Abstract > Title | co-author review |
| 10 Submission | 22-24 | - | GitHub+Zenodo, preprint, cover letter | journal chosen after results |
