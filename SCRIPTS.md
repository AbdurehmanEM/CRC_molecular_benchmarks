# Scripts

Generated from the first comment line of each script.

| Script | Purpose |
|---|---|
| `scripts/00_run_all.R` | ============================================================================= |
| `scripts/00_setup_and_versions.R` | 00_setup_and_versions.R |
| `scripts/01_download_GSE132465.R` | 01_download_GSE132465.R |
| `scripts/02_pilot_load_GSE132465.R` | 02_pilot_load_GSE132465.R |
| `scripts/03_pilot_qc_and_counts.R` | 03_pilot_qc_and_counts.R |
| `scripts/04_download_GSE178341_GSE144735.R` | 04_download_GSE178341_GSE144735.R |
| `scripts/05_download_spatial_and_GSE200997.R` | 05_download_spatial_and_GSE200997.R |
| `scripts/06_audit_GSE178341.R` | 06_audit_GSE178341.R |
| `scripts/07_GSE178341_compartment_counts.R` | 07_GSE178341_compartment_counts.R |
| `scripts/08_GSE178341_pseudobulk.R` | 08_GSE178341_pseudobulk.R |
| `scripts/09_programs_and_gene_partition.R` | 09_programs_and_gene_partition.R |
| `scripts/10_panel_development.R` | 10_panel_development.R |
| `scripts/11_validation_GSE132465.R` | 11_validation_GSE132465.R |
| `scripts/12_validation_GSE144735_pooled.R` | 12_validation_GSE144735_pooled.R |
| `scripts/12b_KUL_from_saved.R` | 12b_KUL_from_saved.R |
| `scripts/13_validation_GSE166555.R` | 13_validation_GSE166555.R  (protocol v0.9; frozen panel v2) - Uhlitz 2021, 12 patients, tumour samples |
| `scripts/14_GSE200997_own_annotation.R` | 14_GSE200997_own_annotation.R  (protocol v1.0, section 10b) - own annotation of GSE200997 |
| `scripts/15_validation_GSE200997.R` | 15_validation_GSE200997.R  (protocol v1.0; frozen panel v2; own annotation per section 10b) |
| `scripts/16_cross_cohort_meta.R` | 16_cross_cohort_meta.R  (protocol v1.0, section 10c) - specified after per-cohort results (disclosed) |
| `scripts/17_spatial_reference_RCTD.R` | 17_spatial_reference_RCTD.R  (protocol v1.2, section 10d) - RCTD reference from GSE178341 (composition only) |
| `scripts/17b_check_reference_donor_balance.R` | 17b_check_reference_donor_balance.R  (diagnostic, read-only; added 19 Oct 2026) |
| `scripts/18_spatial_RCTD_deconvolution.R` | 18_spatial_RCTD_deconvolution.R  (protocol v1.2, 10d) - spot QC + RCTD full mode (spacexr 1.4 Bioconductor API) |
| `scripts/18b_RCTD_label_check.R` | 18b_RCTD_label_check.R - restore spot/cell-type labels on RCTD weights and verify them with marker genes |
| `scripts/19_spatial_support.R` | 19_spatial_support.R  (protocol v1.2, section 10d) - criteria A (composition-adjusted concordance) and B (anatomy) |
| `scripts/20_evidence_and_chip_table.R` | 20_evidence_and_chip_table.R  (protocol v1.3, section 11) - evidence categories, chip assay table, spatial LOPO sensitivity |
| `scripts/21_figures.R` | 21_figures.R - publication Figures 1-6 from saved results (no new analysis) |
| `scripts/22_spatial_extra_cohorts.R` | 22_spatial_extra_cohorts.R  (protocol v1.5) - extra spatial cohorts: GSE280318 (Visium HD, 56-um pseudo-spots) [+ GSE226997 later] |
| `scripts/23_spatial_combined.R` | 23_spatial_combined.R  (protocol v1.5) - pooled spatial support over all eligible patients from all spatial cohorts |
| `scripts/24_discovery_sensitivity.R` | 24_discovery_sensitivity.R  (protocol sections 6, 8, 9 - pre-specified sensitivity analyses; DISCOVERY only, read-only) |
| `scripts/25_audit_manuscript_numbers.R` | 25_audit_manuscript_numbers.R - traces every numerical claim in the manuscript to its source file |
| `scripts/26_compare_outputs.R` | 26_compare_outputs.R - compare two runs of the pipeline file by file |
| `scripts/27_build_submission_files.R` | 27_build_submission_files.R - assemble everything needed for submission |
| `scripts/28_CD8_feasibility.R` | 28_CD8_feasibility.R - phase 1: can a lymphoid compartment be added? |
| `scripts/29_lymphoid_pseudobulk.R` | 29_lymphoid_pseudobulk.R  (protocol v1.8, section 14) - CD8 compartment pseudobulk |
| `scripts/30_lymphoid_coverage.R` | 30_lymphoid_coverage.R  (protocol v1.8, section 14.4) - coverage gate for the lymphoid programs |
| `scripts/31_lymphoid_partition.R` | 31_lymphoid_partition.R  (protocol v1.8, section 14.4) - additive gene partition for the lymphoid compartment |
| `scripts/32_lymphoid_panel.R` | 32_lymphoid_panel.R  (protocol v1.8, section 14.4) - compact panel for LYM_IL2_STAT5, DISCOVERY only. |
| `scripts/33_lymphoid_validation.R` | 33_lymphoid_validation.R  (protocol v1.8, section 14.5) - frozen lymphoid panel in GSE132465. |
| `scripts/34_lymphoid_extra_cohorts.R` | 34_lymphoid_extra_cohorts.R  (protocol v1.9, section 14.10) - frozen lymphoid panel in GSE166555 and GSE144735. |
| `scripts/35_lymphoid_meta.R` | 35_lymphoid_meta.R  (protocol v1.9, section 14.10) - pooled lymphoid validation. |
| `scripts/36_lymphoid_spatial_valdeolivas.R` | 36_lymphoid_spatial_valdeolivas.R  (protocol v1.9, section 14.11) - criteria A and B for LYM_IL2_STAT5. |
| `scripts/37_lymphoid_spatial_extra.R` | 37_lymphoid_spatial_extra.R  (protocol v1.9, section 14.11) - LYM_IL2_STAT5 in GSE280318 + GSE226997. |
| `scripts/38_lymphoid_spatial_combined.R` | 38_lymphoid_spatial_combined.R  (protocol v1.9, section 14.11) - pooled lymphoid spatial support, all cohorts. |
| `scripts/39_evidence_and_chip_table_with_lymphoid.R` | 39_evidence_and_chip_table_with_lymphoid.R  (protocol v1.9) - 8-program evidence + chip tables. |
| `scripts/40_figures_with_lymphoid.R` | 21_figures.R - publication Figures 1-6 from saved results (no new analysis) |
| `scripts/41_build_submission_with_lymphoid.R` | 41_build_submission_with_lymphoid.R - rebuild figures, draw Figures S3-S4, copy final figure names (tables: script 44) |
| `scripts/42_spatial_sensitivity_random_reference.R` | 42_spatial_sensitivity_random_reference.R  (sensitivity analysis, added 20 Oct 2026) |
| `scripts/43_publication_figures.R` | 43_publication_figures.R - publication figures (main 1-4, supplementary S1, S2, S5) |
| `scripts/44_build_supplementary_tables.R` | 44_build_supplementary_tables.R - submission-ready supplementary workbook |
| `scripts/50_revision_helpers.R` | 50_revision_helpers.R - shared functions for revision scripts 51-64. |
| `scripts/51_inspect_inputs.R` | 51_inspect_inputs.R - STEP 1: check that every input exists and looks right. |
| `scripts/52_reproduce_published.R` | 52_reproduce_published.R - STEP 2: REPRODUCTION GATE |
| `scripts/53_null_random_panels.R` | 53_null_random_panels.R - Reviewer point 7 (predictor-target dependence) |
| `scripts/54_leave_one_donor_out.R` | 54_leave_one_donor_out.R - Reviewer point 4 (precision with n = 11) |
| `scripts/55_meta_analysis_sensitivity.R` | 55_meta_analysis_sensitivity.R - Reviewer point 6 (heterogeneity, k = 2-3) |
| `scripts/56_calibration_transport.R` | 56_calibration_transport.R - Reviewer point 8 (external platform transport) |
| `scripts/57_multiplicity_fdr.R` | 57_multiplicity_fdr.R - Reviewer point 15 (no multiplicity correction) |
| `scripts/58_spatial_table_sensitivity.R` | 58_spatial_table_sensitivity.R - Reviewer points 11 and 14 (and code-review B2) |
| `scripts/59_spatial_rerun.R` | 59_spatial_rerun.R - Reviewer points 9 and 10 (pseudoreplication, circularity) |
| `scripts/60_independent_rctd_reference.R` | 60_independent_rctd_reference.R - Reviewer point 10 (reference circularity) |
| `scripts/61_alternative_panels.R` | 61_alternative_panels.R - Reviewer point 13 (panel gene stability) |
| `scripts/62_mmr_site_strata.R` | 62_mmr_site_strata.R - Reviewer point 26 (clinical generalisation) |
| `scripts/63_reproducibility_clean_room.R` | 63_reproducibility_clean_room.R - Reviewer point 27 (reproducibility) |
| `scripts/64_collect_revision_tables.R` | 64_collect_revision_tables.R - STEP 14: assemble outputs for the manuscript |
| `scripts/install_packages_new_pc.R` | install_packages_new_pc.R |
| `scripts/revision/00_run_revision.R` | 00_run_revision.R - runs the revision analyses in order, each in a fresh R |
| `scripts/revision/45_revision_figures.R` | 45_revision_figures.R - revised publication figures (main 1-4, supplementary S1-S5) |
| `scripts/revision/50_revision_config.R` | 50_revision_config.R  -  THE ONLY FILE YOU NEED TO EDIT |
| `scripts/revision/50_revision_helpers.R` | 50_revision_helpers.R - shared functions for revision scripts 51-64. |
| `scripts/revision/50b_make_feature_map.R` | 50b_make_feature_map.R - writes the discovery feature map (row_id -> symbol) |
| `scripts/revision/50c_save_lymphoid_validation.R` | 50c_save_lymphoid_validation.R - saves the CD8 (lymphoid) validation log-CPM |
| `scripts/revision/51_inspect_inputs.R` | 51_inspect_inputs.R - STEP 1: check that every input exists and looks right. |
| `scripts/revision/52_reproduce_published.R` | 52_reproduce_published.R - STEP 2: REPRODUCTION GATE |
| `scripts/revision/53_null_random_panels.R` | 53_null_random_panels.R - Reviewer point 7 (predictor-target dependence) |
| `scripts/revision/54_leave_one_donor_out.R` | 54_leave_one_donor_out.R - Reviewer point 4 (precision with n = 11) |
| `scripts/revision/55_meta_analysis_sensitivity.R` | 55_meta_analysis_sensitivity.R - Reviewer point 6 (heterogeneity, k = 2-3) |
| `scripts/revision/56_calibration_transport.R` | 56_calibration_transport.R - Reviewer point 8 (external platform transport) |
| `scripts/revision/57_multiplicity_fdr.R` | 57_multiplicity_fdr.R - Reviewer point 15 (no multiplicity correction) |
| `scripts/revision/58_spatial_table_sensitivity.R` | 58_spatial_table_sensitivity.R - Reviewer points 11 and 14 (and code-review B2) |
| `scripts/revision/59_spatial_rerun.R` | 59_spatial_rerun.R - Reviewer points 9 and 10 (pseudoreplication, circularity) |
| `scripts/revision/59a_export_spatial_sections.R` | 59a_export_spatial_sections.R - writes one file per spatial section in the |
| `scripts/revision/60_independent_rctd_reference.R` | 60_independent_rctd_reference.R - Reviewer point 10 (reference circularity) |
| `scripts/revision/61_alternative_panels.R` | 61_alternative_panels.R - Reviewer point 13 (panel gene stability) |
| `scripts/revision/62_mmr_site_strata.R` | 62_mmr_site_strata.R - Reviewer point 26 (clinical generalisation) |
| `scripts/revision/63_reproducibility_clean_room.R` | 63_reproducibility_clean_room.R - Reviewer point 27 (reproducibility) |
| `scripts/revision/64_collect_revision_tables.R` | 64_collect_revision_tables.R - STEP 14: assemble outputs for the manuscript |
| `scripts/revision/65_merge_supplementary_workbook.py` |  |
| `scripts/revision/66_pubmed_search.R` | 66_pubmed_search.R - documented PubMed search for the prior-art comparison (reviewer point 1) |
| `scripts/revision/67_scopus_screening_decisions.py` | 67_scopus_screening_decisions.py |
