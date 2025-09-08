# A functional map of CDK-drug interactions at single amino acid resolution

## Author: Sam Gould

## Repository Structure:

```
├── cluster_scripts #-------> contains .py and .sh scripts for NGS analysis on computing cluster
├── figures #-------> contains figure files (png and pdf images)
├── notebooks #jupyter notebooks for analysis of screening data
│   ├── 01_library_generation
│   ├── 02_library_qc_and_counts
│   ├── 03_editing_analysis
│   ├── 04_LFC_FDR_tables
│   ├── 05_RNA_seq
│   ├── 06_splicing
│   ├── 07_KLIFS
│   ├── 08_off-targets
│   ├── 09_analysis_and_viz
│   └── 10_PyMOL
├── screening_data #-------> data generated from base editing screens
│   ├── 01_count_tables
│   ├── 02_RPM_tables
│   ├── 03_LFC_FDR_tables
│   ├── 04_editing
│   ├── ABE_subpool_1_screen_data
│   ├── CBE_subpool_1_screen_data
│   ├── CDK12_13_screen_data
│   ├── CDK2_4_6_screen_data
│   ├── KB_compound_mut_screen_data
│   ├── SY_5609_screen_data
│   └── sample_key.xlsx
└── source_data #-------> other data generated for project & external data sources used for analysis
    ├── 01_library_generation_files
    ├── 02_library
    ├── 03_NGS_library_validation
    ├── 04_MTT_data
    ├── 05_flow_data
    ├── 06_snapgene_oligo_samples
    ├── 07_primer_design_endogenous_seq
    ├── 08_screening_data_cell_counts
    ├── 09_sequencing_forms
    ├── 10_validation_guides
    ├── 11_rnaseq_data_deseq2
    ├── 12_alphamissense_scores
    ├── 13_cavity_search_pockets
    ├── 14_consurf
    ├── 15_splicing
    ├── 16_KLIFS
    ├── 17_off-targets
    ├── 18_boltz-2
    ├── 19_gnomAD
    ├── 20_metapredict_disorder
    └── 21_dbNSFP
```
