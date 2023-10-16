# mam-reproducibility
This repository contains the code for reproducing the single-cell analysis in the paper "SYSTEMS LEVEL IDENTIFICATION OF A MATRISOME-ASSOCIATED MACROPHAGE POLARIZATION STATE IN MULTI-ORGAN FIBROSIS". 

Reference: [Systems level identification of a matrisome-associated macrophage polarisation state in multi-organ fibrosis](https://elifesciences.org/articles/85530)

The code is arranged according to the order of the figures in the paper:
- fig0_prepData1a.R: Script for tidying Liver / Lung / Heart datasets
- fig0_prepData1b.R: Script for tidying Skin / Endometrium / Kidney datasets
- fig0_prepData2.R: Script for extracting monocytes / macrophages from each dataset

- fig1_inteDatasets.R: Script for integrating single cells across datasets for each of the six tissues (Liver / Lung / Heart / Skin / Endometrium / Kidney)
- fig2_inteMAM.R: Script for integrating SPP1 macrophages across all six tissues
- fig2_inteMAMother.R: Scripts for incorporating other macrophage populations i.e. TAM / LAM / DAM / SAMs

- fig3_slingshot_analysis.R: Script for slingshot analysis for each of the six tissues
- fig3_monocle_analysis.R: Script for monocle3 analysis for each of the six tissues
- fig3_traj_cell_propensities.R: Script for calculating cell propensities scores

- fig4_exec_scenic.Rmd: Script to run SCENIC pipeline
- fig4_regulon_analysis.R: Script to extract regulon activities
- fig4_tabulate_all_object.R: Script to tabulate regulon activities for each of the six tissues

- fig5_hlca_signatures.R: Script to calculate signature scores for human lung cell atlas data
- fig5_tabulamurissenis_signatures.R: Script to calculate signature scores for tabula murris senis data

The Seurat objects for each of the six tissues and SPP1 macrophages can be downloaded at: https://zenodo.org/record/8266711
