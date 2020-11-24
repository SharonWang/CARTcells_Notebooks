# CARTcells_Code

## Xiaonan Wang
## 25Jun2020
## Summary: Re-generate plots used in the Autolus Car T-cells paper

## Introduction
Chimeric antigen receptor (CAR) T-cell adoptive therapy is set to transform the treatment of a rapidly expanding range of malignancies. Although the activation process of normal T cells is well characterised, comparatively little is known about the activation of cells via the CAR. Here we have used flow cytometry together with single cell transcriptome profiling to characterise the starting material (peripheral blood mononuclear cells) and CAR therapeutic products of 3 healthy donors in the presence and absence of antigen specific stimulation. Analysis of 57,676 single cell transcriptomes showed CAR products to contain several subpopulations of cells, with cellular composition reproducible from donor to donor, and all major cellular subsets compatible with CAR expression. Only 50% of CAR-expressing cells displayed transcriptional changes upon CAR-specific antigen exposure. The resulting molecular signature for CAR T-cell activation provides a rich resource for future dissection of underlying mechanisms. Targeted data interrogation also revealed that a small proportion of antigen-responding CAR-expressing cells displayed an exhaustion signature, with both known markers and genes not previously associated with T-cell exhaustion. Comprehensive single cell transcriptomic analysis thus represents a powerful way to guide the assessment and optimization of clinical-grade CAR-T-cells, and inform future research into the underlying molecular processes.

## Notebooks
This folder contains all jupyter notebooks to regenerate figures used in the Autolus CAR T-cells paper.
  - <ins>**[Paper_Figures.ipynb](https://github.com/SharonWang/CARTcells_MT_Code/blob/master/Notebooks/Paper_figures.ipynb)**</ins>: 
    - Main notebook to renerate all figures used in the paper
  - <ins>**[TenX_preanalysis.ipynb](https://github.com/SharonWang/CARTcells_MT_Code/blob/master/Notebooks/TenX_preanalysis.ipynb)**</ins>:
    - Preanalysis of all TenX samples, including concatenating UMI count files and doublet detection
  - <ins>**[PBMC_b4BC_analysis.ipynb](https://github.com/SharonWang/CARTcells_MT_Code/blob/master/Notebooks/PBMC_b4BC_analysis.ipynb)**</ins>:
    - Analysis of PBMC samples using SCANPY, patient variance was observed so that batch correction was performed
  - <ins>**[PBMC_BatchCorrection_R.R](https://github.com/SharonWang/CARTcells_MT_Code/blob/master/Notebooks/PBMC_BatchCorrection_R.ipynb)**</ins>:
    - Batch correction of PBMC samples using MNNcorrect in R
  - <ins>**[PBMC_BC_Analysis.ipynb](https://github.com/SharonWang/CARTcells_MT_Code/blob/master/Notebooks/PBMC_BC_Analysis.ipynb)**</ins>:
    - Analysis of batch correction PBMC samples using SCANPY
  - <ins>**[Project_PBMC_onto_PBMC68Kzheng17.ipynb](https://github.com/SharonWang/CARTcells_MT_Code/blob/master/Notebooks/Project_PBMC_onto_PBMC68Kzheng17.ipynb)**</ins>:
    - Projection of PBMC samples onto [Zheng2017](https://www.nature.com/articles/ncomms14049) landscape
  - <ins>**[Tcells_in_PBMC.ipynb](https://github.com/SharonWang/CARTcells_MT_Code/blob/master/Notebooks/Tcells_in_PBMC.ipynb)**</ins>:
    - Analysis of T cells in batch corrected PBMC samples using SCANPY
  - <ins>**[Tcells_analysis.ipynb](https://github.com/SharonWang/CARTcells_MT_Code/blob/master/Notebooks/Tcells_analysis.ipynb)**</ins>:
    - Analysis of *in vitro* T cell samples using SCANPY
  - <ins>**[Tcells_conST_sep.ipynb](https://github.com/SharonWang/CARTcells_MT_Code/blob/master/Notebooks/Tcells_conST_sep.ipynb)**</ins>:
    - Analysis of *in vitro* T cell samples separated by condition (+Antigen/-Antigen) using SCANPY
  - <ins>**[Tcells_vecPB_sep.ipynb](https://github.com/SharonWang/CARTcells_MT_Code/blob/master/Notebooks/Tcells_vecPN_sep.ipynb)**</ins>:
    - Analysis of *in vitro* T cell samples separated by vector expression (Vector positive/Vector negative) using SCANPY
  - <ins>**[Project_Tcells_onto_Bulk_EMTAB2319.ipynb](https://github.com/SharonWang/CARTcells_MT_Code/blob/master/Notebooks/Project_Tcells_onto_BulkEMTAB2319.ipynb)**</ins>:
    - Projection of T cell samples onto [Bonnal2015](https://www.nature.com/articles/sdata201551) landscape
  - <ins>**[Exhaustion_Sig_calculation.ipynb](https://github.com/SharonWang/CARTcells_MT_Code/blob/master/Notebooks/Exhaustion_Sig_calculation.ipynb)**</ins>:
    - Calculation of T cell Exhaustion signature
