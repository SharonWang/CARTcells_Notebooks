# CARTcells_Notebooks

## Xiaonan Wang
## 24Nov2020
## Single Cell Transcriptome Analysis of CAR T-Cell Products Reveals Subpopulations, Stimulation and Exhaustion Signatures

## Abstract
Chimeric antigen receptor (CAR) T-cell adoptive therapy is set to transform the treatment of a rapidly expanding range of malignancies. Although the activation process of normal T cells is well characterised, comparatively little is known about the activation of cells via the CAR. Here we have used flow cytometry together with single cell transcriptome profiling to characterise the starting material (peripheral blood mononuclear cells) and CAR therapeuticproductsof 3healthy donors in the presence and absence of antigen specific stimulation. Analysis of 53,191 single cell transcriptomes showed APRIL-based CAR products to contain several subpopulations of cells, withcellular composition reproducible from donor to donor, and all major cellular subsets compatible with CAR expression.Only 50% of CAR-expressing cellsdisplayed transcriptional changes upon CAR-specific antigen exposure. The resultingmolecular signature for CAR T-cell activation provides a rich resource for future dissection of underlying mechanisms.Targeted data interrogation also revealed that a small proportion of antigen-responding CAR-expressing cells displayed an exhaustion signature, with both known markers and genes not previously associated with T-cell exhaustion.Comprehensive single cell transcriptomic analysis thus represents a powerful way to guide the assessment and optimization of clinical-grade CAR-T-cells, and inform future research into the underlying molecular processes.

## Notebooks
This folder contains all jupyter notebooks to regenerate figures used in the Autolus CAR T-cells paper.
 
  - <ins>**[PBMC_analysis.ipynb](https://github.com/SharonWang/CARTcells_Notebooks/blob/master/Notebooks/PBMC_analysis.ipynb)**</ins>:
    - Analysis of 3 PBMC samples using SCANPY
  - <ins>**[PBMC_Tcells_analysis.ipynb](https://github.com/SharonWang/CARTcells_Notebooks/blob/master/Notebooks/PBMC_Tcells_analysis.ipynb)**</ins>:
    - Analysis of T cells in PBMC samples using SCANPY
  - <ins>**[PBMC_project_onto_PBMC68K.ipynb](https://github.com/SharonWang/CARTcells_Notebooks/blob/master/Notebooks/PBMC_project_onto_PBMC68K.ipynb)**</ins>:
    - Projection of PBMC samples onto [Zheng2017](https://www.nature.com/articles/ncomms14049) landscape
  - <ins>**[Tcells_analysis.ipynb](https://github.com/SharonWang/CARTcells_Notebooks/blob/master/Notebooks/Tcells_analysis.ipynb)**</ins>:
    - Analysis of *in vitro* T cell samples using SCANPY
  - <ins>**[Tcells_Ex_Sig.ipynb](https://github.com/SharonWang/CARTcells_Notebooks/blob/master/Notebooks/Tcells_Ex_Sig.ipynb)**</ins>:
    - Calculation of T cell Exhaustion signature
  - <ins>**[PBMC_Tcells_and_Tcells.ipynb](https://github.com/SharonWang/CARTcells_Notebooks/blob/master/Notebooks/PBMC_Tcells_and_Tcells.ipynb)**</ins>:
    - Integration of PBMC T cells and *in vitro* T cells
