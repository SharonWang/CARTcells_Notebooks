# CARTcells_Notebooks

## Xiaonan Wang
## 08Dec2020
## Single Cell Transcriptome Analysis of CAR T-Cell Products Reveals Subpopulations, Stimulation and Exhaustion Signatures

## Abstract
Chimeric antigen receptor (CAR) T-cell adoptive therapy is set to transform the treatment of a rapidly expanding range of malignancies. Although the activation process of normal T cells is well characterised, comparatively little is known about the activation of cells via the CAR. Here we have used flow cytometry together with single cell transcriptome profiling to characterise the starting material (peripheral blood mononuclear cells) and CAR therapeutic products of 3 healthy donors in the presence and absence of antigen specific stimulation. Analysis of 53,191 single cell transcriptomes showed APRIL-based CAR products to contain several subpopulations of cells, withcellular composition reproducible from donor to donor, and all major cellular subsets compatible with CAR expression. Only 50% of CAR-expressing cells displayed transcriptional changes upon CAR-specific antigen exposure. The resulting molecular signature for CAR T-cell activation provides a rich resource for future dissection of underlying mechanisms. Targeted data interrogation also revealed that a small proportion of antigen-responding CAR-expressing cells displayed an exhaustion signature, with both known markers and genes not previously associated with T-cell exhaustion.Comprehensive single cell transcriptomic analysis thus represents a powerful way to guide the assessment and optimization of clinical-grade CAR-T-cells, and inform future research into the underlying molecular processes.

## Schematic view
The experiments were done following the strategy below:

<p align="center"><img src="figures/Wang_Figure 1_Resubmission-1.png" alt="schematic" width="50%"></p>

## Notebooks
This folder contains all jupyter notebooks to regenerate figures used in the paper.
 
  - <ins>**[PBMC_analysis.ipynb](https://github.com/SharonWang/CARTcells_Notebooks/blob/master/Notebooks/PBMC_analysis.ipynb)**</ins>:
    - Analysis of 3 PBMC samples using SCANPY
  - <ins>**[PBMC_Tcells_analysis.ipynb](https://github.com/SharonWang/CARTcells_Notebooks/blob/master/Notebooks/PBMC_Tcells_analysis.ipynb)**</ins>:
    - Analysis of T cells in PBMC samples using SCANPY
  - <ins>**[Tcells_analysis.ipynb](https://github.com/SharonWang/CARTcells_Notebooks/blob/master/Notebooks/Tcells_analysis.ipynb)**</ins>:
    - Analysis of *in vitro* T cell samples using SCANPY
  - <ins>**[Tcells_Ex_Sig.ipynb](https://github.com/SharonWang/CARTcells_Notebooks/blob/master/Notebooks/Tcells_Ex_Sig.ipynb)**</ins>:
    - Calculation of T cell Exhaustion signature
  - <ins>**[PBMC_Tcells_and_Tcells.ipynb](https://github.com/SharonWang/CARTcells_Notebooks/blob/master/Notebooks/PBMC_Tcells_and_Tcells.ipynb)**</ins>:
    - Integration of PBMC T cells and *in vitro* T cells
