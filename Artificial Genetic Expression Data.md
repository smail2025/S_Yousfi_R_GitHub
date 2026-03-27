# Artificial Genetic Expression Data Analysis

This project contains R scripts for **genetic expression data analysis**.  
The method used is a **density-based dimensionality reduction**, allowing the comparison and visualization of gene expression profiles.

## Project Folder
All files are located in the folder: `S_Yousfi_R_GitHub`

## Included Files
- `Suplementary_materiel_code.r` : contains all the functions necessary for the analysis
- `Suplementary_materiel_execution_code.r` : contains examples of how to run the functions with test data
- `README.md` : this file

## How to Use
| Package | Purpose / Usage |
|---------|----------------|
| `moments` | Calculation of statistical moments, such as skewness and kurtosis for each batch. |
| `MASS` | Generation of multivariate distributions and handling of covariance structures. |
| `copula` | Construction of Gaussian copulas (`normalCopula`) and multivariate distributions (`mvdc`) for simulating dependent Gamma data. |
| `mvtnorm` | Generation of multivariate normal or related distributions, often used with copulas. |
| `cluster` | Hierarchical clustering and cluster analysis. |
| `mclust` | Computation of clustering metrics such as Adjusted Rand Index (ARI). |
| `kernlab` | Kernel-based methods for computing similarity matrices and functional PCA. |
| `graphics` / `grDevices` | Plotting functions, including `plot()` and annotation of FPCA embeddings. |
| `parallel` | Parallelization of heavy computations across batches and dimensions. |
| `base` | Input/output functions such as `saveRDS()` and `readRDS()`. |

-  Place all files in the folder `S_Yousfi_R_GitHub`.
 - In R or RStudio, load the functions first with:

```r
source("Suplementary_materiel_code.r")
