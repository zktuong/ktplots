[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![codecov](https://codecov.io/gh/zktuong/ktplots/branch/master/graph/badge.svg)](https://codecov.io/gh/zktuong/ktplots)
[![R](https://github.com/zktuong/ktplots/actions/workflows/r.yml/badge.svg)](https://github.com/zktuong/ktplots/actions/workflows/r.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6728350.svg)](https://doi.org/10.5281/zenodo.5717922)

# ktplots
R plotting functions to plot gene expression data of single-cell data.

For a python port of `ktplots`, please check out my other [repository](https://www.github.com/zktuong/ktplotspy).

## Installation instructions
You can install the package via `devtools::install_github()` function in R
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
devtools::install_github('zktuong/ktplots', dependencies = TRUE)
```
## Usage instructions
```R
library(ktplots)
```
There is a test dataset in `SingleCellExperiment` format to test the functions.
```R
library(SingleCellExperiment)
data(kidneyimmune)
#Some functions accept Seurat objects too.
```
The data is downsampled from the [kidney cell atlas](https://kidneycellatlas.org).

For more info, please see [Stewart et al. kidney single cell data set published in Science 2019](https://science.sciencemag.org/content/365/6460/1461).

## Vignette

Please go to the official [vignette](https://zktuong.github.io/ktplots/articles/vignette.html) for more information.

For the legacy version of the `README.md` file, please go [here](https://github.com/zktuong/ktplots/blob/master/_legacy_README.md).

## Citation
If you find these functions useful, please consider leaving a star, citing this repository, and/or citing the following [DOI](https://doi.org/10.5281/zenodo.5717922):

To cite a specific version of `ktplots`, please follow the links on the zenodo repository. e.g. v2.0.0:
```
ZK Tuong. (2021). zktuong/ktplots: 2.0.0 (v2.0.0). Zenodo. https://doi.org/10.5281/zenodo.5717922
```

Thank you!
