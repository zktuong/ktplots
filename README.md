# ktplots
R plotting functions to plot gene expression data of single-cell data.


## Installation instructions
You can install the package via ```devtools::install_github()``` function in R
```R
library(devtools)
devtools::install_github('zktuong/ktplots', dependencies = TRUE)

# one function requires SummarizedExperiment from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
```
## Usage instructions
```R
library(ktplots)
```

### geneDotPlot
plotting genexpression dotplots heatmaps
```R
geneDotPlot(seurat_object, as.factor(seurat_object$split), genes = c("CD68","HLA-DRA"), save.plot = FALSE, groups = c("Tissue","PBMC"))
```

### plot_cpdb
Generates the dot plot for cpdb output
```R
pvals <- read.delim("pvalues.txt", check.names = FALSE)
means <- read.delim("means.txt", check.names = FALSE) 
plot_cpdb("Bcell", "Tcell", means, pvals, groups = c("normal", "tumor"), genes = c("CXCL13", "CD274", "CXCR5"))
```
