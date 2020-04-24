[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![codecov](https://codecov.io/gh/zktuong/ktplots/branch/master/graph/badge.svg)](https://codecov.io/gh/zktuong/ktplots)
[![Build Status](https://travis-ci.com/zktuong/ktplots.svg?branch=master)](https://travis-ci.com/zktuong/ktplots)
# ktplots
R plotting functions to plot gene expression data of single-cell data.


## Installation instructions
You can install the package via ```devtools::install_github()``` function in R
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
devtools::install_github('zktuong/ktplots', dependencies = TRUE)
# also maybe install Seurat if you haven't already
install.packages('Seurat')
```
## Usage instructions
```R
library(ktplots)
```
There is a test dataset for most of these functions that can be loaded via
```R
# note, you need to load Seurat to interact with it
library(Seurat)
data(kidneyimmune)
```
It is a seurat object which i downsampled from the [kidney cell atlas](https://kidneycellatlas.org). For more info, please see [Stewart et al. kidney single cell data set published in Science 2019](https://science.sciencemag.org/content/365/6460/1461).

### geneDotPlot
plotting gene expression dot plots heatmaps
```R
# Note, this conflicts with tidyr devel version
geneDotPlot(scdata = kidneyimmune, # object 
	genes = c("CD68", "CD80", "CD86", "CD74", "CD2", "CD5"), # genes to plot
	idents = "celltype", # column name in meta data that holds the cell-cluster ID/assignment
	split.by = 'Project', # column name in the meta data that you want to split the plotting by. If not provided, it will just plot according to idents
	standard_scale = TRUE) + # whether to scale expression values from 0 to 1. See ?geneDotPlot for other options
theme(strip.text.x = element_text(angle=0, hjust = 0, size =7)) + small_guide() + small_legend()
```
hopefully you end up with something like this
![geneDotPlot](exampleImages/geneDotPlot_example.png)

### plot_cpdb
Generates the dot plot for CellPhoneDB analysis via specifying the cell types and the genes. Importantly, the plotting is largely controlled by what is the meta file provided to CellPhoneDB analysis i.e. the annotations must match if you specify split.by or not.
```R
# read in pvalues.txt and means.txt that are in the out folder after CellPhoneDB analysis
# pvals <- read.delim("pvalues.txt", check.names = FALSE)
# means <- read.delim("means.txt", check.names = FALSE)

# I've provided an example dataset
data(cdpb_outout) 
plot_cpdb(cell_type1 = 'B cell', cell_type2 = 'CD4T cell',
	scdata = kidneyimmune,
	idents = 'celltype', # column name where the cell ids are located in the metadata
	means = means,
	pvals = pvals,
	split.by = 'Experiment', # column name where the grouping column is. Optional.
	genes = c("XCR1", "CXCL10", "CCL5")) + # some helper functions included in ktplots to help with the plotting
small_axis(fontsize = 3) + small_grid() + small_guide() + small_legend(fontsize = 2) 
```
![plot_cpdb](exampleImages/plot_cpdb_example.png)

You can try by a crude grep via the 'gene.family'
```R
some examples
```R
plot_cpdb(cell_type1 = 'B cell', cell_type2 = 'CD4T cell', scdata = kidneyimmune,
	idents = 'celltype', means = means, pvals = pvals, split.by = 'Experiment',
	gene.family = 'chemokines') + small_guide() + small_axis() + small_legend(keysize=.5)
```
![plot_cpdb](exampleImages/plot_cpdb_example1.png)
```R
plot_cpdb(cell_type1 = 'B cell', cell_type2 = 'CD4T cell', scdata = kidneyimmune,
	idents = 'celltype', means = means, pvals = pvals, split.by = 'Experiment',
	gene.family = 'chemokines', col_option = "maroon", highlight = "blue") + small_guide() + small_axis() + small_legend(keysize=.5)
```
![plot_cpdb](exampleImages/plot_cpdb_example2.png)
```R
plot_cpdb(cell_type1 = 'B cell', cell_type2 = 'CD4T cell', scdata = kidneyimmune,
	idents = 'celltype', means = means, pvals = pvals, split.by = 'Experiment',
	gene.family = 'chemokines', col_option = viridis::cividis(50)) + small_guide() + small_axis() + small_legend(keysize=.5)
```
![plot_cpdb](exampleImages/plot_cpdb_example3.png)
```R
plot_cpdb(cell_type1 = 'B cell', cell_type2 = 'CD4T cell', scdata = kidneyimmune,
	idents = 'celltype', means = means, pvals = pvals, split.by = 'Experiment',
	gene.family = 'chemokines', noir = TRUE) + small_guide() + small_axis() + small_legend(keysize=.5)
```
![plot_cpdb](exampleImages/plot_cpdb_example4.png)

if ```genes``` and ```gene.family``` are both not specified, by default it will plot everything.
Specifying ```keep_significant_only``` will only keep those that are p<0.05 (with or without adjustment via ```p.adjust.method```).

### StackedVlnPlot
Generates a stacked violinplot like in scanpy's ```sc.pl.stacked_violin```. Credits to [@tangming2005](https://twitter.com/tangming2005)
```R
features <- c("CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "LGALS3", "S100A8", "GNLY", "NKG7", "KLRB1", "FCGR3A", "FCER1A", "CST3")
StackedVlnPlot(kidneyimmune, features = features) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
```
![StackedVlnPlot](exampleImages/StackedVlnPlot_example.png)

### rainCloudPlot
generates a raincloudplot to use boxplot, scatterplot and violin all at once!
kudos to [https://wellcomeopenresearch.org/articles/4-63](https://wellcomeopenresearch.org/articles/4-63)
```R
rainCloudPlot(data = kidneyimmune@meta.data, groupby = "celltype", parameter = "n_counts") + coord_flip()
```
![rainCloudPlot](exampleImages/rainCloudPlot_example.png)

### small_legend/small_guide/small_axis/small_grid/topright_legend/topleft_legend/bottomleft_legend/bottomright_legend
As shown above, some functions to quickly adjust the size and position of ggplots.
```R
# for example
g <- Seurat::DimPlot(kidneyimmune, color = "celltype")
g1 <- g + small_legend() + small_guide() + small_axis() + bottomleft_legend() 
library(patchwork)
g + g1
```
![gghelperfunctions](exampleImages/gghelperfunctions_example.png)