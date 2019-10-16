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
geneDotPlot(scdata = seurat_object, # object 
	idents = Idents(seurat_object), # a vector holding the cell-cluster ID/assignment or some other vector such as those found in the metadata seurat_object$split
	genes = c("CD68", "CD80", "CD86", "CD74", "CD2", "CD5"), # genes to plot
	split.by = "group", # column name in the meta data that you want to split the plotting by. If not provided, it will just plot according to idents
	save.plot = FALSE) # If TRUE, it will save to a location that you can specify via filepath and filename
```
hopefully you end up with something like this
![heatmap](exampleImages/geneDotPlot_example.png)

### plot_cpdb
Generates the dot plot for cpdb output vai specify the cell types and the genes
```R
pvals <- read.delim("pvalues.txt", check.names = FALSE)
means <- read.delim("means.txt", check.names = FALSE) 
plot_cpdb(cell_type1 = "Bcell", # cell_type1 and cell_type2 will call grep, so this will accept regex arguments
	cell_type2 = "Tcell",
	means,
	pvals,
	groups = c("normal", "tumor"),
	genes = c("CXCL13", "CD274", "CXCR5"))
```

or, you can try by a crude grep via the 'gene.family'
```R
pvals <- read.delim("pvalues.txt", check.names = FALSE)
means <- read.delim("means.txt", check.names = FALSE) 
plot_cpdb(cell_type1 = "Bcell",
	cell_type2 = "Tcell",
	means,
	pvals,
	groups = c("normal", "tumor"),
	gene.family = "chemokines") # can also try Th1, Th2, Th17, Treg, costimulatory, coinhibitory, niche, 
```
example of what appears
![heatmap](exampleImages/plot_cpdb_example.png)
