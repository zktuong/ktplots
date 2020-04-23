#' Plotting stacked violin plot
#'
## main function
#' @name StackedVlnPlot
#' @param obj single-cell data. can be seurat/summarizedexperiment object
#' @param features genes/features to plot
#' @param pt.size size of dots
#' @param plot.margin to adjust the white space between each plot.
#' @param ... passed to Seurat::VlnPlot
#' @return StackedVlnPlot
#' @examples
#' features<- c("CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "LGALS3", "S100A8", "GNLY", "NKG7", "KLRB1", "FCGR3A", "FCER1A", "CST3")
#' StackedVlnPlot(obj = pbmc, features = features)
#' @import ggplot2 
#' @export
StackedVlnPlot <- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

## remove the x-axis text and tick
#' @name misc
#' @param obj single-cell data. can be seurat/summarizedexperiment object
#' @param feature genes/features to plot
#' @param pt.size size of dots
#' @param plot.margin to adjust the white space between each plot.
#' @param ... passed to Seurat::VlnPlot
#' @export
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- Seurat::VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
#' @name misc
#' @param p value
#' @export
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

