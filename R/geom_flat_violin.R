#' geom_flat_violin
#' 
#' @param mapping see ggplot2::layer
#' @param data see ggplot2::layer
#' @param stat see ggplot2::layer
#' @param position see ggplot2::layer
#' @param trim see ggplot2::layer
#' @param scale see ggplot2::layer
#' @param show.legend see ggplot2::layer
#' @param inherit.aes see ggplot2::layer
#' @param ... passed to ggplot2::layer
#' @return geom_flat_violin
#' @examples
#' sourced from https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R
#' somewhat hackish solution to:
#' https://twitter.com/EamonCaddigan/status/646759751242620928
#' based mostly on copy/pasting from ggplot2 geom_violin source:
#' https://github.com/hadley/ggplot2/blob/master/R/geom-violin.r
#' 
#' ggplot(diamonds, aes(cut, carat)) +
#'     geom_flat_violin() +
#'     coord_flip()
#' @import dplyr
#' @import ggplot2 
#' @export
geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                        position = "dodge", trim = TRUE, scale = "area",
                        show.legend = NA, inherit.aes = TRUE, ...) {
  
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}
#' @name misc
#' @export
GeomFlatViolin <- ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %nin%
              params$width %nin% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
)

