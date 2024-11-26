
#' Visualize the spatial distribution of cell type proportion
#'
#' @param proportion Data frame of cell type proportion, with each row representing
#'                   a spot and each column representing a cell type
#' @param location A matrix or data frame of coordinates, with each row representing
#'                 a point and two columns for x and y coordinates.
#' @param ct.visualize vector containing the cell types of interest you want to plot
#' @param colors Vector of color names, default color scale is c("lightblue","lightyellow","red")
#' @param NumCols Number of columns for display when combining plots
#' @param px.size Adjust point size for plotting
#' @param ratio Numeric. Aspect ratio for the plot. Default is 1.
#' @param max.cf Vector of minimum and maximum cutoff values for each cell type
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom scales squish
#'
#' @return Returns a ggplot object
#' @export
#'

plot_ct <- function (proportion, location, ct.visualize, colors = c("lightblue", "lightyellow", "red"), NumCols = 4, px.size = 1, ratio = 1, max.cf = 1)
{
  if (is.null(colors)) {
    colors = c("lightblue", "lightyellow", "red")
  }
  else {
    colors = colors
  }
  prop = as.data.frame(proportion)
  prop = prop[, order(colnames(prop))]
  position = as.data.frame(location)
  if (sum(rownames(prop) == rownames(position)) != nrow(prop)) {
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  ct.select = ct.visualize
  prop = prop[, ct.select]
  prop_scale = as.data.frame(apply(prop, 2, function(x) {
    (x - min(x))/(max(x) - min(x))
  }))
  prop_scale$x = as.numeric(position$x)
  prop_scale$y = as.numeric(position$y)
  mData = melt(prop_scale, id.vars = c("x", "y"))
  colnames(mData)[3] <- "Cell_Type"
  b = c(0, 1)
  p = suppressMessages(ggplot(mData, aes(x, y)) + geom_point(aes(colour = value),size = px.size) +
                         scale_color_gradientn(colours = colors,limits=c(0, max.cf),oob = squish) +
                         ##scale_x_discrete(expand = c(0.03, 1)) + scale_y_discrete(expand = c(0.03, 1)) +
                         scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0, 1)) +
                         facet_wrap(~Cell_Type, ncol = NumCols) + coord_fixed(ratio=ratio) +
                         theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
                               panel.background = element_blank(), plot.background = element_blank(),
                               panel.border = element_rect(colour = "grey89", fill = NA, size = 1),
                               axis.text = element_blank(), axis.ticks = element_blank(),axis.title = element_blank(),
                               legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 9),
                               legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.4, "cm"),
                               strip.text = element_text(size = 12, face = "bold"),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.3,"cm")))
  return(p)
}



#' Visualize the spatial distribution of cell type proportion in a geom scatterpie plot
#'
#' @param proportion Data frame of cell type proportion, with each row representing
#'                   a spot and each column representing a cell type
#' @param location A matrix or data frame of coordinates, with each row representing
#'                 a point and two columns for x and y coordinates.
#' @param colors Vector of color names
#' @param title Text for the title
#' @param legend Whether to remove the legend
#' @param legend.title Text for the legend title
#' @param px.size Adjust point size for plotting
#' @param ratio Numeric. Aspect ratio for the plot. Default is 1.
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scatterpie geom_scatterpie
#' @importFrom grDevices colorRampPalette
#' @importFrom gtools mixedsort
#' @return Returns a ggplot object
#' @export
#'

plot_pie <- function (proportion, location, colors = NULL, px.size = 0.52, title = NULL, legend = T, legend.title = "Cell Type", ratio = 1)
{
  prop = as.data.frame(proportion)
  prop = prop[, mixedsort(colnames(prop))]
  position = as.data.frame(location)
  if (sum(rownames(prop) == rownames(position)) != nrow(prop)) {
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  colorCandidate = c( "#d52828", "#da8ec8", "#ceaaa3", "#2c9f2c",
                      "#babc22", "#9267bc", "#8b544c", "#e277c1","#59BE86",
                      "#d42728", "#adc6e8", "#97df89", "#fe9795", "#4381bd",
                      "#f2941f", "#5aa43a", "#cc4d2e", "#9f83c8", "#91675a",
                      "#ff7d0b", "#929292", "#c3c237", "#b4e0ea", "#bacceb",
                      "#f7c685", "#dcf0d0", "#f4a99f", "#c8bad8", "#F56867",
                      "#DADADA","#FEB915", "#C798EE", "#1e77b4", "#7495D3",
                      "#6D1A9C", "#15821E", "#3A84E6", "#997273", "#787878",
                      "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C", "#93796C",
                      "#F9BD3F", "#DAB370", "#877F6C", "#268785")
  if (is.null(colors)) {
    if (ncol(prop) > length(colorCandidate)) {
      colors = colorRampPalette(colorCandidate)(ncol(prop))
    }
    else {
      colors = colorCandidate[sample(1:length(colorCandidate),
                                     ncol(prop))]
    }
  }
  else {
    colors = colors
  }
  data = cbind(prop, position)
  ct.select = colnames(prop)
  if(legend == F){
    p = suppressMessages(
      ggplot() + geom_scatterpie(aes(x = x, y = y, r = px.size), data = data, cols = ct.select,color = NA) +
        coord_fixed(ratio = ratio) + ggtitle(title) + scale_fill_manual(values = colors) +
        scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0, 1)) +
        theme(title = element_text(size = 18, face = "bold"),
              plot.title = element_text(hjust = 0.5),
              plot.margin = margin(0.1, 0.1, 1, 0.1, "cm"),
              panel.background = element_blank(), plot.background = element_blank(),
              panel.border = element_rect(colour = "grey89", fill = NA, size = 0.5),
              axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
              strip.text = element_text(size = 16, face = "bold")) + NoLegend())
  }
  else{
    p = suppressMessages(
      ggplot() + geom_scatterpie(aes(x = x, y = y, r = px.size), data = data, cols = ct.select, color = NA) +
        coord_fixed(ratio = ratio) + ggtitle(title) + scale_fill_manual(values = colors) +
        scale_x_discrete(expand = c(0.008, 1)) + scale_y_discrete(expand = c(0.008, 1)) +
        theme(title = element_text(size = 18, face = "bold"),
              plot.title = element_text(hjust = 0.5),
              plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
              panel.background = element_blank(), plot.background = element_blank(),
              panel.border = element_rect(colour = "grey89", fill = NA, size = 0.5),
              axis.text = element_blank(), axis.ticks = element_blank(),
              axis.title = element_blank(), legend.title = element_text(size = 16, face = "bold"),
              legend.key = element_rect(colour = "transparent",fill = "white"),
              legend.key.size = unit(0.5, "cm"),legend.text = element_text(size = 12),
              legend.spacing.y = unit(1, 'cm'), legend.position = "bottom",
              strip.text = element_text(size = 16,face = "bold")) +
        guides(fill = guide_legend(title = legend.title)))
  }
  return(p)
}
