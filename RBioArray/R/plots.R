#' @title rbioarray_hcluster
#'
#' @description Wrapper for hierarchical clustering analysis and heatmap visualization.
#' @param plotName File name for the export \code{pdf} plot file. Default is \code{"data"}.
#' @param normlist Input normalized data, either a list, \code{EList} or \code{MAList} object.
#' @param n Number of genes to be clustered, numeric input or \code{"all"}. Default is \code{"all"}.
#' @param fct Input \code{factor} object for samples.
#' @param colGroup Colour group, numeric or dependent on \code{fct}.
#' @param distance Distance calculation method. Default is \code{"euclidean"}. See \code{\link{dist}} for more.
#' @param clust Clustering method. Default is \code{"complete"}. See \code{\link{hclust}} for more.
#' @param colColour Column group colour. Default is \code{"Paired"}. See \code{RColorBrewer} package for more.
#' @param mapColour Heat map colour. Default is \code{"PRGn"}. See \code{RColorBrewer} package for more.
#' @param n_mapColour Number of colours displayed. Default is \code{11}. See \code{RColorBrewer} package for more.
#' @param ... Additional arguments for \code{heatmap.2} function from \code{gplots} package.
#' @param plotWidth Width of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param plotHeight Height of the plot. Unit is \code{inch}. Default is \code{7}.
#' @return A heatmap based on hierarchical clustering analysis in \code{pdf} format.
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' rbioarray_hcluster(normlist = normdata, n = 500, fct = conSum, trace = "none", srtCol = 45, offsetCol = 0, adjCol = c(1, 0), labRow = FALSE, key.title = "", keysize = 1.5, key.xlab = "Normalized expression value", key.ylab = "Probe count")
#' }
#' @export
rbioarray_hcluster <- function(plotName = "data", normlist, n = "all",
                 fct, colGroup = ifelse(length(levels(fct)) < 19, length(levels(fct)), 19),
                 distance = "euclidean", clust = "complete",
                 colColour = "Paired", mapColour = "PRGn", n_mapColour = 11, ...,
                 plotWidth = 7, plotHeight = 7){

  ## prepare matrix for plotting
  dfm <- data.frame(normlist$genes, normlist$E)

  if (n == "all"){
    dfm <- dfm[dfm$ControlType == 0, ] # remove control probes
  } else {
    dfm <- dfm[dfm$ControlType == 0, ] # remove control probes
    dfm <- dfm[1:n, ] # subset
  }

  mtx <- as.matrix(dfm[, -c(1:2)])
  rownames(mtx) <- dfm[, 1]


  ## heatmap
  # set up dis and cluster functions
  distfunc <- function(x)dist(x, method = distance)
  clustfunc <- function(x)hclust(x, method = clust)

  # set ColSideColors

  col_cluster <- clustfunc(distfunc(t(dfm[, -c(1:2)])))
  colG <- cutree(col_cluster, colGroup) # column group
  colC <- brewer.pal(ifelse(colGroup < 3, 3, colGroup), colColour) # column colour

  # draw heatmap
  pdf(file = paste(plotName, "_heatmap.pdf", sep = ""), width = plotWidth, height = plotHeight)
  heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
            col = brewer.pal(n_mapColour, mapColour), ColSideColors = colC[colG], ...)
  dev.off()
}
