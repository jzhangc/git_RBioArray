#' @title rbioarray_hcluster
#'
#' @description Wrapper for hierarchical clustering analysis and heatmap visualization.
#' @param plotName File name for the export \code{pdf} plot file. Default is \code{"data"}.
#' @param fltlist Input filtered data, either a list, \code{EList} or \code{MAList} object.
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
#' rbioarray_hcluster(fltlist = normdata, n = 500, fct = conSum, trace = "none", srtCol = 45, offsetCol = 0, adjCol = c(1, 0), labRow = FALSE, key.title = "", keysize = 1.5, key.xlab = "Normalized expression value", key.ylab = "Probe count")
#' }
#' @export
rbioarray_hcluster <- function(plotName = "data", fltlist, n = "all",
                 fct, colGroup = ifelse(length(levels(fct)) < 19, length(levels(fct)), 19),
                 distance = "euclidean", clust = "complete",
                 colColour = "Paired", mapColour = "PRGn", n_mapColour = 11, ...,
                 plotWidth = 7, plotHeight = 7){

  ## prepare matrix for plotting
  dfm <- data.frame(fltlist$genes, fltlist$E)

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



#' @title rbioarray_hcluster_super
#'
#' @description Wrapper for supervised hierarchical clustering analysis and heatmap visualization.
#' @param plotName File name for the export \code{pdf} plot file. Default is \code{"data"}.
#' @param fltDOI Based on filtered data, a subset corresponding to the comparasion, either a list, \code{EList} or \code{MAList} object.
#' @param dfmDE A subset of the DE list, i.e. a \code{topTable} dataframe, corresponding to the comparasion (i.e., contrast).
#' @param pcutoff P value cut off. Default is \code{NULL}.
#' @param method Thresholding method, "fdr" or "spikein". Default is \code{"spikein"}.
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
#' @return A supervised heatmap based on hierarchical clustering analysis in \code{pdf} format.
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' rbioarray_hcluster_super(normlist = normdata, n = 500, fct = conSum, trace = "none", srtCol = 45, offsetCol = 0, adjCol = c(1, 0), labRow = FALSE, key.title = "", keysize = 1.5, key.xlab = "Normalized expression value", key.ylab = "Probe count")
#' }
#' @export
rbioarray_hcluster_super <- function(plotName = "data", fltDOI, dfmDE, pcutoff = NULL, method = "spikein",
                                     fct, colGroup = ifelse(length(levels(fct)) < 19, length(levels(fct)), 19),
                                     distance = "euclidean", clust = "complete",
                                     colColour = "Paired", mapColour = "PRGn", n_mapColour = 11, ...,
                                     plotWidth = 7, plotHeight = 7){ #DOI: fltered subset data of interest


  ## prepare matrix for plotting
  dfm <- data.frame(fltDOI$genes, fltDOI$E)
  dfm <- dfm[dfm$ControlType == 0, ] # remove control probes

  if (is.null(pcutoff)){
    stop("Please set p value threshold.")
  } else {
    ifelse(method == "fdr", pb_name <- dfmDE[dfmDE$adj.P.Val <= pcutoff, "ProbeName"], pb_name <- dfmDE[dfmDE$P.Value <= pcutoff, "ProbeName"])
    dfm <- dfm[dfm$ProbeName %in% pb_name, ]
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
  pdf(file = paste(plotName, "_heatmap.supervised.pdf", sep = ""), width = plotWidth, height = plotHeight)
  heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
            col = brewer.pal(n_mapColour, mapColour), ColSideColors = colC[colG], ...)
  dev.off()
}
