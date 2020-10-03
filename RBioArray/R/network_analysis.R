#' @title rbio_tom
#'
#' @description A plot function for K means cluster results, with or without PCA.
#' @param mtx Matrix. TBC.
#' @param diag Boolean. TBC.
#' @param power Integer. TBC.
#' @param tom_type String. TBC.
#' @param ... Additional arguments passed to \code{TOMdist()} function.
#' @param hclust.method String.
#' @param k Integer. TBC.
#' @param h Numeric. TBC.
#' @param plot.dendro Boolean. Whether to plot a dendrogram. Default is \code{TRUE}.
#' @param plot.export.name String. The prefix for the exported figure file name. Default is \code{NULL}.
#' @param plot.title Boolean. The dendrogram plot title. Default is \code{NULL}.
#' @param plot.title.size Numeric. The dendrogram plot title size. Default is \code{16}.
#' @param plot.dendroline.size Numeric. TBC.
#' @param plot.dendrolabel Boolean. TBC.
#' @param plot.dendrolabel.size Numeric. TBC.
#' @param plot.dendrolabel.space Numeric. TBC.
#' @param plot.ylabel.size Numeric. TBC.
#' @param plot.width The dendrogram plot width. Default is \code{170}.
#' @param plot.height The dendrogram plot height. Default is \code{150}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return
#'         TBC
#' @details
#'         TBC
#' @import ggplot2
#' @import igraph
#' @importFrom grid grid.draw
#' @importFrom stringr str_pad
#' @importFrom WGCNA TOMdist
#' @importFrom dendextend color_branches as.ggdend
#' @examples
#' \dontrun{
#'          tom_g <- rbio_tom(mtx = t(mtx), k = 10, plot.export.name = "ex1",
#'                            plot.dendroline.size = 0.5,
#'                            plot.dendrolabel = TRUE,
#'                            plot.dendrolabel.size = 0.2,
#'                            plot.dendrolabelspace = 2,
#'                            plot.ylabel.size = 1.5)
#' }
#' @export
rbio_tom <- function(mtx,
                     diag = FALSE,
                     power = 6, tom_type = c("unsigned", "signed"), ...,
                     hclust.method = c("complete", "ward.D", "ward.D2", "single",  "average", "mcquitty", "median", "centroid"),
                     k = NULL, h = NULL,
                     plot.dendro = TRUE,
                     plot.export.name = NULL,
                     plot.title = NULL, plot.title.size = 16,
                     plot.dendroline.size = 0.8,
                     plot.dendrolabel = TRUE, plot.dendrolabel.size = 0.8,
                     plot.dendrolabel.space = 2,
                     plot.ylabel.size = 1.5,
                     plot.width = 150, plot.height = 150,
                     verbose = TRUE) {
  # - argument check -
  if (!any(class(mtx) %in% "matrix")) stop("mtx has to be a matrix.")
  if (is.null(rownames(mtx))) rownames(mtx) <- seq(nrow(mtx))
  if (is.null(colnames(mtx))) rownames(mtx) <- paste0("f_", seq(ncol(mtx)))
  tom_type <- match.arg(tom_type, c("unsigned", "signed"))
  hclust.method <- match.arg(hclust.method, c("complete", "ward.D", "ward.D2", "single",  "average", "mcquitty", "median", "centroid"))
  if (is.null(plot.export.name)){
    plot.export.name <- deparse(substitute(mtx))
  } else {
    plot.export.name <- plot.export.name
  }
  if (!is.null(k) && k > ncol(mtx)) stop(paste0("k cannot be greater than the number of items for correlating/clustering, i.e. ncol(mtx) = ", ncol(mtx)))

  # - TOM calculation -
  adjmat <- cor(mtx)^power
  # tom_dist <- TOMdist(adjmat, TOMType = tom_type, ...)  # matrix, array class, here we use dist
  tom_dist <- TOMdist(adjmat, TOMType = tom_type)  # matrix, array class, here we use dist
  rownames(tom_dist) <- rownames(adjmat)
  colnames(tom_dist) <- colnames(adjmat)

  # - TOM hclust -
  tom_dist <- as.dist(tom_dist)  # convert to an R distance object
  tom_dist_hclust <- hclust(tom_dist, method = hclust.method)

  # tom_membership and tom similarity
  tom_membership <- cutree(tom_dist_hclust, h = h, k = k)
  membersihp_for_dendro <- tom_membership
  tom_similarity <- 1 - tom_dist # edge always uses similarity
  g_adjmat <- as.matrix(tom_similarity)
  g_adjmat <- g_adjmat[order(tom_membership), order(tom_membership)]  # reorder it
  tom_membership <- tom_membership[order(tom_membership)] # update tom_membership

  if(plot.dendro) {
    if (verbose) cat("Saving hcluster dendrogram...")
    dendr <- color_branches(tom_dist_hclust, h = h, k = k)
    dendr <- as.ggdend(dendr)
    dendr$segments$lwd <- rep(plot.dendroline.size, times = length(dendr$segments$lwd))  # dendro line size
    dendr$labels$cex <- rep(plot.dendrolabel.size, times = length(dendr$labels$cex))  # label sizes
    dendr$labels$label <- foreach(i = as.character(dendr$labels$label), .combine = "c", .export = "dendr") %do% str_pad(i, width = nchar(i) + plot.dendrolabel.space, side = "right", pad = " ")
    dendr$labels$label <- factor(dendr$labels$label, levels = unique(dendr$labels$label))

    if (plot.dendrolabel) {
      p <- ggplot(dendr)
    } else {
      p <- ggplot(dendr, labels = FALSE)
    }
    p <- p +
      ggtitle(plot.title) +
      scale_x_continuous(labels = NULL, position = "top") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = plot.title.size, face = "bold"),
            axis.title = element_blank(),
            axis.text.y = element_text(size = rel(plot.ylabel.size)),
            panel.grid = element_blank())
    ggsave(filename = paste0(plot.export.name, "_tom_hclust.pdf"), plot = p,
           width = plot.width, height = plot.height, units = "mm", dpi = 600)
    grid.draw(p)
    if (verbose) cat("Done!")
  }

  # - igraph -
  out <- list(g = graph.adjacency(
    g_adjmat,
    mode = "undirected",
    weighted = TRUE,
    diag = diag
  ),
  hclust = tom_dist_hclust,
  k = k,
  h = h,
  tom_membership = tom_membership)
  class(out) <- "rbio_tom_graph"
  return(out)
}
