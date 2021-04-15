#' @title rbio_tom
#'
#' @description TOM (topological overlap measure) analysis.
#' @param mtx matrix. Input data matrix. The function calculates TOM distance/affinity for the columns items (or features).
#' @param diag Boolean. If to include diagonal for igraph construction. Default is \code{FALSE}.
#' @param cor_method String. Correlation method. Default is \code{"pearson"}.
#' @param power integer. The power to the correlation coefficients. Default is \code{6}.
#' @param tom_type string. Directionality of the TOM analysis. Default is \code{"unsigned"}.
#' @param ... Additional arguments passed to \code{TOMdist()} function.
#' @param manual_membership Boolean. If to use manual membership. Default is \code{FALSE}.
#' @param tom_membership Boolean. Set when \code{manual_membership = TRUE}, a \code{named integer vector} containing membership for each features. Default is \code{NULL}.
#' @param hclust.method String.
#' @param cutree.method String. Method to find the optimal k. Default is \code{"silhouette"}.
#' @param dynamictree.min.size Integer. When \code{cutree.method = "dynamic"}, the minimum cluster size. Default is \code{20}.
#' @param k integer. Set when \code{cutree.method = "manual"}, manually set number of groups to cut to.
#' @param h numeric. Set when \code{cutree.method = "manual"}, cut tree height.
#' @param plot.dendro Boolean. Whether to plot a dendrogram. Default is \code{TRUE}.
#' @param plot.export.name string. The prefix for the exported figure file name. Default is \code{NULL}.
#' @param plot.margins numeric vector. Plot margins, unit is "cm", order: \code{t, r, b, l}. Default is \code{c(2, 2, 2, 2)}.
#' @param plot.title Boolean. The dendrogram plot title. Default is \code{NULL}.
#' @param plot.title.size numeric. The dendrogram plot title size. Default is \code{16}.
#' @param plot.dendroline.size numeric. The dendrogram plot line size. Default is \code{1.5}.
#' @param plot.dendrolabel Boolean. If to display dendrogram label. Default is \code{TRUE}.
#' @param plot.dendrolabel.size numeric. Set when \code{plot.dendrolabel = TRUE}, dendrogram label size. Default is \code{1.5}.
#' @param plot.dendrolabel.space numeric. Set when \code{plot.dendrolabel = TRUE}, dendrogram label space. Default is \code{2}.
#' @param plot.ylabel.size numeric. Dendrogram y label size. Default is \code{1.5}.
#' @param plot.width numeric. The dendrogram plot width. Default is \code{170}.
#' @param plot.height numeric. The dendrogram plot height. Default is \code{150}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return
#'         An \code{rbio_tom_graph} object with the following items:
#'
#'         \code{g}: an igraph object for network analysis and visualization.
#'
#'         \code{tom_distance}: a \code{dist} object for TOM distance.
#'
#'         \code{tom_affinity}: a \code{matrix} object for TOM affinity.
#'
#'         \code{tom_dist_hclust}: a \code{hclust} object for TOM distance.
#'
#'         \code{cutree_method}: method used to cut hclust tree, to generate g membership.
#'
#'         \code{dynamictree_min_size}: when {cutree.method = "dynamic"}, the minimum tree size for dynamic tree cutting.
#'
#'         \code{silhouette_score_mean}: when {cutree.method = "silhouette"}, the mean silhouette score.
#'
#'         \code{k}: when \code{cutree.method = "manual"}, manually set number of groups to cut to.
#'
#'         \code{h}: when \code{cutree.method = "manual"}, cut tree height.
#'
#'         \code{manual_membership}: if manually set membership is used.
#'
#'         \code{tom_membership}: a \code{named integer vector} containing membership for each features clustered.
#'
#' @details
#'         When \code{cutree.method = "manual"}, \code{h} and \code{k} are mutually exclusive. Set one, but not both.
#'
#'         When {manual _membership = TRUE}, both hclust and tree cutting processes are ignored.
#'         When {manual _membership = TRUE}, the tom_membership should be a integer vector whose item names are feature names.
#'
#'         The function uses the \code{"tree"} method from the \code{\link{dynamicTreeCut::cutreeDynamic()}} to cut the tree when
#'         \code{cutree.method = "dynamic"}.
#'
#'         The \code{plot.margins} follow the base R setting in \code{\link{par}} for the positioning:
#'         b: mar[1], l: mar[2], t: mar[3], r: mar[4]
#'
#' @import ggplot2
#' @import igraph
#' @importFrom grid grid.draw
#' @importFrom cluster silhouette
#' @importFrom WGCNA TOMdist
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom dendextend color_branches as.ggdend
#' @examples
#' \dontrun{
#'            tom_g <- rbio_tom(mtx = t(mtx), k = 10,
#'                              plot.export.name = "ex1",
#'                              plot.title = "TOM distance hcluster", plot.title.size = 16,
#'                              plot.margins = c(0.5 ,0.5, 0.5, 0.5),
#'                              plot.dendroline.size = 0.2,
#'                              plot.dendrolabel = TRUE,
#'                              plot.dendrolabel.size = 1.5,
#'                              plot.dendrolabel.space = 2,
#'                              plot.ylabel.size = 1.5)
#' }
#' @export
rbio_tom <- function(mtx, diag = FALSE,
                     cor_method = c("pearson", "kendall", "spearman"),
                     power = 6, tom_type = c("unsigned", "signed"), ...,
                     manual_membership = FALSE, tom_membership = NULL,
                     hclust.method = c("complete", "ward.D", "ward.D2", "single",  "average", "mcquitty", "median", "centroid"),
                     cutree.method = c("manual", "silhouette", "dynamic"),
                     dynamictree.min.size = 20,
                     k = NULL, h = NULL,
                     plot.dendro = TRUE,
                     plot.export.name = NULL,
                     plot.margins = c(0.5, 0.5, 0.5, 0.5),
                     plot.title = NULL, plot.title.size = 16,
                     plot.dendroline.size = 0.8,
                     plot.dendrolabel = TRUE,
                     plot.dendrolabel.size = 1.5,
                     plot.dendrolabel.space = 2,
                     plot.ylabel.size = 1.5,
                     plot.width = 150, plot.height = 150,
                     verbose = TRUE) {
  # - argument check -
  if (!any(class(mtx) %in% "matrix")) stop("mtx has to be a matrix.")
  if (is.null(rownames(mtx))) rownames(mtx) <- seq(nrow(mtx))
  if (is.null(colnames(mtx))) rownames(mtx) <- paste0("f_", seq(ncol(mtx)))
  cor_method <- match.arg(cor_method, c("pearson", "kendall", "spearman"))
  tom_type <- match.arg(tom_type, c("unsigned", "signed"))
  if (manual_membership) {
    if (is.null(tom_membership)) {
      warning("no tom_membership is set when manual_membership = TRUE. Proceed with manual_membership = FALSE.\n")
      manual_membership <- FALSE
    } else if (!any(class(tom_membership) %in% c("integer", "numeric", "vector"))) {
      warning("tom_membership is not a named integer vector. Proceed with manual_membership = FALSE.\n")
      manual_membership <- FALSE
    } else if (ncol(mtx) != length(tom_membership)) {
      warning("tom_membership has different length as the input feature number. Proceed with manual_membership = FALSE.\n")
      manual_membership <- FALSE
    }
  }

  if (manual_membership) {
    tom_membership <- tom_membership
    hclust.method <- NA
    cutree.method <- NA
  } else {
    hclust.method <- match.arg(hclust.method, c("complete", "ward.D", "ward.D2", "single",  "average", "mcquitty", "median", "centroid"))
    cutree.method <- match.arg(cutree.method, c("manual", "silhouette", "dynamic"))
    if (!is.null(k) && k > ncol(mtx)) stop(paste0("k cannot be greater than the number of items for correlating/clustering, i.e. ncol(mtx) = ", ncol(mtx)))
  }

  if (is.null(plot.export.name)){
    plot.export.name <- deparse(substitute(mtx))
  } else {
    plot.export.name <- plot.export.name
  }


  # - TOM calculation -
  if (verbose) cat("TOM calulation...")
  adjmat <- cor(mtx, method = cor_method)^power
  tom_dist <- TOMdist(adjmat, TOMType = tom_type, verbose = FALSE, ...)  # matrix, array class, here we use dist
  # tom_dist <- TOMdist(adjmat, TOMType = tom_type, verbose = FALSE)  # matrix, array class, here we use dist
  rownames(tom_dist) <- rownames(adjmat)
  colnames(tom_dist) <- colnames(adjmat)
  tom_similarity <- 1 - tom_dist # convert to similarity matrix
  if (verbose) cat("Done!\n")

  if (manual_membership) {
    tom_dist_hclust = NA
    cutree.method = NA
    dynamictree.min.size = NA
    ss_mean = NA
    k = max(tom_membership)
    h = NA
  } else {  # hclus and tree cutting
    # - TOM hclust -
    tom_dist <- as.dist(tom_dist)  # convert to an R distance object
    tom_dist_hclust <- hclust(tom_dist, method = hclust.method)

    # decide k or h
    if (cutree.method == "dynamic") {
      if (verbose) cat("Dynamic tree cutting...")
      tom_membership <- cutreeDynamic(tom_dist_hclust, method = "tree", deepSplit = TRUE, minClusterSize = dynamictree.min.size, verbose = FALSE)
      if (all(tom_membership == 0)) stop("Dynamic tree cut failed to identify any functional clusters. Try changing the hcluster.method, minClusterSize, or cutree.method. ")
      names(tom_membership) <- tom_dist_hclust$labels
      k <- max(tom_membership)
      if (any(tom_membership == 0)) {
        n <- length(tom_membership[tom_membership == 0])
        tom_membership[tom_membership == 0] <- seq(from = k+1, to = k+n)
      } else {
        n <- 0
      }
      if (verbose) cat(paste0(k, " clusters, ", n, " items unassigned.\n\n"))
      ss_mean <- NULL
    } else if (cutree.method == "silhouette") {
      if (verbose) cat("Silhouette tree cutting...")
      k_range <- 2:(ncol(adjmat)-1)
      ss_mean <- foreach(i = k_range, .combine = "c") %do% {
        sil_m <- cutree(tom_dist_hclust, k = i)
        ss <- cluster::silhouette(sil_m, tom_dist)
        mean(ss[, 3])
      }
      names(ss_mean) <- k_range
      k <- as.integer(names(ss_mean)[ss_mean == max(ss_mean)])
      if (verbose) cat(paste0(k, " clusters..."))
      tom_membership <- stats::cutree(tom_dist_hclust, k = k)
      if (is.null(names(tom_membership))) names(tom_membership) <- as.character(seq(length(tom_membership)))
      dynamictree.min.size <- NULL
    } else {
      if (is.null(h) && is.null(k)) stop("Manually set h or k when cutree.method = \"manual\".")
      h <- h
      k <- k
      tom_membership <- stats::cutree(tom_dist_hclust, h = h, k = k)
      if (is.null(names(tom_membership))) names(tom_membership) <- as.character(seq(length(tom_membership)))
      ss_mean <- NULL
      dynamictree.min.size <- NULL
    }
    if (verbose) cat("Done!\n")

    if(plot.dendro) {
      if (verbose) cat("Saving hcluster dendrogram...")
      dendr <- color_branches(tom_dist_hclust, h = h, k = k)
      dendr <- as.ggdend(dendr)
      dendr$segments$lwd <- rep(plot.dendroline.size, times = length(dendr$segments$lwd))  # dendro line size
      # dendr$labels$cex <- rep(plot.dendrolabel.size, times = length(dendr$labels$cex))  # label sizes
      # dendr$labels$label <- foreach(i = as.character(dendr$labels$label), .combine = "c", .export = "dendr") %do% str_pad(i, width = nchar(i) + plot.dendrolabel.space, side = "right", pad = " ")
      # dendr$labels$label <- factor(dendr$labels$label, levels = unique(dendr$labels$label))
      # if (plot.dendrolabel) {
      #   p <- ggplot(dendr)
      # } else {
      #   p <- ggplot(dendr, labels = FALSE)
      # }
      p <- ggplot(dendr, labels = FALSE)
      p <- p +
        ggtitle(plot.title) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, size = plot.title.size, face = "bold"),
              axis.title = element_blank(),
              axis.text.y = element_text(size = rel(plot.ylabel.size)),
              panel.grid = element_blank(),
              plot.margin = margin(t = plot.margins[1], r = plot.margins[2],
                                   b = plot.margins[3], l = plot.margins[4],
                                   unit = "cm"))
      if (plot.dendrolabel){
        p <- p +
          scale_x_continuous(breaks = seq(length(dendr$labels$label)), labels = dendr$labels$label) +
          theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = plot.dendrolabel.size, margin = margin(t = plot.dendrolabel.space)))
      } else {
        p <- p +
          scale_x_continuous(breaks = NULL, labels = NULL)
      }
      ggsave(filename = paste0(plot.export.name, "_tom_hclust.pdf"), plot = p,
             width = plot.width, height = plot.height, units = "mm", dpi = 600)
      grid.draw(p)
      if (verbose) cat("Done!\n\n")
    } else {
      if (verbose) cat("\n")
    }
  }

  # igraph and final membership construction
  g_adjmat <- as.matrix(tom_similarity)  # igraph edge always uses similarity matrix
  g_adjmat <- g_adjmat[order(tom_membership), order(tom_membership)]  # reorder it
  tom_membership <- tom_membership[order(tom_membership)] # update tom_membership

  # - igraph -
  out <- list(g = graph.adjacency(
    g_adjmat,
    mode = "undirected",
    weighted = TRUE,
    diag = diag
  ),
  tom_distance = tom_dist,
  tom_affinity = tom_similarity,
  tom_dist_hclust = tom_dist_hclust,
  cutree_method = cutree.method,
  dynamictree_min_size = dynamictree.min.size,
  silhouette_score_mean = ss_mean,
  k = k,
  h = h,
  manual_membership = manual_membership,
  tom_membership = tom_membership)
  class(out) <- "rbio_tom_graph"
  return(out)
}


#' @title circle_text_func
#'
#' @description Companion function for the network visualization function: used for displaying node labels around
#'              the circular network figures.
#' @param g igraph object. Input network object.
#' @param circ_layout layout object. The circular layout object for \code{g}, usually derived from function \code{\link{igraph::layout.circle()}}.
#' @param text.label string vector. The label to display. Default is \code{NULL}.
#' @param text.size numeric or numeric vector. The text size for the labels. Default is \code{0.8}.
#' @param text.distance numeric. The distance multiplier between label and nodes. Default is \code{1.5}.
#' @param text.colour string. Label colour. Default is \code{"black"}.
#' @return Added text on the igraph plot.
#' @details The \code{text.size} argument also accepts a vector of sizes with a length equal to the number of vertices.
#'         When different unequal length is detect, the function uses the first number for a universal text size.
#' @examples
#' \dontrun{
#'      # plot the original circular plot
#'      g_layout <- layout.circle(g)
#'      plot(
#'        g,
#'        layout=g_layout,
#'        edge.curved=TRUE,
#'        vertex.size=vSizes,
#'        vertex.label = NA,
#'        vertex.label.dist=-1.5,
#'        vertex.label.color="black",
#'        vertex.label.cex=0.8,
#'        asp=FALSE,
#'        vertex.label.cex=0.6,
#'        edge.width=edgeweights,
#'        edge.arrow.mode=0,
#'        main="")
#'
#'     # add labels to the plot
#'      circle_text_func(g = g, circ_layout = g_layout)
#' }
#' @export
circle_text_func <- function(g, circ_layout,
                             text.label = NULL,
                             text.size = 0.8, text.distance = 1.5, text.colour = "black", ...){
  # - argument check -
  if (!any(class(g) %in% "igraph")) stop("g needs to be an igraph object.")

  # - text label and size processing -
  if (is.null(text.label)) { # text labels
    tLabels <- V(g)$name
  } else if (length(text.label) != length(V(g))) {
    warning("text labels not equal length with the verticee, proceeding with the internal vertex names.\n")
    tLabels <- V(g)$name
  } else {
    tLabels <- text.label
  }

  if (length(text.size) == 1) {  # text sizes
    tSize <- rep(text.size, times = length(V(g)$name))
  } else if (length(text.size) != length(V(g)$name)) {
    warning("text size not equal length with the verticee, proceeding with the first number for text size.\n")
    tSize <- rep(text.size[1], times = length(V(g)$names))
  } else {
    tSize <- text.size
  }

  # - Apply labels manually -
  # Specify x and y coordinates of labels, adjust outward as desired
  x = circ_layout[,1]*text.distance
  y = circ_layout[,2]*text.distance

  # - create vector of angles and justification for text based on number of nodes  -
  # (flipping the orientation of the words half way around so none appear upside down)
  angle = ifelse(atan(-(circ_layout[,1]/circ_layout[,2]))*(180/pi) < 0,  90 + atan(-(circ_layout[,1]/circ_layout[,2]))*(180/pi), 270 + atan(-circ_layout[,1]/circ_layout[,2])*(180/pi))
  # justification
  adj_val = ifelse(circ_layout[,1] > 0 & circ_layout[,2] > -1,  0, 1)

  #Apply the text labels with a loop with angle as srt
  for (i in 1:length(x)) {
    text(x=x[i], y=y[i], labels=tLabels[i], adj = adj_val[i], pos = NULL, cex = tSize[i], col = text.colour, srt = angle[i], xpd = T, ...)
  }
}


#' @title rbio_network
#'
#' @description (UNDER CONSTRUCTION: NOT FUNCTIONAL) Network construction and visualization function
#' @param object object containing a membership and igraph information.
#' @param ... Additional arguments for corresponding S3 class methods.
#' @return TBC
#' @examples
#'
#' \dontrun{
#' TBC
#' }
#'
#' @export
rbio_network <- function(object, ...){
  # - check object -
  if (!any(class(object) %in% c("rbio_tom_graph"))) stop("object needs to be in correct classes, e.g. rbio_tom_graph.")

  # - use method -
  UseMethod("rbio_network", object)
}


#' @title rbio_network.rbio_tom_graph
#' @rdname rbio_network
#' @method rbio_network rbio_tom_graph
#' @param object An \code{rbio_tom_graph} object.
#' @param export.name string. Optional prefix for export file name. Default is \code{NULL}.
#' @param ... Additional arguments for corresponding S3 class methods.
#' @return TBC
#' @examples
#'
#' \dontrun{
#' TBC
#' }
#'
#' @export
rbio_network.rbio_tom_graph <- function(object, export.name = NULL, ...){
  # - argument check -
  if (is.null(export.name)){
    export.name = deparse(substitute(object))
  }

  # - feed to the default method -
  # <TBC: under construction>
  rbio_network.default(g = object$g,
                       g_membership = object$tom_membership,
                       export.name = export.name, ...)
}


#' @title rbio_network.default
#'
#' @rdname rbio_network
#' @method rbio_network default
#' @param g  Input list cantaining DE dataframes for each comparison.
#' @param g_membership Input dataframe containing F stats.
#' @param export.name string. Optional prefix for export file name. Default is \code{NULL}.
#' @param colour_scheme string. The colour set from \code{\link{RColorBrewer}} function. Default is \code{"Accent"}.
#' @param initial_colour_number int. The number of starting colours to use for the clusters. See details. Default is \code{8}.
#' @param plot.title. string. <TBC: under construction>
#' @param plot.community_separated boolean. If to separate community vertices in graph. Default is \code{FALSE}.
#' @param plot.community_separated.dist Numeric. Set when \code{plot.community_separated}, distance between communities. Default is \code{3}.
#' @param plot.ellipse boolean. If to show member groups with ellipse. Default is \code{FALSE}.
#' @param plot.margins. numeric four-vector. Plot margins. See details.
#' @param plot.font.family string. The font family of the labels in the plot. Default is \code{"sans"}.
#' @param plot.highlight_membership boolean. If to highlight the membership, by making all edges transparent. Default is \code{TRUE}.
#' @param plot.layout_type string. Layout type, same as described in \code{igraph} pacakges.
#' @param plot.vertex.remove.zerodegree boolean. If to remove zero degree edges from the graph. Default is \code{TRUE}.
#' @param plot.vertex.size numeric vector. Vertex size.
#' @param plot.vertex.size.scale numeric length two-vector. A length two vector to scale vertex size.
#' @param plot.vertex.label string vector. Optional custom vertex label. Default is \code{NULL}, meaning V(g)$name.
#' @param plot.vertex.topvsize.filter numeric: 0-1. Set when \code{plot.vertex.label.topvsize = TRUE}, top percetage size to display the vertex labels. Default is \code{0.05}.
#' @param plot.vertex.color.highlighttopvsize boolean. When \code{plot.vertex.label.topvsize = TRUE}, if to make non-top vertices transparent and frameless. Default is \code{TRUE}.
#' @param plot.vertex.label.display boolean. If to display vertex labels. Default is \code{TRUE}.
#' @param plot.vertex.label.topvsize boolean. If to display labels with a threshold on vertex size. \code{default is TRUE}.
#' @param plot.vertex.label.size numeric vector. Vertex label size.
#' @param plot.vertex.label.colour string vector. Vertex label colour.
#' @param plot.vertex.label.dist numeric. Vertex label distance to vertex.
#' @param plot.vertex.label.loc radian numeric. The positioning of the vertex labels. See details.
#' @param plot.edge.type string. The type of edges to display in the network. Default is \code{"all"}.
#' @param plot.edge.filter numeric: 0-1. Percentage edges to keep. Default is \code{0.05}.
#' @param plot.edge.filter.direction string. Edge filter direction, top or bottom. Default is \code{"top"}.
#' @param plot.edge.weight numeric vector. Optional custom edge weight. Make sure it has the same length as the number of edges.
#' @param plot.edge.weight.scale numeric two-vector. A length two vector to scale edge weight.
#' @param plot.edge.arrow.mode Boolean. Edge arrow mode. Default is \code{FASLE}.
#' @param plot.edge.curved boolean. If to display curved edges. Default is \code{TRUE}.
#' @param plot.edge.color.highlighttopvsize boolean. When \code{plot.vertex.label.topvsize = TRUE}, if to make all edges transparent. Default is \code{TRUE}.
#' @param plot.width numeric. Plot width. Default is \code{7}.
#' @param plot.height numeric. Plot height. Default is \code{7}.
#' @param random_state integer. Random state for randomly generated network plot arrangement. Default is \code{1}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @details <TBC: under construction>
#'          For \code{colour_scheme}, use the following as a guide (name, maximum number of colours):
#'            Accent	8
#'            Dark2	8
#'            Paired	12
#'            Pastel1	9
#'            Pastel2	8
#'            Set1	9
#'            Set2	8
#'            Set3	12
#'          NOTE: The maximum number of colours does not reflect the number of clusters - it is simply what \code{\link{RColorBrewer}} requires.
#'                The recommended approach is to set \code{initial_colour_number} to this number.
#'
#'          The \code{plot.margins} follow the base R setting in \code{\link{par}} for the positioning:
#'          b: mar[1], l: mar[2], t: mar[3], r: mar[4]
#'
#'          For \code{plot.vertex.label.topvsize.filter}, the functional will apply the filter per group if g_membership is provided.
#'
#'          When \code{random_state = 0}, no random state is set.
#'
#'          For \code{plot.edge.type}:
#'
#'             \code{"all"}: display all the edges, i.e. both within memberships and across memberships
#'
#'             \code{"across"}: display only the cross-membership edges
#'
#'             \code{"within"}: display only the within membership edges
#'
#'          \code{plot.community_separated} and \code{plot.ellipse} are mutually exclusive.
#'
#'
#'          For \code{plot.vertex.label.loc}, from \code{igraph}:
#'          It defines the position of the vertex labels, relative to the center of the vertices.
#'          It is interpreted as an angle in radian, zero means ‘to the right’, and ‘pi’ means to the left, up is -pi/2 and down is pi/2.
#'          The default value is -pi/4.
#'
#' @import ggplot2
#' @import igraph
#' @importFrom scales alpha rescale
#' @importFrom grid grid.newpage grid.draw grid.grab
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
rbio_network.default <- function(g,
                                 export.name = NULL,
                                 g_membership = NULL,
                                 colour_scheme = c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3"),
                                 initial_colour_number = 8,
                                 plot.community_separated = FALSE, plot.community_separated.commu_dist = 3,
                                 plot.ellipse = FALSE,
                                 plot.title = "Network",
                                 plot.margins = c(5, 5, 5, 5),
                                 plot.font.family = "sans",
                                 plot.highlight_membership = TRUE,
                                 plot.layout_type = c("circular", "fr", "tree", "nicely", "sphere"),
                                 plot.vertex.remove.zerodegree = TRUE,
                                 plot.vertex.size = NULL,
                                 plot.vertex.size.scale = c(1, 4),
                                 plot.vertex.topvsize.filter = 0.05,
                                 plot.vertex.color.highlighttopvsize = TRUE,
                                 plot.vertex.label.display = TRUE,
                                 plot.vertex.label = NULL,
                                 plot.vertex.label.topvsize = TRUE,
                                 plot.vertex.label.size = 0.5,
                                 plot.vertex.label.color = "black",
                                 plot.vertex.label.dist = 0,
                                 plot.vertex.label.loc = -pi/4,
                                 plot.edge.type = c("all", "across", "within"),
                                 plot.edge.filter = 0.05,
                                 plot.edge.filter.direction = c("top", "bottom"),
                                 plot.edge.weight = NULL,
                                 plot.edge.weight.scale = c(1, 4),
                                 plot.edge.arrow.mode = FALSE,
                                 plot.edge.curved = FALSE,
                                 plot.edge.color.highlighttopvsize = TRUE,
                                 plot.height = 7, plot.width = 7,
                                 random_state = 1, verbose = TRUE){
  # - set random state -
  if (random_state != 0) {
    set.seed(random_state)
  }

  # - argument check -
  if (plot.ellipse && plot.community_separated) stop("plot.ellipse and plot.community_separated are mutually exclusive.")

  if (!is.null(g_membership)) {
    if (length(g_membership) != length(V(g))) stop("membership length not equal to number of vertecies. \n")
    if (is.null(names(g_membership))) {
      warning("membership has no vertex indices. Proceed with vertex names.\n")
      names(g_membership) <- V(g)$name
    }
  }

  plot.edge.type <- match.arg(plot.edge.type)
  colour_scheme <- match.arg(colour_scheme)
  if (!any(class(g) %in% "igraph")) stop("Input g needs to be an igraph.")
  if (is.null(export.name)){
    export.name <- deparse(substitute(object))
  } else {
    export.name <- export.name
  }
  plot.layout_type <- match.arg(plot.layout_type)
  plot.edge.filter.direction <- match.arg(plot.edge.filter.direction)

  # - initial network  -
  if (!is.null(g_membership)){
    n_colours <- length(unique(g_membership))
    get_colour_func <- colorRampPalette(brewer.pal(initial_colour_number, colour_scheme))
    colours <- get_colour_func(n_colours)
    membership_colours <- vector(length = length(g_membership))
    for (i in 1:n_colours){
      membership_colours[g_membership == i] <- colours[i]
    }
    names(membership_colours) <- names(g_membership)
    V(g)$color <- membership_colours
    V(g)$membership <- as_membership(g_membership)
  } else {
    V(g)$color <- rep("blue", times = length(V(g)))
  }

  # - edge weight and vertices size rescaling -
  # edge size
  if (!is.null(plot.edge.weight)) {
    if (length(plot.edge.weight) == length(E(g)$weight)) {
      E(g)$weight <- plot.edge.weight
    } else {
      warning("edge weight vector not equal length with the edges, proceeding with the internal edge weights.\n")
    }
  }

  # vertex size
  if (is.null(plot.vertex.size)) {
    g <- set_vertex_attr(g, name = "vsize", value = degree(g))
    vsize_is_degree <- TRUE
  } else {
    if (length(plot.vertex.size) == length(V(g))) {
      g <- set_vertex_attr(g, name = "vsize", value = plot.vertex.size)
      vsize_is_degree <- FALSE
    } else {
      warning("vertex size vector not equal length with the vertices, proceeding with the degree centrality.\n")
      g <- set_vertex_attr(g, name = "vsize", value = degree(g))
      vsize_is_degree <- TRUE
    }
  }

  # vertex labels


  if (is.null(plot.vertex.label)) { # text labels
    g <- set_vertex_attr(g, name = "vlabel", value = V(g)$name)
  } else if (length(plot.vertex.label) != length(V(g))) {
    warning("text labels not equal length with the vertices, proceeding with the internal vertex names.\n")
    g <- set_vertex_attr(g, name = "vlabel", value = V(g)$name)
  } else {
    g <- set_vertex_attr(g, name = "vlabel", value = plot.vertex.label)
  }

  # vertex label size
  if (length(plot.vertex.label.size) == 1) {
    g <- set_vertex_attr(g, name = "vlabelsize", value = rep(plot.vertex.label.size, times = length(V(g))))
  } else if (length(plot.vertex.label.size) != length(V(g))) {
    warning("vertex label size vector not equal length with the vertices, proceeding with the first value in the vector.\n")
    g <- set_vertex_attr(g, name = "vlabelsize", value = rep(plot.vertex.label.size[1], times = length(V(g))))
  } else {
    g <- set_vertex_attr(g, name = "vlabelsize", value = plot.vertex.label.size)
  }

  # vertex frame (outline) colour
  g <- set_vertex_attr(g, name = "vframecolour", value = rep("black", times = length(V(g))))

  # - filter edges -
  # based on membership: experimental
  if (plot.edge.type == "across") {
    for (i in 1:length(unique(g_membership))) {
      tmp_member <- unique(g_membership)[i]
      g <- delete.edges(g, E(g)[V(g)[membership == tmp_member] %--% V(g)[membership == tmp_member]])
    }
  } else if (plot.edge.type == "within") {
    comb_m <- combn(unique(g_membership), 2, simplify = FALSE)  # find all cross-membership combination in a list. 2 means combination of 2 elements
    for (i in 1:length(comb_m)) {  # i the combination index in the comb_m list
      temp_idx <- comb_m[[i]]
      g <- delete.edges(g, E(g)[V(g)[membership == temp_idx[[1]]] %--% V(g)[membership == temp_idx[[2]]]])
    }
  }

  # based on edge weight
  if (plot.edge.filter.direction == "top") {
    g <- delete_edges(g, E(g)[E(g)$weight < quantile(E(g)$weight, p = 1-plot.edge.filter)])
  } else {
    g <- delete_edges(g, E(g)[E(g)$weight > quantile(E(g)$weight, p = plot.edge.filter)])
  }

  edge_df <- as.data.frame(get.edgelist(g))
  if (!is.null(g_membership)) {
    if (plot.highlight_membership) {
      E(g)$color <- foreach(i = seq(nrow(edge_df)), .combine = "c") %do% {
        if (g_membership[edge_df[i, 1]] == g_membership[edge_df[i, 2]]){
          colours[g_membership[edge_df[i, 1]]]
        } else {
          scales::alpha("#EBECF0", alpha = 0.5)
          # "#EBECF0"
        }
      }
    } else {
      E(g)$color <- foreach(i = seq(nrow(edge_df)), .combine = "c") %do% {
        if (g_membership[edge_df[i, 1]] == g_membership[edge_df[i, 2]]){
          scales::alpha(colours[g_membership[edge_df[i, 1]]], alpha = 0.2)
          # colours[membership[edge_df[i, 1]]]
        } else {
          # scales::alpha("#EBECF0", alpha=0.2)
          "#3C3C3C"
        }
      }
    }
  } else {
    E(g)$color <-rep("#EBECF0", times = nrow(edge_df))
  }

  # - filer vertices -
  if (plot.vertex.remove.zerodegree) {
    g <- delete.vertices(g, degree(g) == 0)
  }

  # - finalize network -
  # refersh vsize if used degrees
  if (vsize_is_degree) {
    V(g)$vsize <- degree(g)
  }

  # finalize labels (if to selectively display labels and colours)
  if (is.null(g_membership)) {
    to_remove <- V(g)$vsize < quantile(V(g)$vsize, p = 1-plot.vertex.topvsize.filter)
    # V(g)$vlabel[to_remove] <- ""

    if (plot.vertex.label.topvsize) { # only to display label for top size vertices
      V(g)$vlabel[to_remove] <- ""
    }

    if (plot.vertex.color.highlighttopvsize) {  # only to display colour for top size vertices
      V(g)$color[to_remove] <- alpha(V(g)$color[to_remove], alpha = 0.2)
      V(g)$vframecolour[to_remove] <- "NA"
    }
  } else {
    for (i in 1:length(unique(V(g)$membership))) {
      is_member <- V(g)$membership == unique(V(g)$membership)[i]
      to_remove <- V(g)$vsize[is_member] < quantile(V(g)$vsize[is_member], p = 1-plot.vertex.topvsize.filter)

      # V(g)$vlabel[is_member][to_remove] <- ""

      if (plot.vertex.label.topvsize) {  # only display label for top size vertices
        V(g)$vlabel[is_member][to_remove] <- ""
      }

      if (plot.vertex.color.highlighttopvsize) {
        V(g)$color[is_member][to_remove] <- alpha(V(g)$color[is_member][to_remove], alpha = 0.2)
        V(g)$vframecolour[is_member][to_remove] <- "NA"
      }
    }
  }

  if (plot.edge.color.highlighttopvsize) {
    E(g)$color <- alpha(E(g)$color , alpha = 0.2)
  }

  # rescale sizes
  edgeweights <- scales::rescale(E(g)$weight, to = plot.edge.weight.scale)
  V(g)$vsize <- scales::rescale(V(g)$vsize, to = plot.vertex.size.scale)

  # layout
  if (plot.layout_type == "circular") {
    g_layout <- layout_in_circle(g)
  } else if (plot.layout_type == "fr") {
    g_layout <- layout_with_fr(g, weights=E(g)$weight)
  } else if (plot.layout_type == "nicely") {
    g_layout <- layout_nicely(g)
  } else if (plot.layout_type == "tree") {
    g_layout <- layout_as_tree(g)
  } else if (plot.layout_type == "sphere") {
    g_layout <- layout_on_sphere(g)
  }

  # - plot and export -
  # finalize vlabel
  if (!plot.vertex.label.display) {
    V(g)$vlabel <- NA
  }

  # plot
  g.cluster <- make_clusters(g, membership = V(g)$membership)
  pdf(file = paste0(export.name, "_network.pdf"),
      width = plot.width, height = plot.height)
  par(mar = plot.margins)
  # tiff(filename = paste0(export.name, "_network.tiff"),
  #      width = plot.width, height = plot.height, units = "mm", pointsize = 12,
  #      compression = "lzw",
  #      bg = "white", res = 600,
  #      type = "quartz")
  if (plot.ellipse) {
    if (plot.layout_type == "circular") {
      plot(
        g.cluster, g,
        col = V(g)$color,
        layout = g_layout,
        vertex.size = V(g)$vsize,
        vertex.label = NA,
        vertex.frame.color = V(g)$vframecolour,
        asp = TRUE,
        edge.width = edgeweights,
        edge.arrow.mode = plot.edge.arrow.mode,
        edge.curved = plot.edge.curved,
        main = plot.title)

      if (!plot.vertex.label.display) {
        circle_text_func(g = g, circ_layout = g_layout,
                         text.label = V(g)$vlabel,
                         text.size = V(g)$vlabelsize,
                         text.colour = plot.vertex.label.color, text.distance = plot.vertex.label.dist,
                         family = plot.font.family)
      }
    } else {
      g_layout <- layout_with_fr(g, weights=E(g)$weight)
      plot(
        g.cluster, g,
        col = V(g)$color,
        layout = g_layout,
        vertex.size = V(g)$vsize,
        vertex.label = V(g)$vlabel,
        vertex.label.family = plot.font.family,
        vertex.label.cex = V(g)$vlabelsize,
        vertex.label.dist = plot.vertex.label.dist,
        vertex.label.color = plot.vertex.label.color,
        vertex.label.cex = plot.vertex.label.size,
        vertex.label.degree = plot.vertex.label.loc,
        vertex.frame.color = V(g)$vframecolour,
        asp = TRUE,
        edge.curved = plot.edge.curved,
        edge.width = edgeweights,
        edge.arrow.mode = plot.edge.arrow.mode,
        main = plot.title)
    }
  } else if (plot.community_separated) {
    commu_weight <- seq(length(unique(V(g)$membership)))
    subg_list <- vector(mode = "list", length = length(unique(V(g)$membership)))
    subg_list[] <- foreach(i = 1:length(unique(V(g)$membership))) %do% {
      l_idx <- as_ids(V(g)[V(g)$membership == i])
      sub_g <- induced_subgraph(g, l_idx)
      if (plot.layout_type == "circular") {
        sub_l <- layout_in_circle(sub_g)
      } else if (plot.layout_type == "fr") {
        sub_l <- layout_with_fr(sub_g, weights=E(sub_g)$weight)
      } else if (plot.layout_type == "nicely") {
        sub_l <- layout_nicely(sub_g)
      } else if (plot.layout_type == "tree") {
        sub_l <- layout_as_tree(sub_g)
      } else if (plot.layout_type == "sphere") {
        sub_l <- layout_on_sphere(sub_g)
      }
      sub_l_x <- sub_l[, 1] + commu_weight[i] * plot.community_separated.commu_dist
      sub_l_y <- sub_l[, 2] + commu_weight[i] * plot.community_separated.commu_dist
      sub_l <- cbind(sub_l_x, sub_l_y)
      outlist <- list(sub_g = sub_g, sub_l = sub_l)
      outlist
    }
    names(subg_list) <- c(paste0("membership_", unique(V(g)$membership)))
    sub_layout <- foreach(i = 1:length(subg_list), .combine = "rbind") %do% {
      out <- subg_list[[i]][[2]]
      out
    }
    plot(
      g,
      layout = sub_layout,
      vertex.size = V(g)$vsize,
      vertex.label = V(g)$vlabel,
      vertex.label.family = plot.font.family,
      vertex.label.cex = V(g)$vlabelsize,
      vertex.label.dist = plot.vertex.label.dist,
      vertex.label.color = plot.vertex.label.color,
      vertex.label.cex = plot.vertex.label.size,
      vertex.frame.color = V(g)$vframecolour,
      vertex.label.degree = plot.vertex.label.loc,
      vertex.frame.color = V(g)$vframecolour,
      asp = TRUE,
      edge.width = edgeweights,
      edge.arrow.mode = plot.edge.arrow.mode,
      edge.curved = plot.edge.curved,
      main = plot.title)
  } else {
    if (plot.layout_type == "circular") {
      plot(
        g,
        layout = g_layout,
        vertex.size = V(g)$vsize,
        vertex.label = NA,
        vertex.frame.color = V(g)$vframecolour,
        asp = TRUE,
        edge.width = edgeweights,
        edge.arrow.mode = plot.edge.arrow.mode,
        edge.curved = plot.edge.curved,
        main = plot.title)

      if (!plot.vertex.label.display) {
        circle_text_func(g = g, circ_layout = g_layout,
                         text.label = V(g)$vlabel,
                         text.size = V(g)$vlabelsize,
                         text.colour = plot.vertex.label.color, text.distance = plot.vertex.label.dist,
                         family = plot.font.family)
      }
    } else {
      plot(
        g,
        layout = g_layout,
        vertex.size = V(g)$vsize,
        vertex.label = V(g)$vlabel,
        vertex.label.family = plot.font.family,
        vertex.label.cex = V(g)$vlabelsize,
        vertex.label.dist = plot.vertex.label.dist,
        vertex.label.color = plot.vertex.label.color,
        vertex.label.cex = plot.vertex.label.size,
        vertex.frame.color = V(g)$vframecolour,
        vertex.label.degree = plot.vertex.label.loc,
        asp = TRUE,
        edge.width = edgeweights,
        edge.arrow.mode = plot.edge.arrow.mode,
        edge.curved = plot.edge.curved,
        main = plot.title)
    }
  }
  # # export
  # grid.echo()
  # p <- grid.grab()
  # ggsave(filename = paste0(export.name, "_network3.pdf"), plot = p,
  #        width = plot.width, height = plot.height, units = "mm", dpi = 600)
  dev.off()
}

