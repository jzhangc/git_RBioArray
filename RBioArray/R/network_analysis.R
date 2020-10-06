#' @title rbio_tom
#'
#' @description TOM (topological overlap measure) analysis.
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
#' @param plot.margins Numeric vector. Plot margins, unit is "cm", order: \code{t, r, b, l}. Default is \code{c(2, 2, 2, 2)}.
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
#' @importFrom WGCNA TOMdist
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
rbio_tom <- function(mtx,
                     diag = FALSE,
                     power = 6, tom_type = c("unsigned", "signed"), ...,
                     hclust.method = c("complete", "ward.D", "ward.D2", "single",  "average", "mcquitty", "median", "centroid"),
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
  tom_dist <- TOMdist(adjmat, TOMType = tom_type, verbose = verbose, ...)  # matrix, array class, here we use dist
  rownames(tom_dist) <- rownames(adjmat)
  colnames(tom_dist) <- colnames(adjmat)

  # - TOM hclust -
  tom_dist <- as.dist(tom_dist)  # convert to an R distance object
  tom_dist_hclust <- hclust(tom_dist, method = hclust.method)

  # tom_membership and tom similarity
  tom_membership <- stats::cutree(tom_dist_hclust, h = h, k = k)
  # ßmembersihp_for_dendro <- tom_membership
  tom_similarity <- 1 - tom_dist # edge always uses similarity
  g_adjmat <- as.matrix(tom_similarity)
  g_adjmat <- g_adjmat[order(tom_membership), order(tom_membership)]  # reorder it
  tom_membership <- tom_membership[order(tom_membership)] # update tom_membership

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
    if (verbose) cat("Done!")
  }

  # - igraph -
  out <- list(g = graph.adjacency(
    g_adjmat,
    mode = "undirected",
    weighted = TRUE,
    diag = diag
  ),
  tom_dist_hclust = tom_dist_hclust,
  k = k,
  h = h,
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
#' @param text.label string vector. TBC
#' @param text.size numeric or numeric vector. The text size for the labels. Default is \code{0.8}.
#' @param text.distance numeric. The distance multiplier between label and nodes. Default is \code{1.5}.
#' @param text.colour string. Label colour. Default is \code{"black"}.
#' @return Added text on the igraph plot.
#' @detals
#'         The \code{text.size} argument also accepts a vector of sizes with a length equal to the number of vertices.
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

  # - create vector of angles for text based on number of nodes  -
  # (flipping the orientation of the words half way around so none appear upside down)
  angle = ifelse(atan(-(circ_layout[,1]/circ_layout[,2]))*(180/pi) < 0,  90 + atan(-(circ_layout[,1]/circ_layout[,2]))*(180/pi), 270 + atan(-circ_layout[,1]/circ_layout[,2])*(180/pi))

  #Apply the text labels with a loop with angle as srt
  for (i in 1:length(x)) {
    text(x=x[i], y=y[i], labels=V(g)$name[i], adj = NULL, pos = NULL, cex = tSize[i], col = text.colour, srt = angle[i], xpd = T, ...)
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
#' @param plot.margins. numeric four-vector. <TBC: under construction>
#' @param plot.font.family string. The font family of the labels in the plot. Default is \code{"sans"}.
#' @param plot.highlight_membership boolean. <TBC: under construction>
#' @param plot.layout_type string. <TBC: under construction>
#' @param plot.vertex.size numeric vector. <TBC: under construction>
#' @param plot.vertex.size.scale numeric two-vector. <TBC: under construction>
#' @param plot.vertex.label.size numeric vector. <TBC: under construction>
#' @param plot.vertex.label.colour string vector. <TBC: under construction>
#' @param plot.vertex.label.dist numeric. <TBC: under construction>
#' @param plot.edge.filter numeric: 0~1. <TBC: under construction>
#' @param plot.edge.weight numeric vector. <TBC: under construction>
#' @param plot.edge.weight.scale numeric two-vector. <TBC: under construction>
#' @param plot.edge.arrow.mode boolean. <TBC: under construction>
#' @param plot.edge.curved boolean. <TBC: under construction>
#' @param plot.ellipse boolean. <TBC: under construction>
#' @param plot.width numeric. <TBC: under construction>
#' @param plot.height numeric. <TBC: under construction>
#' @param random_state integer. <TBC: under construction>
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
                                 plot.title = "Network",
                                 plot.margins = c(5, 5, 5, 5),
                                 plot.font.family = "sans",
                                 plot.highlight_membership = TRUE,
                                 plot.layout_type = c("circular", "fr", "tree", "nicely", "sphere"),
                                 plot.vertex.size = NULL,
                                 plot.vertex.size.scale = c(1, 4),
                                 plot.vertex.label = NULL,
                                 plot.vertex.label.size = 0.5,
                                 plot.vertex.label.color = "black",
                                 plot.vertex.label.dist = 0,
                                 plot.edge.filter = 0.95,
                                 plot.edge.weight = NULL,
                                 plot.edge.weight.scale = c(1, 4),
                                 plot.edge.arrow.mode = FALSE,
                                 plot.edge.curved = FALSE,
                                 plot.ellipse = FALSE,
                                 plot.height = 150, plot.width = 150,
                                 random_state = 1, verbose = TRUE){
  # - set random state -
  set.seed(random_state)

  # - argument check -
  g_membership <- g_membership
  colour_scheme <- match.arg(colour_scheme, c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3"))
  if (!any(class(g) %in% "igraph")) stop("Input g needs to be an igraph.")
  if (is.null(export.name)){
    export.name <- deparse(substitute(object))
  } else {
    export.name <- export.name
  }
  plot.layout_type <- match.arg(plot.layout_type, c("circular", "fr", "tree", "nicely", "sphere"))

  # - initial network  -
  if (!is.null(g_membership) %% length(g_membership) != length(V(g))){
    warning("membership length not equal to number of vertecies. \n")
  } else {
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
  }

  # - edge weight and vertices size rescaling -
  # edge size
  if (is.null(plot.edge.weight)) {
    edgeweights <- scales::rescale(E(g)$weight, to = plot.edge.weight.scale)
  } else if (length(plot.edge.weight) != length(E(g)$weight)) {
    warning("edge weight vector not equal length with the edges, proceeding with the internal edge weights.\n")
    edgeweights <- scales::rescale(E(g)$weight, to = plot.edge.weight.scale)
  } else {
    edgeweights <- scales::rescale(plot.edge.weight, to = plot.edge.weight.scale)
  }

  # vertex size
  if (is.null(plot.vertex.size)) {
    vSizes <- scales::rescale(degree(g), to = plot.vertex.size.scale)
  } else if (length(plot.vertex.size) != length(V(g))) {
    warning("vertex size vector not equal length with the vertices, proceeding with the degree centrality.\n")
    vSizes <- scales::rescale(degree(g), to = plot.vertex.size.scale)
  } else {
    vSizes <- scales::rescale(plot.vertex.size, to = plot.vertex.size.scale)
  }

  # vertex labels
  if (is.null(plot.vertex.label)) { # text labels
    vLabel <- V(g)$name
  } else if (length(text.labels) != length(V(g))) {
    warning("text labels not equal length with the vertices, proceeding with the internal vertex names.\n")
    vLabel <- V(g)$name
  } else {
    vLabel <- text.labels
  }

  # vertex label size
  if (length(plot.vertex.label.size) == 1) {
    vLabelSize <- plot.vertex.label.size
  } else if (length(plot.vertex.label.size) != length(V(g))) {
    warning("vertex label size vector not equal length with the vertices, proceeding with the first value in the vector.\n")
    vLabelSize <- plot.vertex.label.size[1]
  }

  # - filter edges -
  g <- delete_edges(g, E(g)[E(g)$weight < quantile(E(g)$weight, p = plot.edge.filter)])
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
  }

  # - filter vertices -
  g <- delete.vertices(g, degree(g) == 0)

  # - plot and export -
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

  # plot and export
  grid.newpage()
  if (plot.ellipse) {
    g.cluster <- make_clusters(g, membership = V(g)$membership)
    if (plot.layout_type == "circular") {
      par(mar=plot.margins)
      plot(
        g.cluster, g,
        layout = g_layout,
        vertex.size = vSizes,
        vertex.label = NA,
        asp = FALSE,
        edge.width = edgeweights,
        edge.arrow.mode = plot.edge.arrow.mode,
        edge.curved = plot.edge.curved,
        main = plot.title)
      circle_text_func(g = g, circ_layout = g_layout,
                       text.label = vLabel,
                       text.size = plot.vertex.label.size,
                       text.colour = vertex.label.color, text.distance = plot.vertex.label.dist,
                       family = plot.font.family)
    } else {
      g_layout <- layout_with_fr(g, weights=E(g)$weight)
      par(mar = plot.margins)
      plot(
        g.cluster, g,
        layout = g_layout,
        vertex.size = vSizes,
        vertex.label = vLabel,
        vertex.label.family = plot.font.family,
        vertex.label.cex = vLabelSize,
        vertex.label.dist = plot.vertex.label.dist,
        vertex.label.color = plot.vertex.label.color,
        vertex.label.cex = plot.vertex.label.size,
        asp = FALSE,
        edge.curved = plot.edge.curved,
        edge.width = edgeweights,
        edge.arrow.mode = plot.edge.arrow.mode,
        main = plot.title)
    }
  } else {
    if (plot.layout_type == "circular") {
      par(mar = plot.margins)
      plot(
        g,
        layout = g_layout,
        vertex.size = vSizes,
        vertex.label = NA,
        asp = FALSE,
        edge.width = edgeweights,
        edge.arrow.mode = plot.edge.arrow.mode,
        edge.curved = plot.edge.curved,
        main = plot.title)
      circle_text_func(g = g, circ_layout = g_layout,
                       text.label = vLabel,
                       text.size = plot.vertex.label.size,
                       text.colour = plot.vertex.label.color, text.distance = 1.21,
                       family = plot.font.family)
    } else {
      par(mar = plot.margins)
      plot(
        g,
        layout = g_layout,
        vertex.size = vSizes,
        vertex.label = vLabel,
        vertex.label.family = plot.font.family,
        vertex.label.cex = vLabelSize,
        vertex.label.dist = plot.vertex.label.dist,
        vertex.label.color = plot.vertex.label.color,
        vertex.label.cex = plot.vertex.label.size,
        asp = FALSE,
        edge.width = edgeweights,
        edge.arrow.mode = plot.edge.arrow.mode,
        edge.curved = plot.edge.curved,
        main = plot.title)
    }
  }
  grid.echo()
  p <- grid.grab()
  ggsave(filename = paste0(export.name, "_network3.pdf"), plot = p,
         width = plot.width, height = plot.height, units = "mm", dpi = 600)
}
