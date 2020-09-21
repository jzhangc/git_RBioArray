#' @title rbio_unsupervised_hcluster
#'
#' @description Generic unsupersived hierarchical clustering function.
#' @param object Input object in either \code{rbioseq_de} or \code{rbioarray_flist} class.
#' @param ... Additional arguments for corresponding S3 class methods.
#' @details Since the function depends on the \code{heatmap.2} function from \code{gplots} package,
#'          arguments can be passed directly, as seen in the examples below.
#'
#' @return A pdf file containing a heatmap for unsupervised hierarchical clustering analysis.
#' @examples
#'
#' \dontrun{
#' # rbioarray_flist class input
#' rbio_unsupervised_hcluster(object = flist, rm.control = TRUE, n = 500, gene_symbol.only = TRUE,
#'                            sample_id.var.name = "SampleName", srtCol = 90, offsetCol = 0, labRow = FALSE,
#'                            key.title = "", cexCol = 0.7, trace = "none",
#'                            keysize = 1.5, key.xlab = "Normalized expression level", key.ylab = "Probe count")
#'
#' # rbioseq_de class input
#' rbio_unsupervised_hcluster(object = mrna_de, n = "all",
#'                            sample_id.var.name = "condition", srtCol = 90, offsetCol = 0,
#'                            key.title = "", cexCol = 0.7, trace = "none",
#'                            keysize = 1.5, key.xlab = "CLR transformed counts", key.ylab = "mRNA count")
#'
#'
#'
#' }
#'
#' @export
rbio_unsupervised_hcluster <- function(object, ...){
  ## check arguments
  if (!class(object) %in% c("rbioarray_flist", "rbioseq_de")) stop("The input object needs to be either \"rbioarray_flist\" or \"rbioseq_de\" class object.")

  ## use methods
  UseMethod("rbio_unsupervised_hcluster", object)
}


#' @title rbio_unsupervised_hcluster.rbioarray_flist
#'
#' @rdname rbio_unsupervised_hcluster
#' @method rbio_unsupervised_hcluster rbioarray_flist
#' @param object Input object in \code{rbioarray_flist} class.
#' @param sample_id.var.name Variable name for sample identification, typically from \code{object$target}.
#' @param export.name Optional user defined export name prefix. Default is \code{NULL}.
#' @param ... Additional arguments for the default method.
#' @return A pdf file containing a heatmap for unsupervised hierarchical clustering analysis.
#' @export
rbio_unsupervised_hcluster.rbioarray_flist <- function(object, sample_id.var.name = NULL, export.name = NULL, ...){
  ## check arguments
  if (!is.null(sample_id.var.name)){
    if (!sample_id.var.name %in% names(object$targets)) {
      cat("The sample_id.var.name not found in targets element of the input object. Proceed without using it.\n")
      sample_id.vector <- seq(ncol(object$E))
    } else {
      sample_id.vector <- object$targets[, sample_id.var.name]
    }
  } else {
    sample_id.vector <- seq(ncol(object$E))
  }

  if (is.null(export.name)){
    export.name <- deparse(substitute(object))
  } else {
    export.name <- export.name
  }

  ## use methods
  rbio_unsupervised_hcluster.default(E = object$E, genes = object$genes,
                                     input.genes_annotation.control_type = object$genes_annotation.control_type,
                                     input.genes_annotation.gene_symbol.var_name = object$genes_annotation.gene_symbol.var_name,
                                     input.genes_annotation.gene_id.var_name = object$genes_annotation.gene_id.var_name,
                                     input.sample_groups = object$sample_groups,
                                     sample_id.vector = sample_id.vector, export.name = export.name, ...)
}


#' @title rbio_unsupervised_hcluster.rbioseq_de
#'
#' @rdname rbio_unsupervised_hcluster
#' @method rbio_unsupervised_hcluster rbioseq_de
#' @param object Input object in \code{rbioseq_de} class.
#' @param sample_id.var.name Variable name for sample identification, typically from \code{object$target}.
#' @param export.name Optional user defined export name prefix. Default is \code{NULL}.
#' @param ... Additional arguments for the default method.
#' @details The function uses filtered count data, as opposed to normalized data.
#'          Due to the compositional nature of NGS data, the count data is transformed using CLR method prior to clustering.
#' @return A pdf file containing a heatmap for unsupervised hierarchical clustering analysis.
#' @export
rbio_unsupervised_hcluster.rbioseq_de <- function(object, sample_id.var.name = NULL, export.name = NULL, ...){
  ## check arguments
  if (!is.null(sample_id.var.name)){
    if (!sample_id.var.name %in% names(object$targets)) {
      cat("The sample_id.var.name not found in targets element of the input object. Proceed without using it.\n")
      sample_id.vector <- seq(ncol(object$filter_results$filtered_counts$counts))
    } else {
      sample_id.vector <- object$targets[, sample_id.var.name]
    }
  } else {
    sample_id.vector <- seq(ncol(object$filter_results$filtered_counts$counts))
  }

  if (is.null(export.name)){
    export.name <- deparse(substitute(object))
  } else {
    export.name <- export.name
  }

  ## transform
  cat("CLR transformation of filtered RNAseq count data...")
  E_transfo <- rbioseq_clr_ilr_transfo(object$filter_results$filtered_counts$counts, offset = 1, mode = "clr")  # clr tranformation
  cat("Done!\n")

  ## use methods
  rbio_unsupervised_hcluster.default(E = E_transfo, genes = object$filter_results$filtered_counts$genes,
                                     input.genes_annotation.control_type = NULL,
                                     input.genes_annotation.gene_symbol.var_name = object$genes_annotation.gene_symbol.var_name,
                                     input.genes_annotation.gene_id.var_name = object$genes_annotation.gene_id.var_name,
                                     input.sample_groups = object$sample_groups,
                                     sample_id.vector = sample_id.vector, export.name = export.name, rm.control = FALSE, ...)
}


#' @title rbio_unsupervised_hcluster.default
#'
#' @rdname rbio_unsupervised_hcluster
#' @method rbio_unsupervised_hcluster default
#' @description Default unsupersived hierarchical clustering function.
#' @param E Expression or count matrix, with rows for genes/probes/genomic features, columns for RNA samples.
#' @param genes Annotation data frame for genes/probes/genomic features.
#' @param input.sample_groups Input \code{factor} object for sample groupping labels.
#' @param n Number of genes/probes/genomic features to be clustered, numeric input or \code{"all"}. Default is \code{"all"}.
#' @param rm.control Whether to remove control probes (Agilent platform) or not. Default is \code{TRUE}.
#' @param input.genes_annotation.control_type Only set when \code{rm.control = TRUE}, input control type variable annotation list.
#' @param gene_symbol.only Whether or not to remove probes without gene symbol. Default is \code{FALSE}.
#' @param input.genes_annotation.gene_symbol.var_name Only set when \code{gene_symbol.only = TRUE}, variable name for gene symbol column in \code{genes} data frame.
#' @param input.genes_annotation.gene_id.var_name Only set when \code{gene_symbol.only = TRUE}, variable name for gene id column in \code{genes} data frame.
#' @param sample_id.vector A \code{vector} containing names to display for each heatmap column. Default is \code{NULL} and the function will use the column name from the input.
#' @param distance Distance calculation method. Default is \code{"euclidean"}. See \code{\link{dist}} for more.
#' @param clust Clustering method. Default is \code{"complete"}. See \code{\link{hclust}} for more.
#' @param col.colour Column group colour. Default is \code{"Paired"}. See \code{RColorBrewer} package for more.
#' @param map.colour Heat map colour. Default is \code{"PRGn"}. See \code{RColorBrewer} package for more.
#' @param n.map.colour Number of colours displayed. Default is \code{11}. See \code{RColorBrewer} package for more.
#' @param ... Additional arguments for \code{heatmap.2} function from \code{gplots} package.
#' @param export.name File name for the export \code{pdf} plot file.
#' @param plot.width Width of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param plot.height Height of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return A heatmap based on hierarchical clustering analysis in \code{pdf} format.
#'         The function also outputs the column and row distance and cluster objects in a list.
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @export
rbio_unsupervised_hcluster.default <- function(E, genes, input.sample_groups, n = "all",
                                               rm.control = FALSE, input.genes_annotation.control_type,
                                               gene_symbol.only = FALSE,
                                               input.genes_annotation.gene_symbol.var_name = NULL,
                                               input.genes_annotation.gene_id.var_name = NULL,
                                               sample_id.vector = NULL,
                                               distance = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                                               clust = c("complete", "ward.D", "ward.D2", "single",  "average", "mcquitty", "median", "centroid"),
                                               col.colour = "Paired", map.colour = "PRGn", n.map.colour = 11, ...,
                                               export.name = NULL, plot.width = 7, plot.height = 7,
                                               verbose = TRUE){
  ## check arguments
  distance <- match.arg(tolower(distance), c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
  clust <- match.arg(clust, c("complete", "ward.D", "ward.D2", "single",  "average", "mcquitty", "median", "centroid"))

  if (n != "all" && n %% 1 != 0) stop("Argument n needs to be either \"all\" or an integer number.")
  if (n.map.colour %% 1 != 0) stop("Argument n.map.colour needs to be an integer number.")
  if (rm.control && is.null(input.genes_annotation.control_type)) {
    cat("Argument input.genes_annotation.control_type is NULL when rm.control = TRUE, automatically set rm.control = FALSE.\n\n")
    rm.control <- FALSE
  }
  if (gene_symbol.only && ! input.genes_annotation.gene_symbol.var_name %in% names(genes)) {
    cat("Argument input.genes_annotation.gene_symbol.var_name not found in genes data frame when gene_symbol.only = TRUE, automatically set gene_symbol.only = FALSE.\n\n")
    gene_symbol.only <- FALSE
  }
  if (missing(export.name) || is.null(export.name)) stop("Please set value for argument export.name.")

  ## variables
  if(length(levels(input.sample_groups)) <= 19) {
    colGroup <- length(levels(input.sample_groups))
  } else {
    cat("The sample groups exceed the maximum allowed colour group number (19). Proceed with 19.\n\n")
    colGroup <- 19
  }
  row.lab.var_name <- input.genes_annotation.gene_id.var_name

  ## set up dis and cluster functions
  distfunc <- function(x)dist(x, method = distance, upper = TRUE, diag = TRUE)
  clustfunc <- function(x)hclust(x, method = clust)

  ## prepare dfm for clustering
  dfm <- data.frame(genes, E, check.names = FALSE)

  if (rm.control){ # remove control
    dfm <- dfm[dfm[, input.genes_annotation.control_type$control_type.var_name] == input.genes_annotation.control_type$exp_type.value, ]
  }
  if (gene_symbol.only) {
    dfm <- dfm[complete.cases(dfm[, input.genes_annotation.gene_symbol.var_name]),]
    row.lab.var_name <- input.genes_annotation.gene_symbol.var_name
  }
  if (n != "all"){ # subset
    dfm <- dfm[1:n, ]
  }

  mtx <- as.matrix(dfm[, !names(dfm) %in% names(genes)])
  colnames(mtx) <- sample_id.vector
  rownames(mtx) <- dfm[, row.lab.var_name]
  row.lab <- dfm[, row.lab.var_name]

  ## calculate and output distance, cluster and set ColSideColors
  # column cluster (samples)
  col_dist <- distfunc(t(mtx))
  col_cluster <- clustfunc(col_dist)
  colG <- cutree(col_cluster, colGroup) # column group
  colC <- brewer.pal(ifelse(colGroup < 3, 3, colGroup), col.colour) # column colour
  col_cluster_list <- list(col_dist = col_dist, col_hclust = col_cluster)

  # row cluster (genes/features)
  row_dist <- distfunc(mtx)
  row_cluster <- clustfunc(row_dist)
  row_cluster_list <- list(row_dist = row_dist, row_hclust = row_cluster)

  # output the distance and cluster
  out <- list(col_dist_hclust = col_cluster_list, row__dist_hclust = row_cluster_list,
              distance_method = distance, cluster_method = clust)
  assign(paste0(export.name, "_dist_clust"), out, envir = .GlobalEnv)

  ## heatmap
  # draw heatmap
  if (verbose) cat(paste0("Unsupervised hierarchical clustering heatmap saved to: ", export.name, "_unsuper_heatmap.pdf..."))
  pdf(file = paste0(export.name, "_unsuper_heatmap.pdf"), width = plot.width, height = plot.height)
  heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
            labRow = row.lab,
            col = brewer.pal(n.map.colour, map.colour), ColSideColors = colC[colG], ...)
  if (verbose) cat("Done!\n")
  dev.off()
}


#' @title rbio_supervised_hcluster
#'
#' @description Supersived hierarchical clustering function.
#' @param object Input \code{sig} class object.
#' @param gene_symbol.only Whether or not to remove probes without gene symbol. Default is \code{FALSE}.
#' @param sample_id.vector A \code{vector} containing names to display for each heatmap column. Default is \code{NULL} and the function will use the column name from the input.
#' @param distance Distance calculation method. Default is \code{"euclidean"}. See \code{\link{dist}} for more.
#' @param clust Clustering method. Default is \code{"complete"}. See \code{\link{hclust}} for more.
#' @param col.colour Column group colour. Default is \code{"Paired"}. See \code{RColorBrewer} package for more.
#' @param map.colour Heat map colour. Default is \code{"PRGn"}. See \code{RColorBrewer} package for more.
#' @param n.map.colour Number of colours displayed. Default is \code{11}. See \code{RColorBrewer} package for more.
#' @param ... Additional arguments for \code{heatmap.2} function from \code{gplots} package.
#' @param plot.width Width of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param plot.height Height of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @details Unlike the unsupervised veresion, the sig data hcluster uses normalized expression data for both RNAseq and microaray.
#'          The column colour group is usually 2, since the function only outputs
#'
#'          NOTE: this function only outputs the pair-wise comparison clusters.
#' @return A heatmap based on hierarchical clustering analysis in \code{pdf} format.
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @export
rbio_supervised_hcluster <- function(object,
                                     gene_symbol.only = FALSE,
                                     sample_id.var.name = NULL,
                                     distance = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                                     clust = c("complete", "ward.D", "ward.D2", "single",  "average", "mcquitty", "median", "centroid"),
                                     col.colour = "Paired", map.colour = "PRGn", n.map.colour = 11, ...,
                                     plot.width = 7, plot.height = 7,
                                     verbose = TRUE){
  ## check arguments
  distance <- match.arg(tolower(distance), c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
  clust <- match.arg(clust, c("complete", "ward.D", "ward.D2", "single",  "average", "mcquitty", "median", "centroid"))

  if (any(class(object) != "sig")) stop("The input object has to be a \"sig\" class.")
  if (!is.null(sample_id.var.name)){
    if (!sample_id.var.name %in% names(object$input_data$targets)) {
      cat("The sample_id.var.name not found in targets element of the input object. Proceed without using it.\n")
      sample_id.vector <- seq(ncol(object$input_data$E))
    } else {
      sample_id.vector <- object$input_data$targets[, sample_id.var.name]
    }
  } else {
    sample_id.vector <- seq(ncol(object$input_data$E))
  }
  if (n.map.colour %% 1 != 0) stop("Argument n.map.colour needs to be an integer number.")
  # if (!distance %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
  #   stop("Argument distance needs to be one of \"euclidean\", \"maximum\", \"manhattan\", \"canberra\", \"binary\"or \"minkowski\".")
  # if (!clust %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))
  #   stop("Argument clust needs to be one of \"ward.D\", \"ward.D2\", \"single\", \"complete\", \"average\", \"mcquitty\", \"median\", \"centroid\".")

  # check contrast levels against sample groups
  contra_levels_all <- unique(foreach(i = 1:length(object$input_data$comparisons$comparison_levels), .combine = "c") %do% {
    object$input_data$comparisons$comparison_levels[[i]]
  })
  if (!all(contra_levels_all %in% unique(levels(object$input_data$sample_groups)))) stop("Contrast levels not matching sample groups. Please check the input.")

  ## variables
  E <- object$input_data$norm_E
  genes <- object$input_data$genes
  input.genes_annotation.gene_symbol.var_name = object$input_data$input.genes_annotation.gene_symbol.var_name
  if (gene_symbol.only && ! input.genes_annotation.gene_symbol.var_name %in% names(object$input_data$genes)) {
    cat("Argument input.genes_annotation.gene_symbol.var_name not found in genes data frame when gene_symbol.only = TRUE, automatically set gene_symbol.only = FALSE.\n\n")
    gene_symbol.only <- FALSE
  }
  input.genes_annotation.gene_id.var_name = object$input_data$input.genes_annotation.gene_id.var_name
  export.name <- deparse(substitute(object))
  input.sample_groups <- object$input_data$sample_groups
  input.genes_annotation.control_type <- object$input_data$input.genes_annotation.control_type
  comparisons <- object$input_data$comparisons$comparisons
  comparison_levels <- object$input_data$comparisons$comparison_levels
  comp_to_remove <- which(as.numeric(object$significant_change_summary[, "True"]) < 2)
  thresholding_summary <- object$thresholding_summary
  row.lab.var_name <- input.genes_annotation.gene_id.var_name

  if (is.null(input.genes_annotation.control_type)) {
    cat("Argument input.genes_annotation.control_type is NULL, no control probes are removed.\n\n")
    rm.control <- FALSE
  } else {
    rm.control <- TRUE
  }

  if (length(comp_to_remove) > 0) {
    cat("Comparisons with less than two significant changes were removed: ", comparisons[comp_to_remove], "\n")
    comparisons <- comparisons[-comp_to_remove]
    comparison_levels <- comparison_levels[-comp_to_remove]
    thresholding_summary <- thresholding_summary[-comp_to_remove]
  }

  ## set up data
  # set up dis and cluster functions
  distfunc <- function(x)dist(x, method = distance)
  clustfunc <- function(x)hclust(x, method = clust)

  # prepare dfm for clustering
  dfm <- data.frame(genes, E, check.names = FALSE)

  if (rm.control){ # remove control
    dfm <- dfm[dfm[, input.genes_annotation.control_type$control_type.var_name] == input.genes_annotation.control_type$exp_type.value, ]
  }

  if (gene_symbol.only) {
    dfm <- dfm[complete.cases(dfm[, input.genes_annotation.gene_symbol.var_name]),]
    row.lab.var_name <- input.genes_annotation.gene_symbol.var_name
  }

  ## subsetting and plotting
  sig_dist_clust_list <- vector(mode = "list", length = length(comparisons))
  for (i in seq(length(comparisons))) {
    # set up plotting matrix
    plt_dfm <- dfm[as.logical(thresholding_summary[[i]]), ]
    plt_mtx <- as.matrix(plt_dfm[, !names(plt_dfm) %in% names(genes)])
    colnames(plt_mtx) <- sample_id.vector
    rownames(plt_mtx) <- plt_dfm[, row.lab.var_name]
    plt_mtx <- plt_mtx[, which(input.sample_groups %in% comparison_levels[[i]])]  # subsetting samples for the comparison levels
    row.lab <- plt_dfm[, row.lab.var_name]

    # set ColSideColors
    colGroup <- length(comparison_levels[[i]])
    col_dist <- distfunc(t(plt_mtx))
    col_cluster <- clustfunc(col_dist)
    colG <- cutree(col_cluster, colGroup) # column group
    colC <- brewer.pal(ifelse(colGroup < 3, 3, colGroup), col.colour) # column colour
    col_cluster_list <- list(col_dist = col_dist, col_hclust = col_cluster)

    # row cluster (genes/features)
    row_dist <- distfunc(plt_mtx)
    row_cluster <- clustfunc(row_dist)
    row_cluster_list <- list(row_dist = row_dist, row_hclust = row_cluster)

    # output
    sig_dist_clust_list[[i]] <- list(sig_col_dist_clust = col_cluster_list, sig_row_dist_clust = row_cluster_list)

    # draw heatmap
    if (verbose) cat(paste0("Sig data hierarchical clustering heatmap saved to: ", comparisons[i], "_sig_heatmap.pdf..."))
    pdf(file = paste0(comparisons[i], "_sig_heatmap.pdf"), width = plot.width, height = plot.height)
    heatmap.2(plt_mtx, distfun = distfunc, hclustfun = clustfunc,
              labRow = row.lab,
              col = brewer.pal(n.map.colour, map.colour), ColSideColors = colC[colG], ...)
    if (verbose) cat("Done!\n")
    dev.off()
  }
  names(sig_dist_clust_list) <- comparisons
  sig_dist_clust_list$distance_method <- distance
  sig_dist_clust_list$cluster_method <- clust
  assign(paste0(export.name, "_sig_dist_clust"), sig_dist_clust_list, envir = .GlobalEnv)
}


#' @title rbio_kmeans
#'
#' @description K means cluster function
#' @param x Matrix or Data frame. The function will cluster row items.
#' @param export.name Optional user defined export name prefix. Default is \code{NULL}.
#' @param k_range Integer array. The range of K to try to find the optimal number of cluster (K). Default is \code{2:15}.
#' @param nstart The number of tries to find the optimal cluster per K. Default is \code{25}.
#' @param plot.title Whether to display plot title on top of the plot. Default is \code{FALSE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
#' @param plot.symbolSize Symbol size. Default is \code{2}.
#' @param plot.vline.label Boolean. If to show the label the on the silhouette
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is available on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.xLabel x axis label. Type with quotation marks. Could be NULL. Default is \code{"Features"}.
#' @param plot.xLabelSize x axis label size. Default is \code{10}.
#' @param plot.xTickLblSize Font size of x axis ticks. Default is \code{10}.
#' @param plot.xTickItalic Set x axis tick font to italic. Default is \code{FALSE}.
#' @param plot.xTickBold Set x axis tick font to bold. Default is \code{FALSE}.
#' @param plot.xAngle The rotation angle (degrees) of the x axis marks. Default is \code{0} - horizontal.
#' @param plot.xhAlign The horizontal alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param plot.xvAlign The vertical alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param plot.yLabel y axis label. Type with quotation marks. Could be NULL. Default is \code{"OOB error rate"}.
#' @param plot.yLabelSize y axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Font size of y axis ticks. Default is \code{10}.
#' @param plot.yTickItalic Set y axis tick font to italic. Default is \code{FALSE}.
#' @param plot.yTickBold Set y axis tick font to bold. Default is \code{FALSE}.
#' @param plot.rightsideY If to display the right side y-axis. Default is \code{TRUE}.
#' @param plot.Width The width of the plot (unit: mm). Default is 170. Default will fit most of the cases.
#' @param plot.Height The height of the plot (unit: mm). Default is 150. Default will fit most of the cases.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @details The function does three things:
#'          1. Conducts a silhouette score test to find out the optimal K (number of clusters)
#'          2. Exports the silhouette score plot
#'          3. Conducts the K means cluster analysis with the optimized K
#'
#'          For the input \code{x}, it typically does not contain annotation columns for microarray or RNA-seq data sets.
#'          The function clusters the row items.
#'          To have proper labels for the cluser plot using \code{rbio_kmeans_plot()}, make sure to have rownames set for \code{x}.
#'
#' @return A silhouette score plot showing the optimized number of clusters (K), and a \code{kmeans} object containing the final K means cluster result.
#' @import ggplot2
#' @import ggrepel
#' @import foreach
#' @importFrom cluster silhouette
#' @importFrom RBioplot rightside_y
#' @importFrom grid grid.newpage grid.draw
#' @export
rbio_kmeans <- function(x, export.name = NULL,
                        k_range = 2:15, nstart = 25,
                        plot.title = NULL, plot.titleSize = 10,
                        plot.symbolSize = 2,
                        plot.vline.label = FALSE,
                        plot.fontType = "sans",
                        plot.xLabel = "Number of clusters (K)", plot.xLabelSize = 10,
                        plot.xTickLblSize = 10, plot.xTickItalic = FALSE, plot.xTickBold = FALSE,
                        plot.xAngle = 0, plot.xhAlign = 0.5, plot.xvAlign = 0.5,
                        plot.yLabel = "Average silhouette score", plot.yLabelSize = 10, plot.yTickLblSize = 10, plot.yTickItalic = FALSE,
                        plot.yTickBold = FALSE,
                        plot.rightsideY = TRUE,
                        plot.Width = 170, plot.Height = 150,
                        verbose = TRUE) {
  # -- argument check --
  if (!any(class(x) %in% c("matrix", "data.frame"))) stop("x needs to be a matrix or data frame.")
  if (any(class(x) %in% "data.frame")){
    x <- as.matrix(x)
  } else {
    x <- x
  }
  if (is.null(export.name)) {
    export.name <- deparse(substitute(x))
  }

  # -- load data --
  mtx <- as.matrix(x)

  # -- silhouette score --
  avg_sil <- foreach(i = k_range, .combine = "c") %do% {
    km <- kmeans(mtx, centers = i, nstart = nstart)
    ss <- cluster::silhouette(km$cluster, dist(mtx))
    mean(ss[, 3])
  }
  names(avg_sil) <- k_range
  optim_k <- as.integer(names(avg_sil)[avg_sil == max(avg_sil)])

  # -- final kmeans cluster --
  final_km <- kmeans(mtx, centers = optim_k, nstart = nstart)

  # -- export results --
  assign(paste0(export.name, "_kmeans"), final_km, envir = .GlobalEnv)

  # -- plotting --
  # - silhouette plot -
  pltdfm <- data.frame(x = k_range, y = avg_sil)
  loclEnv <- environment()
  baseplt <- ggplot(pltdfm, aes(x = k_range, y = avg_sil), environment = loclEnv) +
    geom_line() +
    geom_point(size = plot.symbolSize) +
    scale_x_continuous(expand = c(0.05, 0.05)) +
    ggtitle(plot.title) +
    xlab(plot.xLabel) +
    ylab(plot.yLabel) +
    geom_vline(xintercept = pltdfm$x[pltdfm$y == max(pltdfm$y)], linetype = "dashed", colour = "red") +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType),
          axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
          axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
          axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle,
                                     hjust = plot.xhAlign, vjust = plot.xvAlign),
          axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5),
          axis.ticks.x = if(plot.xTickLblSize == 0) element_blank())

  if (plot.vline.label) {
    baseplt <- baseplt +
      geom_text_repel(data = pltdfm[pltdfm$y == max(pltdfm$y), ],
                      aes(x = pltdfm$x[pltdfm$y == max(pltdfm$y)], y = max(pltdfm$y)*0.4, label = paste0(pltdfm$x[pltdfm$y == max(pltdfm$y)], " cluster(s)")),
                      arrow = arrow(length = unit(0.015, "npc")))
  }

  if (plot.xTickItalic & plot.xTickBold){
    baseplt <- baseplt +
      theme(axis.text.x = element_text(face = "bold.italic"))
  } else if (plot.xTickItalic & !plot.xTickBold){
    baseplt <- baseplt +
      theme(axis.text.x = element_text(face = "italic"))
  } else if (plot.xTickBold & !plot.xTickItalic){
    baseplt <- baseplt +
      theme(axis.text.x = element_text(face = "bold"))
  }

  if (plot.yTickItalic & plot.yTickBold){
    baseplt <- baseplt +
      theme(axis.text.y  = element_text(face = "bold.italic"))
  } else if (plot.yTickItalic & !plot.yTickBold){
    baseplt <- baseplt +
      theme(axis.text.y = element_text(face = "italic"))
  } else if (plot.yTickBold & !plot.yTickItalic){
    baseplt <- baseplt +
      theme(axis.text.y = element_text(face = "bold"))
  }

  # add the right-side y axis
  pltgtb <- RBioplot::rightside_y(baseplt)   # add the right-side y axis
  ggsave(filename = paste0(export.name,".kmeans_silouette.pdf"), plot = pltgtb,
         width = plot.Width, height = plot.Height, units = "mm",dpi = 600) # export the plot
  grid.newpage()
  grid.draw(pltgtb)
}
