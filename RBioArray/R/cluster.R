#' @title rbio_unsupervised_hcluster
#'
#' @description Generic unsupersived hierarchical clustering function.
#' @param object Input object in either \code{rbioseq_de} or \code{rbioarray_flist} class.
#' @param ... Additional arguments for corresponding S3 class methods.
#' @details Since the function depends on the \code{heatmap.2} function from \code{gplots} package,
#'          arguments can be passed directly, as seen in the examples below.
#'
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
  if (!any(class(object) %in% c("rbioarray_flist", "rbioarray_de", "rbioseq_de"))) stop("The input object needs to be either \"rbioarray_flist\" or \"rbioseq_de\" class object.")

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


#' @title rbio_unsupervised_hcluster.rbioarray_de
#'
#' @rdname rbio_unsupervised_hcluster
#' @method rbio_unsupervised_hcluster rbioarray_de
#' @param object Input object in \code{rbioarray_de} class.
#' @param sample_id.var.name Variable name for sample identification, typically from \code{object$targets}.
#' @param export.name Optional user defined export name prefix. Default is \code{NULL}.
#' @param ... Additional arguments for the default method.
#' @export
rbio_unsupervised_hcluster.rbioarray_de <- function(object, sample_id.var.name = NULL, export.name = NULL, ...){
  ## check arguments
  if (!is.null(sample_id.var.name)){
    if (!sample_id.var.name %in% names(object$targets)) {
      cat("The sample_id.var.name not found in targets element of the input object. Proceed without using it.\n")
      sample_id.vector <- seq(ncol(object$input_data$E))
    } else {
      sample_id.vector <- object$targets[, sample_id.var.name]
    }
  } else {
    sample_id.vector <- seq(ncol(object$input_data$E))
  }

  if (is.null(export.name)){
    export.name <- deparse(substitute(object))
  } else {
    export.name <- export.name
  }

  ## use methods
  rbio_unsupervised_hcluster.default(E = object$input_data$E, genes = object$input_data$genes,
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
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
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
#' @param distance Distance calculation method. Default is \code{"euclidean"}. See \code{\link{dist}} for more.
#' @param clust Clustering method. Default is \code{"complete"}. See \code{\link{hclust}} for more.
#' @param col.colour Column group colour. Default is \code{"Paired"}. See \code{RColorBrewer} package for more.
#' @param map.colour Heat map colour. Default is \code{"PRGn"}. See \code{RColorBrewer} package for more.
#' @param n.map.colour Number of colours displayed. Default is \code{11}. See \code{RColorBrewer} package for more.
#' @param ... Additional arguments for \code{heatmap.2} function from \code{gplots} package.
#' @param plot.width Width of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param plot.height Height of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @details Unlike the unsupervised veresion, the sig data hcluster uses normalized expression data for both RNAseq and microarray.
#'          The column colour group is usually 2, since the function only outputs pair-wise comparison clusters
#'
#' @return A heatmap based on hierarchical clustering analysis in \code{pdf} format.
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @export
rbio_supervised_hcluster <- function(object,
                                     export.name = NULL,
                                     gene_symbol.only = FALSE,
                                     sample_id.var.name = NULL,
                                     distance = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                                     clust = c("complete", "ward.D", "ward.D2", "single",  "average", "mcquitty", "median", "centroid"),
                                     col.colour = "Paired", map.colour = "PRGn", n.map.colour = 11, ...,
                                     plot.width = 7, plot.height = 7,
                                     verbose = TRUE){
  ## check arguments
  if (is.null(export.name)){
    export.name <- deparse(substitute(object))
  } else {
    export.name <- export.name
  }

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
  input.genes_annotation.gene_symbol.var_name <- object$input_data$input.genes_annotation.gene_symbol.var_name
  if (gene_symbol.only && ! input.genes_annotation.gene_symbol.var_name %in% names(object$input_data$genes)) {
    cat("Argument input.genes_annotation.gene_symbol.var_name not found in genes data frame when gene_symbol.only = TRUE, automatically set gene_symbol.only = FALSE.\n\n")
    gene_symbol.only <- FALSE
  } else if (is.null(input.genes_annotation.gene_symbol.var_name)) {
    cat("input.genes_annotation.gene_symbol.var_name is NULL, automatically set gene_symbol.only = FALSE.\n\n")
    gene_symbol.only <- FALSE
  }
  input.genes_annotation.gene_id.var_name = object$input_data$input.genes_annotation.gene_id.var_name
  if (is.null(export.name)) {
    cat("Object/file export name prefix automatically set to the input object name. \n")
    export.name <- deparse(substitute(object))
  } else {
    export.name <- export.name
  }
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
    is.control <- dfm[, input.genes_annotation.control_type$control_type.var_name] == input.genes_annotation.control_type$exp_type.value
    dfm <- dfm[is.control, ]
    cat(paste0("Number of control targets removed: ", length(which(is.control)), "\n"))
  }

  if (gene_symbol.only) {
    cat("gene_symbol.only = TRUE, any genes without a gene symbol are removed from clustering.\n\n")
    row.lab.var_name <- input.genes_annotation.gene_symbol.var_name
  }

  ## subsetting and plotting
  sig_dist_clust_list <- vector(mode = "list", length = length(comparisons))
  # pair-wise comparison
  for (i in seq(length(comparisons))) {
    # set up plotting matrix
    is_test <- thresholding_summary[[i]]

    if (is.null(object$input_data$full_de_results[[comparisons[i]]]$PROBE_ID)) {  # for RNAseq data
      de_probes <- object$input_data$full_de_results[[comparisons[i]]]$gene_id[is_test]
      plt_dfm <- dfm[dfm$gene_id %in% de_probes, ]
      if (gene_symbol.only) {
        plt_dfm <- plt_dfm[complete.cases(plt_dfm$gene_name), ]
      }
    } else {
      de_probes <- object$input_data$full_de_results[[comparisons[i]]]$PROBE_ID[is_test]
      plt_dfm <- dfm[dfm$PROBE_ID %in% de_probes, ]
      if (gene_symbol.only) {
        plt_dfm <- plt_dfm[complete.cases(plt_dfm$SYMBOL), ]
      }
    }

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
    pdf(file = paste0(export.name, "_", comparisons[i], "s_sig_heatmap.pdf"), width = plot.width, height = plot.height)
    heatmap.2(plt_mtx, distfun = distfunc, hclustfun = clustfunc,
              labRow = row.lab,
              col = brewer.pal(n.map.colour, map.colour), ColSideColors = colC[colG], ...)
    if (verbose) cat("Done!\n")
    dev.off()
  }
  names(sig_dist_clust_list) <- comparisons

  # F stats clustering
  f_is_test <- thresholding_summary[["F_stats"]]

  if (is.null(object$input_data$full_de_results[["F_stats"]]$PROBE_ID)){
    f_de_probes <- object$input_data$full_de_results[["F_stats"]]$gene_id[f_is_test]
    f_plt_dfm <- dfm[dfm$gene_id %in% f_de_probes, ]
  } else {
    f_de_probes <- object$input_data$full_de_results[["F_stats"]]$PROBE_ID[f_is_test]
    f_plt_dfm <- dfm[dfm$PROBE_ID %in% f_de_probes, ]
  }

  if (is.null(f_plt_dfm$SYMBOL)){  # for RNAseq data
    if (gene_symbol.only) {
      f_plt_dfm <- f_plt_dfm[complete.cases(f_plt_dfm$gene_name), ]
    }
  } else {
    if (gene_symbol.only) {
      f_plt_dfm <- f_plt_dfm[complete.cases(f_plt_dfm$SYMBOL), ]
    }
  }


  f_plt_mtx <- as.matrix(f_plt_dfm[, !names(f_plt_dfm) %in% names(genes)])
  colnames(f_plt_mtx) <- sample_id.vector
  rownames(f_plt_mtx) <- f_plt_dfm[, row.lab.var_name]
  f_row.lab <- f_plt_dfm[, row.lab.var_name]

  # set ColSideColors
  f_colGroup <- length(unique(object$input_data$sample_groups))
  f_col_dist <- distfunc(t(f_plt_mtx))
  f_col_cluster <- clustfunc(f_col_dist)
  f_colG <- cutree(f_col_cluster, f_colGroup) # column group
  f_colC <- brewer.pal(ifelse(f_colGroup < 3, 3, f_colGroup), col.colour) # column colour
  f_col_cluster_list <- list(col_dist = f_col_dist, col_hclust = f_col_cluster)

  # row cluster (genes/features)
  f_row_dist <- distfunc(f_plt_mtx)
  f_row_cluster <- clustfunc(f_row_dist)
  f_row_cluster_list <- list(row_dist = f_row_dist, row_hclust = f_row_cluster)

  # output
  sig_dist_clust_list[["F_stats"]] <- list(sig_col_dist_clust = f_col_cluster_list, sig_row_dist_clust = f_row_cluster_list)

  # draw heatmap
  if (verbose) cat(paste0("Sig data hierarchical clustering heatmap saved to: ", "fstats_sig_heatmap.pdf..."))
  pdf(file = paste0(export.name, "_fstats_sig_heatmap.pdf"), width = plot.width, height = plot.height)
  heatmap.2(f_plt_mtx, distfun = distfunc, hclustfun = clustfunc,
            labRow = f_row.lab,
            col = brewer.pal(n.map.colour, map.colour), ColSideColors = f_colC[f_colG], ...)
  # heatmap.2(f_plt_mtx, distfun = distfunc, hclustfun = clustfunc,
  #           labRow = f_row.lab,
  #           col = brewer.pal(n.map.colour, map.colour), ColSideColors = f_colC[f_colG])
  if (verbose) cat("Done!\n")
  dev.off()

  ## save and export
  sig_dist_clust_list$distance_method <- distance
  sig_dist_clust_list$cluster_method <- clust
  assign(paste0(export.name, "_sig_dist_clust"), sig_dist_clust_list, envir = .GlobalEnv)
}


#' @title rbio_kmeans
#'
#' @description K means cluster function
#' @param x Matrix or Data frame. Input data matrix. The function will cluster row items.
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
#' @param plot.width The width of the plot (unit: mm). Default is 170. Default will fit most of the cases.
#' @param plot.height The height of the plot (unit: mm). Default is 150. Default will fit most of the cases.
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
#' @examples
#' \dontrun{
#' mtx <- group_de_list$ex1$E
#' rownames(mtx) <- group_de_list$ex1$genes$SYMBOL
#' rbio_kmeans(x = mtx, export.name = "ex1", plot.vline.label = TRUE)
#' }
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
                        plot.width = 170, plot.height = 150,
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
         width = plot.width, height = plot.height, units = "mm",dpi = 600) # export the plot
  grid.newpage()
  grid.draw(pltgtb)
}


#' @title rbio_kmeans_plot
#'
#' @description A plot function for K means cluster results, with or without PCA.
#' @param km.object A \code{kmeans} object, from \code{\link{rbio_kmeans()}} function.
#' @param x Matrix or data.frame. The same data used to generate the input \code{km.object}.
#' @param centerScale Boolean. If to center and scale the x data fro clustering. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{prcomp} function.
#' @param plot.comps Integer or vector of integers. Index number(s) for principal component(s) to plot. Default is \code{c(1, 2)}.
#' @param plot.title The biplot title. Default is \code{NULL}.
#' @param plot.sampleLabel.type  If to show the sample labels on the graph. Options are \code{"none"}, \code{"direct"} and \code{"indirect"}. Default is \code{"none"}.
#' @param plot.sampleLabelSize Only set when \code{biplot.sampleLabel.type} is not \code{"none"}, The size of the sample label. Default is \code{2}.
#' @param plot.sampleLabel.padding Set only when \code{biplot.sampleLabel.type = "indirect"}, the padding between sample symbol and the label. Default is \code{0.5}.
#' @param plot.SymbolSize The symbol size for the scatter plot portion of the biplot. Default is \code{2}.
#' @param plot.ellipse If to draw ellipses. Default is \code{FALSE}.
#' @param plot.ellipse_conf The confidence value for the ellipses. Default is \code{0.93}.
#' @param plot.xAngle The rotation angle (degrees) of the biplot x axis marks. Default is \code{0} - horizontal.
#' @param plot.xhAlign The horizontal alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param plot.xvAlign The vertical alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param plot.loadingplot If to superimpose loading plot. Default is \code{FALSE}.
#' @param plot.loadingplot.textsize The font size of the loading plot labels. Default is \code{3}.
#' @param plot.mtx.densityplot If to display a density plot on the diagonal for the correlation scoreplot matrix. Default is \code{FALSE}.
#' @param plot.mtx.stripLblSize The label font size for the correlation scoreplot matrix strips. Default is \code{10}.
#' @param plot.width The biplot width. Default is \code{170}.
#' @param plot.height The biplot height. Default is \code{150}.
#' @param plot.rightsideY If to show the right side y-axis for both boxplot and biplot. For biplot, only applicable when the length of \code{comps} is less than 2, inclusive. Default is \code{FALSE}.
#' @param plot.fontType Font for the figure texts. Default is \code{"sans"}.
#' @param plot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param plot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return A PCA/scatter plot from K means cluster. The format is \code{pdf}.
#' @details To have proper labels for the cluser plot using \code{rbio_kmeans_plot()}, make sure to have rownames set for \code{x}.
#'
#'          The function automatically detects if PCA is needed:
#'          when \code{x} has more than two variables, the function will run PCA first.
#'          \code{plot.comps} not needed and will not be used when \code{x} has less than three variables, i.e. no PCA needed for plotting.
#' @import ggplot2
#' @import ggrepel
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y
#' @examples
#' \dontrun{
#' rbio_kmeans_plot(km.object = ex1_kmeans, x = mtx, export.name = "ex1",
#'                  plot.ellipse = TRUE,
#'                  plot.sampleLabel.type = "direct",
#'                  plot.SymbolSize = 4,
#'                  plot.xLabelSize = 14, plot.yLabelSize = 14,
#'                  plot.width = 140, plot.height = 150)
#' }
#' @export
rbio_kmeans_plot <- function(km.object,
                             x, export.name = NULL,
                             centerScale = TRUE,...,
                             plot.comps = c(1:2),
                             plot.title = NULL,
                             plot.sampleLabel.type = "none", plot.sampleLabelSize = 2,
                             plot.sampleLabel.padding = 0.5,
                             plot.SymbolSize = 2,
                             plot.ellipse = FALSE, plot.ellipse_conf = 0.95,
                             plot.xLabelSize = 10, plot.yLabelSize = 10,
                             plot.xTickLblSize = 10, plot.yTickLblSize = 10,
                             plot.xAngle = 0, plot.xhAlign = 0.5, plot.xvAlign = 0.5,
                             plot.mtx.densityplot = FALSE, plot.mtx.stripLblSize = 10,
                             plot.rightsideY = FALSE,
                             plot.fontType = "sans",
                             plot.width = 170, plot.height = 150,
                             verbose = TRUE){
  # -- check argument --
  if (!any(class(x) %in% c("matrix", "data.frame"))) stop("x needs to be a matrix or data frame.")
  if (any(class(x) %in% "data.frame")){
    x <- as.matrix(x)
  } else {
    x <- x
  }
  if (is.null(rownames(x))) {
    if (verbose) cat("No rownames detected for X. Proceed with row numbers as rownames. \n")
    rownames(x) <- seq(nrow(x))
    }
  if (is.null(export.name)) {
    export.name <- deparse(substitute(x))
  }
  if (!any(class(km.object) %in% "kmeans")) stop("km.object needs to be a kmeans class.")
  if (length(km.object$cluster) != nrow(x)) stop("km.object does not mean the data x. ")

  # cluster information
  clusters <- factor(km.object$cluster, levels = unique(km.object$cluster))

  # -- PCA or not --
  if (ncol(x) <= 2) {
    pca_on <- FALSE
  } else {
    pca_on <- TRUE
    if (!all(plot.comps %in% seq(ncol(x)))) stop("plot.comps contain non-existant PC.")
  }
  if (verbose) cat("PCA:", pca_on, "\n")
  if (pca_on){
    PCA <- prcomp(x, scale. = centerScale, center = centerScale, ...)
    # PCA <- prcomp(x, scale. = centerScale, center = centerScale)
    varpp_x <- 100 * summary(PCA)$importance[2, ] # extract and calculate the proportion of variance
    score_x <- data.frame(PCA$x[, plot.comps, drop = FALSE], check.names = FALSE) # extract rotated sample scores
    score_x$group <- clusters
    score_x$sample.label <- rownames(x)
    var_percentage_x <- varpp_x[paste0("PC", plot.comps)] # extract the proportion of variance for the selected PCs
    pc_axis_lbl <- paste("PC ", plot.comps, " (", round(var_percentage_x, digits = 2), "%)", sep = "")
  } else {
    if (centerScale){
      x <- RBioFS::center_scale(x)$centerX
    }
    score_x <- data.frame(x, check.names = FALSE)
    pc_axis_lbl <- colnames(x)
    score_x$group <- clusters
    score_x$sample.label <- rownames(x)
  }

  # -- plot --
  if (pca_on) {
    if (length(plot.comps) == 1){
      if (plot.ellipse){
        cat("ellipse is not applicable when length(plot.comps) == 1. Proceed without one.\n")
      }
      score_x$sample <- seq(nrow(score_x))
      names(score_x)[1] <- "axis1"

      if (verbose) cat(paste("Single PC plot being saved to file: ", export.name, ".kmeans_cluster.pdf...", sep = ""))  # initial message
      biplt <- ggplot(score_x, aes(x = sample, y = axis1)) +
        geom_line(aes(colour = group, linetype = group))

      # sample labels
      if (plot.sampleLabel.type == "direct"){
        biplt <- biplt +
          geom_text(aes(colour = group, label = sample.label), size = plot.SymbolSize)
      } else if (plot.sampleLabel.type == "indirect"){
        biplt <- biplt +
          geom_point(aes(shape = group, colour = group), size = plot.SymbolSize) +
          scale_shape_manual(values=1:nlevels(score_x$group)) +
          geom_text_repel(aes(label = sample.label), point.padding = unit(plot.sampleLabel.padding, "lines"),
                          show.legend = FALSE, size = plot.sampleLabelSize)
      } else {
        biplt <- biplt +
          geom_point(aes(shape = group, colour = group), size = plot.SymbolSize) +
          scale_shape_manual(values=1:nlevels(score_x$group))
      }

      biplt <- biplt +
        ggtitle(plot.title) +
        ylab(pc_axis_lbl[1]) +
        theme_bw() +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", family = plot.fontType, hjust = 0.5),
              axis.title = element_text(face = "bold", family = plot.fontType),
              legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
              axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle, hjust = plot.xhAlign, vjust = plot.xhAlign),
              axis.text.y = element_text(size = plot.xTickLblSize, family = plot.fontType, hjust = 0.5),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType))

      # grid.newpage()
      if (plot.rightsideY){ # add the right-side y axis
        biplt <- RBioplot::rightside_y(biplt)
      }
    } else if (length(plot.comps) == 2){
      names(score_x)[1:2] <- c("axis1", "axis2")

      if (verbose) cat(paste("plot being saved to file: ", export.name, ".kmeans_cluster.pdf...", sep = ""))  # initial message
      biplt <- ggplot(score_x, aes(x = axis1, y = axis2))

      # sample labels
      if (plot.sampleLabel.type == "direct"){
        biplt <- biplt +
          geom_text(aes(colour = group, label = sample.label), size = plot.SymbolSize)
      } else if (plot.sampleLabel.type == "indirect"){
        biplt <- biplt +
          geom_point(aes(shape = group, colour = group), size = plot.SymbolSize) +
          scale_shape_manual(values=1:nlevels(score_x$group)) +
          geom_text_repel(aes(label = sample.label), point.padding = unit(plot.sampleLabel.padding, "lines"),
                          show.legend = FALSE, size = plot.sampleLabelSize)
      } else {
        biplt <- biplt +
          geom_point(aes(shape = group, colour = group), size = plot.SymbolSize) +
          scale_shape_manual(values=1:nlevels(score_x$group))
      }

      biplt <- biplt +
        ggtitle(plot.title) +
        xlab(pc_axis_lbl[1]) +
        ylab(pc_axis_lbl[2]) +
        theme_bw() +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", family = plot.fontType, hjust = 0.5),
              axis.title = element_text(face = "bold", family = plot.fontType),
              legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
              axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle, hjust = plot.xhAlign, vjust = plot.xvAlign),
              axis.text.y = element_text(size = plot.xTickLblSize, family = plot.fontType, hjust = 0.5),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType))

      if (plot.ellipse){ # circles
        biplt <- biplt +
          stat_ellipse(aes(colour = group, group = group), type = "norm", level = plot.ellipse_conf)
      }

      # grid.newpage()
      if (plot.rightsideY){ # add the right-side y axis
        biplt <- RBioplot::rightside_y(biplt)
      }

    } else if (length(plot.comps) > 2) {
      if (plot.rightsideY){ # add the right-side y axis
        cat("Rightside y-axis is not applicable when length(plot.comps) > 2. Proceed without one.\n")
      }

      # custom functions for the paired scoreplot
      ellipsefunc <- function(data = score_x, mapping, label.method = plot.sampleLabel.type,
                              ellipse = plot.ellipse, ellipse_conf = plot.ellipse_conf,
                              ...){
        g <- ggplot(data = data, mapping = mapping)

        if (label.method == "direct"){
          g <- g + geom_text(aes(colour = group, label = sample.label), ...)
        } else if (label.method == "indirect"){
          g <- g +  geom_point(...) +
            scale_shape_manual(values = 1:nlevels(data$group)) +
            geom_text_repel(aes(label = sample.label), point.padding = unit(plot.sampleLabel.padding, "lines"),
                            show.legend = FALSE, size = plot.sampleLabelSize)
        } else {
          g <- g + geom_point(...) +
            scale_shape_manual(values = 1:nlevels(data$group))
        }
        if (ellipse){
          g <- g +
            stat_ellipse(aes(colour = group, group = group), type = "norm", level = ellipse_conf)
        }
        return(g)
      }

      densityfunc <- function(data = score_x, mapping, alpha = 0.1, densityplot = plot.mtx.densityplot, ...){
        g <- ggplot(data = data, mapping = mapping)
        if (densityplot){
          g <- g + geom_density(alpha = alpha, aes(colour = group, linetype = group, ...))
        }
        return(g)
      }

      # matrix scoreplot
      if (verbose) cat(paste("plot matrix being saved to file: ", export.name,".kmeans_cluster.pdf...", sep = ""))  # initial message
      biplt <- ggpairs(score_x, columns = plot.comps, aes(colour = group, shape = group),
                       axisLabels = "show", columnLabels = pc_axis_lbl,
                       showStrips = NULL,
                       lower = list(continuous = ellipsefunc),
                       upper = list(continuous = ellipsefunc),
                       diag = list(continuous = densityfunc),
                       legend = 2)
      biplt <- biplt +
        ggtitle(plot.title) +
        theme(plot.title = element_text(face = "bold", family = plot.fontType, hjust = 0.5),
              axis.title = element_text(face = "bold", family = plot.fontType),
              strip.background = element_blank(),  # no strip background colour
              strip.text = element_text(face = "bold", size = plot.mtx.stripLblSize),
              panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
              axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle, hjust = plot.xhAlign, vjust = plot.xvAlign),
              axis.text.y = element_text(size = plot.xTickLblSize, family = plot.fontType),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType))
      # grid.newpage()
    }
    ggsave(filename = paste(export.name,".kmeans_cluster.pdf", sep = ""), plot = biplt,
           width = plot.width, height = plot.height, units = "mm",dpi = 600)
    if (verbose) cat("Done!\n") # final message
    grid.draw(biplt)
  } else {  # without PCA
    if (ncol(x) == 1){
      if (plot.ellipse){
        cat("ellipse is not applicable when length(plot.comps) == 1. Proceed without one.\n")
      }
      score_x$sample <- seq(nrow(score_x))
      names(score_x)[1] <- "axis1"

      if (verbose) cat(paste("Single PC plot being saved to file: ", export.name, ".kmeans_cluster.pdf...", sep = ""))  # initial message
      biplt <- ggplot(score_x, aes(x = sample, y = axis1)) +
        geom_line(aes(colour = group, linetype = group))

      # sample labels
      if (plot.sampleLabel.type == "direct"){
        biplt <- biplt +
          geom_text(aes(colour = group, label = sample.label), size = plot.SymbolSize)
      } else if (plot.sampleLabel.type == "indirect"){
        biplt <- biplt +
          geom_point(aes(shape = group, colour = group), size = plot.SymbolSize) +
          scale_shape_manual(values=1:nlevels(score_x$group)) +
          geom_text_repel(aes(label = sample.label), point.padding = unit(plot.sampleLabel.padding, "lines"),
                          show.legend = FALSE, size = plot.sampleLabelSize)
      } else {
        biplt <- biplt +
          geom_point(aes(shape = group, colour = group), size = plot.SymbolSize) +
          scale_shape_manual(values=1:nlevels(score_x$group))
      }

      biplt <- biplt +
        ggtitle(plot.title) +
        ylab(pc_axis_lbl[1]) +
        theme_bw() +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", family = plot.fontType, hjust = 0.5),
              axis.title = element_text(face = "bold", family = plot.fontType),
              legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
              axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle, hjust = plot.xhAlign, vjust = plot.xhAlign),
              axis.text.y = element_text(size = plot.xTickLblSize, family = plot.fontType, hjust = 0.5),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType))

      # grid.newpage()
      if (plot.rightsideY){ # add the right-side y axis
        biplt <- RBioplot::rightside_y(biplt)
      }
    } else if (ncol(x) == 2){
      names(score_x)[1:2] <- c("axis1", "axis2")

      if (verbose) cat(paste("plot being saved to file: ", export.name, ".kmeans_cluster.pdf...", sep = ""))  # initial message
      biplt <- ggplot(score_x, aes(x = axis1, y = axis2))

      # sample labels
      if (plot.sampleLabel.type == "direct"){
        biplt <- biplt +
          geom_text(aes(colour = group, label = sample.label), size = plot.SymbolSize)
      } else if (plot.sampleLabel.type == "indirect"){
        biplt <- biplt +
          geom_point(aes(shape = group, colour = group), size = plot.SymbolSize) +
          scale_shape_manual(values=1:nlevels(score_x$group))
        geom_text_repel(aes(label = sample.label), point.padding = unit(plot.sampleLabel.padding, "lines"),
                        show.legend = FALSE, size = plot.sampleLabelSize)
      } else {
        biplt <- biplt +
          geom_point(aes(shape = group, colour = group), size = plot.SymbolSize) +
          scale_shape_manual(values=1:nlevels(score_x$group))
      }

      biplt <- biplt +
        ggtitle(plot.title) +
        xlab(pc_axis_lbl[1]) +
        ylab(pc_axis_lbl[2]) +
        theme_bw() +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", family = plot.fontType, hjust = 0.5),
              axis.title = element_text(face = "bold", family = plot.fontType),
              legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
              axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle, hjust = plot.xhAlign, vjust = plot.xvAlign),
              axis.text.y = element_text(size = plot.xTickLblSize, family = plot.fontType, hjust = 0.5),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType))

      if (plot.ellipse){ # circles
        biplt <- biplt +
          stat_ellipse(aes(colour = group, group = group), type = "norm", level = plot.ellipse_conf)
      }

      # grid.newpage()
      if (plot.rightsideY){ # add the right-side y axis
        biplt <- RBioplot::rightside_y(biplt)
      }
    }
    ggsave(filename = paste(export.name,".kmeans_cluster.pdf", sep = ""), plot = biplt,
           width = plot.width, height = plot.height, units = "mm",dpi = 600)
    if (verbose) cat("Done!\n") # final message
    grid.draw(biplt)
  }
}
