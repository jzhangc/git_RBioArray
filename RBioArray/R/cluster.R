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
#' @description \code{rbio_unsupervised_hcluster} function for \code{rbioarray_flist} class object.
#' @param object Input object in \code{rbioarray_flist} class.
#' @param sample_id.var.name Variable name for sample identification, typically from \code{object$target}.
#' @param ... Additional arguments for the default method.
#' @return A pdf file containing a heatmap for unsupervised hierarchical clustering analysis.
#' @export
rbio_unsupervised_hcluster.rbioarray_flist <- function(object, sample_id.var.name = NULL, ...){
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

  export.name <- deparse(substitute(object))

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
#' @description \code{rbio_unsupervised_hcluster} function for \code{rbioseq_de} class object.
#' @param object Input object in \code{rbioseq_de} class.
#' @param sample_id.var.name Variable name for sample identification, typically from \code{object$target}.
#' @param ... Additional arguments for the default method.
#' @details The function uses filtered count data, as opposed to normalized data.
#'          Due to the compositional nature of NGS data, the count data is transformed using CLR method prior to clustering.
#' @return A pdf file containing a heatmap for unsupervised hierarchical clustering analysis.
#' @export
rbio_unsupervised_hcluster.rbioseq_de <- function(object, sample_id.var.name = NULL, ...){
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

  export.name <- deparse(substitute(object))

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
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @export
rbio_unsupervised_hcluster.default <- function(E, genes, input.sample_groups, n = "all",
                                               rm.control = FALSE, input.genes_annotation.control_type,
                                               gene_symbol.only = FALSE,
                                               input.genes_annotation.gene_symbol.var_name = NULL,
                                               input.genes_annotation.gene_id.var_name = NULL,
                                               sample_id.vector = NULL,
                                               distance = "euclidean", clust = "complete",
                                               col.colour = "Paired", map.colour = "PRGn", n.map.colour = 11, ...,
                                               export.name = NULL, plot.width = 7, plot.height = 7,
                                               verbose = TRUE){
  ## check arguments
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
  if (!distance %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
    stop("Argument distance needs to be one of \"euclidean\", \"maximum\", \"manhattan\", \"canberra\", \"binary\"or \"minkowski\".")
  if (!clust %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))
    stop("Argument clust needs to be one of \"ward.D\", \"ward.D2\", \"single\", \"complete\", \"average\", \"mcquitty\", \"median\", \"centroid\".")

  ## variables
  if(length(levels(input.sample_groups)) <= 19) {
    colGroup <- length(levels(input.sample_groups))
  } else {
    cat("The sample groups exceed the maximum allowed colour group number (19). Proceed with 19.\n\n")
    colGroup <- 19
  }
  row.lab.var_name <- input.genes_annotation.gene_id.var_name

  ## set up dis and cluster functions
  distfunc <- function(x)dist(x, method = distance)
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
  row.lab <- dfm[, row.lab.var_name]

  ## set ColSideColors
  col_cluster <- clustfunc(distfunc(t(mtx)))
  colG <- cutree(col_cluster, colGroup) # column group
  colC <- brewer.pal(ifelse(colGroup < 3, 3, colGroup), col.colour) # column colour

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
#' @param input.sample_groups Input \code{factor} object for sample groupping labels.
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
#' @param plot.width Width of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param plot.height Height of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @details Unlike the unsupervised veresion, the supervised hcluster uses normalized expression data for both RNAseq and microaray.
#'          The column colour group is usually 2.
#'
#' @return A heatmap based on hierarchical clustering analysis in \code{pdf} format.
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @export
rbio_supervised_hcluster <- function(object,
                                     gene_symbol.only = FALSE,
                                     sample_id.var.name = NULL,
                                     distance = "euclidean", clust = "complete",
                                     col.colour = "Paired", map.colour = "PRGn", n.map.colour = 11, ...,
                                     plot.width = 7, plot.height = 7,
                                     verbose = TRUE){
  ## check arguments
  if (class(object) != "sig") stop("The input object has to be a \"sig\" class.")
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
  if (!distance %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
    stop("Argument distance needs to be one of \"euclidean\", \"maximum\", \"manhattan\", \"canberra\", \"binary\"or \"minkowski\".")
  if (!clust %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))
    stop("Argument clust needs to be one of \"ward.D\", \"ward.D2\", \"single\", \"complete\", \"average\", \"mcquitty\", \"median\", \"centroid\".")

  # check contrast levels against sample groups
  contra_levels_all <- unique(foreach(i = 1:length(object$input_data$comparisons$comparison_levels), .combine = "c") %do% {
    object$input_data$comparisons$comparison_levels[[i]]
  })
  if (!all(contra_levels_all %in% unique(levels(object$input_data$sample_groups)))) stop("Contrast levels not matching sample groups. Please check the input.")

  ## variables
  E <- object$input_data$norm_E
  genes <- object$input_data$genes
  input.genes_annotation.gene_symbol.var_name = object$input_data$input.genes_annotation.gene_symbol.var_name
  input.genes_annotation.gene_id.var_name = object$input_data$input.genes_annotation.gene_id.var_name
  export.name <- deparse(substitute(object))
  input.sample_groups <- object$input_data$sample_groups
  input.genes_annotation.control_type <- object$input_data$input.genes_annotation.control_type
  comparisons <- object$input_data$comparisons$comparisons
  comparison_levels <- object$input_data$comparisons$comparison_levels
  thresholding_summary <- object$thresholding_summary
  row.lab.var_name <- input.genes_annotation.gene_id.var_name

  if (gene_symbol.only && ! input.genes_annotation.gene_symbol.var_name %in% names(genes)) {
    cat("Argument input.genes_annotation.gene_symbol.var_name not found in genes data frame when gene_symbol.only = TRUE, automatically set gene_symbol.only = FALSE.\n\n")
    gene_symbol.only <- FALSE
  }
  if (is.null(input.genes_annotation.control_type)) {
    cat("Argument input.genes_annotation.control_type is NULL, no control probes are removed.\n\n")
    rm.control <- FALSE
  } else {
    rm.control <- TRUE
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
  for (i in seq(length(comparisons))) {
    # set up plotting matrix
    plt_dfm <- dfm[as.logical(thresholding_summary[[i]]), ]
    plt_mtx <- as.matrix(plt_dfm[, !names(plt_dfm) %in% names(genes)])
    colnames(plt_mtx) <- sample_id.vector
    plt_mtx <- plt_mtx[, which(input.sample_groups %in% comparison_levels[[i]])]  # subsetting samples for the comparison levels
    row.lab <- plt_dfm[, row.lab.var_name]

    ## set ColSideColors
    colGroup <- length(comparison_levels[[i]])
    col_cluster <- clustfunc(distfunc(t(plt_mtx)))
    colG <- cutree(col_cluster, colGroup) # column group
    colC <- brewer.pal(ifelse(colGroup < 3, 3, colGroup), col.colour) # column colour

    # draw heatmap
    if (verbose) cat(paste0("Supervised hierarchical clustering heatmap saved to: ", comparisons[i], "_super_heatmap.pdf..."))
    pdf(file = paste0(comparisons[i], "_super_heatmap.pdf"), width = plot.width, height = plot.height)
    heatmap.2(plt_mtx, distfun = distfunc, hclustfun = clustfunc,
              labRow = row.lab,
              col = brewer.pal(n.map.colour, map.colour), ColSideColors = colC[colG], ...)
    if (verbose) cat("Done!\n")
    dev.off()
  }
}
