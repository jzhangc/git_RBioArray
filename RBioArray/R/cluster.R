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
#'                            sample_id.var.name = "condition", srtCol = 90, offsetCol = 0, labRow = FALSE,
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
#' @param ... Additional arguments for corresponding S3 class methods.
#' @return A pdf file containing a heatmap for unsupervised hierarchical clustering analysis.
#' @export
rbio_unsupervised_hcluster.rbioarray_flist <- function(object, sample_id.var.name = NULL, ...){
  ## check arguments
  if (!is.null(sample_id.var.name)){
    if (!sample_id.var.name %in% names(object$targets)) {
      cat("The sample_id.var.name not found in targets element of the input object. Proceed without using it.\n")
      sample_id.vector <- NULL
    } else {
      sample_id.vector <- object$targets[, sample_id.var.name]
    }
  } else {
    sample_id.vector <- NULL
  }

  export.name <- deparse(substitute(object))

  ## use methods
  rbio_unsupervised_hcluster.default(E = object$E, genes = object$genes,
                                     input.gene_annotation.control_type = object$gene_annotation.control_type,
                                     input.genes_annotation.gene_symbol.var_name = object$genes_annotation.gene_symbol.var_name,
                                     input.sample_gourps = object$sample_groups,
                                     sample_id.vector = sample_id.vector, export.name = export.name, ...)
}


#' @title rbio_unsupervised_hcluster.rbioseq_de
#'
#' @description \code{rbio_unsupervised_hcluster} function for \code{rbioseq_de} class object.
#' @param object Input object in \code{rbioseq_de} class.
#' @param sample_id.var.name Variable name for sample identification, typically from \code{object$target}.
#' @param ... Additional arguments for corresponding S3 class methods.
#' @details Due to the compositional nature of NGS data, the count data is transformed using CLR method prior to clustering.
#' @return A pdf file containing a heatmap for unsupervised hierarchical clustering analysis.
#' @export
rbio_unsupervised_hcluster.rbioseq_de <- function(object, sample_id.var.name = NULL, ...){
  ## check arguments
  if (!is.null(sample_id.var.name)){
    if (!sample_id.var.name %in% names(object$targets)) {
      cat("The sample_id.var.name not found in targets element of the input object. Proceed without using it.\n")
      sample_id.vector <- NULL
    } else {
      sample_id.vector <- object$targets[, sample_id.var.name]
    }
  } else {
    sample_id.vector <- NULL
  }

  export.name <- deparse(substitute(object))

  ## transform
  cat("CLR transformation of filtered RNAseq count data...")
  E_transfo <- rbioseq_clr_ilr_transfo(object$filter_results$filtered_counts$counts, offset = 1, mode = "clr")  # clr tranformation
  cat("Done!\n")

  ## use methods
  rbio_unsupervised_hcluster.default(E = E_transfo, genes = object$filter_results$filtered_counts$genes,
                                     input.gene_annotation.control_type = NULL,
                                     input.genes_annotation.gene_symbol.var_name = object$gene_symbol_var_name,
                                     input.sample_gourps = object$sample_groups,
                                     sample_id.vector = sample_id.vector, export.name = export.name, rm.control = FALSE, ...)
}


#' @title rbio_unsupervised_hcluster.default
#'
#' @description Default unsupersived hierarchical clustering function.
#' @param E Expression or count matrix, with rows for genes/probes/genomic features, columns for RNA samples.
#' @param genes Annotation data frame for genes/probes/genomic features.
#' @param input.sample_gourps Input \code{factor} object for sample groupping labels.
#' @param n Number of genes/probes/genomic features to be clustered, numeric input or \code{"all"}. Default is \code{"all"}.
#' @param rm.control Whether to remove control probes (Agilent platform) or not. Default is \code{TRUE}.
#' @param input.gene_annotation.control_type Only set when \code{rm.control = TRUE}, input control type variable annotation list.
#' @param gene_symbol.only Whether or not to remove probes without gene symbol. Default is \code{FALSE}.
#' @param input.genes_annotation.gene_symbol.var_name Only set when \code{gene_symbol.only = TRUE}, variable name for gene symbol column in \code{genes} data frame.
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
#' @return A heatmap based on hierarchical clustering analysis in \code{pdf} format.
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @export
rbio_unsupervised_hcluster.default <- function(E, genes, input.sample_gourps, n = "all",
                                               rm.control = FALSE, input.gene_annotation.control_type,
                                               gene_symbol.only = FALSE, input.genes_annotation.gene_symbol.var_name,
                                               sample_id.vector = NULL,
                                               distance = "euclidean", clust = "complete",
                                               col.colour = "Paired", map.colour = "PRGn", n.map.colour = 11, ...,
                                               export.name = NULL, plot.width = 7, plot.height = 7){
  ## check arguments
  if (n != "all" && n %% 1 != 0) stop("Argument n needs to be either \"all\" or an integer number.")
  if (n.map.colour %% 1 != 0) stop("Argument n.map.colour needs to be an integer number.")
  if (rm.control && is.null(input.gene_annotation.control_type)) {
    cat("Argument input.gene_annotation.control_type is NULL when rm.control = TRUE, automatically set rm.control = FALSE.\n\n")
    rm.control <- FALSE
  }
  if (gene_symbol.only && ! input.genes_annotation.gene_symbol.var_name %in% names(genes)) {
    cat("Argument input.genes_annotation.gene_symbol.var_name not found in genes data frame when gene_symbol.only = TRUE, automatically set gene_symbol.only = FALSE.\n\n")
    gene_symbol.only <- FALSE
  }
  if (missing(export.name) || is.null(export.name)) stop("Please set value for argument export.name.")
  if (!tolower(distance) %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
    stop("Argument distance needs to be one of \"euclidean\", \"maximum\", \"manhattan\", \"canberra\", \"binary\"or \"minkowski\".")
  if (!tolower(clust) %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))
    stop("Argument clust needs to be one of \"ward.D\", \"ward.D2\", \"single\", \"complete\", \"average\", \"mcquitty\", \"median\", \"centroid\".")

  ## variables
  if(length(levels(input.sample_gourps)) <= 19) {
    colGroup <- length(levels(input.sample_gourps))
  } else {
    cat("The sample groups exceed the maximum allowed colour group number (19). Proceed with 19.\n\n")
    colGroup <- 19
  }

  ## set up dis and cluster functions
  distfunc <- function(x)dist(x, method = distance)
  clustfunc <- function(x)hclust(x, method = clust)

  ## prepare dfm for clustering
  dfm <- data.frame(genes, E, check.names = FALSE)

  if (rm.control){ # remove control
    dfm <- dfm[dfm[, input.gene_annotation.control_type$control_type.var_name] == input.gene_annotation.control_type$exp_type.value, ]
  }

  if (gene_symbol.only) {
    dfm <- dfm[complete.cases(dfm[, input.genes_annotation.gene_symbol.var_name]),]
  }

  if (n != "all"){ # subset
    dfm <- dfm[1:n, ]
  }

  mtx <- as.matrix(dfm[, !names(dfm) %in% names(genes)])

  ## set ColSideColors
  col_cluster <- clustfunc(distfunc(t(mtx)))
  colG <- cutree(col_cluster, colGroup) # column group
  colC <- brewer.pal(ifelse(colGroup < 3, 3, colGroup), col.colour) # column colour

  ## prepare mtx for plotting
  if (!is.null(sample_id.vector)){
    colnames(mtx) <- sample_id.vector
  }

  ## heatmap
  # draw heatmap
  cat(paste0("Unsupervised hierarchical clustering heatmap saved to: ", export.name, "_unsuper_heatmap.pdf..."))
  pdf(file = paste0(export.name, "_unsuper_heatmap.pdf"), width = plot.width, height = plot.height)
  heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
            col = brewer.pal(n.map.colour, map.colour), ColSideColors = colC[colG], ...)
  cat("Done!\n")
  dev.off()
}


#' @title rbioarray_hcluster_super
#'
#' @description Legacy function. Wrapper for supervised hierarchical clustering analysis and heatmap visualization for both microarray and RNAseq.
#' @param plotName File name for the export \code{pdf} plot file. Default is \code{"data"}.
#' @param fltDOI Based on filtered data, a subset corresponding to the comparasion, either a list, \code{EList} or \code{MAList} object.
#' @param dfmDE A subset of the DE list, i.e. a \code{topTable} dataframe, corresponding to the comparasion (i.e., contrast).
#' @param ctrlProbe Wether or not the data set has control type variable in \code{fltDOI}, with values \code{-1 (negative control)}, \code{0 (gene probes)} and \code{1 (positive control)}. Default is \code{TRUE}.
#' @param ctrlTypeVar Set only when \code{ctrlProbe = TRUE}, the control type variable. Default is the \code{Agilent} variable name \code{"ControlType"}.
#' @param dataProbeVar \code{dfmDE} variable name for probe name. Default is \code{"ProbeName"}.
#' @param DE.sig.method Thresholding method, "fdr", "spikein" or "none. Default is \code{"none"}.
#' @param DE.sig.p P value cut off. Default is \code{0.05}.
#' @param FC Fold change (FC) filter for the heatmap. Default is \code{NULL}.
#' @param fct Input \code{factor} object for samples.
#' @param sampleName A \code{vector} containing names for column. Default is \code{NULL} and the function will use the column name from the input.
#' @param colGroup Colour group, numeric or dependent on \code{fct}.
#' @param distance Distance calculation method. Default is \code{"euclidean"}. See \code{\link{dist}} for more.
#' @param clust Clustering method. Default is \code{"complete"}. See \code{\link{hclust}} for more.
#' @param rowLabel Whether to display row label or not. Default is \code{FALSE}.
#' @param annot Annotation object, usually a \code{dataframe}. Make sure to name the probe ID variable \code{ProbeName}. Only set this argument when \code{rowLabel = TRUE}. Default is \code{NULL}.
#' @param annotProbeVar \code{annot} variable name for probe name. Default is \code{"ProbeName"}.
#' @param genesymbolVar The name of the variable for gene symbols from the \code{annot} object. Only set this argument when \code{rowLabel = TRUE}. Default is \code{NULL}.
#' @param colColour Column group colour. Default is \code{"Paired"}. See \code{RColorBrewer} package for more.
#' @param mapColour Heat map colour. Default is \code{"PRGn"}. See \code{RColorBrewer} package for more.
#' @param n_mapColour Number of colours displayed. Default is \code{11}. See \code{RColorBrewer} package for more.
#' @param ... Additional arguments for \code{heatmap.2} function from \code{gplots} package.
#' @param plotWidth Width of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param plotHeight Height of the plot. Unit is \code{inch}. Default is \code{7}.
#' @details Note that both \code{annot} and \code{genesymbolVar} need to be set to display gene sysmbols as row labels. Otherwise, the row labels will be probe names. Also note that when set to display gene symbols, all the probes without a gene symbol will be removed.
#' @return A supervised heatmap based on hierarchical clustering analysis in \code{pdf} format.
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' rbioarray_hcluster_super(plotName = "pre_experiVnaive", fltDOI = pre_experiVnaive_super, dfmDE = fltdata_DE$pre_experiVnaive, dataProbeVar = "ProbeName",
#'                          FC = 1.5, DE.sig.p = 0.0003461,
#'                          clust = "ward.D2",
#'                          fct = factor(c("naivepre", "naivepre", "naivepre", "naivepre", "exppre", "exppre", "exppre", "exppre", "exppre"),
#'                          levels = c("naivepre", "exppre")), trace = "none", srtCol = 30, offsetCol = 0.5, adjCol = c(1, 0),
#'                          rowLabel = TRUE, annot = annot, annotProbeVar = "ProbeName", genesymbolVar = "GeneSymbol",
#'                          offsetRow = 0, adjRow = c(0, 0.5), cexRow = 0.6,
#'                          key.title = "", keysize = 1.5,
#'                          key.xlab = "Normalized expression value", key.ylab = "Gene count")
#' }
#' @export
rbioarray_hcluster_super <- function(plotName = "data", fltDOI, dfmDE,
                                     ctrlProbe = TRUE, ctrlTypeVar = "ControlType",
                                     dataProbeVar = "ProbeName",
                                     DE.sig.method = "none",
                                     DE.sig.p = 0.05, FC = 1.5,
                                     fct, sampleName = NULL,
                                     colGroup = ifelse(length(levels(fct)) < 19, length(levels(fct)), 19),
                                     distance = "euclidean", clust = "complete",
                                     rowLabel = FALSE, annot = NULL, annotProbeVar = "ProbeName", genesymbolVar = NULL,
                                     colColour = "Paired", mapColour = "PRGn", n_mapColour = 11, ...,
                                     plotWidth = 7, plotHeight = 7){ #DOI: fltered subset data of interest
  ## argument check
  if (!dataProbeVar %in% names(dfmDE) | !dataProbeVar %in% names(fltDOI$genes)){
    stop(cat("data probe variable not found in dfmDE or fltDOI"))
  }

  ## prepare matrix for plotting
  dfm <- data.frame(fltDOI$genes, fltDOI$E)

  if (ctrlProbe){
    if (ctrlTypeVar %in% names(dfm)){
      dfm <- dfm[dfm[, ctrlTypeVar] == 0, ] # remove control probes
    } else {
      print(cat("Control probe setting on, but invalid control type variable detected, proceed without removing control probes."))
    }
  }


  if (tolower(DE.sig.method) == "fdr"){
    if (length(which(dfmDE$adj.P.Val < DE.sig.p)) == 0){
      warning("No significant results found using FDR correction. Please consider using another thresholding method. For now, alpha is applied on raw p.values.")
      pcutoff <- DE.sig.p
      pb_name <- dfmDE[dfmDE$P.Value < pcutoff, dataProbeVar]
      dfm <- dfm[dfm[, dataProbeVar] %in% pb_name, ]
    } else {
      pcutoff <- max(dfmDE[dfmDE$adj.P.Val < DE.sig.p, ]$P.Value)
      pb_name <- dfmDE[dfmDE$P.Value < pcutoff, dataProbeVar]
      dfm <- dfm[dfm[, dataProbeVar] %in% pb_name, ]
    }
  } else {
    pcutoff <- DE.sig.p
    pb_name <- dfmDE[dfmDE$P.Value < pcutoff, dataProbeVar]
    dfm <- dfm[dfm[, dataProbeVar] %in% pb_name, ]
  }

  ## set FC filter
  if (!is.null(FC)){
    pb_name_fc <- dfmDE[abs(dfmDE$logFC) >= log2(FC), dataProbeVar]
    dfm <- dfm[dfm[, dataProbeVar] %in% pb_name_fc, ]
  }

  ## check the dim
  if (dim(dfm)[1] < 2){
    stop("Only one row left after DE filtering. Nothing to cluster.")
  }

  ## heatmap
  # set up dis and cluster functions
  distfunc <- function(x)dist(x, method = distance)
  clustfunc <- function(x)hclust(x, method = clust)

  # set ColSideColors
  col_cluster <- clustfunc(distfunc(t(dfm[, -c(1:(ncol(dfm) - ncol(fltDOI$E)))])))
  colG <- cutree(col_cluster, colGroup) # column group
  colC <- brewer.pal(ifelse(colGroup < 3, 3, colGroup), colColour) # column colour

  # draw heatmap
  ogNcol <- dim(fltDOI$E)[2] # original numbers of col
  annoNcol <- dim(dfm)[2] # numbers of col with annotation
  s <- (annoNcol - ogNcol + 1):annoNcol # extract only the data by removing the annotation columns

  if (rowLabel){
    if (is.null(annot) | is.null(genesymbolVar)){
      warning("No annotation object or gene sybmol variable detected. Row labels will be the default probe names.")

      # retrive correct data matrix without annoation
      mtx <- as.matrix(dfm[, s])
      rownames(mtx) <- dfm[, dataProbeVar]

      if (!is.null(sampleName)){
        if (length(sampleName) != ncol(mtx)){
          stop("sampleName variable isn't the same length as the samples.")
        } else {
          colnames(mtx) <- sampleName
        }
      }

      pdf(file = paste(plotName, "_heatmap.supervised.pdf", sep = ""), width = plotWidth, height = plotHeight)
      heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
                col = brewer.pal(n_mapColour, mapColour), ColSideColors = colC[colG], ...)
      dev.off()

    } else {
      geneSymbl <- annot[annot[, annotProbeVar] %in% dfm[, dataProbeVar], ][, genesymbolVar]

      # process dfm further to only retain probes with gene symbol
      dfm$geneSymbol <- geneSymbl
      dfm <- dfm[complete.cases(dfm), ] # remove probes withnout a gene symbol
      labrow <- dfm$geneSymbol

      # retrive correct data matrix without annoation
      mtx <- as.matrix(dfm[, s])
      rownames(mtx) <- dfm[, dataProbeVar]

      if (!is.null(sampleName)){
        if (length(sampleName) != ncol(mtx)){
          stop("sampleName variable isn't the same length as the samples.")
        } else {
          colnames(mtx) <- sampleName
        }
      }

      pdf(file = paste(plotName, "_heatmap.supervised.pdf", sep = ""), width = plotWidth, height = plotHeight)
      heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
                col = brewer.pal(n_mapColour, mapColour), ColSideColors = colC[colG], labRow = labrow, ...)
      dev.off()
      cat("Probes with no gene names are removed.")
    }

  } else {
    mtx <- as.matrix(dfm[, s])
    rownames(mtx) <- dfm[, dataProbeVar]

    if (!is.null(sampleName)){
      if (length(sampleName) != ncol(mtx)){
        stop("sampleName variable isn't the same length as the samples.")
      } else {
        colnames(mtx) <- sampleName
      }
    }

    pdf(file = paste(plotName, "_heatmap.supervised.pdf", sep = ""), width = plotWidth, height = plotHeight)
    heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
              col = brewer.pal(n_mapColour, mapColour), ColSideColors = colC[colG], labRow = FALSE,...)
    dev.off()
  }
}
