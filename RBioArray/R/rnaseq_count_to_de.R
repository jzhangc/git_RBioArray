#' @title rnaseq_de
#'
#' @description Generic RNAseq statistical analysis function
#' @param object Object containing expression values. Currently the function supports \code{rbioseq_count} and \code{mir_count} from \code{RBioMIR} pacakge.
#' @param ... Additional arguments for corresponding S3 class methods.
#' @return A differential expression result object.
#' @examples
#'
#' \dontrun{
#' rnaseq_de(object = mrnaseq, design = design, contra = contra)
#' }
#'
#' @export
rnaseq_de <- function(object, ...){
  ## check object
  if (!class(object) %in% c("rbioseq_count", "mir_count")) stop("object needs to be either a \"rbioseq_count\" or \"mir_count\" object")

  ## use methods
  UseMethod("rnaseq_de", object)
}


#' @title rnaseq_de.mir_count
#'
#' @rdname rnaseq_de
#' @method rnaseq_de mir_count
#' @param object A \code{mir_count} object from the \code{mirProcess} function of \code{RBioMIR} package.
#' @param filter.threshold.min.count Minimum count for the smallest library for filter thresholding. Default is \code{10}.
#' @param filter.threshold.min.sample Minimum number of samples meeting the count threshold. Default is \code{NULL}.
#' @param ... Additional arguments for \code{\link{rnaseq_de.defuault}}.
#'
#' @export
rnaseq_de.mir_count <- function(object, filter.threshold.min.count = 10, filter.threshold.min.sample = NULL, ...){
  ## for setting up comparison groups info and minimum sample number
  annot.group <- object$sample_groups

  ## contruct rbioseq_de object
  out <- rnaseq_de.default(x = object$raw_read_count, y = object$genes,
                           filter.threshold.cpm = filter.threshold.min.count * min(object$sample_library_sizes) / 1000000,
                           filter.threshold.min.sample = filter.threshold.min.sample,
                           annot.group = annot.group, ...)
  out <- append(out, object[c("targets", "sample_groups")])
  class(out) <- "rbioseq_de"
  return(out)
}


#' @title rnaseq_de.rbioseq_count
#'
#' @rdname rnaseq_de
#' @method rnaseq_de rbioseq_count
#' @param object A \code{rbioseq_count} object from \code{\code{rbioseq_import_count}} function.
#' @param filter.threshold.min.count Minimum count for the smallest library for filter thresholding. Default is \code{10}.
#' @param filter.threshold.min.sample Minimum number of samples meeting the count threshold. Default is \code{NULL}.
#' @param ... Additional arguments for \code{\link{rnaseq_de.defuault}}.
#'
#' @export
rnaseq_de.rbioseq_count <- function(object, filter.threshold.min.count = 10, filter.threshold.min.sample = NULL, ...){
  ## for setting up comparison groups info and minimum sample number
  annot.group <- object$sample_groups

  ## contruct rbioseq_de object
  out <- rnaseq_de.default(x = object$raw_read_count, y = object$genes,
                           y.gene_id.var.name = "gene_id", y.gene_symbol.var.name = "gene_name",
                           filter.threshold.cpm = filter.threshold.min.count * min(object$sample_library_sizes) / 1000000,
                           filter.threshold.min.sample = filter.threshold.min.sample,
                           annot.group = annot.group, ...)
  out <- append(out, object[c("targets", "sample_groups")])
  class(out) <- "rbioseq_de"
  return(out)
}


#' @title rnaseq_de.default
#'
#' @rdname rnaseq_de
#' @method rnaseq_de default
#' @param x Input read count matrix, with rows for genes and columns for RNA samples.
#' @param y Target (e.g. genes) annotation matrix or vector.
#' @param y.gene_id.var.name Variable name for gene (i.e. target) identification. Default is \code{"genes"}.
#' @param y.gene_symbol.var.name Variable name for gene (i.e. target) symbols. Default is \code{"genes"}.
#' @param filter.threshold.cpm Filtering threshold for counts based on CPM (counts per million). Default is \code{"none"}.
#' @param filter.threshold.min.sample Minimum number of samples meeting the count threshold. Default is \code{NULL}.
#' @param annot.group Ssample group annotation object. Can be a \code{factor} or \code{vector} object.
#' @param library.size.scale.method Library size scaling method. Options are: \code{"none"}, \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"}. Default is \code{"TMM"}.
#' @param design Design matrix.
#' @param contra Contrast matrix.
#' @param qc.plot QC plot for the input read counts
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @details \code{filter.threshold.cpm.count} uses CPM (counts per million) as the basis for filtering. The rule of thumb is to filter reads independently from the groupping information.
#'
#'          The default is based on the paper by Chen et al (2016):
#'          10~15 counts per library, which translate to 10~15/L with L being the smallest library size in millions.
#'          Therefore, if using 10 as an example, the above becomes 10 / min(library_sizes) / 1000000
#'
#'          Chen W, Lun ATL, Smyth GK. 2016. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments
#'          using Rsubread and the edgeR quasi-likelihood pipeline. F1000Research. 5:1438.
#'
#' @return A list containing core elements of a \code{rbioseq_de} class. The items of a \code{rbioseq_de} class are following:
#'
#'         \code{filter_results}: A list containing \code{filter_threshold_cpm}, \code{filter_threshold_min_sample}, \code{filter_summary} and \code{filtered_counts}.
#'                                Note: \code{filtered_counts} is a \code{DGEList} class from \code{edgeR} package.
#'
#'         \code{normalization_results}
#'
#'         \code{normalized_data}: A \code{EList} generated from \code{voomWithQualityWeights} function from \code{limma} package.
#'                                 Note: The \code{targets} is NOT the same as the \code{targets} outside.
#'
#'         \code{genes_annotation.gene_id.var_name}
#'
#'         \code{genes_annotation.gene_symbol.var_name}
#'
#'         \code{F_stats}
#'
#'         \code{DE_results}
#'
#'         \code{comparisons}: a list with comparisons and comparison levles
#'
#'         \code{targets}: sample annotation data frame
#'
#'         \code{sample_groups}
#'
#' @importFrom limma lmFit eBayes topTable contrasts.fit voomWithQualityWeights
#' @importFrom edgeR DGEList calcNormFactors
#'
#' @export
rnaseq_de.default <- function(x, y = NULL,
                              y.gene_id.var.name = "genes",
                              y.gene_symbol.var.name = "genes",
                              filter.threshold.cpm = "none",
                              filter.threshold.min.sample = NULL, annot.group = NULL,
                              library.size.scale.method = "TMM",
                              design, contra, qc.plot = TRUE, verbose = TRUE){
  ## check the key arguments
  if (!class(x) %in% c("data.frame", "matrix")){
    stop("x has to be either a data.frame or matrix object")
  } else {
    x <- as.matrix(x)
  }

  if (is.null(y)){
    cat("y not supplied, using x row numbers as gene identifiers")
    y <- seq(nrow(x))
  } else {
    if (is.null(dim(y))){
      if (length(y) != nrow(x)){
        stop("Read count matrix x doesn't have the same row number as the target annotation vector y.")
      }
    } else {
      if (nrow(y) != nrow(x)){
        stop("Read count matrix x doesn't have the same row number as the target annotation matrix y.")
      }
    }
  }

  if (!is.null(dim(y))) {  # check if the variables are in the gene annotation matix
    if (!y.gene_id.var.name %in% colnames(y) || !y.gene_symbol.var.name %in% colnames(y)){
      stop("y.gene_id.var.name or y.gene_symbol.var.name can't be found in target annotation matrix y")
    }
  }

  if (is.null(filter.threshold.min.sample)) {
    if (is.null(annot.group)) {
      stop("When filter.threshold.min.sample = NULL, annot.group needs to be set and can only be a vector or factor object.")
    } else if (!is.factor(annot.group)) {
      if (!is.null(dim(annot.group))) {
        stop("annot.group needs to be a vector or factor object.")
      } else {
        annot.group <- factor(annot.group, levels = unique(annot.group))
      }
    }
    filter.threshold.min.sample <- min(table(annot.group))
  }

  if (is.null(design)){
    stop("Please set design matrix.")
  }

  if (is.null(contra)){
    stop("Please set contrast object.")
  }

  if (!library.size.scale.method %in% c("none", "TMM","RLE", "upperquartile")){
    stop("Please set the library.size.scale.method with exactly one of the following: \"none\", \"TMM\", \"RLE\", or \"upperquartile\".")
  }

  ## extract coefficients
  cf <- colnames(contra) # extract coefficient
  contra_levels <- vector(mode = "list", length = length(cf))
  contra_levels[] <- foreach(i = seq(length(cf))) %do% {
    rownames(contra)[which(contra[, i] != 0)]
  }
  names(contra_levels) <- cf

  ## DE
  if (verbose) cat("Data filtering and normalization...") # message
  dge <- DGEList(counts = x, genes = y)

  if (filter.threshold.cpm != "none"){ # set the count threshold for filtering
    isexpr <- rowSums(cpm(dge$counts) > filter.threshold.cpm) >= filter.threshold.min.sample  # cpm threshold, cite the paper
    dge <- dge[isexpr, , keep.lib.size = FALSE] # filtering
    flt_summary <- as.numeric(table(isexpr))
    names(flt_summary) <- c("filtered", "remaning")
    filter_results <- list(filter_threshold_cpm = filter.threshold.cpm,
                           filter_threshold_min_sample = filter.threshold.min.sample,
                           filter_summary = flt_summary,
                           filtered_counts = dge)
  } else {
    filter_results <- NULL
  }

  # library size scaling normalization
  dgenormf <- calcNormFactors(dge, method = library.size.scale.method)  # between-genes
  # between-genes: Voom normalization with quality weights
  # vmwt <- voomWithQualityWeights(dgenormf, design = design, plot = qc.plot, normalization = "quantile") # old
  vmwt <- voomWithQualityWeights(dgenormf, design = design, plot = qc.plot, normalize.method = "quantile")
  if (verbose) cat("DONE!\n") # message

  # fitting
  if (verbose) cat("Linear fitting...") # message
  fit <- lmFit(vmwt, design = design) # linear fitting
  fit <- contrasts.fit(fit, contrasts = contra)
  fit <- eBayes(fit)
  if (verbose) cat("DONE!\n") # message

  ## output
  # below: set sort.by = "none" for both f_stats and de_list to preserve the order, for supervsied clustering analysis.
  # specifically, it is for constructing threshoding vector used for subsetting E matrix,
  f_stats <- topTable(fit, number = Inf, sort.by = "none")
  de_list <- vector(mode = "list", length(cf))
  de_list[] <- foreach(i = seq(length(cf))) %do% topTable(fit = fit, coef = cf[i], number = Inf, sort.by = "none")
  names(de_list) <- cf
  comparisons <- list(comparisons = cf, comparison_levels = contra_levels)

  out <- list(filter_results = filter_results,
              normalization_method = list(library_scalling = library.size.scale.method, between_genes = "voom process with quantile normalization"),
              normalized_data = vmwt,
              genes_annotation.gene_id.var_name = y.gene_id.var.name,
              genes_annotation.gene_symbol.var_name = y.gene_symbol.var.name,
              F_stats = f_stats,
              DE_results = de_list,
              comparisons = comparisons)
  return(out)
}

#' @export
print.rbioseq_de <- function(x, ...){
  cat("\n")
  cat("--- RNA-seq gene differential expression analysis ---\n")
  cat("\n")
  cat("Gene filtering summary: \n")
  cat(paste0("\tCPM threshold: ", x$filter_results$filter_threshold_cpm, "\n"))
  cat(paste0("\tMinimum sample group threshold: ", x$filter_results$filter_threshold_min_sample, "\n"))
  cat(paste0("\tFiltered genes: ", x$filter_results$filter_summary[1], "\n"))
  cat(paste0("\tRemaining genes: ", x$filter_results$filter_summary[2], "\n"))
  cat("\n")
  cat("Reads normalization methods: \n")
  cat(paste0("\tLibrary size-scaling: ", x$normalization_method$library_scalling, "\n"))
  cat(paste0("\tBetween-genes: ", x$normalization_method$between_genes, "\n"))
  cat("\n")
  cat("Comparisons assessed: \n")
  cat(paste0("\t", x$comparisons$comparisons, "\n"))
  cat("\n")
}
