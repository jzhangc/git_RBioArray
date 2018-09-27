#' Title rnaseq_de
#'
#' @description Generic RNAseq differential expression analysis function
#' @param object Object containing expression values. Currently the function supports \code{rbioseq_count} and \code{mir_count} from \code{RBioMIR} pacakge.
#' @param ... Additional arguments for corresponding S3 class methods.
#' @return A differential expression result object.
#' @examples
#'
#' \dontrun{
#' rnaseq_de(object = mrnaseq)
#' }
#'
#' @export
rnaseq_de <- function(object, ...){
  ## check object
  if (!class(object) %in% c("rbioseq_count", "mir_count")) stop("object needs to be either a \"rbioseq_count\" or \"mir_count\" object")

  ## use methods
  UseMethod("rnaseq_de", object)
}


#' Title rnaseq_de.mir_count
#'
#' @description The \code{rnaseq_de} function for \code{mir_count} object from \code{RBioMIR} object.
#' @param object A \code{mir_count} object from the \code{mirProcess} function of \code{RBioMIR} package.
#' @param filter.threshold.min.count Minimum count for the smallest library for filter thresholding. Default is \code{10}.
#' @param ... Additional arguments for \code{\link{rnaseq_de.defuault}}.
#'
#' @export
rnaseq_de.mir_count <- function(object, filter.threshold.min.count = 10, ...){
  out <- rnaseq_de.default(x = object$raw_read_count, y = object$genes,
                           filter.threshold.cpm = filter.threshold.min.count * min(object$sample_library_sizes) / 1000000,
                           ...)
  return(out)
}


#' Title rnaseq_de.rbioseq_count
#'
#' @description The \code{rnaseq_de} function for \code{rbioseq_count} object from \code{\code{rbioseq_import_count}} function.
#' @param object A \code{rbioseq_count} object from \code{\code{rbioseq_import_count}} function.
#' @param filter.threshold.min.count Minimum count for the smallest library for filter thresholding. Default is \code{10}.
#' @param ... Additional arguments for \code{\link{rnaseq_de.defuault}}.
#'
#' @export
rnaseq_de.rbioseq_count <- function(object, filter.threshold.min.count = 10, ...){
  out <- rnaseq_de.default(x = object$raw_read_count, y = object$genes,
                           y.gene_id.var.name = "gene_id", y.gene_symbol.var.name = "gene_name",
                           filter.threshold.cpm = filter.threshold.min.count * min(object$sample_library_sizes) / 1000000,
                           ...)
  return(out)
}


#' Title rnaseq_de.default
#'
#' @description The default \code{rnaseq_de} function
#' @param x Input read count matrix, with rows for genes and columns for RNA samples
#' @param y Target (e.g. genes) annotation matrix or vector.
#' @param y.gene_id.var.name Variable name for gene (i.e. target) identification. Default is \code{"genes"}.
#' @param y.gene_symbol.var.name Variable name for gene (i.e. target) symbols. Default is \code{"genes"}.
#' @param filter.threshold.cpm Filtering threshold for counts based on CPM (counts per million). Default is \code{"none"}.
#' @param filter.threshold.min.sample Minimum number of samples meeting the count threshold. Default is \code{NULL}.
#' @param annot.group Set only when \code{filter.threshold.min.sample} is not \code{NULL}, sample group annotation object. Can be a \code{factor} or \code{vector} object.
#' @param between.genes.norm.method Between gene normalization method. Options are: \code{"none"}, \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"}. Default is \code{"TMM"}.
#' @param design Design matrix.
#' @param contra Contrast matrix.
#' @param qc.plot QC plot for the input read counts
#' @details \code{filter.threshold.cpm.count} uses CPM (counts per million) as the basis for filtering. The rule of thumb is to filter reads independently from the groupping information.
#'
#'          The default is based on the paper by Chen et al (2016):
#'          10~15 counts per library, which translate to 10~15/L with L being the smallest library size in millions.
#'          Therefore, if using 10 as an example, the above becomes 10 / min(library_sizes) / 1000000
#'
#'          Chen W, Lun ATL, Smyth GK. 2016. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments
#'          using Rsubread and the edgeR quasi-likelihood pipeline. F1000Research. 5:1438.
#'
#' @return A \code{rbioseq_de} object. The items of the object are following:
#'
#'         \code{filter_results}
#'
#'         \code{normalization_results}
#'
#'         \code{gene_id_var_name}
#'
#'         \code{gene_symbol_var_name}
#'
#'         \code{F_stats}
#'
#'         \code{DE_results}
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
                              between.genes.norm.method = "TMM",
                              design, contra, qc.plot = TRUE){
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

  if (!between.genes.norm.method %in% c("none", "TMM","RLE","upperquartile")){
    stop("Please set the norm.method with exactly one of the following: \"none\", \"TMM\", \"RLE\", or \"upperquartile\".")
  }

  ## extract coefficients
  cf <- colnames(contra) # extract coefficient

  ## DE
  cat("Data filtering and normalization...") # message
  dge <- DGEList(counts = x, genes = y)

  if (filter.threshold.cpm != "none"){ # set the count threshold for filtering
    isexpr <- rowSums(cpm(dge$counts) > filter.threshold.cpm) >= filter.threshold.min.sample  # cpm threshold, cite the paper
    dge <- dge[isexpr, , keep.lib.size = FALSE] # filtering
    flt_summary <- as.numeric(table(isexpr))
    names(flt_summary) <- c("filtered", "remaning")
    filter_results <- list(filter_threshold_cpm = filter.threshold.cpm,
                           filter_threshold_min_sample = filter.threshold.min.sample,
                           filter_summary = flt_summary)
  } else {
    filter_results <- NULL
  }

  # between-genes normalization
  dgenormf <- calcNormFactors(dge, method = between.genes.norm.method)  # between-genes
  # between-samples: Voom normalization with quality weights
  vmwt <- voomWithQualityWeights(dgenormf, design = design, plot = qc.plot, normalization = "quantile")
  cat("DONE!\n") # message

  # fitting
  cat("Linear fitting...") # message
  fit <- lmFit(vmwt, design = design) # linear fitting
  fit <- contrasts.fit(fit, contrasts = contra)
  fit <- eBayes(fit)
  cat("DONE!\n") # message

  ## output
  f_stats <- topTable(fit, number = Inf)
  de_list <- vector(mode = "list", length(cf))
  de_list[] <- foreach(i = seq(length(cf))) %do% topTable(fit = fit, coef = cf[i], number = Inf)
  names(de_list) <- cf

  out <- list(filter_results = filter_results,
              normalization_method = list(between_genes = between.genes.norm.method, between_samples = "voom"),
              gene_id_var_name = y.gene_id.var.name,
              gene_symbol_var_name = y.gene_symbol.var.name,
              F_stats = f_stats,
              DE_results = de_list)
  class(out) <- "rbioseq_de"
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
  cat(paste0("\tBetween-genes: ", x$normalization_method$between_genes, "\n"))
  cat(paste0("\tBetween-samples: ", x$normalization_method$between_samples, "\n"))
  cat("\n")
  cat("Comparisons assessed: \n")
  cat(paste0("\t", names(x$DE_results), "\n"))
  cat("\n")
}
