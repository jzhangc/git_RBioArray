#' @title rbioseq_de_analysis
#'
#' @description All-in-one wrapper for RNAseq differential expression (DE) analysis, i.e. from read count files to volcano plots.
#' @param raw.file.path Path to raw files. Default is the system working directory.
#' @param species Optional species code, following the traditional abbreviated naming convention, e.g. "hsa", "mmu".
#' @param target.annot.file Annotation file describing filenames and targets, and should be in \code{csv} format.
#' @param sample_groups.var.name Sample group annotation variable name in the \code{target.annot.file}.
#' @param gtf Parsed gtf/gff annotation matirx. Can be obtained by function \code{\link{rbioseq_import_gtf}}.
#' @param raw.file.ext Raw file extention. Default is \code{".txt"}.
#' @param raw.file.sep Raw read count file separators. Default is \code{""\"\"}, i.e. white space.
#' @param raw.file.source Raw file source, i.e. program used to generate read counts. Currently only supports \code{"htseq-count"}.
#' @param filter.threshold.cpm Filtering threshold for counts based on CPM (counts per million). Default is \code{"none"}.
#' @param filter.threshold.min.sample Minimum number of samples meeting the count threshold. Default is \code{NULL}.
#' @param design Design matrix.
#' @param contra Contrast matrix.
#' @param FC Threshold for fold change. Default is \code{1.5}.
#' @param alpha Threshold for p values. Default is \code{0.05}.
#' @param FDR If to use FDR correction for p values. Default is \code{TRUE}.
#' @param export.name Name used for output objects to the environment and directory. Not optional. Default is \code{NULL}.
#' @param export.mode Mode used to export results to the directory. Options are \code{"all"}, \code{"all.gene_symbol"}, \code{"sig"}. Default is \code{"all"}. See details.
#' @param ... Additional arguments for \code{\link{sig}} function.
#' @param parallelComputing Wether to use parallel computing or not. Default is \code{TRUE}.
#' @param cluterType clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @details When \code{raw.file.source = "htseq-count"}, the function will cut off the last five summary raws.
#'
#'          For more on the classes \code{rbioseq_count}, \code{rbioseq_de} and \code{sig}, see the help pages for functions \code{\link{rbioseq_import_count}}, \code{\link{rnaseq_de}} and \code{\link{sig}}, respectively.
#'
#'          The function is not suitable for small RNA RNAseq analysis, e.g. miRNA analysis. For miRNA analysis,
#'          it is advised to use \code{RBioMIR} pacakge in conjunction with functions \code{\link{rnaseq_de}} and \code{\link{sig}} from the current pacakge.
#'
#' @return Outputs the following objects:
#'
#'         \code{rbioseq_count}: non small RNA RNA-seq read count object
#'
#'         \code{rbioseq_de}: RNAseq DE results object
#'
#'         \code{sig}: RNAseq DE significant analysis results object.
#'
#'         The function also outputs \code{csv} files containing
#'
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @examples
#' \dontrun{
#' mrna_count <- rbioseq_import_count(path = "~/dataset/",
#'                                    species = "hsa",
#'                                    target.annot.file = "target.csv", sample_groups.var.name = "condition",
#'                                    gtf.matrix = gtf,
#'                                    raw.file.ext = ".out", raw.file.sep = "",
#'                                    raw.file.source = "htseq-count",
#'                                    parallelComputing = TRUE, clusterType = "FORK")
#' }
#' @export
rbioseq_de_analysis <- function(raw.file.path, raw.file.ext = ".txt", raw.file.sep = "", raw.file.source = "htseq-count",
                                species,
                                target.annot.file, sample_groups.var.name = NULL,
                                gtf,
                                filter.threshold.min.count = 10,
                                filter.threshold.min.sample = NULL,
                                design, contra,
                                alpha = 0.05, FC = 1.5, FDR = TRUE,
                                export.name = "data", export.mode = "all",
                                ...,
                                parallelComputing = FALSE, clusterType = "FORK", verbose = TRUE){
  ## import raw reads
  # import
  count <- rbioseq_import_count(htseq_file_dir = raw.file.path,
                                htseq_file.ext = raw.file.ext,
                                htseq_file.sep = raw.file.sep,
                                count_data_type = raw.file.source,
                                species = species,
                                htseq_sample.annot.file = target.annot.file,
                                htseq_sample.annot.group.var = sample_groups.var.name,
                                gtf = gtf,
                                parallelComputing = parallelComputing, clusterType = clusterType,
                                verbose = verbose)
  # export
  assign(paste0(export.name, "_count"), count, envir = .GlobalEnv)
  if (verbose) cat("\n\n")

  ## DE analysis
  # DE
  if (is.null(filter.threshold.min.sample)) {
    annot.group <- count$sample_groups
  } else {
    annot.group = NULL
  }
  de <- rnaseq_de(object = count, design = design, contra = contra,
                  filter.threshold.min.count = 5,
                  filter.threshold.min.sample = filter.threshold.min.sample,
                  annot.group = annot.group, verbose = verbose)
  # export
  assign(paste0(export.name, "_de"), de, envir = .GlobalEnv)
  if (verbose) cat("\n\n")

  ## sig analysis
  sig <- sig(object = de, alpha = alpha, FC = FC, FDR = FDR, export.name = export.name, export.mode = export.mode, ...,
             verbose = verbose)

  # export
  assign(paste0(export.name, "_sig"), sig, envir = .GlobalEnv)
}
