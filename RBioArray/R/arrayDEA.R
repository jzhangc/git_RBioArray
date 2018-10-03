#' @title rbioarray_rlist
#'
#' @description Function to constuct \code{rbioarray_rlist} class object from microarary and annotation data. The \code{rbioarray_rlist} object is the starting point for all microarray data analysis.
#' @param raw.dataframe Input data frame containing microarray hybridization signals with rows as probe/gene/genomic features and columns as samples. Note: the data frame should contain at least one annotation column.
#' @param raw.annot.var.name A string vector containing variable (i.e. column) name(s) for all the annotation columns in \code{raw.dataframe}.
#' @param raw.gene_id.var.name Variable (i.e. column) name for gene/probe/genomic feature identification from \code{raw.dataframe}.
#' @param gene.annot.dataframe Optional annotation data frame for gene/probe/genomic feature annotation.
#' @param gene.annot.gene_id.var.name Set only when \code{gene.annot.dataframe} is provided, variable name for probe/gene/genomic features identification from \code{gene.annot.dataframe}.
#' @param gene.annot.gene_symbol.var.name Set only when \code{gene.annot.dataframe} is provided, variable name for probe/gene/genomic features display name from \code{gene.annot.dataframe}, e.g. gene symbols.
#' @param gene.annot.control_type.var.name Optional and set only when \code{gene.annot.dataframe} is provided, name for the variable that denotes if the gene is a control and its control type. Default is \code{NULL}.
#' @param gene.annot.control_type.val.pos Set only when \code{gene.annot.control_type.var.name} is provided, value indicating a positive control probe. Default is \code{1}.
#' @param gene.annot.control_type.val.neg Set only when \code{gene.annot.control_type.var.name} is provided, value indicating a negative control probe. Default is \code{-1}.
#' @param gene.annot.control_type.val.exp Set only when \code{gene.annot.control_type.var.name} is provided, value indicating a none control probe. Default is \code{0}.
#' @param target.annot.file File name for the target (i.e. sample) annotation \code{.csv} file.
#' @param target.annot.file.path The directory for \code{target.annot.file}. Default is \code{getwd()}.
#' @param sample_groups.var.name The variable name for sample groupping information from \code{target.annot.file}.
#' @details The word "gene" used in argument names and output item names is in its broader meaning of gene/probe/genomic feature.
#'
#'          The optional annotation data frame \code{gene.annot.dataframe} should contain any additional annotation information in addition to the annotation column(s) from \code{raw.dataframe}.
#'          It is noted that \code{gene.annot.dataframe} should contain at least one column for the sample type of gene/probe/genomic feature identification as the identification variable from  \code{raw.dataframe}.
#'          Such column is set via argument \code{raw_file.gene_id.var_name} and \code{genes_annotation.gene_id.var_name}.
#'          The gene display name will only be used if \code{gene.annot.dataframe} is used. Otherwise, it wil use \code{raw.gene_id.var.name} from \code{raw.dataframe}.
#'
#'          To keep things consistent with the \code{Elist} from the dependent \code{limma} package. The \code{rbioarray_rlist} class contains many common elements from \code{Elist} class.
#'
#' @return A \code{rbioarray_rlist} object. The \code{rbioarray_rlist} class contains the following items:
#'
#'         \code{E}: raw expression (i.e. hybridization signal)
#'
#'         \code{raw_file.gene_annotation.var_name}
#'
#'         \code{raw_file.gene_id.var_name}
#'
#'         \code{gene_display_name_used}: if the gene display name is extracted from \code{gene.annot.dataframe}.
#'
#'         \code{genes_annotation.gene_id.var_name}
#'
#'         \code{genes_annotation.gene_symbol.var_name}
#'
#'         \code{gene_annotation.control_type}: if \code{gene.annot.control_type.var.name} is set, a list containing all the control type information
#'
#'         \code{targets}: the sample annotation data frame.
#'
#'         \code{sample_groups}
#'
#' @examples
#' \dontrun{
#' raw_list <- rbioarray_rlist(raw.dataframe = raw, raw.annot.var.name = c("PROBE_ID","ILMN_GENE"), raw.gene_id.var.name = "PROBE_ID",
#'                             gene.annot.dataframe = annot, gene.annot.gene_id.var.name = "PROBE_ID", gene.annot.gene_symbol.var.name = "SYMBOL",
#'                             target.annot.file = "sample_index_end6.csv", sample_groups.var.name = "time_point")
#' }
#' @export
rbioarray_rlist <- function(raw.dataframe, raw.annot.var.name = NULL, raw.gene_id.var.name = NULL,
                            gene.annot.dataframe = NULL, gene.annot.gene_id.var.name = NULL, gene.annot.gene_symbol.var.name = NULL,
                            gene.annot.control_type.var.name = NULL,
                            gene.annot.control_type.val.pos = 1, gene.annot.control_type.val.neg = -1, gene.annot.control_type.val.exp = 0,
                            target.annot.file = NULL, target.annot.file.path = getwd(), sample_groups.var.name = NULL){
  ## check arguments and set up variables
  # raw data
  if (is.null(dim(raw.dataframe)))
    stop("raw.matrix has to be a data.frame, with rows for genes/probes/genomic features, columns for samples.")
  if (is.null(raw.annot.var.name) || is.null(raw.gene_id.var.name))
    stop("Please set raw.annot.var.name AND raw.gene_id.var.name arguments.")
  if (!raw.annot.var.name %in% names(raw.dataframe)  || !raw.gene_id.var.name %in% names(raw.dataframe))
    stop("Annoation variables (i.e. columns) and/or the gene_id variable (i.e. column) not found in raw.dataframe")

  # gene annotation
  if (is.null(gene.annot.dataframe)){  # check and load gene annoation
    cat("Note: gene.annot.dataframe not provided. Proceed with raw.dataframe annoation information.\n")
    gene.annot.dataframe <- raw.dataframe[, raw.annot.var.name]
    gene.annot.gene_id.var.name <- raw.gene_id.var.name
    gene.annot.gene_symbol.var.name <- raw.gene_id.var.name
    gene.symbol <- FALSE
  } else {
    if (!is.data.frame(gene.annot.dataframe)) stop("gene.annot needs to be a dataframe")
    if (nrow(gene.annot.dataframe) < nrow(raw.dataframe)) stop("gene.annot.dataframe has less record than the raw.dataframe") # check size

    if (is.null(gene.annot.gene_id.var.name) || is.null(gene.annot.gene_symbol.var.name)) { # gene id and symbol variables
      stop("Please set gene.annot.gene_id.var.name AND gene.annot.gene_symbol.var.name arguments according to gene.annot")
    } else if (!gene.annot.gene_id.var.name %in% names(gene.annot.dataframe) || !gene.annot.gene_symbol.var.name %in% names(gene.annot.dataframe)) {
      stop("Gene ID variable and/or gene symbol variable not found in gene.annot.dataframe")
    } else {
      gene.symbol <- TRUE
    }

    if (!is.null(gene.annot.control_type.var.name)) {  # optional control_type variable
      if (!gene.annot.control_type.var.name %in% names(gene.annot.dataframe)) {
        cat("The set control type variable not found in gene.annot.dataframe. Proceed without using it.\n")
        gene.annot.control_type = NULL
      } else {
        control_type.values <- unique(gene.annot.dataframe[, gene.annot.control_type.var.name])
        input.control_type.values <- c(gene.annot.control_type.val.pos, gene.annot.control_type.val.neg, gene.annot.control_type.val.exp)
        if (all(input.control_type.values %in% control_type.values)){
          gene.annot.control_type = list(control_type.var_name = gene.annot.control_type.var.name,
                                         pos_type.value = gene.annot.control_type.val.pos,
                                         neg_type.value = gene.annot.control_type.val.neg,
                                         exp_type.value = gene.annot.control_type.val.exp)
        } else {
          cat("Input probe control type values not found in the gene control type variable from gene.annot.dataframe. Proceed without using them.\n")
          gene.annot.control_type = NULL
        }
      }
    } else {
      gene.annot.control_type = NULL
    }
  }

  # target annotation
  if (is.null(target.annot.file)){  # check and load target (sample) annotation
    stop("Please provide a target annotation file for target.annot.file arugment.")
  } else {
    target.annot_name_length <- length(unlist(strsplit(target.annot.file, "\\.")))
    target.annot_ext <- unlist(strsplit(target.annot.file, "\\."))[target.annot_name_length]
    if (target.annot_ext != "csv") {
      stop("target.annot.file is not in csv format.")
    } else {
      cat("Loading target annotation file...")
      tgt <- read.csv(file = paste0(target.annot.file.path, "/",target.annot.file), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
      cat("Done!\n")
    }
  }
  if (is.null(sample_groups.var.name)) stop("Please provide sample_groups.var.name.")
  if (!sample_groups.var.name %in% names(tgt)){
    stop("Sample group annotation variable not found in the target annotation file.")
  } else {
    sample.groups <- factor(tgt[, sample_groups.var.name], levels = unique(tgt[, sample_groups.var.name]))
  }

  # additional variables
  raw_dfm <- raw.dataframe
  genes_annot_dfm <- gene.annot.dataframe

  ## Set up the information
  cat("Constructing rlist...")
  # merge gene annotation with raw dataframe
  raw_dfm$merge_id <- raw_dfm[, raw.gene_id.var.name]
  genes_annot_dfm$merge_id <- genes_annot_dfm[, gene.annot.gene_id.var.name]
  merged_raw_gene_annot_dfm <- merge(raw_dfm, genes_annot_dfm)  # this merge will extract annotation info from gene_annot_dfm and merge to the smaller data dataframe.
  all_annot_var_names <- unique(c(raw.annot.var.name, names(genes_annot_dfm)))  # all annotation variable names

  # Set up output E matrix
  E <- as.matrix(merged_raw_gene_annot_dfm[, !names(merged_raw_gene_annot_dfm) %in% all_annot_var_names]) # remove annotation columns

  # set up output annotation information
  genes <- merged_raw_gene_annot_dfm[, names(merged_raw_gene_annot_dfm) %in% all_annot_var_names]

  ## output
  out <- list(E = E,
              raw_file.gene_annotation.var_name = raw.annot.var.name,
              raw_file.gene_id.var_name = raw.gene_id.var.name,
              genes = genes,
              gene_display_name_used = gene.symbol,
              genes_annotation.gene_id.var_name = gene.annot.gene_id.var.name,
              genes_annotation.gene_symbol.var_name = gene.annot.gene_symbol.var.name,
              gene_annotation.control_type = gene.annot.control_type,
              targets = tgt,
              sample_groups = sample.groups)
  class(out) <- "rbioarray_rlist"
  cat("Done!\n")
  cat("\n")
  cat(paste0("The resulted rbioarray_rlist object contains ", nrow(E), " genes/probes/genomic features, ", nrow(tgt), " samples for ", length(unique(sample.groups)), " groups."))
  return(out)
}


#' @export
print.rbioarray_rlist <- function(x, ...){
  cat("---- rbioarray_rlist information ----\n")
  cat(paste0("Number of genes/probes/genomic features: ", nrow(x$E), "\n"))
  cat(paste0("Number of samples: ", nrow(x$targets), "\n"))
  cat(paste0("Groups: "))
  cat(paste0(levels(x$sample_groups)))
  cat("\n\n")
  cat(paste0("Gene display name from annotation information: ", ifelse(x$gene_display_name_used, "Available\n", "Unavailable\n")))
  cat("\n")
}


#' Title rbioarray_transfo_normalize
#'
#' @description Generic data log transformation and nomalization function for microarray.
#' @param object Input obejct with raw data and annotation information. Should be a \code{rbioarray_rlist} class.
#' @param ... Additional arguments for corresponding S3 class methods.
#' @details The \code{rbioarray_rlist} object can be obtained from \code{\link{rbioarray_rlist}} function.
#'          The \code{Elist} and \code{MAList} classes are from the dependent \code{limma} package.
#' @return A \code{rbioarray_plist} class object with the following items:
#'
#'         \code{E}: Processed expression matrix
#'
#'         \code{design}: Microrray experiment sample design matrix
#'
#'         \code{ArrayWeight}
#'
#'         \code{extra_E_data}: A list with original (i.e. raw) expression matrix and, if applicable, log transformed expression matrix
#'
#'         Additionally, \code{rbioarray_plist} object will include additional items from the input \code{rbioarray_rlist} object.
#'
#' @examples
#'
#' \dontrun{
#' # for \code{rbioarray_rlist} object:
#' microarray_plist <- rbioarray_transfo_normalize(object = microarray_rlist, logTransfo = TRUE, logTransfo.method = "log2",
#'                                                 logTransfo.parallelComputing = TRUE, logTransfo.cluterType = "FORK",
#'                                                 bgc.method = "auto", between.sample.norm.method = "quantile")
#' }
#'
#' @export
rbioarray_transfo_normalize <- function(object, ...){
  ## check arguments
  if (class(object) != "rbioarray_rlist") stop("The input object needs to be either a \"rlist\", \"EList\" or \"MAList\"")

  ## use method
  UseMethod("rbioarray_transfo_normalize", object)
}


#' Title rbioarray_transfo_normalize.rbioarray_rlist
#'
#' @description \code{\link{rbioarray_transfo_normalize}} for \code{rbioarray_rlist} class object.
#' @param object Input obejct with raw data and annotation information. Could be \code{rbioarray_rlist}, \code{Elist} or \code{MAList} classes.
#' @param design Microarray experiment sample design matrix.
#' @param ... Additional arguments the default method \code{\link{rbioarray_transfo_normalize.default}}.
#' @details The \code{rbioarray_rlist} object can be obtained from \code{\link{rbioarray_rlist}} function.
#' @return A \code{rbioarray_plist} class object.
#'
#' @export
rbioarray_transfo_normalize.rbioarray_rlist <- function(object, design, ...){
  ## processing
  default_out <- rbioarray_transfo_normalize.default(E = object$E,
                                                     logTransfo.gene.annot = object$genes[, object$raw_file.gene_annotation.var_name],
                                                     logTransfo.export.name = deparse(substitute(object)),
                                                     between.sample.weight.design = design, ...)
  ## output
  cat("\n")
  cat("Constructing rbioarray_plist...")
  out <- append(default_out, object[!names(object) %in% "E"])
  class(out) <- "rbioarray_plist"
  cat("Done!\n\n")
  return(out)
}


#' Title rbioarray_transfo_normalize.default
#'
#' @description \code{\link{rbioarray_transfo_normalize}} for \code{rbioarray_rlist} class object.
#' @param E Input raw expression value matrix with columns for samples, rows for genes/probes/genomic features.
#' @param logTransfo If to perfom a log tranformation on the expresson data prior to background correction and beetween-sample normalization. Default is \code{FALSE}.
#' @param logTransfo.method Set only when \code{logTransfo = TRUE},
#' @param logTransfo.gene.annot Set only when \code{logTransfo = TRUE},
#' @param logTransfo.export.name Set only when \code{logTransfo = TRUE}, file name prefix for  the \code{csv} file containing log transformed data.
#' @param logTransfo.parallelComputing Set only when \code{logTransfo = TRUE}, if to use parallel computing for the transformation or not. Default is \code{FALSE}.
#' @param logTransfo.cluterType Set only when \code{logTransfo = TRUE} and \code{logTransfo.parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param bgc.method Background correction method. Default is \code{"auto"}. See \code{backgroundCorrect()} function from \code{limma} package for details.
#' @param between.sample.norm.method Normalization method. Default is \code{"quantile"}. See \code{normalizeBetweenArrays()} function from \code{limma} package for details.
#' @param between.sample.weight.design Microarray experiment sample design matrix for beteween sample weight calculation.
#' @param ... Additional arguments the default method \code{\link{rbioarray_transfo_normalize.default}}.
#' @details The \code{rbioarray_rlist} object can be obtained from \code{\link{rbioarray_rlist}} function.
#'          The word "gene" used in argument names and output item names is in its broader meaning of gene/probe/genomic feature.
#' @return A list with the core items for a \code{rbioarray_plist} class.
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom limma backgroundCorrect normalizeBetweenArrays backgroundCorrect.matrix arrayWeights
#' @export
rbioarray_transfo_normalize.default <- function(E,
                                                logTransfo = FALSE, logTransfo.method = "log2", logTransfo.gene.annot = NULL,
                                                logTransfo.export.name = "data",
                                                logTransfo.parallelComputing = FALSE, logTransfo.cluterType = "FORK",
                                                bgc.method = "auto", between.sample.norm.method = "quantile",
                                                between.sample.weight.design = NULL, ...){
  ## check arguments
  if (!bgc.method %in% c("auto", "auto", "none", "subtract", "half", "minimum", "movingmin", "edwards", "normexp"))
    stop("Agument bgc.method needs to be set with one of \"auto\", \"auto\", \"none\", \"subtract\", \"half\", \"minimum\", \"movingmin\", \"edwards\", and \"normexp\" exactly.")
  if (!between.sample.norm.method %in% c("quantile", "Aquantile")) stop("Argument between.sample.norm.method needs to be set with ")
  if (is.null(between.sample.weight.design)) stop("Please provide the microarray design matrix for between.sample.weight.design argument.")

  ## pre-processing
  ## log transform  or not
  transfo_E_mtx <- NULL
  if (logTransfo){  # log transform
    cat(paste0("Data ", logTransfo.method, " tranformation..."))
    # additional argument check
    if (!logTransfo.method %in% c("log2", "log10")) stop("Argument logTransfo.method needs to be set with either \"log2\" or \"log10\" exactly.")

    if (is.null(logTransfo.gene.annot)) {
      stop("Please provide logTransfo.gene.annot when logTransfo = TRUE.")
    } else {
      if (is.null(dim(logTransfo.gene.annot))){
        if (length(logTransfo.gene.annot) != nrow(E)) stop("logTransfo.gene.annot vector length not same as row number of E")
      } else {
        if (nrow(logTransfo.gene.annot) != nrow(E)) stop("Row number of logTransfo.gene.annot matrix/data frame not the same length as the row number of E")
      }
    }

    # log2 transformation
    if (!logTransfo.parallelComputing){
      transfo_E_mtx <- apply(E, c(1,2), FUN = ifelse(logTransfo.method == "log10", log10, log2))
    } else {  # parallel computing
      # set up cpu cluster
      n_cores <- detectCores() - 1
      cl <- makeCluster(n_cores, type = logTransfo.cluterType)
      registerDoParallel(cl) # part of doParallel package
      on.exit(stopCluster(cl)) # close connect when exiting the function
      # log transform
      transfo_E_mtx <- foreach(i = E, .combine = "c") %dopar% {
        out <- apply(i, c(1,2), FUN = ifelse(logTransfo.method == "log10", log10, log2))
        out
      }
      transfo_E_mtx <- matrix(data = transfo_E_mtx, nrow = nrow(E), ncol = ncol(E))
      colnames(transfo_E_mtx) <- colnames(E)
    }
    E_mtx <- transfo_E_mtx
    cat("Done!\n")

    # store and export log transformed data into a csv file if applicable
    cat(paste0("The log transformed data saved to file: ", logTransfo.export.name, "_log_transformed.csv\n"))
    logTransfo_out <- data.frame(logTransfo.gene.annot, transfo_E_mtx)
    write.csv(logTransfo_out, file = paste0(logTransfo.export.name, "_log_transformed.csv"), row.names = FALSE)
  } else {
    E_mtx <- E
  }

  ## normalization
  cat("\n")
  cat("Background correction: \n")
  BgC <- backgroundCorrect.matrix(E_mtx, method = bgc.method, ...) #background correction
  cat("\n")
  cat(paste0("Data normalization using ", between.sample.norm.method, " method..."))
  Norm <- suppressMessages(normalizeBetweenArrays(BgC, between.sample.norm.method)) # quantile normalization
  Wgt <- arrayWeights(Norm, design = between.sample.weight.design) # array weight
  cat("Done!\n")

  ## output
  default_out <- list(E = Norm,
                      design = between.sample.weight.design,
                      ArrayWeight = Wgt,
                      background_correction_method = bgc.method,
                      between_sample_normalization_method = between.sample.norm.method,
                      extra_E_data = list(log_transformation = logTransfo,
                                          original_E = E,
                                          log_transfo_E = transfo_E_mtx))
  return(default_out)
}


#' @export
print.rbioarray_plist <- function(x, ...){
  cat("---- rbioarray_plist information ----\n")
  cat(paste0("Log transformation: ", ifelse(x$extra_E_data$log_transformation, TRUE, FALSE), "\n"))
  cat(paste0("Background correction method: ", x$background_correction_method, "\n"))
  cat(paste0("Between-sample normalization method: ", x$between_sample_normalization_method))
  cat("\n\n")
  cat(paste0("Number of genes/probes/genomic features: ", nrow(x$E), "\n"))
  cat(paste0("Number of samples: ", nrow(x$targets), "\n"))
  cat(paste0("Groups: "))
  cat(paste0(levels(x$sample_groups)))
  cat("\n")
}


#' Title rbioarray_filter_combine
#'
#' @description Function to filter, averaging and (if set) combine genes/probes/genomic features from the \code{rbioarray_plist} objects.
#' @param object Input \code{rbioarray_plist} object from function \code{\link{rbioarray_transfo_normalize}}.
#' @param filter.percentile The percentile threshold for filtering. Default is \code{0.05}. See details for more.
#' @param filter.threshold.min.sample Minimum number of samples meeting the filtering threshold. Default is \code{NULL}. See details for more.
#' @param combine.gene.duplicate If to combine different transcripts from the same gene/genomic feature. Default is \code{FALSE}. See details for more.
#' @param parallelComputing Wether to use parallel computing or not. Default is \code{TRUE}.
#' @param cluterType clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @details For \code{filter.percentile},
#'
#'          For \code{filter.threshold.min.sample}
#'
#'          When \code{combine.gene.duplicate = TRUE}
#'
#' @return A \code{rbioarray_flist} class, including the following core items:
#'
#'         \code{E}
#'
#'         \code{genes}
#'
#'         \code{gene_duplicates_combined}: either \code{TRUE} or \code{FALSE}
#'
#'         \code{filter_results}: a list containing \code{neg_control_used}, \code{filter_percentile}, \code{filter_threshold_min_sample}, and \code{filter_summary}
#'
#'         Additionally, \code{rbioarray_flist} object will include additional items from the input \code{rbioarray_plist} object.
#'
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom limma avereps
#' @export
rbioarray_filter_combine <- function(object,
                                     filter.percentile = 0.05,
                                     filter.threshold.min.sample = NULL,
                                     combine.gene.duplicate = FALSE,
                                     parallelComputing = FALSE, clusterType = "PSOCK"){
  ## check arguments
  if (!"rbioarray_plist" %in% class(object)) stop("The input object needs to be, but not exclusive to, rbioarray_plist class.")

  ## filter
  cat("Filtering low expresson genes/probes/genomic features...")
  # set negative/low expression threshold
  if (!is.null(object$gene_annotation.control_type)){ # 0.95 percentile of all negative control values
    cat("Gene control type variable detected in the input object. The filter.percentile argument value reset to 0.95.\n")
    filter.percentile = 0.95
    control_type.var.name = object$gene_annotation.control_type$genes_annotation.control_type.var_name
    control_type.neg.value = object$gene_annotation.control_type$gene.annot.control_type.val.neg

    if (class(object$E[object$genes[, control_type.var.name] == control_type.neg.value, ]) == "numeric"){ # if there is only one entry in the neg values
      neg <- object$E[object$genes[, control_type.var.name] == control_type.neg.value, ] # no 95% percentile required as only one neg entry
    } else {
      neg <- apply(object$E[object$genes[, control_type.var.name] == control_type.neg.value, ], 2, function(x)quantile(x, p = filter.percentile)) # neg95
    }
    neg_control_used = TRUE
  } else {
    neg <- apply(object$E, 2, function(x)quantile(x, p = filter.percentile)) # 5% percentile of all the data
    neg_control_used = FALSE
  }

  # low expression cuttoff set at at least 10% hihger than the neg
  LE_cutoff <- matrix(1.1 * neg, nrow(object$E), ncol(object$E), byrow = TRUE)
  # set filtering matrix
  if (is.null(filter.threshold.min.sample)) {
    filter.threshold.min.sample <- min(table(object$sample_groups))
  }
  isexpr <- rowSums(object$E > LE_cutoff) >= filter.threshold.min.sample
  flt_summary <- as.numeric(table(isexpr))
  names(flt_summary) <- c("filtered", "remaning")

  # filter
  flt_E <- object$E[isexpr, ] # this is a way of extracting samples logically considered TRUE by certain tests
  flt_genes <- object$genes[isexpr, ]
  cat("Done!\n")

  ## averaging technical replicates
  cat("Averaging technical replicates...")
  flt_E_avg <- avereps(flt_E, ID = flt_genes[, object$genes_annotation.gene_id.var_name])
  flt_genes_avg <- unique(flt_genes[flt_genes[, object$genes_annotation.gene_id.var_name] %in% rownames(flt_E_avg), ])
  cat("Done!\n")

  ## combine duplicate genes if set
  if (combine.gene.duplicate) {
    if (length(flt_genes_avg[, object$genes_annotation.gene_id.var_name]) != length(unique(flt_genes_avg[, object$genes_annotation.gene_symbol.var_name]))) {
      if (object$gene_display_name_used){
        cat("Combining gene duplicates (i.e. different transcripts belonging to the same gene)...")
        if (!parallelComputing) {
          combGeneProbe <- foreach(i = unique(flt_genes_avg[, object$genes_annotation.gene_symbol.var_name]), .combine = "c") %do% {
            tmp <- flt_E_avg[which(flt_genes_avg[,object$genes_annotation.gene_symbol.var_name] %in% i), ]
            if (is.null(dim(tmp))){
              out <- var(tmp)
            } else {
              out <- foreach(j = seq(nrow(tmp)), .combine = "c") %do% {
                var(as.numeric(tmp[j, ]))
              }
            }
            out <- data.frame(flt_genes_avg[which(flt_genes_avg[, object$genes_annotation.gene_symbol.var_name] %in% i), ], var = out)
            out[which(out$var %in% max(out$var)), object$genes_annotation.gene_id.var_name]
          }
        } else {
          n_cores <- detectCores() - 1
          cl <- makeCluster(n_cores, type = clusterType, outfile = "")
          registerDoParallel(cl) # part of doParallel package
          on.exit(stopCluster(cl)) # close connect when exiting the function

          combGeneProbe <- foreach(i = unique(flt_genes_avg[, object$genes_annotation.gene_symbol.var_name]), .combine = "c", .packages = "foreach") %dopar% {
            tmp <- flt_E_avg[which(flt_genes_avg[,object$genes_annotation.gene_symbol.var_name] %in% i), ]
            if (is.null(dim(tmp))){
              out <- var(tmp)
            } else {
              out <- foreach(j = seq(nrow(tmp)), .combine = "c") %do% {
                var(as.numeric(tmp[j, ]))
              }
            }
            out <- data.frame(flt_genes_avg[which(flt_genes_avg[, object$genes_annotation.gene_symbol.var_name] %in% i), ], var = out)
            out[which(out$var %in% max(out$var)), object$genes_annotation.gene_id.var_name]
          }
        }
        out_E <- flt_E_avg[flt_genes_avg[, object$genes_annotation.gene_id.var_name] %in% combGeneProbe, ]
        out_genes <- flt_genes_avg[flt_genes_avg[, object$genes_annotation.gene_id.var_name] %in% combGeneProbe, ]
        cat("Done!\n")
      } else {
        cat("\n")
        cat("Gene symbol not used in the input object. Proceed without combining gene duplicates.\n")
        out_E <- flt_E_avg
        out_genes <- flt_genes_avg
      }
    } else {
      cat("Number of transcripts same as number of genes, no need to combine.\n")
      out_E <- flt_E_avg
      out_genes <- flt_genes_avg
    }
  } else {
    out_E <- flt_E_avg
    out_genes <- flt_genes_avg
  }

  ## output
  cat("Constucting rbioarray_flist...")
  filter.results <- list(neg_control_used = neg_control_used,
                         filter_percentile = filter.percentile,
                         filter_threshold_min_sample = filter.threshold.min.sample,
                         filter_summary = flt_summary)
  out <- list(E = out_E,
              genes = out_genes,
              gene_duplicates_combined = combine.gene.duplicate,
              filter_results = filter.results)
  out <- append(out, object[!names(object) %in% c("E", "genes")])
  class(out) <- "rbioarray_flist"
  cat("Done!\n")
  return(out)
}


#' @export
print.rbioarray_flist <- function(x, ...){
  cat("---- rbioarray_flist information ----\n")
  cat(paste0("Number of genes/probes/genomic features upon filtering: ", nrow(x$E), "\n"))
  cat(paste0("Number of genes/probes/genomic features before filtering: ", nrow(x$extra_E_data$original_E), "\n"))
  cat(paste0("Gene duplicates combined: ", ifelse(x$gene_duplicates_combined, TRUE, FALSE), "\n"))
  cat("\n")
  cat(paste0("Number of samples: ", nrow(x$targets), "\n"))
  cat(paste0("Groups: "))
  cat(paste0(levels(x$sample_groups)))
  cat("\n")
}
