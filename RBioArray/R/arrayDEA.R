#' @title rbioarray_rlist
#'
#' @description Function to construct \code{rbioarray_rlist} class object from microarary and annotation data. The \code{rbioarray_rlist} object is the starting point for all microarray data analysis.
#' @param object Input object.
#' @details When \code{object} is missing, the default function is used - make sure to pass all the arguments.
#'
#'          To keep things consistent with the \code{Elist} from the dependent \code{limma} package. The \code{rbioarray_rlist} class contains many common elements from \code{Elist} class.
#'
#' @return A \code{rbioarray_rlist} object. The \code{rbioarray_rlist} class contains the following items:
#'
#'         \code{E}: raw expression matrix (i.e. hybridization signal)
#'
#'         \code{E_background}: background signal matrix if applicable
#'
#'         \code{raw_file.gene_annotation.var_name}
#'
#'         \code{raw_file.gene_id.var_name}
#'
#'         \code{gene_display_name_used}: if the gene display name is extracted from \code{extra.gene.annot.dataframe}.
#'
#'         \code{genes_annotation.gene_id.var_name}
#'
#'         \code{genes_annotation.gene_symbol.var_name}
#'
#'         \code{genes_annotation.control_type}: if \code{gene.annot.control_type.var.name} is set, a list containing all the control type information
#'
#'         \code{genes_annotation.to_remove.var.name}: variable names (string or string vector) to remove from the final output \code{gene} element.
#'                                                    This is important if starting with \code{EListRaw} object, since array positional variables need to be remove for \code{\link{rbioarray_filter_combine}} function.
#'
#'         \code{targets}: the sample annotation data frame.
#'
#'         \code{sample_groups_var_name}
#'
#'         \code{sample_groups}
#'
#' @examples
#' \dontrun{
#' \code{# ElistRaw}
#' rlist <- rbioarray_rlist(object = ElistRaw_list, raw.gene_id.var.name = "ProbeName", extra.gene.annot.dataframe = annot,
#'                          gene.annot.gene_id.var.name = "ProbeName", gene.annot.gene_symbol.var.name = "Gene.Symbol",
#'                          gene.annot.control_type.var.name = "ControlType", gene.annot.control_type.val.pos = 1,
#'                          gene.annot.control_type.val.neg = -1, gene.annot.control_type.val.exp = 0,
#'                          target.annot.file = "Targets.csv", target.annot.file.path = getwd(),
#'                          sample_groups.var.name = "Group")
#'
#' \code{# non-EListRaw}
#' raw_list <- rbioarray_rlist(raw.dataframe = raw, raw.annot.var.name = c("PROBE_ID","ILMN_GENE"), raw.gene_id.var.name = "PROBE_ID",
#'                             extra.gene.annot.dataframe = annot, gene.annot.gene_id.var.name = "PROBE_ID", gene.annot.gene_symbol.var.name = "SYMBOL",
#'                             target.annot.file = "sample_index_end6.csv", sample_groups.var.name = "time_point")
#' }
#' @export
rbioarray_rlist <- function(object, ...){
  ## processing
  if (missing(object)) {
    rlist <- rbioarray_rlist.default(...)
    return(rlist)
  } else {
    UseMethod("rbioarray_rlist", object)
  }
}

#' @title rbioarray_rlist.EListRaw
#'
#' @rdname rbioarray_rlist
#' @method rbioarray_rlist EListRaw
#' @param object A input \code{EListRaw} class object from \code{limma} package
#' @param ... Additional argument for the default function.
#' @export
rbioarray_rlist.EListRaw <- function(object, ...){
  ## set up variables
  raw.dataframe <- data.frame(object$genes, object$E, stringsAsFactors = FALSE, check.names = FALSE)
  raw.background.signal.matrix <- object$Eb
  raw.annot.var.name <- names(object$genes)

  ## processing
  rlist <- rbioarray_rlist.default(raw.dataframe = raw.dataframe, raw.background.signal.matrix = raw.background.signal.matrix,
                                   raw.annot.var.name = raw.annot.var.name, ...)

  return(rlist)
}

#' @title rbioarray_rlist.default
#'
#' @rdname rbioarray_rlist
#' @method rbioarray_rlist default
#' @param raw.dataframe Input data frame containing microarray hybridization signals with rows as probe/gene/genomic features and columns as samples. Note: the data frame should contain at least one annotation column, for the probe/gene/genomic features.
#' @param raw.background.signal.matrix A optional matrix containing background signals. The dimension should be the same as the input expression data without annotation columns.
#' @param raw.annot.var.name A string vector containing variable (i.e. column) name(s) for all the annotation columns in \code{raw.dataframe}.
#' @param raw.gene_id.var.name Variable (i.e. column) name for gene/probe/genomic feature identification from \code{raw.dataframe}.
#' @param extra.gene.annot.dataframe Optional annotation data frame for gene/probe/genomic feature annotation.
#' @param gene.annot.gene_id.var.name Set only when \code{extra.gene.annot.dataframe} is provided, variable name for probe/gene/genomic features identification from \code{extra.gene.annot.dataframe}.
#' @param gene.annot.gene_symbol.var.name Set only when \code{extra.gene.annot.dataframe} is provided, variable name for probe/gene/genomic features display name from \code{extra.gene.annot.dataframe}, e.g. gene symbols.
#' @param gene.annot.control_type.var.name Optional, from either \code{raw.annot} or \code{extra.gene.annot.dataframe} is provided, name for the variable that denotes if the gene is a control and its control type. Default is \code{NULL}.
#' @param gene.annot.control_type.val.pos Set only when \code{gene.annot.control_type.var.name} is provided, value indicating a positive control probe. Default is \code{1}.
#' @param gene.annot.control_type.val.neg Set only when \code{gene.annot.control_type.var.name} is provided, value indicating a negative control probe. Default is \code{-1}.
#' @param gene.annot.control_type.val.exp Set only when \code{gene.annot.control_type.var.name} is provided, value indicating a none control probe. Default is \code{0}.
#' @param gene.annot.rm.var.name Optional variable names for the columns to remove from the gene annotation \code{genes} in the output.
#' @param raw.gene.annot.merge_mode How to merge the raw dataframe and the extra gene annotation dataframe. Options are \code{by_raw}, \code{intersect} and {no_merge}.
#' @param target.annot.file File name for target (i.e. sample/subject) annotation \code{.csv} file.
#' @param target.annot.dataframe Name for a dataframe for target (i.e. sample/subject) annotation.
#' @param target.annot.subject_id.var.name Variable name for subject id from \code{target.annot.file} or \code{target.annot.dataframe}
#' @param target.sample_group.var.name The variable name for sample grouping information from \code{target.annot.file}.
#' @param sample_groups.var.name The same as the \code{target.sample_group.var.name}, legacy argument for compatiblity.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @details The \code{raw.background.signal.matrix} is useful when processing a \code{EListRaw} class object from \code{limma} package,
#'          and using \code{bgc.method = "subtract"} for the \code{\link{rbioarray_transfo_normalize()}} function.
#'
#'          The word "gene" used in argument names and output item names is in its broader meaning of gene/probe/genomic feature.
#'
#'          The function uses \code{data.table} to speed up merging operations.
#'
#'          The optional annotation data frame \code{extra.gene.annot.dataframe} should contain any additional annotation information in addition to the annotation column(s) from \code{raw.dataframe}.
#'          It is noted that \code{extra.gene.annot.dataframe} should contain at least one column for the sample type of gene/probe/genomic feature identification as the identification variable from  \code{raw.dataframe}.
#'          Such column is set via argument \code{raw_file.gene_id.var_name} and \code{genes_annotation.gene_id.var_name}.
#'          The gene display name will only be used if \code{extra.gene.annot.dataframe} is used. Otherwise, it will use \code{raw.gene_id.var.name} from \code{raw.dataframe}.
#'
#'          To merge the raw dataframe and extra gene annotation, the \code{by_raw} mode make sure the raw dataframe is all preserved, which may lead to NAs when the extra annotation dataframe does not include all probes.
#'          For the \code{intersect} mode, the raw dataframe might be subset and only contrain probes with annotation information from the extra annotation dataframe.
#'          With \code{no_merge} mode, the function acts as if no extra gene annotation dataframe is provided.
#'
#'          The function automatically compare, process and re-order the \code{E} matrix according to the values of \code{target.annot.file}, \code{target.annot.dataframe}, and \code{target.annot.subject_id.var.name}.
#'          Behavior:
#'          1. \code{E} matrix is subset by \code{targets$target.annot.subject_id.var.name} if \code{colnames(E)} is a subset of \code{targets$target.annot.subject_id.var.name}, vice versa.
#'          2. Both \code{E} and \code{targets} dataframe are subset if there is intersect between \code{colnames(E)} and \code{targets$target.annot.subject_id.var.name}, and the intersect order is set by \code{targets$target.annot.subject_id.var.name}.
#'          3. \code{E} matrix column will be re-ordered by \code{targets$target.annot.subject_id.var.name} if the orders of both are different.
#'
#' @importFrom data.table setDT
#' @export
rbioarray_rlist.default <- function(raw.dataframe,
                                    raw.background.signal.matrix = NULL,
                                    raw.annot.var.name = NULL, raw.gene_id.var.name = NULL,

                                    extra.gene.annot.dataframe = NULL,
                                    gene.annot.gene_id.var.name = NULL, gene.annot.gene_symbol.var.name = NULL,
                                    gene.annot.control_type.var.name = NULL,
                                    gene.annot.control_type.val.pos = 1, gene.annot.control_type.val.neg = -1, gene.annot.control_type.val.exp = 0,
                                    gene.annot.rm.var.name = NULL,
                                    raw.gene.annot.merge_mode = c("by_raw", "intersect", "no_merge"),

                                    target.annot.file = NULL, target.annot.file.path = getwd(),
                                    target.annot.dataframe = NULL,
                                    target.annot.subject_id.var.name = NULL,
                                    target.sample_group.var.name = NULL,
                                    sample_groups.var.name = NULL,
                                    verbose = TRUE){
  # ------ check arguments and set up variables ------
  # ---- raw data ----
  if (verbose) cat("Loading raw measurement data...")
  if (is.null(dim(raw.dataframe)))
    stop("raw.matrix has to be a dataframe, with rows for genes/probes/genomic features, columns for samples.")
  if (is.null(raw.annot.var.name) || is.null(raw.gene_id.var.name))
    stop("Please set raw.annot.var.name AND raw.gene_id.var.name arguments.")
  if (!all(raw.annot.var.name %in% names(raw.dataframe)) || !raw.gene_id.var.name %in% names(raw.dataframe))
    stop("Annoation variables (i.e. columns) and/or the gene_id variable (i.e. column) not found in raw.dataframe")
  raw_dfm <- setDT(raw.dataframe)
  if (verbose) cat("Done!\n")

  # ---- gene annotation ----
  if (verbose) cat("Loading gene/probe/feature annotation data and contructing \"E\" matrix and \"genes\" dataframe...\n")
  raw.gene.annot.merge_mode <- match.arg(raw.gene.annot.merge_mode)

  if (is.null(extra.gene.annot.dataframe) || raw.gene.annot.merge_mode == "no_merge"){  # check and load gene annotation
    if (verbose) cat("Either no extra gene annotation provided or  \"raw.gene.annot.merge_mode\" set to \"no_merge\"...\n")
    # E <- as.matrix(raw_dfm[, !names(raw_dfm) %in% raw.annot.var.name, drop = FALSE]) # remove annotation columns # vanila R
    E <- as.matrix(raw_dfm[, !names(raw_dfm) %in% raw.annot.var.name, with = FALSE, drop = FALSE]) # remove annotation columns, data.table syntax
    rownames(E) <- raw_dfm$probe_id

    cat("No extra gene annotation used. Proceed with raw.dataframe annoation information.\n")
    gene.annot.gene_id.var.name <- raw.gene_id.var.name
    gene.annot.gene_symbol.var.name <- raw.gene_id.var.name
    genes <- raw.dataframe[, ..raw.annot.var.name, drop = FALSE]
    gene.symbol <- FALSE
  } else {
    if (verbose) cat(paste0("\"raw.gene.annot.merge_mode\" set to ", raw.gene.annot.merge_mode, "..."))
    if (!is.data.frame(extra.gene.annot.dataframe)) stop("gene.annot needs to be a dataframe")
    # if (nrow(extra.gene.annot.dataframe) < nrow(raw.dataframe)) stop("extra.gene.annot.dataframe has less record than the raw.dataframe") # check size

    if (is.null(gene.annot.gene_id.var.name) || is.null(gene.annot.gene_symbol.var.name)) { # gene id and symbol variables
      stop("Please set gene.annot.gene_id.var.name AND gene.annot.gene_symbol.var.name arguments according to gene.annot")
    } else if (!gene.annot.gene_id.var.name %in% names(extra.gene.annot.dataframe) || !gene.annot.gene_symbol.var.name %in% names(extra.gene.annot.dataframe)) {
      stop("Gene ID variable and/or gene symbol variable not found in extra.gene.annot.dataframe")
    } else {
      gene.symbol <- TRUE
    }
    # additional variables
    genes_annot_dfm <- setDT(extra.gene.annot.dataframe)

    ## Set up the information
    # merge gene annotation with raw dataframe
    raw_dfm$row_id <- as.integer(rownames(raw_dfm))  # used to make sure the same order as the background matrix for EListRaw
    # raw_dfm$merge_id <- raw_dfm[, raw.gene_id.var.name]  # vanila R syntax
    # genes_annot_dfm$merge_id <- genes_annot_dfm[, gene.annot.gene_id.var.name]  # vanila R syntax
    raw_dfm$merge_id <- raw_dfm[, ..raw.gene_id.var.name]  # DT[, ..var_name] is a data.table syntax
    genes_annot_dfm$merge_id <- genes_annot_dfm[, ..gene.annot.gene_id.var.name]  # DT[, ..var_name] is a data.table syntax

    if (raw.gene.annot.merge_mode == "by_raw") {
      merged_raw_gene_annot_dfm <- merge(raw_dfm, genes_annot_dfm, all.x = TRUE)  # this merge will extract annotation info from gene_annot_dfm and merge to the smaller data dataframe.
    } else if (raw.gene.annot.merge_mode == "intersect") {
      merged_raw_gene_annot_dfm <- merge(raw_dfm, genes_annot_dfm, by = "merge_id")  # this merge will extract annotation info from gene_annot_dfm and merge to the smaller data dataframe.
    }

    if (any((is.na(merged_raw_gene_annot_dfm[, ..gene.annot.gene_id.var.name])))){
      warning(paste0("The annotation \"", gene.annot.gene_id.var.name, "\" column of the output rbioarray_rlist element \"gene\" dataframe contains NA. Recommend trying another merge mode for the raw dataframe and extra gene annotation dataframe."))
    }

    merged_raw_gene_annot_dfm <- merged_raw_gene_annot_dfm[order(merged_raw_gene_annot_dfm$row_id), ]  # restore the original order
    all_annot_var_names <- unique(c(raw.annot.var.name, names(genes_annot_dfm), "row_id"))  # all annotation variable names

    # Set up output E matrix
    # E <- as.matrix(merged_raw_gene_annot_dfm[, !names(merged_raw_gene_annot_dfm) %in% all_annot_var_names]) # remove annotation columns: vanila R syntax
    E <- as.matrix(merged_raw_gene_annot_dfm[, !names(merged_raw_gene_annot_dfm) %in% all_annot_var_names, with = FALSE]) # data.table syntax
    rownames(E) <- merged_raw_gene_annot_dfm$merge_id
    if (!is.null(raw.background.signal.matrix) && dim(raw.background.signal.matrix) != dim(E)) {
      cat("The dimension of raw.background.signal.matrix not the same as the expressoin matrix. Proceed without using it. ")
      raw.background.signal.matrix <- NULL
    }

    # set up output annotation information
    # genes <- merged_raw_gene_annot_dfm[, names(merged_raw_gene_annot_dfm) %in% all_annot_var_names, drop = FALSE]  # vanila R syntax
    genes <- merged_raw_gene_annot_dfm[, names(merged_raw_gene_annot_dfm) %in% all_annot_var_names, drop = FALSE, with = FALSE]
    # genes <- genes[, !names(genes) %in% c("merge_id", "row_id"), drop = FALSE]  # row_id has to be removed for the filtering step (aka to avoid same gene symbol but different row id): vanila R syntax
    genes <- genes[, !names(genes) %in% c("merge_id", "row_id"), drop = FALSE, with = FALSE]  # row_id has to be removed for the filtering step (aka to avoid same gene symbol but different row id): data.table syntax
    if (verbose) cat("Done!\n")
  }

  # ---- target annotation ----
  if (verbose) cat("Loading target annotation data and contructing \"targets\" dataframe...")
  if (is.null(target.annot.file) && is.null(target.annot.dataframe)){  # check and load target (sample) annotation
    stop("Please provide a target annotation file or dataframe for target.annot.file arugment.")
  } else if (!is.null(target.annot.file) && !is.null(target.annot.dataframe)){
    stop("Please only set one or the same value: \"target.annot.file\" or \"target.annot.dataframe\".")
  } else if (!is.null(target.annot.file)){
    target.annot_name_length <- length(unlist(strsplit(target.annot.file, "\\.")))
    target.annot_ext <- unlist(strsplit(target.annot.file, "\\."))[target.annot_name_length]
    if (target.annot_ext != "csv") {
      stop("target.annot.file is not in csv format.")
    } else {
      if(verbose) cat("Loading target annotation file...")
      tgt <- read.csv(file = paste0(target.annot.file.path, "/",target.annot.file), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
      if(verbose) cat("Done!\n")
    }
  } else {
    tgt <- target.annot.dataframe
  }

  if (is.null(target.sample_group.var.name) && is.null(sample_groups.var.name)){  # check and load target (sample) annotation
    stop("Please provide \"sample_groups.var.name\" or \"target.sample_group.var.name\" from the target annotation dataframe or file.")
  } else if (!is.null(target.sample_group.var.name) && !is.null(sample_groups.var.name)){
    if (!(all(target.sample_group.var.name %in% sample_groups.var.name) && length(target.sample_group.var.name) == length(sample_groups.var.name))) stop("Please only set one or with the same value: \"target.sample_group.var.name\" or \"sample_groups.var.name\".")
  } else if (!is.null(target.sample_group.var.name)){
    sample_groups.var.name <- target.sample_group.var.name
  }

  if (!all(sample_groups.var.name %in% names(tgt))){
    stop("Sample group annotation variable not found in the target annotation file.")
  } else if (length(sample_groups.var.name) == 1) {
    sample.groups <- factor(tgt[, sample_groups.var.name], levels = unique(tgt[, sample_groups.var.name]))
  } else {
    sample.groups <- vector(mode = "list", length = length(sample_groups.var.name))
    names(sample.groups) <- sample_groups.var.name
    for (var in sample_groups.var.name) {
      sample.groups[[var]] <- factor(tgt[, var], levels = unique(tgt[, var]))
    }
  }
  if (verbose) cat("Done!\n")

  # ---- optional control_type variable ----
  if (!is.null(gene.annot.control_type.var.name)) {
    if (verbose) cat("Loading control type information and contructing \"genes_annotation.control_type\" list...")
    if (!gene.annot.control_type.var.name %in% names(genes)) {
      cat("The set control type variable not found in extra.gene.annot.dataframe. Proceed without using it.\n")
      gene.annot.control_type = NULL
    } else {
      control_type.values <- unique(genes[, gene.annot.control_type.var.name])
      input.control_type.values <- c(gene.annot.control_type.val.pos, gene.annot.control_type.val.neg, gene.annot.control_type.val.exp)
      if (all(input.control_type.values %in% control_type.values)){
        gene.annot.control_type = list(control_type.var_name = gene.annot.control_type.var.name,
                                       pos_type.value = gene.annot.control_type.val.pos,
                                       neg_type.value = gene.annot.control_type.val.neg,
                                       exp_type.value = gene.annot.control_type.val.exp)
      } else {
        if (verbose) cat("Input probe control type values not found in the gene control type variable from gene annotation. Proceed without using them.\n")
        gene.annot.control_type = NULL
      }
    }
    if (verbose) cat("Done!\n")
  } else {
    gene.annot.control_type = NULL
  }

  # ---- variable names to remove ----
  if (!is.null(gene.annot.rm.var.name) && !gene.annot.rm.var.name %in% names(genes)) {
    cat("gene.annot.rm.var.name not found in gene annation. Proceed without using it")
    gene.annot.rm.var.name <- NULL
  }

  # ------ output ------
  # -- subset E or tgt, and rearrange E columns (if necessary) according to the tgt subject var --
  if (verbose) cat("Processing \"E\" matrix and \"targets\" dataframe...\n")
  i = 1
  while (i < length(names(tgt)[!names(tgt) %in% sample_groups.var.name])){
    if (!is.null(target.annot.subject_id.var.name) && target.annot.subject_id.var.name %in% names(tgt) ){
      tgt_var <- target.annot.subject_id.var.name
      i <- length(names(tgt)[!names(tgt) %in% sample_groups.var.name]) + 1
      if (verbose) cat("no loop\n")
    } else {
      tgt_var <- names(tgt)[!names(tgt) %in% sample_groups.var.name][i]
    }

    if (length(tgt[[tgt_var]]) == length(colnames(E)) && all(tgt[[tgt_var]] %in% colnames(E))) {
      if (verbose) cat(paste0("No subsetting of E or targets dataframe required", tgt_var, "\n"))
      break
    } else if (all(tgt[[tgt_var]] %in% colnames(E))) {
      if (verbose) cat(paste0("E subset by: ", tgt_var, " variable from the targets dataframe\n"))
      break
    } else if (all(colnames(E) %in% tgt[[tgt_var]])) {
      tgt <- tgt[tgt[tgt_var] %in% colnames(E), , drop = FALSE]
      if (verbose) cat(paste0("tgt subset by colnames(E) according with the ", tgt_var, " variable from the targets dataframe\n"))
      break
    } else if (length(intersect(tgt[tgt_var], colnames(E))) >= 1) {
      common_subject_id <- intersect(tgt[tgt_var], colnames(E))
      E <- E[, common_subject_id, drop = FALSE]
      tgt <-  tgt[tgt[tgt_var] %in% common_subject_id, , drop = FALSE]
      if (verbose) cat(paste0("tgt and E subset by intersect between colnames(E) and the values of ",  tgt_var, " variable from the targets dataframe\n"))
      break
    } else {
      i <- i + 1
      next
    }
    stop(paste0("\"raw.dataframe\" column name as subject id not found in the ", tgt_var, " variable from the \"target.annot.file\" or \" target.annot.dataframe\""))
  }
  E <- E[, tgt[[tgt_var]], drop = FALSE]
  if (verbose) cat(paste0("E columns might be re-arranged according to the target dataframe variabl: ", tgt_var, "\n"))

  # -- output --
  if (verbose) cat("Constructing rlist...")
  out <- list(E = E,
              E_background = raw.background.signal.matrix,
              raw_file.gene_annotation.var_name = raw.annot.var.name,
              raw_file.gene_id.var_name = raw.gene_id.var.name,
              genes = genes,
              gene_display_name_used = gene.symbol,
              genes_annotation.gene_id.var_name = gene.annot.gene_id.var.name,
              genes_annotation.gene_symbol.var_name = gene.annot.gene_symbol.var.name,
              genes_annotation.control_type = gene.annot.control_type,
              genes_annotation.to_remove.var.name = gene.annot.rm.var.name,
              targets = tgt,
              sample_groups_var_name = sample_groups.var.name,
              sample_groups = sample.groups)
  class(out) <- "rbioarray_rlist"
  if (verbose) {
    cat("Done!\n")
    cat("\n")
    cat(paste0("The resulted rbioarray_rlist object contains ", nrow(E), " genes/probes/genomic features, ", nrow(tgt), "\n"))
    if ("list" %in% class(sample.groups)) {
      for (l in names(sample.groups)) {
        cat(paste0("Group variable: ", l, "\n   Levels: ", paste0(as.character(unique(sample.groups[[l]])), collapse = ", "), "\n"))
      }
    } else {
      cat(paste0("Groups: ", paste0(as.character(unique(sample.groups)), collapse = ", "), "\n"))
    }
  }

  return(out)
}


#' @export
print.rbioarray_rlist <- function(x, ...){
  cat("\n")
  cat("---- rbioarray_rlist information ----\n")
  cat(paste0("Number of genes/probes/genomic features: ", nrow(x$E), "\n"))
  cat(paste0("Number of samples: ", nrow(x$targets), "\n"))
  cat(paste0("Groups: "))
  if ("list" %in% class(x$sample_groups)){
    for (l in names(x$sample_groups)) {
      cat(paste0("Group variable: ", l, "\n   Levels: ", paste0(as.character(unique(x$sample_groups[[l]])), collapse = ", "), "\n"))
    }
  } else {
    cat(paste0(levels(x$sample_groups)))
  }
  cat("\n\n")
  cat(paste0("Gene display name from annotation information: ", ifelse(x$gene_display_name_used, "Available\n", "Unavailable\n")))
  cat("\n")
}


#' @title rbioarray_transfo_normalize
#'
#' @description Data log transformation and normalization function for microarray.
#' @param object Input object with raw data and annotation information. Should be a \code{rbioarray_rlist} class.
#' @param ... Additional arguments for corresponding S3 class methods.
#' @details The expression matrix will be normalized and then log2 transformed for output.
#'
#'          A note on the \code{bgc.method}: "subtract" or "normexp"?
#'          In most cases, we use "normexp" (or "auto" without providing a background matrix).
#'          However the \code{limma} manual provided a table of background matrix information for the popular microarray platforms.
#'          If available, we can use "subtract" (or "auto" with the correct background matrix).
#'
#'          The function currently only supports \code{"quantile"} method for between array normalization.
#'          A note to the \code{normalizeBetweenArrays} function from \code{limma} package:
#'          The function normalizes data BEFORE log2 transformation when the input is \code{EListRaw} object.
#'             However, when input is \code{matrix}, it assumes the data has already been log2 transformed, meaning this function will not
#'          log2 transform the data after normalization.
#'             After comparison, these two methods will NOT produce the same results for the same expression data:
#'          In other words, applying \code{normalizeBetweenArrays} directly to log2 transformed E matrix is NOT the same as apply the function
#'          to the \code{EListRaw} that contains E.
#'
#'             In fact, when using \code{"quantile"} method, applying \code{normalizeBetweenArrays} to log2 transformed E matrix is the
#'          same as applying \code{normalizeBetweenArrays} to the \code{EListRaw} that has object$E using \code{"cyclicloess"} method.
#'          Indeed, the source code of \code{normalizeBetweenArrays} suggests that's the case since log2 transformation is applied BEFORE
#'          \code{"cyclicloess"} normalization. However, transformation happens AFTER normalization for \code{"quantile"} and other methods.
#'
#'          In \code{limma}, the \code{EListRaw} is data before log2 transformation, whereas \code{EList} contains log2 transformed data.
#'
#'          To avoid such confusion and provide a unified experience, the current \code{rbioarray_transfo_normalize.rbioarray_rlist} method only
#'          treats input E data as a matrix, which is ensured by the \code{rbioarray_rlist} function. The \code{rbioarray_rlist} accepts
#'          unlogged data or converts \code{ElistRaw} objects into to a \code{matrix}.
#'
#'          The \code{limma} author suggests doing quantile on logged or unlogged data remains to be debatable, but "slowly leaning towards"
#'          quantile on raw and then log transform. As such, the \code{rbioarray_transfo_normalize}will first conduct normalization then log2
#'          transform the data.
#'
#'
#' @return A \code{rbioarray_plist} class object with the following items:
#'
#'         \code{E}: Normalized and then log2 transformed expression matrix
#'
#'         \code{design}: Microarray experiment sample design matrix
#'
#'         \code{ArrayWeight}
#'
#'         \code{raw_data}: A list including original E matrix and background matrix from the input \code{rbioarray_rlist}.
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
  if (any(class(object) != "rbioarray_rlist")) stop("The input object needs to be \"rbioarray_rlist\" class.\"")

  ## use method
  UseMethod("rbioarray_transfo_normalize", object)
}


#' Title rbioarray_transfo_normalize.rbioarray_rlist
#'
#' @rdname rbioarray_transfo_normalize
#' @method rbioarray_transfo_normalize rbioarray_rlist
#' @param object Input object with raw data and annotation information. Could be \code{rbioarray_rlist}, \code{Elist} or \code{MAList} classes.
#' @param design Microarray experiment sample design matrix. Make sure the design colnames are the same as the levels of \code{object$sample_groups}.
#' @param ... Additional arguments the default method \code{\link{rbioarray_transfo_normalize.default}}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @details The \code{rbioarray_rlist} object can be obtained from \code{\link{rbioarray_rlist}} function.
#' @return A \code{rbioarray_plist} class object.
#'
#' @export
rbioarray_transfo_normalize.rbioarray_rlist <- function(object, design, ..., verbose = TRUE){
  ## processing
  default_out <- rbioarray_transfo_normalize.default(E = object$E, E.background = object$E_background,
                                                     between.sample.weight.design = design, ..., verbose = verbose)
  ## output
  if (verbose) cat("\n")
  if (verbose) cat("Constructing rbioarray_plist...")
  out <- append(default_out, object[!names(object) %in% c("E", "E_background")])
  class(out) <- "rbioarray_plist"
  if (verbose) cat("Done!\n\n")
  return(out)
}


#' @title rbioarray_transfo_normalize.default
#'
#' @rdname rbioarray_transfo_normalize
#' @method rbioarray_transfo_normalize default
#' @param E Input raw expression value matrix with columns for samples, rows for genes/probes/genome features.
#' @param E.background A optional matrix containing background signals. The dimension should be the same as the input expression data without annotation columns.
#' @param bgc.method Background correction method. Default is \code{"auto"}. See \code{backgroundCorrect()} function from \code{limma} package for details.
#' @param between.sample.norm.method Normalization method. Default is \code{"quantile"}. See \code{normalizeBetweenArrays()} function from \code{limma} package for details.
#' @param between.sample.weight.design Microarray experiment sample design matrix for between sample weight calculation.
#' @param ... Additional arguments the default method \code{backgroundCorrect.matrix} function from \code{limma} package.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @details The word "gene" used in argument names and output item names is in its broader meaning of gene/probe/genome feature.
#'          ADD the fact that \code{normalizeBetweenArrays} function is where \code{limma} log transforms data.
#' @return A list with the core items for a \code{rbioarray_plist} class.
#' @importFrom limma backgroundCorrect normalizeBetweenArrays backgroundCorrect.matrix arrayWeights
#' @export
rbioarray_transfo_normalize.default <- function(E, E.background = NULL,
                                                bgc.method = c("auto", "none", "subtract", "half", "minimum", "movingmin", "edwards", "normexp"),
                                                between.sample.norm.method = c("quantile", "Aquantile"),
                                                between.sample.weight.design = NULL, ...,
                                                verbose = TRUE){
  ## check arguments
  bgc.methd <- match.arg(tolower(bgc.method), c("auto", "auto", "none", "subtract", "half", "minimum", "movingmin", "edwards", "normexp"))
  between.sample.norm.method <- match.arg(between.sample.norm.method, c("quantile", "Aquantile"))
  if (between.sample.norm.method != "quantile") {
    stop("Currently the function only supports quantile method")
  }
  if (is.null(between.sample.weight.design)) stop("Please provide the microarray design matrix for between.sample.weight.design argument.")

  ## pre-processing
  if (verbose) cat("Background correction: \n")
  # E_gbc is a "matrix" "array" object
  E_bgc <- backgroundCorrect.matrix(E = E, Eb = E.background, method = bgc.method, ...) #background correction
  if (verbose) cat("Done!\n")

  ## normalization
  if (verbose) cat("\n")
  if (verbose) cat(paste0("Data normalization using ", between.sample.norm.method, " method..."))
  Norm <- log2(normalizeBetweenArrays(E_bgc, between.sample.norm.method)) # quantile normalization then log2 transformation
  Wgt <- arrayWeights(Norm, design = between.sample.weight.design) # array weight
  if (verbose) cat("Done!\n")

  ## output
  raw.data <- list(original_E = E, original_background = E.background)
  default_out <- list(E = Norm,
                      design = between.sample.weight.design,
                      ArrayWeight = Wgt,
                      background_correction_method = bgc.method,
                      between_sample_normalization_method = between.sample.norm.method,
                      raw_data = raw.data)
  return(default_out)
}


#' @export
print.rbioarray_plist <- function(x, ...){
  cat("\n")
  cat("---- rbioarray_plist information ----\n")
  cat(paste0("Background correction method: ", x$background_correction_method, "\n"))
  cat(paste0("Between-sample normalization method: ", x$between_sample_normalization_method))
  cat("\n\n")
  cat(paste0("Number of genes/probes/genomic features: ", nrow(x$E), "\n"))
  cat(paste0("Number of samples: ", nrow(x$targets), "\n"))
  cat(paste0("Groups: "))
  cat(paste0(levels(x$sample_groups)))
  cat("\n")
}


#' @title rbioarray_filter_combine
#'
#' @description Function to filter, averaging and (if set) combine genes/probes/genome features from the \code{rbioarray_plist} objects.
#' @param object Input \code{rbioarray_plist} object from function \code{\link{rbioarray_transfo_normalize}}.
#' @param filter.bg If to filter according to the threshold set by \code{filter.percentile} and \code{filter.threshold.min.sample}.
#' @param filter.percentile The percentile threshold for filtering. Default is \code{0.05}. See details for more.
#' @param filter.threshold.min.sample Minimum number of samples meeting the filtering threshold. Default is \code{NULL}. See details for more.
#' @param combine.gene.duplicate If to combine different transcripts from the same gene/genome feature. Default is \code{FALSE}. See details for more.
#' @param parallelComputing Whether to use parallel computing or not. Default is \code{TRUE}.
#' @param cluterType clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @details For \code{filter.percentile}, when the input data has negative probes, the value is set to \code{0.95} so that the 95 percentile of the negative values is the considered the threhold
#'          When there is no negative control probes, the value is set to \code{0.05}, so that the 5 percentile of the entire data is considered the lower threshold
#'
#'          For \code{filter.threshold.min.sample}, usually make sure to ensure the target gene has at least three samples, so that stats can be done.
#'
#'          When \code{combine.gene.duplicate = TRUE},the function combines the expression values from different probes from the same gene.
#'          This depends on if the input object has a valid gene symbol variable. The function retains probes of a gene with highest variation
#'          across groups.
#'
#' @return A \code{rbioarray_flist} class, including the following core items:
#'
#'         \code{E}: filtered (and if set, combined) and normalized expression matrix
#'
#'         \code{genes}: gene annotation with same row number as E
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
                                     filter.bg = TRUE,
                                     filter.percentile = 0.05,
                                     filter.threshold.min.sample = NULL,
                                     combine.gene.duplicate = FALSE,
                                     parallelComputing = FALSE, clusterType = "PSOCK",
                                     verbose = TRUE){
  ## check arguments
  if (!"rbioarray_plist" %in% class(object)) stop("The input object needs to be, but not exclusive to, rbioarray_plist class.")

  ## filter
  if (filter.bg) {
    # set negative/low expression threshold
    if (!is.null(object$genes_annotation.control_type)){ # 0.95 percentile of all negative control values
      cat("Gene control type variable detected in the input object. The filter.percentile argument value reset to 0.95.\n")
      filter.percentile = 0.95
      control_type.var.name = object$genes_annotation.control_type$control_type.var_name
      control_type.neg.value = object$genes_annotation.control_type$neg_type.value

      if (any(class(object$E[object$genes[, control_type.var.name] == control_type.neg.value, ]) == "numeric")){ # if there is only one entry in the neg values
        neg <- object$E[object$genes[, control_type.var.name] == control_type.neg.value, ] # no 95% percentile required as only one neg entry
      } else {
        neg <- apply(object$E[object$genes[, control_type.var.name] == control_type.neg.value, ], 2, function(x)quantile(x, p = filter.percentile)) # neg95
      }
      neg_control_used = TRUE
    } else {
      neg <- apply(object$E, 2, function(x)quantile(x, p = filter.percentile)) # 5% percentile of all the data
      neg_control_used = FALSE
    }

    if (verbose) cat("Filtering low expresson genes/probes/genomic features...")
    # low expression cuttoff set at at least 10% hihger than the neg
    LE_cutoff <- matrix(1.1 * neg, nrow(object$E), ncol(object$E), byrow = TRUE)
    # set filtering matrix
    if (is.null(filter.threshold.min.sample)) {
      filter.threshold.min.sample <- min(table(object$sample_groups))
    }
    isexpr <- rowSums(object$E > LE_cutoff) >= filter.threshold.min.sample

    flt_summary <- as.numeric(table(isexpr))
    isexpr_fact <- factor(isexpr, levels = unique(isexpr))
    if (length(levels(isexpr_fact)) == 1){
      if (levels(isexpr_fact) == "TRUE"){
        flt_summary <- c(0, flt_summary)
      } else {
        flt_summary <- c(flt_summary, 0)
      }
    }
    names(flt_summary) <- c("filtered", "remaning")

    # filter
    flt_E <- object$E[isexpr, ] # this is a way of extracting samples logically considered TRUE by certain tests
    flt_genes <- object$genes[isexpr, !names(object$genes) %in% object$genes_annotation.to_remove.var.name, drop = FALSE]
    if (verbose) cat("Done!\n")
  } else {
    flt_summary <- NULL
    flt_E <- object$E  # no filter
    flt_genes <- object$genes  # no filter
  }


  ## averaging technical replicates
  if (verbose) cat("Averaging technical replicates...")
  if ("data.table" %in% class(flt_genes)) {
    average_id <- as.matrix(flt_genes[, object$genes_annotation.gene_id.var_name, with = FALSE])[, 1]
    flt_E_avg <- avereps(flt_E, ID = average_id)
    flt_genes_avg <- unique(flt_genes[average_id %in% rownames(flt_E_avg), , drop = FALSE])
  } else {
    flt_E_avg <- avereps(flt_E, ID = flt_genes[, object$genes_annotation.gene_id.var_name])
    flt_genes_avg <- unique(flt_genes[flt_genes[, object$genes_annotation.gene_id.var_name] %in% rownames(flt_E_avg), , drop = FALSE])
  }
  if (verbose) cat("Done!\n")

  ## combine duplicate genes if set
  # this part needs to be updated to make compatible with data.table
  if (combine.gene.duplicate) {
    if (length(flt_genes_avg[, object$genes_annotation.gene_id.var_name]) != length(unique(flt_genes_avg[, object$genes_annotation.gene_symbol.var_name]))) {
      if (object$gene_display_name_used){
        if (verbose) cat("Combining gene duplicates (i.e. different transcripts belonging to the same gene/genomic feature)...")
        if (!parallelComputing) {
          combGeneProbe <- foreach(i = unique(flt_genes_avg[, object$genes_annotation.gene_symbol.var_name]), .combine = "c") %do% {
            tmp <- flt_E_avg[which(flt_genes_avg[, object$genes_annotation.gene_symbol.var_name] %in% i), ]
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
        if (verbose) cat("Done! Records without a gene symbol have been automatically removed. \n")
      } else {
        if (verbose) cat("\n")
        if (verbose) cat("Gene symbol not used in the input object. Proceed without combining gene duplicates.\n")
        out_E <- flt_E_avg
        out_genes <- flt_genes_avg
      }
    } else {
      if (verbose) cat("Number of transcripts same as number of genes, no need to combine.\n")
      out_E <- flt_E_avg
      out_genes <- flt_genes_avg
    }
  } else {
    out_E <- flt_E_avg
    out_genes <- flt_genes_avg
  }
  ## output
  if (verbose) cat("Constucting rbioarray_flist...")
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
  if (verbose) cat("Done!\n")
  return(out)
}

#' @export
print.rbioarray_flist <- function(x, ...){
  cat("\n")
  cat("---- rbioarray_flist information ----\n")
  cat(paste0("Number of genes/probes/genomic features upon filtering/averaging/combing: ", nrow(x$E), "\n"))
  cat(paste0("Original number of genes/probes/genomic features: ", nrow(x$raw_data$original_E), "\n"))
  cat(paste0("Gene duplicates combined: ", ifelse(x$gene_duplicates_combined, TRUE, FALSE), "\n"))
  cat("\n")
  cat(paste0("Number of samples: ", nrow(x$targets), "\n"))
  cat(paste0("Groups: "))
  cat(paste0(levels(x$sample_groups)))
  cat("\n")
}

#' @title microarray_de
#'
#' @description Function that performs statistical analysis for microarray data from \code{rbioarray_flist} class object.
#' @param object The input \code{rbioarray_flist} object.
#' @param contra contra Contrast matrix.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return A \code{rbioarray_de} class object with DE results, containing the following core items:
#'
#'         \code{F_stats}
#'
#'         \code{DE_results}: a list containing the DE results
#'
#'         \code{comparisons}: a list with comparisons and comparison levels
#'
#'         \code{fit}: the limma fitted DE object as a reference
#'
#'         \code{input_data}: input E matrix and genes data frame from rbioarray_flist object
#'                            note the E matrix is \code{E} from the \code{rbioarray_flist} object, meaning the filtered and normalized data
#'
#'         Additionally, the \code{rbioarray_de} includes additional items from the input \code{rbioarray_flist} object.
#'
#' @importFrom limma lmFit eBayes topTable contrasts.fit
#' @examples
#' \dontrun{
#'
#' de_list <- microarray_de(object = filtered_list, contra = contra)
#'
#' }
#' @export
microarray_de <- function(object, contra, verbose = TRUE){
  if (verbose) cat("Constructing rbioarray_de object...")
  ## argument check
  if (!"rbioarray_flist" %in% class(object)) stop("The input object needs to be a \"rbioarray_flist\" class.")

  ## variable initation
  cf <- colnames(contra)
  contra_levels <- vector(mode = "list", length = length(cf))
  contra_levels[] <- foreach(i = seq(length(cf))) %do% {
    rownames(contra)[which(contra[, i] != 0)]
  }
  names(contra_levels) <- cf

  ## fitting
  fit <- lmFit(object$E, design = object$design, weights = object$ArrayWeight)
  fit <- contrasts.fit(fit, contrasts = contra)
  fit <- eBayes(fit)
  fit$genes <- object$genes

  ## output
  f_stats <- topTable(fit, number = Inf, sort.by = "none")
  f_stats[, object$genes_annotation.gene_id.var_name] <- rownames(f_stats)
  f_stats <- merge(f_stats, object$genes)

  de_list <- vector(mode = "list", length(cf))
  de_list[] <- foreach(i = seq(length(cf))) %do% {
    # below: set sort.by = "none" to preserve the order, for supervsied clustering analysis.
    # specifically, it is for constructing threshoding vector used for subsetting E matrix,
    de_dfm <- topTable(fit = fit, coef = cf[i], number = Inf, sort.by = "none")
    de_dfm[, object$genes_annotation.gene_id.var_name] <- rownames(de_dfm)
    de_dfm <- merge(de_dfm, object$genes)
    de_dfm
  }
  names(de_list) <- cf

  flist_data <- list(E = object$E, genes = object$genes)
  comparisons <- list(comparisons = cf, comparison_levels = contra_levels)

  out <- list(F_stats = f_stats,
              DE_results = de_list,
              comparisons = comparisons,
              fit = fit,
              input_data = flist_data)
  out <- append(out, object[!names(object) %in% c("E", "genes", "raw_data", "raw_file.genes_annotation.var_name", "raw_file.gene_id.var_name")])
  class(out) <- "rbioarray_de"
  if (verbose) cat("Done!\n")

  return(out)
}


#' @export
print.rbioarray_de <- function(x, ...){
  cat("\n")
  cat("--- Microarray gene differential expression analysis ---\n")
  cat("\n")
  cat("Comparisons assessed: \n")
  cat(paste0("\t", x$comparisons$comparisons, "\n"))
  cat("\n")
}
