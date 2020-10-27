#' @title rbioseq_clr_ilr_transfo
#'
#' @description Log ratio transformation function for read count data. Row: sample, column: features.
#' @param x Input read count data matrix.
#' @param offset Read count offset value added to avoid zero. Default is \code{1}.
#' @param mode Log ratio transformation method. Options are "clr" (centered log transformation) and "ilr" (isometric log transformation). Default is \code{"clr"}.
#' @param ilr.method.fast Useful only when \code{mode = "ilr"}. Default is \code{TRUE}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return A data matrix with log ratio transformed values.
#' @details This function is needed as part of the data pre-processing procedure to run multivariate and machine learning analysis featured in \code{RBioFS} package.
#'
#'          As per Quinn et al. (2018), NGS data can be considered as compositional data. As such, data must be transformed for usual statistical analysis and visualization.
#'          Log ratio transformation serves such purpose. Note that the number of features will be one less when using "ilr" method.
#'
#'          It is not to be combined with the other normalization methods featured in the \code{\link{rbioseq_DE}}.
#'
#'          Ref: Quinn TP, et al. 2018. Understanding sequencing data as compositions: an outlook and review. Bioinformatics. 2018: 1 - 9.
#' @examples
#' \dontrun{
#' tstX <- rbioseq_clr_ilr_transfo(tstdata, offset = 1, mode = "clr")
#' }
#' @export
rbioseq_clr_ilr_transfo <- function(x, offset = 1, mode = c("clr", "ilr"), ilr.method.fast = TRUE,
                                    verbose = TRUE){
  # data and arguments check
  mode <- match.arg(tolower(mode), c("clr", "ilr"))
  if (!is.matrix(x))stop("x needs to b e a matrix")
  if (any(x == 0) & offset == 0)stop("zero detected in x. set offset to avoid it for ratio transformation")
  # if (!tolower(mode) %in% c("clr", "ilr"))stop("choose the proper transformation mode: \"clr\" or \"ilr\"")

  # log ratio transformation
  if (tolower(mode) == "clr"){  # clr calculation
    if (dim(x)[2] == 1){
      out <- list(x.clr = x, gm = rep(1, dim(x)[1]))
    } else {
      gm <- apply(x, 1, function(x) exp(mean(log(x + offset))))  # geometric mean = exp(mean(log(X)))
      clrX <- log((x + offset) / (gm))  # clr formula
      out <- clrX
    }
  } else if (tolower(mode) == "ilr"){  # ilr calculation, modified from ilr.transfo function from mixOmics package
    if (verbose) cat(paste0("ilr mode result in one less variable: ", ncol(x) - 1, " variables left for this dataset upon transformation.\n"))
    ilrX = matrix(NA, nrow = nrow(x), ncol = ncol(x) - 1)
    D = ncol(x)
    if (ilr.method.fast) {
      for (i in 1:ncol(ilrX)) {
        ilrX[, i] = sqrt((D - i)/(D - i + 1)) * log(((apply(as.matrix(x[,
                                                                        (i + 1):D, drop = FALSE]), 1, prod) + offset)^(1/(D -
                                                                                                                            i)))/(x[, i] + offset))
      }
    } else {
      for (i in 1:ncol(ilrX)) {
        ilrX[, i] = sqrt((D - i)/(D - i + 1)) * log(apply(as.matrix(x[,
                                                                      (i + 1):D]), 1, function(x) {
                                                                        exp(log(x))
                                                                      })/(x[, i] + offset) + offset)
      }
    }
    out <- as.matrix(ilrX)
  }

  return(out)
}


#' Title rbioseq_import_gtf
#'
#' @description Import GTF/GFF files
#' @param file GTF/GFF file. Add path if needed.
#' @param parallelComputing Wether to use parallel computing or not. Default is \code{TRUE}.
#' @param cluterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return A data frame with items from GTF/GFF file.
#' @details The following items are extracted from GTF/GFF file: \code{chromosome}, \code{gene_id}, \code{gene_type}, and \code{gene_name}.
#' @examples
#'
#' \dontrun{
#' gtf <- rbioseq_import_gtf(file = "gencode_v19.gtf")
#' }
#'
#' @export
rbioseq_import_gtf <- function(file,
                               verbose = TRUE){
  # open
  gtf_gff <- file(file)
  # gtf_gff <- file("/Users/jingzhang/OneDrive/my papers/my papers/(published)human_RNAseq_paper/dataset/gtf/gencode_v19.gtf")
  # check if the file can be read
  filecheck <- try(suppressWarnings(open(gtf_gff)), silent = TRUE)

  # load file
  if (any(class(filecheck) == "try-error")) {
    stop("Bad gtf/gff file.")
  } else {
    if (verbose) cat("Loading GTF/GFF file (speed depending on the hardware configurations)...")
    tmpfile <- scan(gtf_gff, what = "", quiet = TRUE, sep = "\n")
    # close connection
    close(gtf_gff)
    # tmpfile <- data.table::fread(file = file, sep = "\n", header = FALSE)
    tmpfile <- tmpfile[-c(1:5)]
    if (verbose) cat("Done!\n")
    if (verbose) cat("Parsing annotation information (speed depending on the hardware configurations)...")

    out_mtx <- matrix(ncol = 8, nrow = length(tmpfile))
    for (i in seq(length(tmpfile))) {
      tmp <- tmpfile[i]
      # below: column 1-8 parsing
      tmp1 <- unlist(strsplit(tmp, split = "\t"))[1]
      tmp4_start <- as.numeric(unlist(strsplit(tmp, split = "\t"))[4])
      tmp5_end <- as.numeric(unlist(strsplit(tmp, split = "\t"))[5])
      tmp_length <- tmp5_end - tmp4_start + 1

      # below: column 9 parsing
      tmp9 <- unlist(strsplit(tmp, split = "\t"))[9]  # column 9, which contains all the annotation info
      tmp9 <- gsub("\"", "", tmp9)  # remove the \" pattern
      tmp9 <- gsub("; ", ";", tmp9)  # remove the first space in each string
      tmp9 <- unlist(strsplit(tmp9, split = ";", fixed = TRUE))
      tmp9 <- tmp9[c(1, 2, 3, 5)]  # only the fisrt 8 are useful
      tmp9 <- strsplit(tmp9, split = " ", fixed = TRUE)
      tmp_colnames9 <- vector(length = 4)  # will be used outside of the loop
      tmpout9 <- vector(length = 4)
      for (j in 1:length(tmpout9)){
        tmp_colnames9[j] <- tmp9[[j]][1]
        tmpout9[j] <- tmp9[[j]][2]
      }
      # output
      tmpout <- c(tmpout9, tmp1, tmp4_start, tmp5_end, tmp_length)
      out_mtx[i, ] <- tmpout
    }
    out_colnames <- c(tmp_colnames9, "chromosome", "start", "end", "length")
    colnames(out_mtx) <- out_colnames
    out_mtx <- unique(out_mtx)
    out_dfm <- as.data.frame(out_mtx)
    if (verbose) cat("Done!\n")
  }

  # output
  if (verbose) cat(paste(nrow(out_dfm), " records sucessfully loaded from the inoput GTF/GFF file.", sep = ""))
  return(out_dfm)
}


#' @title rbioseq_import_count
#'
#' @description Data pre-processing for RNA-seq read count files.
#' @param path Path to raw files. Default is the system working directory.
#' @param species Optional species code, following the traditional abbreviated naming convention, e.g. "hsa", "mmu".
#' @param target.annot.file Annotation file describing filenames and targets, and should be in \code{csv} format.
#' @param sample_id.var.name Sample id variable name in the \code{target.annot.file}.
#' @param sample_groups.var.name Sample group annotation variable name in the \code{target.annot.file}.
#' @param gtf matrix or data.frame. Parsed gtf/gff annotation. Can be obtained by function \code{\link{rbioseq_import_gtf}}.
#' @param raw.file.ext Raw file extention. Default is \code{".txt"}.
#' @param raw.file.sep Raw read count file separators. Default is \code{""\"\"}, i.e. white space.
#' @param raw.file.source Raw file source, i.e. program used to generate read counts. Currently only supports \code{"htseq-count"}.
#' @param parallelComputing Wether to use parallel computing or not. Default is \code{TRUE}.
#' @param cluterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @details When \code{raw.file.source = "htseq-count"}, the function will cut off the last five summary raws.
#'
#'          For \code{target.annot.file}, the argument doesn't accept full file path.
#'          The function will only seek the file under working directory. So, the file should be placed under working directory.
#'
#'          Since the HTSeq-count program uses GTF/GFF annotation file for read couting,
#'          the results will always contain "\code{gene_id}" as the gene identification item.
#'          Therefore, when and \code{count_source = "htseq-count"} and \code{gtf.matrix} is set,
#'          the rest of the GTF/GFF information is merged into the \code{genes} item in the resulting \code{rbioseq_count} class object.
#'
#'          The items from GTF/GFF information are as following:
#'
#'          \code{gene_name}
#'
#'          \code{gene_type}
#'
#'          \code{chromosome}
#'
#'          \code{start}
#'
#'          \code{end}
#'
#'          \code{length}
#'
#'          Since the current HTSeq-count setting is to examine genes, NOT transcript.
#'          The \code{transcript_id} item is used to find the gene length.
#'          Specifically, the gene length is the length from the record where \code{transcript_id == gene_id}
#'
#'
#'          Transcript assessment will be added through future updates.
#'
#' @return Outputs a \code{rbioseq_count} object with merged read counts from mutliple files, with annotation. The \code{rbioseq_count} object contains the following:
#'
#'          \code{raw_read_count}
#'
#'          \code{sample_library_sizes}
#'
#'          \code{targets}: Sample annotation matrix
#'
#'          \code{sample_groups}: Factor object for sample group annotation
#'
#'          \code{genes}: The associated feature names. The use of "gene" here is in a generic sense.
#'
#'          \code{count_source}: program used to generate the reads
#'
#'          \code{GTF_annotation}: if GTF annoation matrix was used
#'
#'          \code{species}
#'
#'          \code{files_processed}
#'
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @examples
#' \dontrun{
#' mrna_count <- rbioseq_import_count(path = "~/dataset/",
#'                                    species = "hsa",
#'                                    target.annot.file = "target.csv", sample_groups.var.name = "condition",
#'                                    gtf = gtf,
#'                                    raw.file.ext = ".out", raw.file.sep = "",
#'                                    raw.file.source = "htseq-count",
#'                                    parallelComputing = TRUE, clusterType = "FORK")
#' }
#' @export
rbioseq_import_count <- function(path = getwd(), species = NULL,
                                 target.annot.file = NULL,
                                 sample_id.var.name = NULL,
                                 sample_groups.var.name = NULL,
                                 gtf = NULL,
                                 raw.file.ext = ".txt", raw.file.sep = "", raw.file.source = "htseq-count",
                                 parallelComputing = FALSE, clusterType = "FORK",
                                 verbose = TRUE){
  ## check argument
  if (is.null(target.annot.file)){  # check and load target (sample) annotation
    stop("Please provide a target annotation file for target.annot.file arugment.")
  } else {
    target.annot_name_length <- length(unlist(strsplit(target.annot.file, "\\.")))
    target.annot_ext <- unlist(strsplit(target.annot.file, "\\."))[target.annot_name_length]
    if (target.annot_ext != "csv") {
      stop("target.annot.file is not in csv format.")
    } else {
      if (verbose) cat("Loading target annotation file...")
      tgt <- read.csv(file = target.annot.file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
      if (verbose) cat("Done!\n")
    }
  }

  if (is.null(sample_groups.var.name)) stop("Please provide sample_groups.var.name.")
  if (is.null(sample_id.var.name)) stop("Please provide sample_id.var.name.")
  if (!all(c(sample_groups.var.name, sample_id.var.name) %in% names(tgt))){
    stop("Sample id or group annotation variables not found in the target annotation file.")
  } else {
    sample.groups <- factor(tgt[, sample_groups.var.name], levels = unique(tgt[, sample_groups.var.name]))
  }

  ## load files
  # set read files
  filename <- list.files(path = path, pattern = raw.file.ext)
  filename_wo_ext <- sub("[.][^.]*$", "", filename)  # general expression to remove extension, i.e. a.b.c becomes a.b

  # load reads
  if (verbose) cat("Processing read count files...")

  raw_list <- vector(mode = "list", length = length(filename))
  if (!parallelComputing){ # single core
    raw_list[] <- foreach(i = seq(length(filename))) %do% {
      tmp <- read.table(file = paste0(path, "/", filename[i]), header = FALSE, sep = raw.file.sep, stringsAsFactors = FALSE,
                        col.names = c("gene_id", filename_wo_ext[i]), row.names = NULL)
      if (raw.file.source == "htseq-count") tmp <- head(tmp, n = -5)  # remove last five rows from htseq-count results
      tmp
    }
  } else {  # parallel
    # set clusters
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores, clusterType = clusterType)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    # file processing
    raw_list[] <- foreach(i = seq(length(filename)), .packages = "foreach") %dopar% {
      tmp <- read.table(file = paste0(path, "/", filename[i]), header = FALSE, sep = raw.file.sep, stringsAsFactors = FALSE,
                        col.names = c("gene_id", filename_wo_ext[i]), row.names = NULL)
      if (raw.file.source == "htseq-count") tmp <- head(tmp, n = -5)  # remove last five rows from htseq-count results
      tmp
    }
  }
  names(raw_list) <- filename_wo_ext
  if (verbose) cat("Done!\n")

  # file loaded message
  if (verbose) cat("\n")
  if (verbose) cat("Files loaded: \n")
  for (i in filename){
    if (verbose) cat(paste0("\t", i, "\n"))
  }

  # merge reads
  out_dfm <- Reduce(function(i, j)merge(i, j, all = TRUE), raw_list)
  out_dfm[is.na(out_dfm) == TRUE] <- 0

  # check and load gtf/gff annotation
  if (is.null(gtf)){
    features <- out_dfm[, 1]
  } else {
    if (any(class(gtf) %in% "matrix")) {
      gtf_dfm <- as.data.frame(gtf, stringsAsFactors = FALSE)
    } else {
      gtf_dfm <- gtf
    }

    gtf_dfm_working <- gtf_dfm[gtf_dfm$transcript_id %in% out_dfm$gene_id, ]
    feature_out_dfm <- merge(gtf_dfm_working, out_dfm)
    features <- feature_out_dfm[, c("gene_id", "gene_type", "gene_name", "chromosome", "start", "end", "length")]
  }

  ## output
  counts <- out_dfm[, -1]
  counts <- as.matrix(counts)
  counts <- counts[, tgt[, sample_id.var.name]]
  lib_size <- colSums(counts)
  out <- list(raw_read_count = counts,
              sample_library_sizes = lib_size,
              targets = tgt,
              sample_groups = sample.groups,
              genes = features,
              count_source = raw.file.source,
              GTF_annotation = ifelse(is.null(gtf), FALSE, TRUE),
              species = species,
              files_processed = filename)
  class(out) <- "rbioseq_count"
  return(out)
}


#' @export
print.rbioseq_count <- function(x, ...){
  cat("RNAseq raw reads processing summary:\n")
  cat("\n")
  cat(paste0(" Total number of genomic features: ", ifelse(x$GTF_annotation, nrow(x$genes), length(x$genes)), "\n"))
  cat("\n")
  cat(paste0(" Files read: ", "\n"))
  cat(paste0(" ", x$files_processed, collapse = "\n"))
  cat("\n\n")
}


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
  count <- rbioseq_import_count(path = raw.file.path,
                                raw.file.ext = raw.file.ext,
                                raw.file.sep = raw.file.sep,
                                raw.file.source = raw.file.source,
                                species = species,
                                target.annot.file = target.annot.file,
                                sample_groups.var.name = sample_groups.var.name,
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
