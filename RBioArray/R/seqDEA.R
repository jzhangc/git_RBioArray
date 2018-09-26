#' @title rbioseq_clr_ilr_transfo
#'
#' @description Log ratio tansformation function for read count data. Row: sample, column: features.
#' @param x Input read count data matrix.
#' @param offset Read count offset value added to avoid zero. Default is \code{1}.
#' @param mode Log ratio transformation method. Options are "clr" (centered log transformation) and "ilr" (isometric log transformation). Default is \code{"clr"}.
#' @param ilr.method.fast Useful only when \code{mode = "ilr"}. Default is \code{TRUE}.
#' @return A data matrix with log ratio transformed values.
#' @details This funciton is needed as part of the data pre-processing procedure to run multivariate and machine learning analysis featured in \code{RBioFS} package.
#'
#'          As per Quinn et al. (2018), NGS data can be considered as compositional data. As such, data must be transformed for usual statistical analysis and visualization.
#'          Log ratio transfromation serves such purpose. Note that the number of features will be one less when using "ilr" method.
#'
#'          It is not to be combined with the other normalization methods featured in the \code{\link{rbioseq_DE}}.
#'
#'          Ref: Quinn TP, et al. 2018. Understanding sequencing data as compositions: an outlook and review. Bioinformatics. 2018: 1 - 9.
#' @examples
#' \dontrun{
#' tstX <- rbioseq_clr_ilr_transfo(tstdata, offset = 1, mode = "clr")
#' }
#' @export
rbioseq_clr_ilr_transfo <- function(x, offset = 1, mode = "clr", ilr.method.fast = TRUE){
  # data and arguments check
  if (!is.matrix(x))stop("x needs to b e a matrix")
  if (any(x == 0) & offset == 0)stop("zero detected in x. set offset to avoid it for ratio transformation")
  if (!tolower(mode) %in% c("clr", "ilr"))stop("choose the proper transformation mode: \"clr\" or \"ilr\"")

  # log ratio transformation
  if (tolower(mode) == "clr"){  # clr calculation
    if (dim(x)[2] == 1){
      cat()
      out <- list(x.clr = x, gm = rep(1, dim(x)[1]))
    } else {
      gm <- apply(x, 1, function(x) exp(mean(log(x + offset))))  # geometric mean = exp(mean(log(X)))
      clrX <- log((x + offset) / (gm))  # clr formula
      out <- clrX
    }
  } else if (tolower(mode) == "ilr"){  # ilr calculation, modified from ilr.transfo function from mixOmics package
    cat(paste0("ilr mode result in one less variable: ", ncol(x) - 1, " variables left for this dataset upon transformation.\n"))
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
#' @return A matrix with items from GTF/GFF file.
#' @details The following items are extracted from GTF/GFF file: \code{chromosome}, \code{gene_id}, \code{gene_type}, and \code{gene_name}.
#' @examples
#'
#' \dontrun{
#' gtf <- rbioseq_import_gtf(file = "gencode_v19.gtf")
#' }
#'
#' @export
rbioseq_import_gtf <- function(file){
  # open connection
  gtf_gff <- file(file)

  # check if the file can be read
  filecheck <- try(suppressWarnings(open(gtf_gff)), silent = TRUE)
  if (class(filecheck) == "try-error") {
    stop("Bad gtf/gff file.")
  } else {
    cat("Loading GTF/GFF file (speed depending on the hardware configurations)...")
    tmpfile <- scan(gtf_gff, what = "", quiet = TRUE, sep = "\n")
    # tmpfile <- data.table::fread(file = file, sep = "\n", header = FALSE)
    tmpfile <- tmpfile[-c(1:5)]
    cat("Done!\n")
    cat("Parsing annotation information (speed depending on the hardware configurations)...")
    out_mtx <- matrix(ncol = 4, nrow = length(tmpfile))
    for (i in seq(length(tmpfile))) {
      tmp <- tmpfile[i]
      # below: column 1-8 parsing
      tmp1 <- unlist(strsplit(tmp, split = "\t"))[1]

      # below: column 9 parsing
      tmp9 <- unlist(strsplit(tmp, split = "\t"))[9]  # column 9, which contains all the annotation info
      tmp9 <- gsub("\"", "", tmp9)  # remove the \" pattern
      tmp9 <- gsub("; ", ";", tmp9)  # remove the first space in each string
      tmp9 <- unlist(strsplit(tmp9, split = ";", fixed = TRUE))
      tmp9 <- tmp9[c(1, 3, 5)]  # only the fisrt 8 are useful
      tmp9 <- strsplit(tmp9, split = " ", fixed = TRUE)
      tmp_colnames9 <- vector(length = 3)  # will be used outside of the loop
      tmpout9 <- vector(length = 3)
      for (j in 1:length(tmpout9)){
        tmp_colnames9[j] <- tmp9[[j]][1]
        tmpout9[j] <- tmp9[[j]][2]
      }
      # output
      tmpout <- c(tmpout9, tmp1)
      out_mtx[i, ] <- tmpout
    }
    out_colnames <- c(tmp_colnames9, "chromosome")
    colnames(out_mtx) <- out_colnames
    out_mtx <- unique(out_mtx)
    cat("Done!\n")
  }
  # close connection
  close(gtf_gff)

  # output
  cat(paste(nrow(out_mtx), " records sucessfully loaded from the inoput GTF/GFF file.", sep = ""))
  return(out_mtx)
}


#' @title rbioseq_import_count
#'
#' @description Data pre-processing for RNA-seq read count files.
#' @param path Path to raw files. Default is the system working directory.
#' @param species Optional species code, following the traditional abbreviated naming convention, e.g. "hsa", "mmu".
#' @param target.annot.file Annotation file describing filenames and targets, and should be in \code{csv} format.
#' @param gtf.matrix Parsed gtf/gff annotation matirx. Can be obtained by function \code{\link{rbioseq_import_gtf}}.
#' @param raw.file.ext Raw file extention. Default is \code{".txt"}.
#' @param raw.file.sep Raw read count file separators. Default is \code{""\"\"}, i.e. white space.
#' @param raw.file.source Raw file source, i.e. program used to generate read counts. Currently only supports \code{"htseq-count"}.
#' @param parallelComputing Wether to use parallel computing or not. Default is \code{TRUE}.
#' @param cluterType clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @details When \code{raw.file.source = "htseq-count"}, the function will cut off the last five summary raws.
#' @return Outputs a \code{rbioseq_count} object with merged read counts from mutliple files, with annotation. The \code{rbioseq_count} object contains the following:
#'
#'          \code{raw_read_count}
#'
#'          \code{sample_library_sizes}
#'
#'          \code{targets}: Sample annotation
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
#'                                    target.annot.file = "target.csv", gtf.matrix = gtf,
#'                                    raw.file.ext = ".out", raw.file.sep = "",
#'                                    raw.file.source = "htseq-count",
#'                                    parallelComputing = TRUE, clusterType = "FORK")
#' }
#' @export
rbioseq_import_count <- function(path = getwd(), species = NULL,
                                 target.annot.file = NULL, gtf.matrix = NULL,
                                 raw.file.ext = ".txt", raw.file.sep = "", raw.file.source = "htseq-count",
                                 parallelComputing = FALSE, clusterType = "FORK"){
  ## check argument
  if (is.null(target.annot.file)){  # check and load target (sample) annotation
    tgt <- NULL
  } else {
    target.annot_name_length <- length(unlist(strsplit(target.annot.file, "\\.")))
    target.annot_ext <- unlist(strsplit(target.annot.file, "\\."))[target.annot_name_length]
    if (target.annot_ext != "csv") {
      cat("target.annot.file is not in csv format. Proceed without using the file.\n")
      tgt <- NULL
    } else {
      cat("Loading target annotation file...")
      tgt <- read.csv(file = target.annot.file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
      cat("Done!\n")
    }
  }

  ## load files
  # set read files
  filename <- list.files(path = path, pattern = raw.file.ext)
  filename_wo_ext <- sub("[.][^.]*$", "", filename)  # general expression to remove extension, i.e. a.b.c becomes a.b

  # load reads
  cat("Processing read count files...")

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
  cat("Done!\n")

  # file loaded message
  cat("\n")
  cat("Files loaded: \n")
  for (i in filename){
    cat(paste0("\t", i, "\n"))
  }

  # merge reads
  out_dfm <- Reduce(function(i, j)merge(i, j, all = TRUE), raw_list)
  out_dfm[is.na(out_dfm) == TRUE] <- 0

  # check and load gtf/gff annotation
  if (is.null(gtf.matrix)){
    features <- out_dfm[, 1]
  } else {
    gtf_dfm <- data.frame(gtf.matrix, stringsAsFactors = FALSE)
    feature_out_dfm <- merge(gtf_dfm, out_dfm)
    features <- feature_out_dfm[, c("gene_id", "gene_type", "gene_name", "chromosome")]
  }
  ## output
  counts <- out_dfm[, -1]
  counts <- as.matrix(counts)
  lib_size <- colSums(counts)
  out <- list(raw_read_count = counts,
              sample_library_sizes = lib_size,
              targets = tgt,
              genes = features,
              count_source = raw.file.source,
              GTF_annotation = ifelse(is.null(gtf.matrix), FALSE, TRUE),
              species = species,
              files_processed = filename)
  class(out) <- "rbioseq_count"
  return(out)
}


#' @export
print.rbioseq_count <- function(x, ...){
  cat("RNAseq raw reads processing summary:\n")
  cat("\n")
  cat(paste0(" Total number of mRNA: ", ifelse(x$GTF_annotation, nrow(x$genes), length(x$genes)), "\n"))
  cat("\n")
  cat(paste0(" Files read: ", "\n"))
  cat(paste0(" ", x$files_processed, collapse = "\n"))
  cat("\n\n")
}
