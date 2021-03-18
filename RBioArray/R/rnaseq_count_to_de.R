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


#' Title rbioseq_gtf
#'
#' @description Import GTF/GFF files
#' @param file GTF/GFF file. Add path if needed.
#' @param parallelComputing Whether to use parallel computing or not. Default is \code{TRUE}.
#' @param cluterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param n_cores Only set when \code{parallelComputing = TRUE}, number of cores to use. When \code{NULL}, the function uses max-1.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return A data frame with items from GTF/GFF file.
#' @details The following items are extracted from GTF/GFF file: \code{gene_id}, \code{gene_name}, \code{transcript_id}, \code{transcript_name}, \code{chr}, \code{type}, .
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom data.table fread
#' @examples
#'
#' \dontrun{
#' gtf <- rbioseq_gtf(file = "gencode_v19.gtf")
#' }
#'
#' @export
rbioseq_gtf <- function(file, verbose = TRUE, parallelComputing = FALSE, clusterType = "FORK", n_cores = NULL){
  # -- check file --
  f_check <- file(file)
  filecheck <- try(suppressWarnings(open(f_check)), silent = TRUE)
  close(f_check)
  if (any(class(filecheck) == "try-error")) {
    stop("Bad gtf file.")
  }

  # -- open file --
  if (verbose) cat("Loading GTF file (speed depending on the hardware configurations)...")
  gtf <- suppressWarnings(fread(file, verbose = FALSE))
  setnames(gtf , names(gtf), c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes") )
  if (verbose) cat("Done!\n")
  if (verbose) cat("Parsing annotation information (speed depending on the hardware configurations)...")
  # -- attributes --
  att <- gtf$attributes
  attr_names <- c("gene_id", "gene_name" , "transcript_id", "transcript_type", "transcript_name")
  if (parallelComputing) {
    # set clusters
    if (!is.null(n_cores)) {
      n_cores = n_cores
    } else {
      n_cores <- detectCores() - 1
    }
    cl <- makeCluster(n_cores, clusterType = clusterType)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    # computing
    att_df <- foreach(i = att, .combine = "rbind") %dopar% {
      input_r <- i
      tmpatt <- strsplit(input_r, "; ")
      tmpatt <- gsub("\"", "", unlist(tmpatt))

      tmp_o <- foreach(j = attr_names, .combine = "cbind") %do% {
        if(!is.null(unlist(strsplit(tmpatt[grep(j, tmpatt)], " ")))){
          return(unlist(strsplit(tmpatt[grep(j, tmpatt)], " "))[2])
        }else{
          return(NA)
        }
      }
      tmp_o
    }
  } else {
    att_df <- foreach(i = att, .combine = "rbind") %do% {
      input_r <- i
      tmpatt <- strsplit(input_r, "; ")
      tmpatt <- gsub("\"", "", unlist(tmpatt))

      tmp_o <- foreach(j = attr_names, .combine = "cbind") %do% {
        if(!is.null(unlist(strsplit(tmpatt[grep(j, tmpatt)], " ")))){
          return(unlist(strsplit(tmpatt[grep(j, tmpatt)], " "))[2])
        }else{
          return(NA)
        }
      }
      tmp_o
    }
  }
  colnames(att_df) <- attr_names
  if (verbose) cat("Done!\n")

  # -- out --
  out <- data.frame(cbind(att_df, gtf[, c("chr", "type")]))
  out <- unique(out)
  return(out)
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
#' @param raw.file.source Raw file source, i.e. program used to generate read counts. Currently only supports \code{"htseq-count"}.
#' @param raw.file.ext Raw file extension. Default is \code{".txt"}.
#' @param raw.file.sep Raw read count file separators. Default is \code{""\"\"}, i.e. white space.
#' @param parallelComputing Whether to use parallel computing or not. Default is \code{TRUE}.
#' @param cluterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @details When \code{raw.file.source = "htseq-count"}, the function will cut off the last five summary rows.
#'
#'          For \code{target.annot.file}, the argument doesn't accept full file path.
#'          The function will only seek the file under working directory. So, the file should be placed under working directory.
#'
#'          Since the HTSeq-count program uses GTF/GFF annotation file for read counting,
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
  if (!any(class(object) %in% c("rbioseq_count", "mir_count"))) stop("object needs to be either a \"rbioseq_count\" or \"mir_count\" object")

  ## use methods
  UseMethod("rnaseq_de", object)
}


#' @title rnaseq_de.mir_count
#'
#' @rdname rnaseq_de
#' @method rnaseq_de mir_count
#' @param object A \code{mir_count} object from the \code{mirProcess} function of \code{RBioMIR} package.
#' @param filter.threshold.min.count Minimum count for the smallest library for filter threshold. Default is \code{10}.
#' @param filter.threshold.min.sample Minimum number of samples meeting the count threshold. Default is \code{NULL}.
#' @param ... Additional arguments for \code{\link{rnaseq_de.default}}.
#'
#' @export
rnaseq_de.mir_count <- function(object, filter.threshold.min.count = 10, filter.threshold.min.sample = NULL, ...){
  ## for setting up comparison groups info and minimum sample number
  annot.group <- object$sample_groups
  if (within.sample.norm.method == "rpkm") stop("within.sample.norm.method = \"rpkm\" is not supported for mir_count object")

  ## construct rbioseq_de object
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
#' @param ... Additional arguments for \code{\link{rnaseq_de.default}}.
#'
#' @export
rnaseq_de.rbioseq_count <- function(object, filter.threshold.min.count = 10, filter.threshold.min.sample = NULL, ...){
  ## for setting up comparison groups info and minimum sample number
  annot.group <- object$sample_groups

  ## construct rbioseq_de object
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
#' @param annot.group Sample group annotation object. Can be a \code{factor} or \code{vector} object.
#' @param library.size.scale.method aka between sample normalization method. Options are: \code{"none"}, \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"}. Default is \code{"TMM"}.
#' @param design Design matrix.
#' @param contra Contrast matrix.
#' @param qc.plot QC plot for the input read counts
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @details \code{filter.threshold.cpm.count} uses CPM (counts per million) as the basis for filtering. The rule of thumb is to filter reads independently from the groupping information.
#'
#'          The default is based on the paper by Chen et al (2016):
#'          10~15 counts per library, which translate to 10~15/L with L being the smallest library size in millions.
#'          Therefore, if using 10 as an example, the above becomes 10 / min(library_sizes) / 1000000
#'
#'          Chen W, Lun ATL, Smyth GK. 2016. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments
#'          using Rsubread and the edgeR quasi-likelihood pipeline. F1000Research. 5:1438.
#'
#'
#'          It absolutely important to note that VOOM (or any other full-fledged RNA-seq methods like DEseq2) DOES NOT use FPKM for DE analysis
#'          VOOM uses logCPM with library size scaling, e.g. TMM.
#'          In fact, FPKM/RPKM/TPM are not suitable for DE analysis as they are designed to compare WITHIN sample genes, i.e. gene A vs gene B in one sample.
#'          Further, users DO NOT need to convert their raw reads beforehand. VOOM will do the conversion for you.
#'
#'          SO, ONLY supply raw read counts to the function.
#'
#' @return A list containing core elements of a \code{rbioseq_de} class. The items of a \code{rbioseq_de} class are following:
#'
#'         \code{filter_results}: A list containing \code{filter_threshold_cpm}, \code{filter_threshold_min_sample}, \code{filter_summary} and \code{filtered_counts}.
#'                                Note: \code{filtered_counts} is a \code{DGEList} class from \code{edgeR} package.
#'
#'         \code{normalization_method}
#'
#'         \code{calcNormFactors_output}: \code{edgeR} function \code{calcNormFactors} output. A DGElist with \code{object$samples$norm.factors}.
#'                                        Good for \code{cpm()} conversion for PCA, clustering, etc. The data is filtered.
#'                                        A "filtered raw count with lib size" object.
#'
#'         \code{voom_output}: A \code{EList} generated from \code{voomWithQualityWeights} function from \code{limma} package.
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
#' @import foreach
#' @importFrom limma lmFit eBayes topTable contrasts.fit voomWithQualityWeights
#' @importFrom edgeR DGEList calcNormFactors cpm rpkm
#'
#' @export
rnaseq_de.default <- function(x, y = NULL,
                              y.gene_id.var.name = "genes",
                              y.gene_symbol.var.name = "genes",
                              filter.threshold.cpm = "none",
                              filter.threshold.min.sample = NULL, annot.group = NULL,
                              library.size.scale.method = c("TMM", "RLE", "upperquartile", "none"),
                              design, contra, qc.plot = TRUE, verbose = TRUE){
  ## check the key arguments
  if (!any(class(x) %in% c("data.frame", "matrix"))){
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
  # no need to put choices argument below as the match.arg grabs from the top
  library.size.scale.method <- match.arg(library.size.scale.method)

  ## extract coefficients
  cf <- colnames(contra) # extract coefficient
  contra_levels <- vector(mode = "list", length = length(cf))
  contra_levels[] <- foreach(i = seq(length(cf))) %do% {
    rownames(contra)[which(contra[, i] != 0)]
  }
  names(contra_levels) <- cf

  ## DE
  if (verbose) cat("Data filtering and normalization...") # message
  dge <- DGEList(counts = x, genes = y, group = annot.group)

  if (filter.threshold.cpm != "none"){ # set the count threshold for filtering
    isexpr <- rowSums(cpm(dge$counts) > filter.threshold.cpm) >= filter.threshold.min.sample  # cpm threshold, cite the paper
    # isexpr <- rep(FALSE, times = 611)
    dge <- dge[isexpr, , keep.lib.size = FALSE] # filtering

    # isexpr <- c(TRUE, TRUE, TRUE)
    # isexpr <- c(FALSE, FALSE, FALSE)
    # isexpr <- c(FALSE, FALSE, TRUE)

    if (all(isexpr)) {  # if all TRUE
      warning('Nothing is filtered out. \n')
      flt_summary <- table(isexpr)
      flt_summary['FALSE'] <- 0
      flt_summary <- sort(flt_summary)  # make sure always this order: FALSE, TRUE
    } else if (all(!isexpr)){ # if all FALSE
      stop('Nothing remained after filtering. Function stopped.')
      flt_summary <- table(isexpr)
      flt_summary['TRUE'] <- 0
    } else {
      flt_summary <- table(isexpr)
    }

    names(flt_summary) <- c("filtered", "remaning")
    filter_results <- list(filter_threshold_cpm = filter.threshold.cpm,
                           filter_threshold_min_sample = filter.threshold.min.sample,
                           filter_summary = flt_summary,
                           filtered_counts = dge)
  } else {
    cat("No filtering applied as filter.threshold.cpm = none. \n")
    filter_results <- NULL
  }

  # library size scaling normalization: within sample
  # new_counts <- switch(within.sample.norm.method,
  #                      cpm = cpm(dge$counts),
  #                      rpkm = rpkm(dge$counts, gene.length = as.numeric(y$length)[isexpr]),
  #                      none = dge$counts)
  # dge$counts <- new_counts
  dgenormf <- calcNormFactors(dge, method = library.size.scale.method)  # between samples

  # between-genes: Voom normalization with quality weights
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
              normalization_method = list(within_sample_method = "(automatic) CPM, logCPM",
                                          library_scalling = library.size.scale.method,
                                          between_genes = "voom process with quantile normalization"),
              calcNormFactors_output = dgenormf,
              voom_output = vmwt,
              genes_annotation.gene_id.var_name = y.gene_id.var.name,
              genes_annotation.gene_symbol.var_name = y.gene_symbol.var.name,
              F_stats = f_stats,
              DE_results = de_list,
              comparisons = comparisons)
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
  cat(paste0("\tWithin sample method (automatic): CPM, logCPM. \n"))
    cat(paste0("\tLibrary size-scaling: ", x$normalization_method$library_scalling, "\n"))
  cat(paste0("\tBetween-genes: ", x$normalization_method$between_genes, "\n"))
  cat("\n")
  cat("Comparisons assessed: \n")
  cat(paste0("\t", x$comparisons$comparisons, "\n"))
  cat("\n")
}
