#' @title rbioarray_rlist
#'
#' @description Function to constuct \code{rlist} class object from microarary and annotation data. The \code{rlist} object is the starting point for all microarray data analysis.
#' @param raw.dataframe Input data frame containing microarray hybridization signals with rows as probe/gene/genomic features and columns as samples. Note: the data frame should contain at least one annotation column.
#' @param raw.annot.var.name A string vector containing variable (i.e. column) name(s) for all the annotation columns in \code{raw.dataframe}.
#' @param raw.gene_id.var.name Variable (i.e. column) name for gene/probe/genomic feature identification from \code{raw.dataframe}.
#' @param gene.annot.dataframe Optional annotation data frame for gene/probe/genomic feature annotation.
#' @param gene.annot.gene_id.var.name Set only when \code{gene.annot.dataframe} is provided, variable name for probe/gene/genomic features identification from \code{gene.annot.dataframe}.
#' @param gene.annot.gene_symbol.var.name Set only when \code{gene.annot.dataframe} is provided, variable name for probe/gene/genomic features display name from \code{gene.annot.dataframe}, e.g. gene symbols.
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
#'          To keep things consistent with the \code{Elist} from the dependent \code{limma} package. The \code{rlist} class contains many common elements from \code{Elist} class.
#'
#' @return A \code{rlist} object. The \code{rlist} class contains the following items:
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
    cat("Note: gene.annot.dataframe not provided. Proceed with raw.dataframe annoation information. \n")
    gene.annot.dataframe <- raw.dataframe[, raw.annot.var.name]
    gene.annot.gene_id.var.name <- raw.gene_id.var.name
    gene.annot.gene_symbol.var.name <- raw.gene_id.var.name
    gene.symbol <- FALSE
  } else {
    if (!is.data.frame(gene.annot.dataframe)) stop("gene.annot needs to be a dataframe")
    if (nrow(gene.annot.dataframe) < nrow(raw.dataframe)) stop("gene.annot.dataframe has less record than the raw.dataframe") # check size
    if (is.null(gene.annot.gene_id.var.name) || is.null(gene.annot.gene_symbol.var.name)) {
      stop("Please set gene.annot.gene_id.var.name AND gene.annot.gene_symbol.var.name arguments according to gene.annot")
    } else if (!gene.annot.gene_id.var.name %in% names(gene.annot.dataframe) || !gene.annot.gene_symbol.var.name %in% names(gene.annot.dataframe)) {
      stop("Gene ID variable and/or gene symbol variable not found in gene.annot")
    } else {
      gene.symbol <- TRUE
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
              targets = tgt,
              sample_groups = sample.groups)
  class(out) <- "rlist"
  cat("Done!\n")
  cat("\n")
  cat(paste0("The resulted rlist contains ", nrow(E), " genes/probes/genomic features, ", nrow(tgt), " samples for ", length(unique(sample.groups)), " groups."))
  return(out)
}


#' @export
print.rlist <- function(x, ...){
  cat("---- rlist information ----\n")
  cat(paste0("Number of genes/probes/genomic features: ", nrow(x$E), "\n"))
  cat(paste0("Number of samples: ", nrow(x$targets), "\n"))
  cat(paste0("Groups: "))
  cat(paste0(levels(x$sample_groups)))
  cat("\n\n")
  cat(paste0("Gene display name from annotation information: ", ifelse(x$gene_display_name_used, "Available\n", "Unavailable\n")))
  cat("\n")
}


#' @title rbioarray_PreProc
#'
#' @description Data pre-processing function for the microarary data.
#' @param rawlist Input data, either a list, \code{EList} or \code{MAList} object.
#' @param logTrans If to perfom a log transformation on the data or not. Default is \code{FALSE}.
#' @param logTransMethod If \code{logTrans = TRUE}, set which method to use for the transformation, \code{"log2"} or \code{"log10"}. Default is \code{"log2"}.
#' @param logTransObjT If \code{logTrans = TRUE}, set the file name for the output \code{csv} file containing the log transformed data.
#' @param logTransParallelComputing If \code{logTrans = TRUE}, set if to use parallel computing for the transformation or not. Default is \code{FALSE}.
#' @param bgMethod Background correction method. Default is \code{"auto"}. See \code{backgroundCorrect()} function from \code{limma} package for details.
#' @param normMethod Normalization method. Default is \code{"quantile"}. See \code{normalizeBetweenArrays()} function from \code{limma} package for details.
#' @param ... arguments for \code{backgroundCorrect.matrix()} or \code{backgroundCorrect()} functions from \code{limma} package.
#' @return Depending on the input type, the function outputs a \code{list}, \code{Elist} or \code{MAList} object with corrected and normalized expression values. If \code{logTrans = TRUE}, the function also outputs a \code{csv} file containing the log transformed data.
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom limma backgroundCorrect normalizeBetweenArrays backgroundCorrect.matrix arrayWeights
#' @examples
#' \dontrun{
#' normdata <- rbioarray_PreProc(mydata)
#' }
#' @export
rbioarray_PreProc <- function(rawlist, logTrans = FALSE, logTransMethod = "log2",
                              logTransObjT = "data", logTransParallelComputing = FALSE,
                              bgMethod = "auto", normMethod = "quantile", ...){
  if (class(rawlist) == "list"){
    ## log transform  or not
    if (logTrans){
      if (!logTransParallelComputing){
        # log transform
        mtx <- apply(rawlist$E, c(1,2), FUN = ifelse(logTransMethod == "log10", log10, log2))
      } else {
        # parallel computing
        # set up cpu cluster
        n_cores <- detectCores() - 1
        cl <- makeCluster(n_cores, type = "PSOCK")
        on.exit(stopCluster(cl)) # close connect when exiting the function
        # log transform
        mtx <- foreach(i = rawlist$E) %dopar% {
          out <- apply(i, c(1,2), FUN = ifelse(logTransMethod == "log10", log10, log2))
        }
      }
      tmpdata <- list(E = mtx, genes = rawlist$genes, targets = rawlist$targets)
      # store and export log transformed data into a csv file
      logTransOut <- data.frame(rawlist$genes, mtx)
      write.csv(logTransOut, file = paste(logTransObjT, "_log_transformed.csv", sep = ""), row.names = FALSE)
    } else {
      tmpdata <- rawlist
    }

    ## normalization
    BgC <- backgroundCorrect.matrix(tmpdata$E, method = bgMethod, ...) #background correction
    Norm <- normalizeBetweenArrays(BgC, normMethod) # quantile normalization
    Wgt <- arrayWeights(Norm) # array weight
    output <- list(E = Norm, genes = rawlist$genes, targets = rawlist$targets, ArrayWeight = Wgt)
  } else {
    ## normalization
    BgC <- backgroundCorrect(rawlist, method = bgMethod, ...) #background correction
    Norm <- normalizeBetweenArrays(BgC, method = normMethod) # quantile normalization
    Wgt <- arrayWeights(Norm)
    Norm$ArrayWeight <- Wgt
    output <- Norm
  }
  return(output)
}

#' @title rbioarray_flt
#'
#' @description data filter function based on spike-in negative control.
#' @param normlst Normalized data, either a list, \code{EList} or \code{MAList} object.
#' @param ctrlProbe Wether or not the data set has control type variable, with values \code{-1 (negative control)}, \code{0 (gene probes)} and \code{1 (positive control)}. Default is \code{TRUE}.
#' @param ctrlTypeVar Set only when \code{ctrlProbe = TRUE}, the control type variable. Default is the \code{Agilent} variable name \code{"ControlType"}.
#' @param percentile The percentile cutoff. When \code{ctrlProbe = TRUE} and muliptle negative control probes are detected, the default is \code{0.95}. When \code{ctrlProbe = FALSE}, default value is \code{0.05}.
#' @param combineGeneDup Wether or not to combine gene duplicates (different probe ID) by probe signal variance. Default is \code{FALSE}.
#' @param geneSymbolVar Set only when \code{combineGeneDup = TRUE}, the name for variables contained in normlst or annotation dataframe. Default is \code{NULL}.
#' @param annot Set only when \code{combineGeneDup = TRUE} and normlst is a \code{EList} object, the probe annotation dataframe.
#' @param parallelComputing Set only when \code{combineGeneDup = TRUE}, if to use parallel computing. Default is \code{FALSE}.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @details When ctrlProbe is present, the retained probes are the ones with expression value above 10 percent of the 95 percentile of the negative control probe signal by default. When \code{ctrlProbe = FALSE}, the probes with a expression value 10 percent higher than the 5 percentile of total expression values are retained by default. When \code{combineGeneDup = TRUE}, the probe with the highest variance in signal will be retained for the gene of interest.
#' @return Depending on the input type, the function outputs a \code{list}, \code{Elist} or \code{MAList} object with filtered expression values.
#' @importFrom limma avereps
#' @examples
#' \dontrun{
#' fltdata <- rbioarray_flt(normdata)
#' }
#' @export
rbioarray_flt <- function(normlst, ctrlProbe = TRUE, ctrlTypeVar = "ControlType",
                          percentile = ifelse(ctrlProbe, 0.95, 0.05),
                          combineGeneDup = FALSE, geneSymbolVar = NULL, annot = NULL,
                          parallelComputing = FALSE, clusterType = "PSOCK"){
  ## check key arguments
  if (combineGeneDup){
    if (is.null(geneSymbolVar)){
      stop(cat("Please set variable name for gene symbol when combineGeneDup = TRUE. Function terminated.\n"))
    }
  }

  if (!"ProbeName" %in% names(normlst$genes)){
    stop(cat("Make sure to name the variable containing probe name \"ProbeName\". Function terminated.\n"))
  }

  if (ctrlProbe){
    if (!ctrlTypeVar %in% names(normlst$genes)){
      stop(cat("Make sure to include the correct variable name for control probes. Function terminated.\n"))
    }
  }

  if (!class(normlst) == "list" & combineGeneDup){
    if(is.null(annot)){
      stop(cat("Since combineGeneDup = TRUE and norlst is an Elist objects, please set annotation dataframe for annot argument.
               Function terminated.\n"))
    }
    if (!"ProbeName" %in% names(annot)){
      stop(cat("Make sure to name the variable containing probe name \"ProbeName\" in annotation dataframe. Function terminated.\n"))
    }
    if (!geneSymbolVar %in% names(annot)){
      stop(cat("Make sure to name the variable containing gene symbols in annotation data.frame. Function terminated.\n"))
    }
  }

  ## extract the 95% percentile of the negative control signals
  if (ctrlProbe){ # if there are neg control probes
    if (class(normlst$E[normlst$genes[, ctrlTypeVar] == -1, ]) == "numeric"){ # if there is only one entry in the neg values
      neg <- normlst$E[normlst$genes[, ctrlTypeVar] == -1, ] # no 95% percentile required as only one neg entry
    } else {
      neg <- apply(normlst$E[normlst$genes[, ctrlTypeVar] == -1, ], 2, function(x)quantile(x, p = percentile)) # neg95
    }
  } else { # no neg control probes, we use the
    neg <- apply(normlst$E, 2, function(x)quantile(x, p = 0.05)) # 5% percentile of all the data
  }

  if (class(normlst) == "list"){
    ## low expression cuttoff set at at least 10% hihger than the neg
    LE_cutoff <- matrix(1.1 * neg, nrow(normlst$E), ncol(normlst$E), byrow = TRUE)

    ## summary after applying LE_cutoff (T/F for each sample)
    # this only compares the element Nrm$E
    isexpr <- rowSums(normlst$E > LE_cutoff) >= 3 # We keep probes that meet the criterion on at least 3 arrays (minimum for stats)

    ## filter out only the low expressed probes (not control) in the dataset
    # LE means low expression removed (only)
    flt_E <- normlst$E[isexpr, ] # this is a way of extracting samples logically considered TRUE by certain tests
    flt_gene <- normlst$gene[isexpr, ]
    fltlst <- list(E = flt_E, genes = flt_gene, targets = normlst$targets)
    flt_E_avg <- avereps(fltlst$E, ID = fltlst$genes$ProbeName)
    genes <- unique(fltlst$genes[fltlst$genes$ProbeName %in% rownames(flt_E_avg), ])

    ## to combine gene repeats by probe signel variance (retain maximum variance probes for the same gene)
    if (combineGeneDup){
      tmplst <- list(E = flt_E_avg, genes = genes)

      if (parallelComputing){  # parallel computing
        # set up clusters for PSOCK
        n_cores <- detectCores() - 1
        cl <- makeCluster(n_cores, type = clusterType, outfile = "")
        registerDoParallel(cl) # part of doParallel package
        on.exit(stopCluster(cl)) # close connect when exiting the function

        cat("Combining gene duplicates by signal variance...") # initiation message
        combGeneProbe <- foreach(i = unique(tmplst$genes[, geneSymbolVar]), .combine = "c", .packages = "foreach") %dopar% {
          tmp <- tmplst$E[which(tmplst$genes[, geneSymbolVar] %in% i), ]
          if (is.null(nrow(tmp))){
            out <- var(tmp)
          } else {
            out <- foreach(j = seq(nrow(tmp)), .combine = "c") %do% {
              var(as.numeric(tmp[j, ]))
            }
          }
          out <- data.frame(tmplst$genes[which(tmplst$genes[, geneSymbolVar] %in% i), ], var = out)
          out[which(out$var %in% max(out$var)), "ProbeName"]
        }
        cat("DONE!\n") # termination message
      } else {  # single core
        cat("Combining gene duplicates by signal variance...") # initiation message
        combGeneProbe <- foreach(i = unique(tmplst$genes[, geneSymbolVar]), .combine = "c", .packages = "foreach") %do% {
          tmp <- tmplst$E[which(tmplst$genes[, geneSymbolVar] %in% i), ]
          if (is.null(nrow(tmp))){
            out <- var(tmp)
          } else {
            out <- foreach(j = seq(nrow(tmp)), .combine = "c") %do% {
              var(as.numeric(tmp[j, ]))
            }
          }
          out <- data.frame(tmplst$genes[which(tmplst$genes[, geneSymbolVar] %in% i), ], var = out)
          out[which(out$var %in% max(out$var)), "ProbeName"]
        }
        cat("DONE!\n") # termination message
      }
      avgProbes <- list(E = flt_E_avg[which(genes$ProbeName %in% combGeneProbe), ],
                        genes = genes[which(genes$ProbeName %in% combGeneProbe), ],
                        targets = normlst$targets, ArrayWeight = normlst$ArrayWeight)
    } else {
      avgProbes <- list(E = flt_E_avg, genes = genes,
                        targets = normlst$targets, ArrayWeight = normlst$ArrayWeight)
    }

  } else {
    ## low expression cuttoff set at at least 10% hihger than the neg95
    LE_cutoff <- matrix(1.1 * neg, nrow(normlst), ncol(normlst), byrow = TRUE)

    ## summary after applying LE_cutoff (T/F for each sample)
    # this only compares the element Nrm$E
    isexpr <- rowSums(normlst$E > LE_cutoff) >= 3 # We keep probes that meet the criterion on at least 3 arrays (minimum for stats)

    ## filter out only the low expressed probes (not control) in the dataset
    # LE means low expression removed (only)
    fltNrm <- normlst[isexpr, ] # this is a way of extracting samples logically considered TRUE by certain tests
    avgProbes <- avereps(fltNrm, ID = fltNrm$genes$ProbeName) #average the probes. note: use ProbeName instead of SystematicName

    ## to combine gene repeats by probe signel variance (retain maximum variance probes for the same gene)
    if (combineGeneDup){
      tmplst <- avgProbes
      tmplst$genes <- merge(tmplst$genes, annot, by = "ProbeName", all.x = TRUE)

      if (parallelComputing){
        # set up clusters for PSOCK
        n_cores <- detectCores() - 1
        cl <- makeCluster(n_cores, type = clusterType, outfile = "")
        registerDoParallel(cl) # part of doParallel package
        on.exit(stopCluster(cl)) # close connect when exiting the function

        cat("Combining gene duplicates by signal variance...") # initiation message
        combGeneProbe <- foreach(i = unique(tmplst$genes[, geneSymbolVar]), .combine = "c", .packages = "foreach") %dopar% {
          tmp <- tmplst$E[which(tmplst$genes[, geneSymbolVar] %in% i), ]
          if (is.null(nrow(tmp))){
            out <- var(tmp)
          } else {
            out <- foreach(j = seq(nrow(tmp)), .combine = "c") %do% {
              var(as.numeric(tmp[j, ]))
            }
          }
          out <- data.frame(tmplst$genes[which(tmplst$genes[, geneSymbolVar] %in% i), ], var = out)
          out[which(out$var %in% max(out$var)), "ProbeName"]
        }
        cat("DONE!\n") # termination message
      } else {
        cat("Combining gene duplicates by signal variance...") # initiation message
        combGeneProbe <- foreach(i = unique(tmplst$genes[, geneSymbolVar]), .combine = "c", .packages = "foreach") %do% {
          tmp <- tmplst$E[which(tmplst$genes[, geneSymbolVar] %in% i), ]
          if (is.null(nrow(tmp))){
            out <- var(tmp)
          } else {
            out <- foreach(j = seq(nrow(tmp)), .combine = "c") %do% {
              var(as.numeric(tmp[j, ]))
            }
          }
          out <- data.frame(tmplst$genes[which(tmplst$genes[, geneSymbolVar] %in% i), ], var = out)
          out[which(out$var %in% max(out$var)), "ProbeName"]
        }
        cat("DONE!\n") # termination message
      }

      avgProbes <- avgProbes[avgProbes$genes[, "ProbeName"] %in% combGeneProbe, ]
    }
  }

  ## output
  cat("Probes retained upon backgroud filtering:\n")
  print(table(isexpr)) # output the isexpr summary
  return(avgProbes)
}

