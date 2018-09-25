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
rbioarray_PreProc <- function(rawlist, logTrans = FALSE, logTransMethod = "log2", logTransObjT = "data", logTransParallelComputing = FALSE,
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
rbioarray_flt <- function(normlst, ctrlProbe = TRUE, ctrlTypeVar = "ControlType", percentile = ifelse(ctrlProbe, 0.95, 0.05),
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

