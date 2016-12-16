#' @title rbioarray_PreProc
#'
#' @description Data pre-processing function for the
#' @param RawData Input data, either a list, \code{EList} or \code{MAList} object.
#' @param BgMtd Background correction method. Default is \code{"auto"}. See \code{backgroundCorrect()} function from \code{limma} package for details.
#' @param BgOffst The offset to add to the background corrected values, as a base line. Default is \code{"50"}.
#' @param NormMtd Normalization method. Default is \code{"quantile"}. See \code{normalizeBetweenArrays()} function from \code{limma} package for details.
#' @param ... arguments for \code{backgroundCorrect.matrix()} or \code{backgroundCorrect()} functions from \code{limma} package.
#' @return Depending on the input type, the function outputs a \code{list}, \code{Elist} or \code{MAList} object with corrected and normalized expression values.
#' @importFrom limma backgroundCorrect normalizeBetweenArrays backgroundCorrect.matrix
#' @examples
#' \dontrun{
#' normdata <- rbioarray_PreProc(mydata)
#' }
#' @export
rbioarray_PreProc <- function(RawData, bgMethod = "auto", normMethod = "quantile", ...){
  if (class(RawData) == "list"){
    BgC <- backgroundCorrect.matrix(RawData$E, method = bgMethod, ...) #background correction
    Norm <- normalizeBetweenArrays(BgC, normMethod) # quantile normalization
    Wgt <- arrayWeights(Norm) # array weight

    output <- list(E = Norm, genes = RawData$gene, target = RawData$target, ArrayWeight = Wgt)
  } else {
    BgC <- backgroundCorrect(RawData, method = bgMethod, ...) #background correction
    Norm <- normalizeBetweenArrays(BgC, method = normMethod) # quantile normalization

    output <- Norm
  }

  return(output)
}



#' @title rbioarray_flt
#'
#' @description data filter function based on spike-in negative control.
#' @param normlst Normalized data, either a list, \code{EList} or \code{MAList} object.
#' @param percentile The percentile cutoff, only used when muliptle negative control probes are detected. Default is \code{0.95}.
#' @param ... arguments for \code{backgroundCorrect.matrix()} or \code{backgroundCorrect()} functions from \code{limma} package.
#' @return Depending on the input type, the function outputs a \code{list}, \code{Elist} or \code{MAList} object with filtered expression values.
#' @examples
#' \dontrun{
#' fltdata <- rbioarray_flt(normdata)
#' }
#' @export
rbioarray_flt <- function(normlst, percentile = 0.95){


  if (class(normlst$E[normlst$genes$ControlType == -1, ]) == "numeric"){
    ### LE keeps expression values more than 10% brighter than NC (dark corner) on at least 3 arrays
    ## extract the 95% quanitle of the negative control signals
    neg <- normlst$E[normlst$genes$ControlType == -1, ] # no 95% percentile as only one entry
  } else {
    ### LE keeps expression values more than 10% brighter than NC (dark corner) on at least 3 arrays
    ## extract the 95% quanitle of the negative control signals
    neg <- apply(normlst$E[normlst$genes$ControlType == -1, ], 2, function(x)quantile(x, p = percentile)) # neg95
  }

  if (class(normlst) == "list"){

    ## low expression cuttoff set at at least 10% hihger than the neg95
    LE_cutoff <- matrix(1.1 * neg, nrow(normlst$E), ncol(normlst$E), byrow = TRUE)


    ## summary after applying LE_cutoff (T/F for each sample)
    # this only compares the element Nrm$E
    isexpr <- rowSums(normlst$E > LE_cutoff) >= 3 # We keep probes that meet the criterion on at least 3 arrays (minimum for stats)

    ## filter out only the low expressed probes (not control) in the dataset
    # LE means low expression removed (only)
    flt_E <- normlst$E[isexpr, ] # this is a way of extracting samples logically considered TRUE by certain tests
    flt_gene <- normlst$gene[isexpr, ]
    fltlst <- list(E = flt_E, genes = flt_gene, target = normlst$target)

    flt_E_avg <- avereps(fltlst$E, ID = fltlst$genes$ProbeName)

    avgProbesLE <- list(E = flt_E_avg, genes = unique(fltlst$genes[fltlst$genes$ProbeName %in% rownames(flt_E_avg), ]),
                        target = normlst$target, ArrayWeight = normlst$ArrayWeight)

  } else {

    ## low expression cuttoff set at at least 10% hihger than the neg95
    LE_cutoff <- matrix(1.1 * neg, nrow(normlst), ncol(normlst), byrow = TRUE)


    ## summary after applying LE_cutoff (T/F for each sample)
    # this only compares the element Nrm$E
    isexpr <- rowSums(normlst$E > LE_cutoff) >= 3 # We keep probes that meet the criterion on at least 3 arrays (minimum for stats)


    ## filter out only the low expressed probes (not control) in the dataset
    # LE means low expression removed (only)
    fltNrmLE <- normlst[isexpr, ] # this is a way of extracting samples logically considered TRUE by certain tests

    avgProbesLE <- avereps(fltNrmLE, ID = fltNrmLE$genes$ProbeName) #average the probes. note: use ProbeName instead of SystematicName

  }


  ## output
  print(table(isexpr)) # output the isexpr summary

  return(avgProbesLE)
}


#' @title rbioarray_DE
#'
#' @description DE analysis function.
#' @param objTitle Name for the output list. Default is \code{"data_filtered"}.
#' @param fltdata filtered data, either a list, \code{EList} or \code{MAList} object.
#' @param anno Annotation object, usually a \code{dataframe}.
#' @param design Design matrix.
#' @param contra Contrast matrix.
#' @param weights Array weights, determined by \code{arrayWeights()} function from \code{limma} package. Default is \code{NULL}.
#' @param multicore If to use parallel computing. Default is \code{FALSE}.
#' @param ... arguments for \code{topTable()} from \code{limma} package.
#' @param DE DE methods set for p value thresholding.
#' @return The function outputs a \code{list} object with DE results, merged with annotation. The function also exports DE reuslts to the working directory in \code{csv} format.
#' @importFrom limma lmFit eBayes topTable
#' @importFrom parallel detectCores makeCluster stopCluster parApply parLapply
#' @examples
#' \dontrun{
#' rbiorrary_DE(objTitle = "data", fltdata, design, anno = Anno)
#' }
#' @export
rbioarray_DE <- function(objTitle = "data_filtered", fltdata, anno,
                         design, contra, weights = NULL,
                         multicore = FALSE, ...,
                         pvalue = 0.05, DE = "fdr"){

  ## extract coefficients
  coef <- colnames(contra) # extract coefficient

  if(!multicore){

    ## prepare input data object according to the DE method, "fdr" or "spikein".
    if (class(fltdata) == "list"){

      ## DE
      fit <- lmFit(fltdata$E, design, weights = weights)
      fit <- contrasts.fit(fit, contrasts = contra)
      fit <- eBayes(fit)
      fit$genes <- fltdata$genes # add genes matrix to the DE results

      out <- fit[fit$genes$ControlType == 0, ] # remove control probes

      outlist <- lapply(coef, function(i){
        tmp <- topTable(out, coef = i, number = Inf, ...)
        tmp$ProbeName <- rownames(tmp)
        tmp <- merge(tmp, anno, by = "ProbeName")
        return(tmp)
      })

      names(outlist) <- coef

    } else {

      ## DE
      fit <- lmFit(fltdata, design, weights = weights)
      fit <- contrasts.fit(fit, contrasts = contra)
      fit <- eBayes(fit)

      out <- fit[fit$genes$ControlType == 0, ] # remove control probes

      outlist <- lapply(coef, function(i){
        tmp <- topTable(out, coef = i, number = Inf, ...)
        tmp$ProbeName <- rownames(tmp)
        tmp <- merge(tmp, anno, by = "ProbeName")
        return(tmp)
      })

      names(outlist) <- coef

    }

    ## write DE results into files
    lapply(1:length(coef), function(j){
      write.csv(outlist[[j]], file = paste(coef[[j]], "_DE.csv", sep = ""), na = "NA", row.names = FALSE)
    })

  } else {

    ## parallel computing
    # set up cpu cluster
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores)
    clusterExport(cl, varlist = c("fltdata", "design", "anno"), envir = environment())
    on.exit(stopCluster(cl)) # close connect when exiting the function

    ## DE

    if (class(fltdata) == "list"){

      fit <- lmFit(fltdata$E, design, weights = weights)
      fit <- contrasts.fit(fit, contrasts = contra)
      fit <- eBayes(fit)
      fit$genes <- fltdata$genes # add genes matrix to the DE results

      out <- fit[fit$genes$ControlType == 0, ] # remove control probes

      outlist <- parLapply(cl, coef, fun = function(i){
        tmp <- limma::topTable(out, coef = i, number = Inf, ...)
        tmp$ProbeName <- rownames(tmp)
        tmp <- merge(tmp, anno, by = "ProbeName")
        return(tmp)
      })

      names(outlist) <- coef

    } else {

      fit <- lmFit(fltdata, design, weights = weights)
      fit <- contrasts.fit(fit, contrasts = contra)
      fit <- eBayes(fit)

      out <- fit[fit$genes$ControlType == 0, ] # remove control probes

      outlist <- parLapply(cl, coef, fun = function(i){
        tmp <- limma::topTable(out, coef = i, number = Inf, ...)
        tmp$ProbeName <- rownames(tmp)
        tmp <- merge(tmp, anno, by = "ProbeName")
        return(tmp)
      })

      names(outlist) <- coef

    }

    ## write DE results into files
    parLapply(cl, 1:length(coef), fun = function(j){
      write.csv(outlist[[j]], file = paste(coef[[j]], "_DE.csv", sep = ""),  na = "NA", row.names = FALSE)
    })

  }


  # output the DE object to the environment
  assign(paste(objTitle, "_DE", sep = ""), outlist, envir = .GlobalEnv)

}
