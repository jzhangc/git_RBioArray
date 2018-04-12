#' @title rbioarray_venn_DE
#'
#' @description Venn diagrame for DE results. Uses \code{vennDiagram()} from \code{limma} package.
#' @param plotName Name for the output venn diagram. Default is \code{"DE"}.
#' @param plotWidth The width of the figure for the final output figure file. Default is \code{170}.
#' @param plotHeight The height of the figure for the final output figure file. Default is \code{150}.
#' @param ... arguments for \code{vennDiagram()} from \code{limma} package.
#' @param annot Annotation data frame. If set, the output csv file will have a gene symbol column. The function will seek \code{genesymbolVar} value as the variable name for the gene symbol information. Default is \code{NULL}.
#' @param DEdata Input DE object from \code{\link{rbioarray_DE}} or \code{\link{rbioseq_DE}}.
#' @param dataProbeVar \code{DEdata} variable name for probe name. Default is \code{"ProbeName"}.
#' @param geneName If to only use probes with a gene name. Default is \code{FALSE}.
#' @param annotProbeVar \code{annot} variable name for probe name. Default is \code{"ProbeName"}.
#' @param genesymbolVar Only needed when \code{geneName = TRUE}. The name of the variable for gene symbols from the \code{annot} object. Only set this argument when \code{geneName = TRUE}. Default is \code{NULL}.
#' @param sig.method DE methods set for p value thresholding. Values are \code{"fdr"}, \code{"spikein"} or \code{"none"}. Default is \code{"fdr"}.
#' @param fltlist Only needed when \code{sig.method = "spikein"}. Filtered data, either a list, \code{EList} or \code{MAList} object. Default is \code{NULL}.
#' @param design Only needed when \code{sig.method = "spikein"}. Design matrix. Default is \code{NULL}.
#' @param contra Only needed when \code{sig.method = "spikein"}. Contrast matrix. Default is \code{NULL}.
#' @param weights Only needed when \code{sig.method = "spikein"}. Array weights, determined by \code{arrayWeights()} function from \code{limma} package. Default is \code{NULL}.
#' @param sig.p P value threshold. Only needed for \code{sig.method = "fdr"} and \code{sig.method = "spikein"} when calculated p value is larger than sig.p. Default is \code{0.05}.
#' @param FC Fold change threshold. Default is \code{1.5}.
#' @param parallelComputing If to use parallel computing. Default is \code{FALSE}.
#' @return The function outputs a \code{pdf} file for venn diagrams (total, up- and down-regulations). The function also exports overlapping gene or probe into a \code{csv} file.
#' @details When \code{"fdr"} set for sig.method, the p value threshold is set as \code{0.05}. When there is no significant genes or probes identified under \code{sig.method = "fdr"}, the threshold is set to \code{1}. If the arugments for \code{sig.method = "spikein"} are not complete, the function will automatically use \code{"fdr"}.
#' @import doParallel
#' @import foreach
#' @importFrom limma lmFit eBayes topTable contrasts.fit vennDiagram
#' @importFrom parallel detectCores makeCluster stopCluster
#' @examples
#' \dontrun{
#' rbioarray_venn_DE(plotName = "DE", cex = c(1, 2, 2), mar = rep(0.5,4), names = c("control", "stress1", "stress2"),
#'                   DEdata = fltdata_DE, geneName = TRUE, genesymbolVar = "GeneSymbol",
#'                   sig.method = "spikein", fltlist = fltdata, annot = annot, design = design, contra = contra, weights = fltdata$ArrayWeight,
#'                   parallelComputing = FALSE)
#' }
#' @export
rbioarray_venn_DE <- function(objTitle = "DE", plotName = "DE", plotWidth = 5, plotHeight = 5, ...,
                              annot = NULL,
                              DEdata = NULL, dataProbeVar = "ProbeName",
                              geneName = FALSE, annotProbeVar = "ProbeName", genesymbolVar = NULL,
                              sig.method = "fdr", fltlist = NULL, design = NULL, contra = NULL, weights = NULL, sig.p = 0.05, FC = 1.5,
                              parallelComputing = FALSE){
  ## check the key arguments
  if (is.null(DEdata)){
    stop(cat("Please set input DE data object. Hint: it is the output list from rbioarray_DE(). Function terminated.\n"))
  }

  if (is.null(design)){
    stop(cat("Please set design matrix. Function terminated.\n"))
  }
  if (is.null(contra)){
    stop(cat("Please set contrast object. Function terminated.\n"))
  }

  ## set up the DE dataframe
  vennDE <- vector(mode = "list", length = length(names(DEdata)))
  names(vennDE) <- names(DEdata)

  if (geneName){
    if (!is.null(genesymbolVar)){
      for (i in 1:length(names(DEdata))){
        vennDE[[i]] <- DEdata[[i]][complete.cases(DEdata[[i]][, genesymbolVar]), ]
      }
    } else {
      warning("No variable name for gene symbol set. Proceed with probe names with no probes removed.")
      vennDE <- DEdata
    }
  } else {
    vennDE <- DEdata
  }

  ## tempfunc
  # m - vennDE; n - individual coef
  cal_pcutoff <- function(m, n){
    # set up tmpdfm
    tmpdfm <- m[[n]]
    # set the cutoff
    if (tolower(sig.method) == "fdr"){
      if (length(which(tmpdfm$adj.P.Val < sig.p)) == 0){
        warning("No significant results found using FDR correction. Please consider using another thresholding method. For now, sig.p is applied on raw p.values.")
        p_threshold <- sig.p
      } else {
        p_threshold <- max(tmpdfm[tmpdfm$adj.P.Val < sig.p, ]$P.Value)
      }

    } else if (tolower(sig.method) == "spikein") {
      # check arugments
      if (is.null(fltlist) | is.null(design) | is.null(contra) | is.null(weights)){
        warning(cat("Arguments not complete for \"spikein\" method. Proceed with \"fdr\" instead."))

        if (length(which(tmpdfm$adj.P.Val < sig.p)) == 0){
          p_threshold <- 1
        } else {
          p_threshold <- max(tmpdfm[tmpdfm$adj.P.Val < sig.p, ]$P.Value)
        }

      } else {
        # DE for extracting positive control
        if (class(fltlist) == "list"){
          fit <- lmFit(fltlist$E, design, weights = weights)
          fit <- contrasts.fit(fit, contrasts = contra)
          fit <- eBayes(fit)
          fit$genes <- fltlist$genes # add genes matrix to the DE results
        } else {
          fit <- lmFit(fltlist, design, weights = weights)
          fit <- contrasts.fit(fit, contrasts = contra)
          fit <- eBayes(fit)
        }
        cf <- names(m)
        PC <- fit[fit$genes$ControlType == 1, ]
        ifelse(min(PC$p.value[, cf[n]]) > sig.p, p_threshold <- sig.p, p_threshold <- min(PC$p.value[, cf[n]]))
      }

    } else if (tolower(sig.method) == "none"){
      p_threshold <- sig.p
    } else {stop(cat("Please set p value thresholding method, \"fdr\", \"spikein\", or \"none\"."))}
    return(p_threshold)
  }

  ## prepare matrices
  # log fold change
  lfc <- matrix(nrow = length(vennDE[[1]][, dataProbeVar]), ncol = length(vennDE))
  rownames(lfc) <- vennDE[[1]][, dataProbeVar]
  colnames(lfc) <- names(vennDE)

  # p value
  p <- array(NA, dim(lfc), dimnames = dimnames(lfc)) # this is another way to create a matrix. Give the value of NA.

  # matrix for plotting
  mtx <- array(NA, dim(lfc), dimnames = dimnames(lfc))

  ## populate the matrices
  if (!parallelComputing){
    for (j in 1:length(names(vennDE))){
      lfc[, j] <- vennDE[[j]]$logFC # extract log fold change
      p[, j] <- vennDE[[j]]$P.Value # extract p value (p) to a matrix
      pcutoff <- cal_pcutoff(m = vennDE, n = j)
      # note we are using factors here e.g. "-1L, 0L, 1L". the starting value is 0L
      mtx[, j] <- ifelse(p[, j] >= pcutoff | abs(lfc[, j]) < log2(FC), 0L, ifelse(lfc[, j] > 0, 1L, -1L))
    }
  } else { # parallel computing
    ## set up cpu cluster
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores, type = "PSOCK")
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    ## populate matrices
    lfc[] <- foreach(j = 1:length(names(vennDE)), .combine = "cbind") %dopar% {
      out <- vennDE[[j]]$logFC
    }
    p[] <- foreach(j = 1:length(names(vennDE)), .combine = "cbind") %dopar% {
      out <- vennDE[[j]]$P.Value
    }
    mtx[] <- foreach(j = 1:length(names(vennDE)), .combine = "cbind", .packages = "limma") %dopar% {
      pcutoff <- cal_pcutoff(m = vennDE, n = j)
      out <- ifelse(p[, j] >= pcutoff | abs(lfc[, j]) < log2(FC), 0L, ifelse(lfc[, j] > 0, 1L, -1L))
    }
  }

  ## venn diagram plotting
  pdf(file = paste(plotName, "_venn_total.pdf", sep = ""), width = plotWidth, height = plotHeight)
  vennDiagram(mtx, circle.col = 1:length(names(DEdata)), ...)
  dev.off()

  pdf(file = paste(plotName, "_venn_up.pdf", sep = ""), width = plotWidth, height = plotHeight)
  vennDiagram(mtx, circle.col = 1:length(names(DEdata)), include = "up", ...)
  dev.off()

  pdf(file = paste(plotName, "_venn_down.pdf", sep = ""), width = plotWidth, height = plotHeight)
  vennDiagram(mtx, circle.col = 1:length(names(DEdata)), include = "down", ...)
  dev.off()

  ## output a csv file with annotation
  if (is.null(annot)){
    outdfm <- data.frame(ProbeName = rownames(mtx), mtx)
    write.csv(outdfm, file = paste(objTitle, "_venn_table.csv", sep = ""), row.names = FALSE)
  } else {
    outdfm <- data.frame(ProbeName = rownames(mtx), mtx)
    outdfm <- merge(annot[, c(annotProbeVar, genesymbolVar)], outdfm, by = dataProbeVar, all.y = TRUE)
    write.csv(outdfm, file = paste(objTitle, "_venn_table.csv", sep = ""), row.names = FALSE)
  }
}
