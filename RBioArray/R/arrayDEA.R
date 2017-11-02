#' @title rbioarray_PreProc
#'
#' @description Data pre-processing function for the
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

      tmpdata <- list(E = mtx, genes = rawlist$genes, target = rawlist$target)

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

    output <- list(E = Norm, genes = rawlist$genes, target = rawlist$target, ArrayWeight = Wgt)

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
#' @param percentile The percentile cutoff, only used when muliptle negative control probes are detected. Default is \code{0.95}.
#' @param ... arguments for \code{backgroundCorrect.matrix()} or \code{backgroundCorrect()} functions from \code{limma} package.
#' @return Depending on the input type, the function outputs a \code{list}, \code{Elist} or \code{MAList} object with filtered expression values.
#' @importFrom limma avereps
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
#' @param fltdata filtered data, either a list, \code{EList} or \code{MAList} object. Default is \code{NULL}.
#' @param anno Annotation object, usually a \code{dataframe}. Make sure to name the probe ID variable \code{ProbeName}. Default is \code{NULL}.
#' @param design Design matrix. Default is \code{NULL}.
#' @param contra Contrast matrix. Default is \code{NULL}.
#' @param weights Array weights, determined by \code{arrayWeights()} function from \code{limma} package. Default is \code{NULL}.
#' @param ... arguments for \code{topTable()} from \code{limma} package.
#' @param plot If to generate volcano plots for the DE results. Defualt is \code{TRUE}. Plots are exported as \code{pdf} files.
#' @param geneName If to only plot probes with a gene name. Default is \code{FALSE}.
#' @param genesymbolVar The name of the variable for gene symbols from the \code{anno} object. Only set this argument when \code{geneName = TRUE}. Default is \code{NULL}.
#' @param topgeneLabel If to display the gene identification, i.e., probem name or gene name, on the plot. Default is \code{FALSE}.
#' @param nGeneSymbol When \code{topgeneLabel = TRUE}, to set how many genes to display. Default is \code{5}.
#' @param padding When \code{topgeneLabel = TRUE}, to set the distance between the dot and the gene symbol. Default is \code{0.5}.
#' @param FC Threshold for fold change (FC) for volcano plot. Default is \code{1.5}.
#' @param DE DE methods set for p value thresholding. Values are \code{"fdr"} and \code{"spikein"}. Default is \code{"fdr"}.
#' @param q.value Only used when DE set as \code{"spikein"}, backup threshold for the p value if spikein p values is larger than \code{0.05}.
#' @param Title Figure title. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param xLabel X-axis label. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param yLabel Y-axis label. Make sure to use quotatio marks. Use \code{NULL} to hide. Default is \code{"Mean Decrease in Accurac"}
#' @param symbolSize Size of the symbol. Default is \code{2}.
#' @param sigColour Colour of the significant genes or probes. Default is \code{"red"}.
#' @param nonsigColour Colour of the non-significant genes or probes. Default is \code{"gray"}.
#' @param xTxtSize Font size for the x-axis text. Default is \code{10}.
#' @param yTxtSize Font size for the y-axis text. Default is \code{10}.
#' @param plotWidth The width of the figure for the final output figure file. Default is \code{170}.
#' @param plotHeight The height of the figure for the final output figure file. Default is \code{150}.
#' @param parallelComputing If to use parallel computing. Default is \code{FALSE}.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @return The function outputs a \code{list} object with DE results, a \code{data frame} object for the F test results, merged with annotation. The function also exports DE reuslts to the working directory in \code{csv} format.
#' @details When \code{"fdr"} set for DE, the p value threshold is set as \code{0.05}. When there is no significant genes or probes identified under \code{DE = "fdr"}, the threshold is set to \code{1}. Also note that both \code{geneName} and \code{genesymbolVar} need to be set to display gene sysmbols on the plot. Otherwise, the labels will be probe names. Additionally, when set to display gene symbols, all the probes without a gene symbol will be removed.
#' @import ggplot2
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @importFrom limma lmFit eBayes topTable contrasts.fit
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom grid grid.newpage grid.draw
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom ggrepel geom_text_repel
#' @examples
#' \dontrun{
#' rbioarray_DE(objTitle = "fltdata2", fltdata, anno = Anno, design, contra = contra,
#'              weights = fltdata$ArrayWeight, parallelComputing = TRUE,
#'              plot = TRUE, geneName = TRUE, genesymbolVar = "GeneSymbol",
#'              DE = "spikein")
#' }
#' @export
rbioarray_DE <- function(objTitle = "data_filtered", fltdata = NULL, anno = NULL,
                         design = NULL, contra = NULL, weights = NULL,
                         ...,
                         plot = TRUE, geneName = FALSE, genesymbolVar = NULL, topgeneLabel = FALSE, nGeneSymbol = 5, padding = 0.5,
                         FC = 1.5, DE = "fdr", q.value = 0.05,
                         Title = NULL, xLabel = "log2(fold change)", yLabel = "-log10(p value)",
                         symbolSize = 2, sigColour = "red", nonsigColour = "gray",
                         xTxtSize = 10, yTxtSize =10,
                         plotWidth = 170, plotHeight = 150,
                         parallelComputing = FALSE, clusterType = "PSOCK"){

  ## check the key arguments
  if (is.null(fltdata)){
    stop("Please set input data object. Hint: either a list, EList or MAList object with pre-processed and flitered expression data.")
  }

  if (is.null(anno)){
    stop(cat("Please set the annotatioin dataframe. Make sure to name the variable containing probe name \"ProbeName\"."))
  }

  if (is.null(design)){
    stop("Please set design matrix.")
  }

  if (is.null(contra)){
    stop("Please set contrast object.")
  }

  ## extract coefficients
  cf <- colnames(contra) # extract coefficient

  # set an empty matrix for exporting the threolding summery
  threshold_summary <- matrix(nrow = length(cf), ncol = 5)
  colnames(threshold_summary) <- c("coeffcient", "p.value.threshold", "fold.change.threshold", "True", "False")
  threshold_summary <- as.matrix(threshold_summary)

  ## temp func for plotting
  # i: outlist (object) listed below
  tmpfunc <- function(i, j, PC = NULL){

    # set the data frame
    if (geneName){
      if (!is.null(genesymbolVar)){
        tmpdfm <- i[[j]][complete.cases(i[[j]][, genesymbolVar]), ]
      } else {
        warning("No variable name for gene symbol set. Proceed with probe names with no probes removed.")
        tmpdfm <- i[[j]]
      }
    } else {
      tmpdfm <- i[[j]]
    }

    # set the cutoff
    if (DE == "fdr"){

      if (length(which(tmpdfm$adj.P.Val < q.value)) == 0){
        pcutoff <- 1
      } else {
        pcutoff <- max(tmpdfm[tmpdfm$adj.P.Val < q.value, ]$P.Value)
      }

      cutoff <- as.factor(abs(tmpdfm$logFC) >= log2(FC) & tmpdfm$P.Value < pcutoff)

    } else if (DE == "spikein") {

      ifelse(min(PC$p.value[, cf[j]]) > q.value, pcutoff <- q.value, pcutoff <- min(PC$p.value[, cf[j]]))
      cutoff <- as.factor(abs(tmpdfm$logFC) >= log2(FC) & tmpdfm$P.Value < pcutoff)

    } else {stop(cat("Please set p value thresholding method, \"fdr\" or \"spikein\"."))}


    # plot
    loclEnv <- environment()
    plt <- ggplot(tmpdfm, aes(x = logFC, y = -log10(P.Value)), environment = loclEnv) +
      geom_point(alpha = 0.4, size = symbolSize, aes(colour = cutoff)) +
      scale_color_manual(values = c(nonsigColour, sigColour)) +
      ggtitle(Title) +
      scale_y_continuous(expand = c(0.02, 0)) +
      xlab(xLabel) +
      ylab(yLabel) +
      geom_vline(xintercept = log2(FC), linetype = "dashed") +
      geom_vline(xintercept = - log2(FC), linetype = "dashed") +
      geom_hline(yintercept = - log10(pcutoff), linetype = "dashed") +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none",
            legend.title = element_blank(),
            axis.text.x = element_text(size = xTxtSize),
            axis.text.y = element_text(size = yTxtSize, hjust = 0.5))

    if (topgeneLabel){
      tmpfltdfm <- tmpdfm[abs(tmpdfm$logFC) >= log2(FC) & tmpdfm$P.Value < pcutoff, ]
      tmpfltdfm <- tmpfltdfm[order(tmpfltdfm$P.Value), ]
      plt <- plt + geom_text_repel(data = head(tmpfltdfm, n = nGeneSymbol),
                                   aes(x = logFC, y = -log10(P.Value), label = head(tmpfltdfm, n = nGeneSymbol)[, genesymbolVar]),
                                   point.padding = unit(padding, "lines"))
    }

    grid.newpage()

    # extract gtable
    pltgtb <- ggplot_gtable(ggplot_build(plt))

    # add the right side y axis
    Aa <- which(pltgtb$layout$name == "axis-l")
    pltgtb_a <- pltgtb$grobs[[Aa]]
    axs <- pltgtb_a$children[[2]]
    axs$widths <- rev(axs$widths)
    axs$grobs <- rev(axs$grobs)
    axs$grobs[[1]]$x <- axs$grobs[[1]]$x - unit(1, "npc") + unit(0.08, "cm")
    Ap <- c(subset(pltgtb$layout, name == "panel", select = t:r))
    pltgtb <- gtable_add_cols(pltgtb, pltgtb$widths[pltgtb$layout[Aa, ]$l], length(pltgtb$widths) - 1)
    pltgtb <- gtable_add_grob(pltgtb, axs, Ap$t, length(pltgtb$widths) - 1, Ap$b)

    # export the file and draw a preview
    ggsave(filename = paste(cf[[j]],".volcano.pdf", sep = ""), plot = pltgtb,
           width = plotWidth, height = plotHeight, units = "mm",dpi = 600)
    grid.draw(pltgtb) # preview

    # dump the info to the threshold dataframe
    if (length(levels(cutoff)) == 1){

      if (levels(cutoff) == "TRUE"){

        tmp <- c(cf[[j]], signif(pcutoff, digits = 4), FC, summary(cutoff)[["TRUE"]], 0)

      } else {
        tmp <- c(cf[[j]], signif(pcutoff, digits = 4), FC, 0, summary(cutoff)[["FALSE"]])
      }

    } else {
      tmp <- c(cf[[j]], signif(pcutoff, digits = 4), FC, summary(cutoff)[["TRUE"]], summary(cutoff)[["FALSE"]])
    }

  }


  ## DE
  cat("Linear fitting...") # message
  if (class(fltdata) == "list"){

    fit <- lmFit(fltdata$E, design, weights = weights)
    fit <- contrasts.fit(fit, contrasts = contra)
    fit <- eBayes(fit)
    fit$genes <- fltdata$genes # add genes matrix to the DE results

  } else {

    fit <- lmFit(fltdata, design, weights = weights)
    fit <- contrasts.fit(fit, contrasts = contra)
    fit <- eBayes(fit)

  }
  cat("DONE!\n") # message

  out <- fit[fit$genes$ControlType == 0, ] # remove control probes

  if (DE == "spikein"){ # extract PC stats for spikein method
    PCntl <- fit[fit$genes$ControlType == 1, ]
  } else {
    PCntl <- NULL
  }

  ## output and plotting
  if(!parallelComputing){

    # compile resutls into a list
    outlist <- lapply(cf, function(i){
      tmp <- topTable(out, coef = i, number = Inf, ...)
      tmp$ProbeName <- rownames(tmp)
      tmp <- merge(tmp, anno, by = "ProbeName")
      return(tmp)
    })

    names(outlist) <- cf

    # write DE results into files
    lapply(1:length(cf), function(j){
      write.csv(outlist[[j]], file = paste(cf[[j]], "_DE.csv", sep = ""), na = "NA", row.names = FALSE)
    })

    # volcano plot and output summary
    if (plot){

      threshold_summary[] <- t(sapply(1: length(cf), function(x)tmpfunc(i = outlist, j = x, PC = PCntl)))

    }

  } else { ## parallel computing

    # check the cluster type
    if (clusterType != "PSOCK" & clusterType != "FORK"){
      stop("Please set the cluter type. Options are \"PSOCK\" (default) and \"FORK\".")
    }

    # set up cpu cores
    n_cores <- detectCores() - 1

    if (clusterType == "PSOCK"){ # all OS types
      # set up cpu cluster
      cl <- makeCluster(n_cores, type = "PSOCK")
      registerDoParallel(cl)
      on.exit(stopCluster(cl)) # close connect when exiting the function

      outlist <- foreach(i = 1: length(cf), .packages = "limma") %dopar% {

        tmp <- limma::topTable(out, coef = cf[i], number = Inf, ...)
        tmp$ProbeName <- rownames(tmp)
        tmp <- merge(tmp, anno, by = "ProbeName")
        return(tmp)

      }

      names(outlist) <- cf

      # write DE results into files
      foreach(j = 1: length(cf)) %dopar% {
        write.csv(outlist[[j]], file = paste(cf[[j]], "_DE.csv", sep = ""),  na = "NA", row.names = FALSE)
      }

      # volcano plot and output summary
      if (plot){

        threshold_summary[] <- foreach(j = 1: length(cf), .combine = "rbind", .packages = c("limma", "ggplot2", "gtable", "grid")) %dopar% {
          tmpfunc(i = outlist, j = j, PC = PCntl)

        }
      }

    } else { # macOS and Unix-like systems

      outlist <- mclapply(cf, FUN = function(i){
        tmp <- topTable(out, coef = i, number = Inf, ...)
        tmp$ProbeName <- rownames(tmp)
        tmp <- merge(tmp, anno, by = "ProbeName")
        return(tmp)
      }, mc.cores = n_cores, mc.preschedule = FALSE)

      names(outlist) <- cf

      # write DE results into files
      mclapply(1:length(cf), FUN = function(j){
        write.csv(outlist[[j]], file = paste(cf[[j]], "_DE.csv", sep = ""), na = "NA", row.names = FALSE)
      }, mc.cores = n_cores, mc.preschedule = FALSE)

      # volcano plot and output summary
      if (plot){

        threshold_summary[] <- t(sapply(1: length(cf), function(x)tmpfunc(i = outlist, j = x, PC = PCntl)))

      }

    }


  }

  ## output the DE/fit objects to the environment, as well as the DE csv files into wd
  fitout <- topTable(out, number = Inf)
  fitout <- merge(fitout, anno, by = "ProbeName")
  assign(paste(objTitle, "_fit", sep = ""), fitout, envir = .GlobalEnv)
  write.csv(fitout, file = paste(objTitle, "_DE_Fstats.csv", sep = ""), row.names = FALSE)

  assign(paste(objTitle, "_DE", sep = ""), outlist, envir = .GlobalEnv)


  write.csv(threshold_summary, file = paste(objTitle, "_thresholding_summary.csv", sep = ""), row.names = FALSE)
  ## message
  if (geneName & !is.null(genesymbolVar)){
    print("Probes without a gene symbol are removed from the volcano plots")
  }

}
