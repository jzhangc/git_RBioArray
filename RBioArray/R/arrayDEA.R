#' @title rbioarray_PreProc
#'
#' @description Data pre-processing function for the
#' @param rawlist Input data, either a list, \code{EList} or \code{MAList} object.
#' @param logTrans If to perfom a log transformation on the data or not. Default is \code{FALSE}.
#' @param logTransMethod If \code{logTrans = TRUE}, set which method to use for the transformation, \code{"log2"} or \code{"log10"}. Default is \code{"log2"}.
#' @param logTransObjT If \code{logTrans = TRUE}, set the file name for the output \code{csv} file containing the log transformed data.
#' @param logTransMulticore If \code{logTrans = TRUE}, set if to use parallel computing for the transformation or not. Default is \code{FALSE}.
#' @param bgMethod Background correction method. Default is \code{"auto"}. See \code{backgroundCorrect()} function from \code{limma} package for details.
#' @param normMethod Normalization method. Default is \code{"quantile"}. See \code{normalizeBetweenArrays()} function from \code{limma} package for details.
#' @param ... arguments for \code{backgroundCorrect.matrix()} or \code{backgroundCorrect()} functions from \code{limma} package.
#' @return Depending on the input type, the function outputs a \code{list}, \code{Elist} or \code{MAList} object with corrected and normalized expression values. If \code{logTrans = TRUE}, the function also outputs a \code{csv} file containing the log transformed data.
#' @importFrom limma backgroundCorrect normalizeBetweenArrays backgroundCorrect.matrix arrayWeights
#' @examples
#' \dontrun{
#' normdata <- rbioarray_PreProc(mydata)
#' }
#' @export
rbioarray_PreProc <- function(rawlist, logTrans = FALSE, logTransMethod = "log2", logTransObjT = "data", logTransMulticore = FALSE,
                              bgMethod = "auto", normMethod = "quantile", ...){

  if (class(rawlist) == "list"){

    ## log transform  or not
    if (logTrans){

      if (!logTransMulticore){

        # log transform
        mtx <- apply(rawlist$E, c(1,2), FUN = ifelse(logTransMethod == "log10", log10, log2))

      } else {

        # parallel computing
        # set up cpu cluster
        n_cores <- detectCores() - 1
        cl <- makeCluster(n_cores)
        clusterExport(cl, varlist = "rawlist", envir = environment())
        on.exit(stopCluster(cl)) # close connect when exiting the function

        # log transform
        mtx <- parApply(cl, rawlist$E, c(1,2), FUN = ifelse(logTransMethod == "log10", log10, log2))

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
#' @param plot If to generate volcano plots for the DE results. Defualt is \code{TRUE}. Plots are exported as \code{pdf} files.
#' @param FC Threshold for fold change (FC) for volcano plot. Default is \code{1.5}.
#' @param DE DE methods set for p value thresholding. Values are \code{"fdr"} and \code{"spikein"}. Default is \code{"fdr"}.
#' @param q.value Only used when DE set as \code{"spikein"}, backup threshold for the p value if spikein p values is larger than \code{0.05}.
#' @param Title Figure title. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param xLabel X-axis label. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param yLabel Y-axis label. Make sure to use quotatio marks. Use \code{NULL} to hide. Default is \code{"Mean Decrease in Accurac"}
#' @param symbolSize Size of the symbol. Default is \code{2}.
#' @param xTxtSize Font size for the x-axis text. Default is \code{10}.
#' @param yTxtSize Font size for the y-axis text. Default is \code{10}.
#' @param plotWidth The width of the figure for the final output figure file. Default is \code{170}.
#' @param plotHeight The height of the figure for the final output figure file. Default is \code{150}.
#' @return The function outputs a \code{list} object with DE results, merged with annotation. The function also exports DE reuslts to the working directory in \code{csv} format.
#' @details When \code{"fdr"} set for DE, the p value threshold is set as \code{0.05}.
#' @import ggplot2
#' @importFrom limma lmFit eBayes topTable
#' @importFrom parallel detectCores makeCluster stopCluster parApply parLapply
#' @importFrom grid grid.newpage grid.draw
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @examples
#' \dontrun{
#' rbiorrary_DE(objTitle = "data", fltdata, design, anno = Anno)
#' }
#' @export
rbioarray_DE <- function(objTitle = "data_filtered", fltdata, anno,
                         design, contra, weights = NULL,
                         multicore = FALSE, ...,
                         plot = TRUE,
                         FC = 1.5, DE = "fdr", q.value = 0.05,
                         Title = NULL, xLabel = "log2(fold change)", yLabel = "-log10(p value)",
                         symbolSize = 2, xTxtSize = 10, yTxtSize =10,
                         plotWidth = 170, plotHeight = 150){

  ## extract coefficients
  cf <- colnames(contra) # extract coefficient

  # set an empty matrix for exporting the threolding summery
  threshold_summary <- matrix(nrow = length(cf), ncol = 5)
  colnames(threshold_summary) <- c("coeffcient", "p.value.threshold", "fold.change.threshold", "True", "False")
  threshold_summary <- as.matrix(threshold_summary)

  if(!multicore){

    ## prepare input data object according to the DE method, "fdr" or "spikein".
    if (class(fltdata) == "list"){

      ## DE
      fit <- lmFit(fltdata$E, design, weights = weights)
      fit <- contrasts.fit(fit, contrasts = contra)
      fit <- eBayes(fit)
      fit$genes <- fltdata$genes # add genes matrix to the DE results

      out <- fit[fit$genes$ControlType == 0, ] # remove control probes

      outlist <- lapply(cf, function(i){
        tmp <- topTable(out, coef = i, number = Inf, ...)
        tmp$ProbeName <- rownames(tmp)
        tmp <- merge(tmp, anno, by = "ProbeName")
        return(tmp)
      })

      names(outlist) <- cf

    } else {

      ## DE
      fit <- lmFit(fltdata, design, weights = weights)
      fit <- contrasts.fit(fit, contrasts = contra)
      fit <- eBayes(fit)

      out <- fit[fit$genes$ControlType == 0, ] # remove control probes

      outlist <- lapply(cf, function(i){
        tmp <- topTable(out, coef = i, number = Inf, ...)
        tmp$ProbeName <- rownames(tmp)
        tmp <- merge(tmp, anno, by = "ProbeName")
        return(tmp)
      })

      names(outlist) <- cf

    }

    ## write DE results into files
    lapply(1:length(cf), function(j){
      write.csv(outlist[[j]], file = paste(cf[[j]], "_DE.csv", sep = ""), na = "NA", row.names = FALSE)
    })

    ## volcano plot
    # plot
    if (plot){

      if (DE == "fdr"){
        n_probe <- length(rownames(topTable(out, n = Inf))) #extract the total probe number

        threshold_summary[] <- t(sapply(1: length(cf), function(j){

          # set the cutoff
          cutoff <- as.factor(abs(outlist[[j]]$logFC) >= log2(FC) & outlist[[j]]$P.Value <= q.value / n_probe) # divide by the probe number

          # plot
          loclEnv <- environment()
          plt <- ggplot(outlist[[j]], aes(x = logFC, y = -log10(P.Value), colour = cutoff), environment = loclEnv) +
            geom_point(alpha = 0.4, size = symbolSize) +
            ggtitle(Title) +
            scale_y_continuous(expand = c(0.02, 0)) +
            xlab(xLabel) +
            ylab(yLabel) +
            geom_vline(xintercept = log2(FC), linetype = "dashed") +
            geom_vline(xintercept = - log2(FC), linetype = "dashed") +
            geom_hline(yintercept = - log10(q.value / length(rownames(topTable(out, n = Inf)))), linetype = "dashed") +
            theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
                  plot.title = element_text(hjust = 0.5),
                  legend.position = "none",
                  legend.title = element_blank(),
                  axis.text.x = element_text(size = xTxtSize),
                  axis.text.y = element_text(size = yTxtSize, hjust = 0.5))


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
                 width = plotWidth, height = plotHeight, units = "mm",dpi = 600) # deparse(substitute(x)) converts object name into a character string
          grid.draw(pltgtb) # preview

          # dump the info to the threshold dataframe
          if (length(levels(cutoff)) == 1){

            if (levels(cutoff) == "TRUE"){

              tmp <- c(cf[[j]], signif(q.value / n_probe, digits = 4), FC, summary(cutoff)[["TRUE"]], 0)

            } else {
              tmp <- c(cf[[j]], signif(q.value / n_probe, digits = 4), FC, 0, summary(cutoff)[["FALSE"]])
            }

          } else {
            tmp <- c(cf[[j]], signif(q.value / n_probe, digits = 4), FC, summary(cutoff)[["TRUE"]], summary(cutoff)[["FALSE"]])
          }

        }))



      } else if (DE == "spikein"){


        PCntl <- fit[fit$genes$ControlType == 1, ] # extract PC stats

        threshold_summary[] <- t(sapply(1: length(cf), function(j){

          # set cutoff
          ifelse(min(PCntl$p.value[, cf[j]]) > 0.05, pcutoff <- q.value, pcutoff <- min(PCntl$p.value[, cf[j]]))

          cutoff <- as.factor(abs(outlist[[j]]$logFC) >= log2(FC) & outlist[[j]]$P.Value <= pcutoff)

          # plot
          loclEnv <- environment()
          plt <- ggplot(outlist[[j]], aes(x = logFC, y = - log10(P.Value), colour = cutoff), environment = loclEnv) +
            geom_point(alpha = 0.4, size = symbolSize) +
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
                 width = plotWidth, height = plotHeight, units = "mm", dpi = 600) # deparse(substitute(x)) converts object name into a character string
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

        }))

      } else {stop("Please choose a proper DE method for p value thresholding")}

    }

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

      outlist <- parLapply(cl, cf, fun = function(i){
        tmp <- limma::topTable(out, coef = i, number = Inf, ...)
        tmp$ProbeName <- rownames(tmp)
        tmp <- merge(tmp, anno, by = "ProbeName")
        return(tmp)
      })

      names(outlist) <- cf

    } else {

      fit <- lmFit(fltdata, design, weights = weights)
      fit <- contrasts.fit(fit, contrasts = contra)
      fit <- eBayes(fit)

      out <- fit[fit$genes$ControlType == 0, ] # remove control probes

      outlist <- parLapply(cl, cf, fun = function(i){
        tmp <- limma::topTable(out, coef = i, number = Inf, ...)
        tmp$ProbeName <- rownames(tmp)
        tmp <- merge(tmp, anno, by = "ProbeName")
        return(tmp)
      })

      names(outlist) <- cf

    }

    ## write DE results into files
    parLapply(cl, 1:length(cf), fun = function(j){
      write.csv(outlist[[j]], file = paste(cf[[j]], "_DE.csv", sep = ""),  na = "NA", row.names = FALSE)
    })

    ## volcano plot
    if (plot){

      if (DE == "fdr"){

        n_probe <- length(rownames(topTable(out, n = Inf))) #extract the total probe number

        threshold_summary[] <- t(parSapply(cl, 1: length(cf), FUN = function(j){

          # set the cutoff
          cutoff <- as.factor(abs(outlist[[j]]$logFC) >= log2(FC) & outlist[[j]]$P.Value <= q.value / n_probe) # divide by the probe number

          # plot
          loclEnv <- environment()
          plt <- ggplot2::ggplot(outlist[[j]], ggplot2::aes(x = logFC, y = -log10(P.Value), colour = cutoff), environment = loclEnv) +
            ggplot2::geom_point(alpha = 0.4, size = symbolSize) +
            ggplot2::ggtitle(Title) +
            ggplot2::scale_y_continuous(expand = c(0.02, 0)) +
            ggplot2::xlab(xLabel) +
            ggplot2::ylab(yLabel) +
            ggplot2::geom_vline(xintercept = log2(FC), linetype = "dashed") +
            ggplot2::geom_vline(xintercept = - log2(FC), linetype = "dashed") +
            ggplot2::geom_hline(yintercept = - log10(q.value / length(rownames(limma::topTable(out, n = Inf)))), linetype = "dashed") +
            ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white', colour = 'black'),
                           panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.5),
                           plot.title = ggplot2::element_text(hjust = 0.5),
                           legend.position = "none",
                           legend.title = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_text(size = xTxtSize),
                           axis.text.y = ggplot2::element_text(size = yTxtSize, hjust = 0.5))


          grid::grid.newpage()

          # extract gtable
          pltgtb <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plt))

          # add the right side y axis
          Aa <- which(pltgtb$layout$name == "axis-l")
          pltgtb_a <- pltgtb$grobs[[Aa]]
          axs <- pltgtb_a$children[[2]]
          axs$widths <- rev(axs$widths)
          axs$grobs <- rev(axs$grobs)
          axs$grobs[[1]]$x <- axs$grobs[[1]]$x - ggplot2::unit(1, "npc") + ggplot2::unit(0.08, "cm")
          Ap <- c(subset(pltgtb$layout, name == "panel", select = t:r))
          pltgtb <- gtable::gtable_add_cols(pltgtb, pltgtb$widths[pltgtb$layout[Aa, ]$l], length(pltgtb$widths) - 1)
          pltgtb <- gtable::gtable_add_grob(pltgtb, axs, Ap$t, length(pltgtb$widths) - 1, Ap$b)

          # export the file and draw a preview
          ggplot2::ggsave(filename = paste(cf[[j]],".volcano.pdf", sep = ""), plot = pltgtb,
                          width = plotWidth, height = plotHeight, units = "mm",dpi = 600) # deparse(substitute(x)) converts object name into a character string
          grid::grid.draw(pltgtb) # preview

          # dump the info to the threshold dataframe
          if (length(levels(cutoff)) == 1){

            if (levels(cutoff) == "TRUE"){

              tmp <- c(cf[[j]], signif(q.value / n_probe, digits = 4), FC, summary(cutoff)[["TRUE"]], 0)

            } else {
              tmp <- c(cf[[j]], signif(q.value / n_probe, digits = 4), FC, 0, summary(cutoff)[["FALSE"]])
            }

          } else {
            tmp <- c(cf[[j]], signif(q.value / n_probe, digits = 4), FC, summary(cutoff)[["TRUE"]], summary(cutoff)[["FALSE"]])
          }
        }))

      } else if (DE == "spikein"){


        PCntl <- fit[fit$genes$ControlType == 1, ] # extract PC stats

        threshold_summary[] <- t(parSapply(cl, 1: length(cf), FUN = function(j){

          # set cutoff
          ifelse(min(PCntl$p.value[, cf[j]]) > 0.05, pcutoff <- q.value, pcutoff <- min(PCntl$p.value[, cf[j]]))

          cutoff <- as.factor(abs(outlist[[j]]$logFC) >= log2(FC) & outlist[[j]]$P.Value <= pcutoff)

          # plot
          loclEnv <- environment()
          plt <- ggplot2::ggplot(outlist[[j]], ggplot2::aes(x = logFC, y = -log10(P.Value), colour = cutoff), environment = loclEnv) +
            ggplot2::geom_point(alpha = 0.4, size = symbolSize) +
            ggplot2::ggtitle(Title) +
            ggplot2::scale_y_continuous(expand = c(0.02, 0)) +
            ggplot2::xlab(xLabel) +
            ggplot2::ylab(yLabel) +
            ggplot2::geom_vline(xintercept = log2(FC), linetype = "dashed") +
            ggplot2::geom_vline(xintercept = - log2(FC), linetype = "dashed") +
            ggplot2::geom_hline(yintercept = - log10(pcutoff), linetype = "dashed") +
            ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white', colour = 'black'),
                           panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.5),
                           plot.title = ggplot2::element_text(hjust = 0.5),
                           legend.position = "none",
                           legend.title = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_text(size = xTxtSize),
                           axis.text.y = ggplot2::element_text(size = yTxtSize, hjust = 0.5))


          grid::grid.newpage()

          # extract gtable
          pltgtb <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plt))

          # add the right side y axis
          Aa <- which(pltgtb$layout$name == "axis-l")
          pltgtb_a <- pltgtb$grobs[[Aa]]
          axs <- pltgtb_a$children[[2]]
          axs$widths <- rev(axs$widths)
          axs$grobs <- rev(axs$grobs)
          axs$grobs[[1]]$x <- axs$grobs[[1]]$x - ggplot2::unit(1, "npc") + ggplot2::unit(0.08, "cm")
          Ap <- c(subset(pltgtb$layout, name == "panel", select = t:r))
          pltgtb <- gtable::gtable_add_cols(pltgtb, pltgtb$widths[pltgtb$layout[Aa, ]$l], length(pltgtb$widths) - 1)
          pltgtb <- gtable::gtable_add_grob(pltgtb, axs, Ap$t, length(pltgtb$widths) - 1, Ap$b)

          # export the file and draw a preview
          ggplot2::ggsave(filename = paste(cf[[j]],".volcano.pdf", sep = ""), plot = pltgtb,
                          width = plotWidth, height = plotHeight, units = "mm", dpi = 600) # deparse(substitute(x)) converts object name into a character string
          grid::grid.draw(pltgtb) # preview

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

        }))

      } else {stop("Please choose a proper DE method for p value thresholding")}

    }

  }

  ## output the DE object to the environment
  assign(paste(objTitle, "_DE", sep = ""), outlist, envir = .GlobalEnv)
  write.csv(threshold_summary, file = paste(objTitle, "_thresholding_summary.csv", sep = ""), row.names = FALSE)

}
