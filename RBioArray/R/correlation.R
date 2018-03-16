#' @title cor_pvalue
#'
#' @description Calculate p value for the correlation analysis
#' @param r Correlation coefficient
#' @param n Sample size
#' @return P value for the correlation analysis
#' @examples
#' \dontrun{
#'
#' p <- cor_pvalue(r = 0.806687547, n = 6)
#'
#' }
#' @export
cor_pvalue <- function(r, n){
  t <- r/sqrt((1 - r^2)/(n - 2))
  p <- 2*pt(-abs(t), df = n - 2)
  return(p)
}

#' @title rbioarray_corcluster_super
#'
#' @description Wrapper for supervised (or unsupervised) Pearson correlation clustering analysis and heatmap visualization for both microarray and RNAseq, for gene co-expression analysis.
#' @param plotName File name for the export \code{pdf} plot file. Default is \code{"data"}.
#' @param fltlist Based on filtered data, a subset corresponding to the comparasion, either a list, \code{EList} or \code{MAList} object.
#' @param rmControl If to remove control probes (Agilent platform). Default is \code{TRUE}.
#' @param n_subgroup A vector of sample index (row number) for phenotype group. Default is \code{NULL}. The setting can be obtained from the corresponding condition summary object.
#' @param dfmDE A subset of the DE list, i.e. a \code{topTable} dataframe, corresponding to the comparasion (i.e., contrast).
#' @param FDR If to use FDR corrcted p value. Default is \code{TRUE}.
#' @param q.value P value cut off. Default is \code{0.05}. For unsupervised clustering, set \code{q.value = 1}.
#' @param FC Fold change (FC) filter for the heatmap. Default is \code{NULL}.
#' @param dataProbeVar \code{dfmDE} variable name for probe name. Default is \code{NULL}.
#' @param method The correlation method, options are "pearson", "spearman" and "pearson". Default is \code{"pearson"}.
#' @param sigPlot If to generate a significance heatmap. Default is \code{FALSE},
#' @param cor.sig Only set when \code{sigPlot = TRUE}, the alpha value for correlation p value. Default is \code{0.05}
#' @param cor.sigLabelColour Only set when \code{sigPlot = TRUE}, the colour for label for the significant pairs. Default is \code{"red"}.
#' @param cor.sigLabelSize Only set when \code{sigPlot = TRUE}, the size for label for the significant pairs. Default is \code{3}.
#' @param cor.labelColour Only set when \code{sigPlot = TRUE}, the significance heatmap axis label colour. Default is \code{"black"}.
#' @param cor.labelSize Only set when \code{sigPlot = TRUE}, the significance heatmap axis label size. Default is \code{1}.
#' @param cor.labelAngle Only set when \code{sigPlot = TRUE}, the significance heatmap axis label angle. Default is \code{90}.
#' @param cor.keySize Only set when \code{sigPlot = TRUE}, the significance heatmap colour key size. Default is \code{1}.
#' @param axisLabel Whether to display label for both x- and y- axes. Default is \code{FALSE}.
#' @param annot The optional annotation matrix. Only needs to be set to display inforamtions for
#' @param genesymbolVar The name of the variable for gene symbols from the \code{annot} object. Only set this argument when \code{rowLabel = TRUE}. Default is \code{NULL}.
#' @param mapColour Heat map colour. Default is \code{"PRGn"}. See \code{RColorBrewer} package for more.
#' @param n_mapColour Number of colours displayed. Default is \code{11}. See \code{RColorBrewer} package for more.
#' @param ... Additional arguments for \code{heatmap.2} function from \code{gplots} package.
#' @param plotWidth Width of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param plotHeight Height of the plot. Unit is \code{inch}. Default is \code{7}.
#' @details Note that both \code{annot} and \code{genesymbolVar} need to be set to display gene sysmbols as row labels. Otherwise, the row labels will be probe names. Also note that when set to display gene symbols, all the probes without a gene symbol will be removed.
#' @return A supervised heatmap based on hierarchical clustering analysis in \code{pdf} format, along with correaltion coefficient and p value matrices. If set, the function also outputs a significant value heatmap.
#' @import corrplot
#' @import foreach
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#'
#' # n_subgroup = c(1:4) means the correlation uses samples from 1 to 4 (control in this case).
#' # The settings can be obtained from the corresponding condition summary object.
#' rbioarray_corcluster_super(fltlist = all_nrm, n_subgroup = c(1:4),
#'                            dataProbeVar = "gene_id", FDR = TRUE, q.value = 0.02,
#'                            dfmDE = all_DE$`conSumPost - conSumPre`,
#'                            axisLabel = TRUE, genesymbolVar = "gene_name",
#'                            key.title = "", cexRow = 0.3, cexCol = 0.3, offsetRow = 0.001,
#'                            offsetCol = 0.001, margins = c(4, 4))
#'
#' }
#' @export
rbioarray_corcluster_super <- function(plotName = "data",
                                       fltlist = NULL, rmControl = TRUE,
                                       n_subgroup = NULL,
                                       dfmDE = NULL, FDR = TURE, q.value = 0.05, FC = NULL,
                                       dataProbeVar = NULL,
                                       method = "pearson",
                                       sigPlot = FALSE, cor.sig = 0.05, cor.sigLabelColour = "red", cor.sigLabelSize = 3,
                                       cor.labelColour = "black", cor.labelSize = 1, cor.labelAngle = 90, cor.keySize = 1,
                                       axisLabel = FALSE, annot = NULL, genesymbolVar = NULL,
                                       mapColour = "PRGn", n_mapColour = 11, ...,
                                       plotWidth = 7, plotHeight = 7){
  #### test variables
  if (is.null(fltlist)){
    stop(cat("Please set processed data object via fltlist. Function terminated.\n"))
  }

  if (is.null(dfmDE)){
    stop(cat("Please set DE object via dfmDE. Function terminated.\n"))
  }

  if (is.null(dataProbeVar)){
    stop(cat("Please set unique genomic feature ID via dataProbeVar. Function terminated.\n"))
  }


  if (is.null(n_subgroup)){
    stop(cat("Please set the index for phenotype group via n_subgroup. Function terminated.\n"))
  }

  if (rmControl){
    if (!"ControlType" %in% names(fltlist$genes)){
      stop(cat("make sure to have/name ControlType variable in the fltlist"))
    }
  }

  #### fiter and normalization
  vmwt <- fltlist
  dfm <- data.frame(vmwt$genes, vmwt$E)

  if (rmControl){ # remove control
    dfm <- dfm[dfm$ControlType == 0, ]
  }

  if (!is.null(annot)){
    dfm <- merge(annot[, c(dataProbeVar, genesymbolVar)], dfm, by = dataProbeVar, all.y = TRUE)
  }

  #### dfm subsetting using DE resutls (dfmDE)
  ## p value filter
  if (FDR){
    if (length(which(dfmDE$adj.P.Val < q.value)) == 0){
      warning("No significant results found using FDR correction. Please consider using another thresholding method. For now, q.value is applied on raw p.values.")
      pcutoff <- q.value
      pb_name <- dfmDE[dfmDE$P.Value < pcutoff, dataProbeVar]
      dfm <- dfm[dfm[, dataProbeVar] %in% pb_name, ]
    } else {
      pcutoff <- max(dfmDE[dfmDE$adj.P.Val < q.value, ]$P.Value)
      pb_name <- dfmDE[dfmDE$P.Value < pcutoff, dataProbeVar]
      dfm <- dfm[dfm[, dataProbeVar] %in% pb_name, ]
    }
  } else {
    pcutoff <- q.value
    pb_name <- dfmDE[dfmDE$P.Value < pcutoff, dataProbeVar]
    dfm <- dfm[dfm[, dataProbeVar] %in% pb_name, ]
  }

  ## set FC filter, if applicable
  if (!is.null(FC)){
    pb_name_fc <- dfmDE[abs(dfmDE$logFC) >= log2(FC), dataProbeVar]
    dfm <- dfm[dfm[, dataProbeVar] %in% pb_name_fc, ]
  }

  #### heatmap
  ogNcol <- dim(vmwt$E)[2] # original numbers of col
  annoNcol <- dim(dfm)[2] # numbers of col with annotation
  s <- (annoNcol - ogNcol + 1):annoNcol # extract only the data by removing the annotation columns

  if (axisLabel){
    if (!is.null(genesymbolVar) & !is.null(annot)){
      dfm <- dfm[complete.cases(dfm[, genesymbolVar]), ]
      axisrow <- dfm[, genesymbolVar]
      print("Probes with no gene names are removed.")

      mtx <- as.matrix(dfm[, s])
      rownames(mtx) <- dfm[, dataProbeVar]
      cormtx <- t(mtx)
      corcoef <- cor(cormtx[n_subgroup, ], method = method)
      rownames(corcoef) <- axisrow
      colnames(corcoef) <- axisrow
      corp <- foreach(i = corcoef, .combine = "cbind") %do% cor_pvalue(i, n = nrow(cormtx[n_subgroup, ])) # p value matrix
      diag(corp) <- NA
      rownames(corp) <- axisrow
      colnames(corp) <- axisrow

      pdf(file = paste(plotName, "_corheatmap.pdf", sep = ""), width = plotWidth, height = plotHeight)
      heatmap.2(corcoef, symm = TRUE, trace = "none",
                col = brewer.pal(n_mapColour, mapColour), labRow = axisrow, labCol = axisrow, ...)
      dev.off()

      if (sigPlot){
        tryCatch(
          {
            pdf(file = paste(plotName, "_corheatmap.sigplot.pdf", sep = ""), width = plotWidth, height = plotHeight)
            corrplot(corr = corcoef, method = "color", type = "upper", p.mat = corp, sig.level = cor.sig,
                     insig = c("label_sig"), pch.col = cor.sigLabelColour, pch.cex = cor.sigLabelSize,
                     col = brewer.pal(n_mapColour, mapColour),
                     tl.col = cor.labelColour, tl.cex = cor.labelSize, tl.srt = cor.labelAngle, cl.length = 3, cl.cex = cor.keySize)
            dev.off()
          },
          error = function(err){
            print("No significance found. Therefore no significance plot generated.")
            dev.off()
          }
        )
      }
    } else {
      print("No gene symbol variable or annotation dataframe detected. Proceed without one.")

      mtx <- as.matrix(dfm[, s])
      rownames(mtx) <- dfm[, dataProbeVar]
      cormtx <- t(mtx)
      corcoef <- cor(cormtx[n_subgroup, ], method = method)
      corp <- foreach(i = corcoef, .combine = "cbind") %do% cor_pvalue(i, n = nrow(cormtx[n_subgroup, ])) # p value matrix
      diag(corp) <- NA

      pdf(file = paste(plotName, "_corheatmap.pdf", sep = ""), width = plotWidth, height = plotHeight)
      heatmap.2(corcoef, symm = TRUE, trace = "none",
                col = brewer.pal(n_mapColour, mapColour), labRow = FALSE, labCol = FALSE,...)
      dev.off()

      if (sigPlot){
        tryCatch(
          {
            pdf(file = paste(plotName, "_corheatmap.sigplot.pdf", sep = ""), width = plotWidth, height = plotHeight)
            corrplot(corr = corcoef, method = "color", type = "upper", p.mat = corp, sig.level = cor.sig,
                     insig = c("label_sig"), pch.col = cor.sigLabelColour, pch.cex = cor.sigLabelSize,
                     col = brewer.pal(n_mapColour, mapColour),
                     tl.col = cor.labelColour, tl.cex = cor.labelSize, tl.srt = cor.labelAngle, cl.length = 3, cl.cex = cor.keySize)
            dev.off()
          },
          error = function(err){
            print("No significant correlation found. Therefore no significance plot generated.")
            dev.off()
          }
        )
      }
    }
  } else {
    mtx <- as.matrix(dfm[, s])
    rownames(mtx) <- dfm[, dataProbeVar]
    cormtx <- t(mtx)

    corcoef <- cor(cormtx[n_subgroup, ], method = method)
    corp <- foreach(i = corcoef, .combine = "cbind") %do% cor_pvalue(i, n = nrow(cormtx[n_subgroup, ])) # p value matrix
    diag(corp) <- NA

    pdf(file = paste(plotName, "_corheatmap.pdf", sep = ""), width = plotWidth, height = plotHeight)
    heatmap.2(corcoef, symm = TRUE, trace = "none",
              col = brewer.pal(n_mapColour, mapColour), labRow = FALSE, labCol = FALSE,...)
    dev.off()

    if (sigPlot){
      tryCatch(
        {
          pdf(file = paste(plotName, "_corheatmap.sigplot.pdf", sep = ""), width = plotWidth, height = plotHeight)
          corrplot(corr = corcoef, method = "color", type = "upper", p.mat = corp, sig.level = cor.sig,
                   insig = c("label_sig"), pch.col = cor.sigLabelColour, pch.cex = cor.sigLabelSize,
                   col = brewer.pal(n_mapColour, mapColour),
                   tl.col = cor.labelColour, tl.cex = cor.labelSize, tl.srt = cor.labelAngle, cl.length = 3, cl.cex = cor.keySize)
          dev.off()
        },
        error = function(err){
          print("No significance found. Therefore no significance plot generated.")
          dev.off()
        }
      )
    }
  }

  # export correlation matrix
  write.csv(corcoef, file = paste(plotName, ".cor.csv", sep = ""))
  write.csv(corp, file = paste(plotName, ".cor.pvalue.csv", sep = ""))
}
