#' @title rbioarray_hcluster
#'
#' @description Wrapper for hierarchical clustering analysis and heatmap visualization.
#' @param plotName File name for the export \code{pdf} plot file. Default is \code{"data"}.
#' @param fltlist Input filtered data, either a list, \code{EList} or \code{MAList} object.
#' @param dataProbeVar \code{fltlist} variable name for probe name. Default is \code{"ProbeName"}.
#' @param genesymbolOnly Whether or not to remove probes without gene symbol. Default is \code{FALSE}.
#' @param annot Annotation object, usually a \code{dataframe}. Make sure to name the probe ID variable \code{ProbeName}. Only set this argument when \code{genesymbolVar = TRUE}. Default is \code{NULL}.
#' @param annotProbeVar \code{annot} variable name for probe name. Default is \code{"ProbeName"}.
#' @param genesymbolVar The name of the variable for gene symbols from the \code{annot} object. Only set this argument when \code{genesymbolVar = TRUE}. Default is \code{NULL}.
#' @param rmControl Set only when \code{ctrlProbe = TRUE} and \code{ctrlTypeVar} is properly set,  whether to remove control probes (Agilent platform) or not. Default is \code{TRUE}.
#' @param ctrlProbe Wether or not the data set has control type variable, with values \code{-1 (negative control)}, \code{0 (gene probes)} and \code{1 (positive control)}. Default is \code{TRUE}.
#' @param ctrlTypeVar Set only when \code{ctrlProbe = TRUE}, the control type variable. Default is the \code{Agilent} variable name \code{"ControlType"}.
#' @param n Number of genes to be clustered, numeric input or \code{"all"}. Default is \code{"all"}.
#' @param fct Input \code{factor} object for sample groupping labels.
#' @param sampleName A \code{vector} containing names to display for each heatmap column. Default is \code{NULL} and the function will use the column name from the input.
#' @param colGroup Colour group, numeric or dependent on \code{fct}.
#' @param distance Distance calculation method. Default is \code{"euclidean"}. See \code{\link{dist}} for more.
#' @param clust Clustering method. Default is \code{"complete"}. See \code{\link{hclust}} for more.
#' @param colColour Column group colour. Default is \code{"Paired"}. See \code{RColorBrewer} package for more.
#' @param mapColour Heat map colour. Default is \code{"PRGn"}. See \code{RColorBrewer} package for more.
#' @param n_mapColour Number of colours displayed. Default is \code{11}. See \code{RColorBrewer} package for more.
#' @param ... Additional arguments for \code{heatmap.2} function from \code{gplots} package.
#' @param plotWidth Width of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param plotHeight Height of the plot. Unit is \code{inch}. Default is \code{7}.
#' @return A heatmap based on hierarchical clustering analysis in \code{pdf} format.
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#'
#' # standard usage
#' rbioarray_hcluster(fltlist = normdata, n = 500, fct = conSum, trace = "none", srtCol = 45,
#' offsetCol = 0, adjCol = c(1, 0), labRow = FALSE, key.title = "", keysize = 1.5,
#' key.xlab = "Normalized expression value", key.ylab = "Probe count")
#'
#' # for non microarray or RNAseq data sets
#' ###### unsupervised heatmap ######
#' ## load the file
#' raw <- read.csv(file = "all_data.csv", na.strings = " ", stringsAsFactors = FALSE, check.names = FALSE)
#'
#' ## build the index
#' idx <- raw[, 1:2] # extract the sample information
#'
#' conSum <- factor(idx$Condition, levels = unique(idx$Condition)) # extract the factor
#'
#' # create the input data for RBioArray pacakge
#' rawT <- t(raw[, -(1:2)])
#' colnames(rawT) <- idx$SampleID
#' rawT <- apply(rawT, c(1,2), FUN = log2) # log2 tranforamtion
#' inputlist <- list(E = as.matrix(rawT),
#'                   genes = data.frame(GeneNames = rownames(rawT), ControlType = rep(0, length(rownames(rawT))),
#'                   stringsAsFactors = FALSE), targets = idx)
#'
#' rbioarray_hcluster(plotName = "all", fltlist = inputlist, dataProbeVar = "GeneNames", n = "all", rmControl = FALSE, fct = conSum,
#'                    trace = "none",
#'                    distance = "euclidean", clust = "complete", colColour = "Paired",
#'                    mapColour = "RdBu", n_mapColour = 11,
#'                    srtCol = 30, offsetCol = 0, labRow = inputlist$genes$GeneNames,
#'                    key.title = "", keysize = 1.5, cexCol = 0.9, cexRow = 0.6,
#'                    key.xlab = "Log expression value", key.ylab = "miRNA count",
#'                    plotWidth = 10, plotHeight = 10)
#'
#' }
#' @export
rbioarray_hcluster <- function(plotName = "data", fltlist = NULL, dataProbeVar = "ProbeName",
                               genesymbolOnly = FALSE, annot = NULL, annotProbeVar = "ProbeName", genesymbolVar = NULL,
                               n = "all",
                               rmControl = TRUE, ctrlProbe = TRUE, ctrlTypeVar = "ControlType",
                               sampleName = NULL,
                               fct = NULL, colGroup = ifelse(length(levels(fct)) < 19, length(levels(fct)), 19),
                               distance = "euclidean", clust = "complete",
                               colColour = "Paired", mapColour = "PRGn", n_mapColour = 11, ...,
                               plotWidth = 7, plotHeight = 7){
  ## chekc arguments
  if (is.null(fltlist)){
    stop(cat("Please provide filtered input data. Function terminated.\n"))
  }
  if (is.null(fct)){
    stop(cat("Please provide smaple index with argument fct. Function terminated.\nFunction terminated.\n"))
  }
  if (rmControl){
    if (!ctrlProbe){
      stop(cat("rmControl can only be set TRUE when ctrlProbe = TRUE. Function terminated.\n"))
    } else {
      if (!ctrlTypeVar %in% names(fltlist$genes)){
        stop(cat("ctrlTypeVar not found. Function terminated.\n"))
      }
    }
  }

  ## set up dis and cluster functions
  distfunc <- function(x)dist(x, method = distance)
  clustfunc <- function(x)hclust(x, method = clust)

  ## prepare dfm for clustering
  dfm <- data.frame(fltlist$genes, fltlist$E, check.names = FALSE)

  if (rmControl){ # remove control
    dfm <- dfm[dfm[, ctrlTypeVar] == 0, ]
  }

  if (n != "all"){ # subset
    dfm <- dfm[1:n, ]
  }

  ## set ColSideColors
  col_cluster <- clustfunc(distfunc(t(dfm[, -c(1:(ncol(dfm) - ncol(fltlist$E)))])))

  colG <- cutree(col_cluster, colGroup) # column group
  colC <- brewer.pal(ifelse(colGroup < 3, 3, colGroup), colColour) # column colour

  ## prepare mtx for plotting
  dfm2 <- dfm
  if (genesymbolOnly){ # remove probes without gene symbol or not
    if (is.null(annot) | is.null(genesymbolVar)){
      warning("No annotation object or gene sybmol variable detected. Cluster will proceed with all probes.")
      dfm2 <- dfm2
      mtx <- as.matrix(dfm2[, -c(1:(ncol(dfm2) - ncol(fltlist$E)))]) # remove annotation info. same as below.
      rownames(mtx) <- dfm2[, dataProbeVar]
    } else {
      geneSymbl <- annot[annot[, annotProbeVar] %in% dfm2[, dataProbeVar], ][, genesymbolVar]
      dfm2$geneSymbol <- geneSymbl
      dfm2 <- dfm2[complete.cases(dfm2), ] # remove probes withnout a gene symbol
      mtx <- as.matrix(dfm2[, -c(1:(ncol(dfm2) - ncol(fltlist$E) -1), ncol(dfm2))])
      rownames(mtx) <- dfm2[, "geneSymbol"] # row names are now gene symbols. note that the variable name is geneSymbol, NOT the argument value.
    }
  } else {
    dfm2 <- dfm2
    mtx <- as.matrix(dfm2[, -c(1:(ncol(dfm2) - ncol(fltlist$E)))])
    rownames(mtx) <- dfm2[, dataProbeVar]
  }

  if (!is.null(sampleName)){
    colnames(mtx) <- sampleName
  }

  ## heatmap
  # draw heatmap
  pdf(file = paste(plotName, "_heatmap.pdf", sep = ""), width = plotWidth, height = plotHeight)
  heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
            col = brewer.pal(n_mapColour, mapColour), ColSideColors = colC[colG], ...)
  dev.off()
}

#' @title rbioseq_hcluster
#'
#' @description Wrapper for hierarchical clustering analysis and heatmap visualization for RNA seq data.
#' @param plotName File name for the export \code{pdf} plot file. Default is \code{"data"}.
#' @param dfm_count Dataframe contains the feature read counts, with rows as genomic featues (or genes) and column as samples. Default is \code{NULL}.
#' @param dfm_annot Dataframe contains the gene annotation information, with rows as genmic features and columns as annotation variables. The row lengths of this dataframe should be the same as \code{dfm_count}.
#' @param count_threshold Read count threshold. No filtering will be applied when set \code{"none"}. Otherwise, a numeric number can be set as the minimum read count for filtering. DDefault is \code{"none"}.
#' @param qc_plot Wether or not to produce a QC plot upon filtering, normalization and weight calculation. Default is \code{FALSE}.
#' @param n Number of genes to be clustered, numeric input or \code{"all"}. Default is \code{"all"}.
#' @param fct Input \code{factor} object for samples. Default is \code{NULL}.
#' @param sampleName A \code{vector} containing names for column. Default is \code{NULL} and the function will use the column name from the input.
#' @param colGroup Colour group, numeric or dependent on \code{fct}.
#' @param distance Distance calculation method. Default is \code{"euclidean"}. See \code{\link{dist}} for more.
#' @param clust Clustering method. Default is \code{"complete"}. See \code{\link{hclust}} for more.
#' @param colColour Column group colour. Default is \code{"Paired"}. See \code{RColorBrewer} package for more.
#' @param mapColour Heat map colour. Default is \code{"PRGn"}. See \code{RColorBrewer} package for more.
#' @param n_mapColour Number of colours displayed. Default is \code{11}. See \code{RColorBrewer} package for more.
#' @param ... Additional arguments for \code{heatmap.2} function from \code{gplots} package.
#' @param plotWidth Width of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param plotHeight Height of the plot. Unit is \code{inch}. Default is \code{7}.
#' @return A heatmap based on hierarchical clustering analysis in \code{pdf} format
#' @details The data filtering and normalization functions are also included in the function.
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#'
#' # standard usage
#'
#' }
#' @export
rbioseq_hcluster <- function(plotName = "data", dfm_count = NULL, dfm_annot = NULL, geneidVar = "gene_id",
                             count_threshold = "none", design = NULL, qc_plot = FALSE,
                             n = "all", sampleName = NULL,
                             fct = NULL, colGroup = ifelse(length(levels(fct)) < 19, length(levels(fct)), 19),
                             distance = "euclidean", clust = "complete",
                             colColour = "Paired", mapColour = "PRGn", n_mapColour = 11, ...,
                             plotWidth = 7, plotHeight = 7){

  ## chekc variables
  if (is.null(dfm_count) | is.null(dfm_annot) | class(dfm_count) != "data.frame" | class(dfm_annot) != "data.frame"){
    stop(cat("Please provide the read count and annotation dataframes. Please also make sure the type as data.frame. Function terminated.\n"))
  }

  if (is.null(design)){
    stop(cat("Please provide design matrix. Function terminated.\n"))
  }

  if (is.null(fct)){
    stop(cat("Please provide smaple index with argument fct. Function terminated.\n"))
  }

  ## set up dis and cluster functions
  distfunc <- function(x)dist(x, method = distance)
  clustfunc <- function(x)hclust(x, method = clust)

  ## clustering analysis and colour setup
  # normalization and filtering
  cat("Data filtering and normalization...") # message
  dge <- DGEList(counts = dfm_count, genes = dfm_annot)

  if (count_threshold != "none"){ # set the count threshold for filtering
    count_s <- rowSums(dge$counts) # thresholdd
    isexpr <- count_s > count_threshold

    dge <- dge[isexpr, , keep.lib.size = FALSE] # filtering
  }

  # for data Voom normalization
  dgenormf <- calcNormFactors(dge)
  vmwt <- voomWithQualityWeights(dgenormf, design = design, plot = qc_plot, normalization = "quantile") # Voom normalization with quality weights
  cat("DONE!\n") # message


  # cluster and colour
  col_cluster <- clustfunc(distfunc(t(vmwt$E)))
  colG <- cutree(col_cluster, colGroup) # column group
  colC <- brewer.pal(ifelse(colGroup < 3, 3, colGroup), colColour) # column colour

  ## prepare mtx for plotting
  ## prepare dfm for clustering
  dfm_E <- vmwt$E
  dfm_A <- vmwt$genes

  if (n != "all"){ # subset
    dfm_E <- dfm_E[1:n, ]
    dfm_A <- dfm_A[1:n, ]
  }

  mtx <- as.matrix(dfm_E)
  rownames(mtx) <- dfm_A[, geneidVar]

  if (!is.null(sampleName)){
    colnames(mtx) <- sampleName
  } else {
    colnames(mtx) <- colnames(vmwt$E)
  }

  ## heatmap
  # draw heatmap
  pdf(file = paste(plotName, "_heatmap.pdf", sep = ""), width = plotWidth, height = plotHeight)
  heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
            col = brewer.pal(n_mapColour, mapColour), ColSideColors = colC[colG], ...)
  dev.off()
}


#' @title rbioarray_hcluster_super
#'
#' @description Wrapper for supervised hierarchical clustering analysis and heatmap visualization for both microarray and RNAseq.
#' @param plotName File name for the export \code{pdf} plot file. Default is \code{"data"}.
#' @param fltDOI Based on filtered data, a subset corresponding to the comparasion, either a list, \code{EList} or \code{MAList} object.
#' @param dfmDE A subset of the DE list, i.e. a \code{topTable} dataframe, corresponding to the comparasion (i.e., contrast).
#' @param ctrlProbe Wether or not the data set has control type variable in \code{fltDOI}, with values \code{-1 (negative control)}, \code{0 (gene probes)} and \code{1 (positive control)}. Default is \code{TRUE}.
#' @param ctrlTypeVar Set only when \code{ctrlProbe = TRUE}, the control type variable. Default is the \code{Agilent} variable name \code{"ControlType"}.
#' @param dataProbeVar \code{dfmDE} variable name for probe name. Default is \code{"ProbeName"}.
#' @param DE.sig.method Thresholding method, "fdr", "spikein" or "none. Default is \code{"none"}.
#' @param DE.sig.p P value cut off. Default is \code{0.05}.
#' @param FC Fold change (FC) filter for the heatmap. Default is \code{NULL}.
#' @param fct Input \code{factor} object for samples.
#' @param sampleName A \code{vector} containing names for column. Default is \code{NULL} and the function will use the column name from the input.
#' @param colGroup Colour group, numeric or dependent on \code{fct}.
#' @param distance Distance calculation method. Default is \code{"euclidean"}. See \code{\link{dist}} for more.
#' @param clust Clustering method. Default is \code{"complete"}. See \code{\link{hclust}} for more.
#' @param rowLabel Whether to display row label or not. Default is \code{FALSE}.
#' @param annot Annotation object, usually a \code{dataframe}. Make sure to name the probe ID variable \code{ProbeName}. Only set this argument when \code{rowLabel = TRUE}. Default is \code{NULL}.
#' @param annotProbeVar \code{annot} variable name for probe name. Default is \code{"ProbeName"}.
#' @param genesymbolVar The name of the variable for gene symbols from the \code{annot} object. Only set this argument when \code{rowLabel = TRUE}. Default is \code{NULL}.
#' @param colColour Column group colour. Default is \code{"Paired"}. See \code{RColorBrewer} package for more.
#' @param mapColour Heat map colour. Default is \code{"PRGn"}. See \code{RColorBrewer} package for more.
#' @param n_mapColour Number of colours displayed. Default is \code{11}. See \code{RColorBrewer} package for more.
#' @param ... Additional arguments for \code{heatmap.2} function from \code{gplots} package.
#' @param plotWidth Width of the plot. Unit is \code{inch}. Default is \code{7}.
#' @param plotHeight Height of the plot. Unit is \code{inch}. Default is \code{7}.
#' @details Note that both \code{annot} and \code{genesymbolVar} need to be set to display gene sysmbols as row labels. Otherwise, the row labels will be probe names. Also note that when set to display gene symbols, all the probes without a gene symbol will be removed.
#' @return A supervised heatmap based on hierarchical clustering analysis in \code{pdf} format.
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' rbioarray_hcluster_super(plotName = "pre_experiVnaive", fltDOI = pre_experiVnaive_super, dfmDE = fltdata_DE$pre_experiVnaive, dataProbeVar = "ProbeName",
#'                          FC = 1.5, DE.sig.p = 0.0003461,
#'                          clust = "ward.D2",
#'                          fct = factor(c("naivepre", "naivepre", "naivepre", "naivepre", "exppre", "exppre", "exppre", "exppre", "exppre"),
#'                          levels = c("naivepre", "exppre")), trace = "none", srtCol = 30, offsetCol = 0.5, adjCol = c(1, 0),
#'                          rowLabel = TRUE, annot = annot, annotProbeVar = "ProbeName", genesymbolVar = "GeneSymbol",
#'                          offsetRow = 0, adjRow = c(0, 0.5), cexRow = 0.6,
#'                          key.title = "", keysize = 1.5,
#'                          key.xlab = "Normalized expression value", key.ylab = "Gene count")
#' }
#' @export
rbioarray_hcluster_super <- function(plotName = "data", fltDOI, dfmDE,
                                     ctrlProbe = TRUE, ctrlTypeVar = "ControlType",
                                     dataProbeVar = "ProbeName",
                                     DE.sig.method = "none",
                                     DE.sig.p = 0.05, FC = 1.5,
                                     fct, sampleName = NULL,
                                     colGroup = ifelse(length(levels(fct)) < 19, length(levels(fct)), 19),
                                     distance = "euclidean", clust = "complete",
                                     rowLabel = FALSE, annot = NULL, annotProbeVar = "ProbeName", genesymbolVar = NULL,
                                     colColour = "Paired", mapColour = "PRGn", n_mapColour = 11, ...,
                                     plotWidth = 7, plotHeight = 7){ #DOI: fltered subset data of interest
  ## argument check
  if (!dataProbeVar %in% names(dfmDE) | !dataProbeVar %in% names(fltDOI$genes)){
    stop(cat("data probe variable not found in dfmDE or fltDOI"))
  }

  ## prepare matrix for plotting
  dfm <- data.frame(fltDOI$genes, fltDOI$E)

  if (ctrlProbe){
    if (ctrlTypeVar %in% names(dfm)){
      dfm <- dfm[dfm[, ctrlTypeVar] == 0, ] # remove control probes
    } else {
      print(cat("Control probe setting on, but invalid control type variable detected, proceed without removing control probes."))
    }
  }


  if (tolower(DE.sig.method) == "fdr"){
    if (length(which(dfmDE$adj.P.Val < DE.sig.p)) == 0){
      warning("No significant results found using FDR correction. Please consider using another thresholding method. For now, alpha is applied on raw p.values.")
      pcutoff <- DE.sig.p
      pb_name <- dfmDE[dfmDE$P.Value < pcutoff, dataProbeVar]
      dfm <- dfm[dfm[, dataProbeVar] %in% pb_name, ]
    } else {
      pcutoff <- max(dfmDE[dfmDE$adj.P.Val < DE.sig.p, ]$P.Value)
      pb_name <- dfmDE[dfmDE$P.Value < pcutoff, dataProbeVar]
      dfm <- dfm[dfm[, dataProbeVar] %in% pb_name, ]
    }
  } else {
    pcutoff <- DE.sig.p
    pb_name <- dfmDE[dfmDE$P.Value < pcutoff, dataProbeVar]
    dfm <- dfm[dfm[, dataProbeVar] %in% pb_name, ]
  }

  ## set FC filter
  if (!is.null(FC)){
    pb_name_fc <- dfmDE[abs(dfmDE$logFC) >= log2(FC), dataProbeVar]
    dfm <- dfm[dfm[, dataProbeVar] %in% pb_name_fc, ]
  }

  ## check the dim
  if (dim(dfm)[1] < 2){
    stop("Only one row left after DE filtering. Nothing to cluster.")
  }

  ## heatmap
  # set up dis and cluster functions
  distfunc <- function(x)dist(x, method = distance)
  clustfunc <- function(x)hclust(x, method = clust)

  # set ColSideColors
  col_cluster <- clustfunc(distfunc(t(dfm[, -c(1:(ncol(dfm) - ncol(fltDOI$E)))])))
  colG <- cutree(col_cluster, colGroup) # column group
  colC <- brewer.pal(ifelse(colGroup < 3, 3, colGroup), colColour) # column colour

  # draw heatmap
  ogNcol <- dim(fltDOI$E)[2] # original numbers of col
  annoNcol <- dim(dfm)[2] # numbers of col with annotation
  s <- (annoNcol - ogNcol + 1):annoNcol # extract only the data by removing the annotation columns

  if (rowLabel){
    if (is.null(annot) | is.null(genesymbolVar)){
      warning("No annotation object or gene sybmol variable detected. Row labels will be the default probe names.")

      # retrive correct data matrix without annoation
      mtx <- as.matrix(dfm[, s])
      rownames(mtx) <- dfm[, dataProbeVar]

      if (!is.null(sampleName)){
        if (length(sampleName) != ncol(mtx)){
          stop("sampleName variable isn't the same length as the samples.")
        } else {
          colnames(mtx) <- sampleName
        }
      }

      pdf(file = paste(plotName, "_heatmap.supervised.pdf", sep = ""), width = plotWidth, height = plotHeight)
      heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
                col = brewer.pal(n_mapColour, mapColour), ColSideColors = colC[colG], ...)
      dev.off()

    } else {
      geneSymbl <- annot[annot[, annotProbeVar] %in% dfm[, dataProbeVar], ][, genesymbolVar]

      # process dfm further to only retain probes with gene symbol
      dfm$geneSymbol <- geneSymbl
      dfm <- dfm[complete.cases(dfm), ] # remove probes withnout a gene symbol
      labrow <- dfm$geneSymbol

      # retrive correct data matrix without annoation
      mtx <- as.matrix(dfm[, s])
      rownames(mtx) <- dfm[, dataProbeVar]

      if (!is.null(sampleName)){
        if (length(sampleName) != ncol(mtx)){
          stop("sampleName variable isn't the same length as the samples.")
        } else {
          colnames(mtx) <- sampleName
        }
      }

      pdf(file = paste(plotName, "_heatmap.supervised.pdf", sep = ""), width = plotWidth, height = plotHeight)
      heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
                col = brewer.pal(n_mapColour, mapColour), ColSideColors = colC[colG], labRow = labrow, ...)
      dev.off()
      cat("Probes with no gene names are removed.")
    }

  } else {
    mtx <- as.matrix(dfm[, s])
    rownames(mtx) <- dfm[, dataProbeVar]

    if (!is.null(sampleName)){
      if (length(sampleName) != ncol(mtx)){
        stop("sampleName variable isn't the same length as the samples.")
      } else {
        colnames(mtx) <- sampleName
      }
    }

    pdf(file = paste(plotName, "_heatmap.supervised.pdf", sep = ""), width = plotWidth, height = plotHeight)
    heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
              col = brewer.pal(n_mapColour, mapColour), ColSideColors = colC[colG], labRow = FALSE,...)
    dev.off()
  }
}
