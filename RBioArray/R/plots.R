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
#' @param fct Input \code{factor} object for samples.
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
#'                   stringsAsFactors = FALSE), target = idx)
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
#' @param dataProbeVar \code{dfmDE} variable name for probe name. Default is \code{"ProbeName"}.
#' @param pcutoff P value cut off. Default is \code{NULL}.
#' @param FC Fold change (FC) filter for the heatmap. Default is \code{NULL}.
#' @param method Thresholding method, "fdr" or "spikein". Default is \code{"spikein"}.
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
#'                          FC = 1.5, pcutoff = 0.0003461,
#'                          clust = "ward.D2",
#'                          fct = factor(c("naivepre", "naivepre", "naivepre", "naivepre", "exppre", "exppre", "exppre", "exppre", "exppre"),
#'                          levels = c("naivepre", "exppre")), trace = "none", srtCol = 30, offsetCol = 0.5, adjCol = c(1, 0),
#'                          rowLabel = TRUE, annot = annot, annotProbeVar = "ProbeName", genesymbolVar = "GeneSymbol",
#'                          offsetRow = 0, adjRow = c(0, 0.5), cexRow = 0.6,
#'                          key.title = "", keysize = 1.5,
#'                          key.xlab = "Normalized expression value", key.ylab = "Gene count")
#' }
#' @export
rbioarray_hcluster_super <- function(plotName = "data", fltDOI, dfmDE, dataProbeVar = "ProbeName",
                                     pcutoff = NULL, FC = NULL,
                                     method = "spikein",
                                     fct, sampleName = NULL,
                                     colGroup = ifelse(length(levels(fct)) < 19, length(levels(fct)), 19),
                                     distance = "euclidean", clust = "complete",
                                     rowLabel = FALSE, annot = NULL, annotProbeVar = "ProbeName", genesymbolVar = NULL,
                                     colColour = "Paired", mapColour = "PRGn", n_mapColour = 11, ...,
                                     plotWidth = 7, plotHeight = 7){ #DOI: fltered subset data of interest


  ## prepare matrix for plotting
  dfm <- data.frame(fltDOI$genes, fltDOI$E)
  dfm <- dfm[dfm$ControlType == 0, ] # remove control probes

  if (is.null(pcutoff)){
    stop("Please set p value threshold.")
  } else {
    ifelse(method == "fdr", pb_name <- dfmDE[dfmDE$adj.P.Val < pcutoff, dataProbeVar], pb_name <- dfmDE[dfmDE$P.Value < pcutoff, dataProbeVar])
    dfm <- dfm[dfm$ProbeName %in% pb_name, ] # the reason to use $ProbeName is because dfm is from fltDOI
  }

  ## set FC filter
  if (!is.null(FC)){
    pb_name_fc <- dfmDE[abs(dfmDE$logFC) >= log2(FC), dataProbeVar]
    dfm <- dfm[dfm$ProbeName %in% pb_name_fc, ]
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
        colnames(mtx) <- sampleName
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
        colnames(mtx) <- sampleName
      }

      pdf(file = paste(plotName, "_heatmap.supervised.pdf", sep = ""), width = plotWidth, height = plotHeight)
      heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
                col = brewer.pal(n_mapColour, mapColour), ColSideColors = colC[colG], labRow = labrow, ...)
      dev.off()

      print("Probes with no gene names are removed.")
    }

  } else {
    mtx <- as.matrix(dfm[, s])
    rownames(mtx) <- dfm[, dataProbeVar]

    if (!is.null(sampleName)){
      colnames(mtx) <- sampleName
    }

    pdf(file = paste(plotName, "_heatmap.supervised.pdf", sep = ""), width = plotWidth, height = plotHeight)
    heatmap.2(mtx, distfun = distfunc, hclustfun = clustfunc,
              col = brewer.pal(n_mapColour, mapColour), ColSideColors = colC[colG], labRow = FALSE,...)
    dev.off()
  }
}


#' @title rbioarray_corcluster_super
#'
#' @description Wrapper for supervised (or unsupervised) Pearson correlation clustering analysis and heatmap visualization for both microarray and RNAseq, for gene co-expression analysis.
#' @param plotName File name for the export \code{pdf} plot file. Default is \code{"data"}.
#' @param fltData Based on filtered data, a subset corresponding to the comparasion, either a list, \code{EList} or \code{MAList} object.
#' @param rmControl If to remove control probes (Agilent platform). Default is \code{TRUE}.
#' @param n_subgroup A vector of sample index (row number) for phenotype group. Default is \code{NULL}. The setting can be obtained from the corresponding condition summary object.
#' @param dfmDE A subset of the DE list, i.e. a \code{topTable} dataframe, corresponding to the comparasion (i.e., contrast).
#' @param dataProbeVar \code{dfmDE} variable name for probe name. Default is \code{NULL}.
#' @param FDR If to use FDR corrcted p value. Default is \code{TRUE}.
#' @param q.value P value cut off. Default is \code{0.05}. For unsupervised clustering, set \code{q.value = 1}.
#' @param FC Fold change (FC) filter for the heatmap. Default is \code{NULL}.
#' @param method Thresholding method, "fdr" or "spikein". Default is \code{"spikein"}.
#' @param axisLabel Whether to display label for both x- and y- axes. Default is \code{FALSE}.
#' @param annot The optional annotation matrix. Only needs to be set to display inforamtions for
#' @param genesymbolVar The name of the variable for gene symbols from the \code{annot} object. Only set this argument when \code{rowLabel = TRUE}. Default is \code{NULL}.
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
#'
#' # n_subgroup = c(1:4) means the correlation uses samples from 1 to 4 (control in this case).
#' # The settings can be obtained from the corresponding condition summary object.
#' rbioarray_corcluster_super(fltData = all_nrm, n_subgroup = c(1:4),
#'                            dataProbeVar = "gene_id", FDR = TRUE, q.value = 0.02,
#'                            dfmDE = all_DE$`conSumPost - conSumPre`,
#'                            axisLabel = TRUE, genesymbolVar = "gene_name",
#'                            key.title = "", cexRow = 0.3, cexCol = 0.3, offsetRow = 0.001,
#'                            offsetCol = 0.001, margins = c(4, 4))
#'
#' }
#' @export
rbioarray_corcluster_super <- function(plotName = "data",
                                       fltData = NULL, rmControl = TRUE,
                                       n_subgroup = NULL,
                                       dfmDE = NULL, FDR = TURE, q.value = 0.05, FC = NULL,
                                       dataProbeVar = NULL,
                                       axisLabel = FALSE,
                                       annot = NULL, genesymbolVar = NULL,
                                       mapColour = "PRGn", n_mapColour = 11, ...,
                                       plotWidth = 7, plotHeight = 7){

  #### test variables
  if (is.null(fltData)){
    stop(cat("Please set processed data object via fltData. Function terminated.\n"))
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

  #### fiter and normalization
  vmwt <- fltData
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
    dfm <- dfm[dfm$ProbeName %in% pb_name_fc, ]
  }

  #### heatmap
  ogNcol <- dim(vmwt$E)[2] # original numbers of col
  annoNcol <- dim(dfm)[2] # numbers of col with annotation
  s <- (annoNcol - ogNcol + 1):annoNcol # extract only the data by removing the annotation columns

  if (axisLabel){
    if (!is.null(genesymbolVar)){
      dfm <- dfm[complete.cases(dfm[, genesymbolVar]), ]
      axisrow <- dfm[, genesymbolVar]
      print("Probes with no gene names are removed.")

      mtx <- as.matrix(dfm[, s])
      rownames(mtx) <- dfm[, dataProbeVar]
      cormtx <- t(mtx)

      pdf(file = paste(plotName, "_corheatmap.supervised.pdf", sep = ""), width = plotWidth, height = plotHeight)
      heatmap.2(cor(cormtx[n_subgroup, ], method = "pearson"), symm = TRUE, trace = "none",
                col = brewer.pal(n_mapColour, mapColour), labRow = axisrow, labCol = axisrow, ...)
      dev.off()
    } else {
      print("No gene symbol variable detected. Proceed without one.")

      mtx <- as.matrix(dfm[, s])
      rownames(mtx) <- dfm[, dataProbeVar]
      cormtx <- t(mtx)

      pdf(file = paste(plotName, "_corheatmap.supervised.pdf", sep = ""), width = plotWidth, height = plotHeight)
      heatmap.2(cor(cormtx[n_subgroup, ], method = "pearson"), symm = TRUE, trace = "none",
                col = brewer.pal(n_mapColour, mapColour), labRow = FALSE, labCol = FALSE,...)
      dev.off()
    }
  } else {
    mtx <- as.matrix(dfm[, s])
    rownames(mtx) <- dfm[, dataProbeVar]
    cormtx <- t(mtx)

    pdf(file = paste(plotName, "_corheatmap.supervised.pdf", sep = ""), width = plotWidth, height = plotHeight)
    heatmap.2(cor(cormtx[n_subgroup, ], method = "pearson"), symm = TRUE, trace = "none",
              col = brewer.pal(n_mapColour, mapColour), labRow = FALSE, labCol = FALSE,...)
    dev.off()
  }
}


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
#' @param DE DE methods set for p value thresholding. Values are \code{"fdr"} and \code{"spikein"}. Default is \code{"fdr"}.
#' @param fltdata Only needed when \code{DE = "spikein"}. Filtered data, either a list, \code{EList} or \code{MAList} object. Default is \code{NULL}.
#' @param design Only needed when \code{DE = "spikein"}. Design matrix. Default is \code{NULL}.
#' @param contra Only needed when \code{DE = "spikein"}. Contrast matrix. Default is \code{NULL}.
#' @param weights Only needed when \code{DE = "spikein"}. Array weights, determined by \code{arrayWeights()} function from \code{limma} package. Default is \code{NULL}.
#' @param q.value P value threshold. Only needed for \code{DE = "fdr"} and \code{DE = "spikein"} when calculated p value is larger than q.value. Default is \code{0.05}.
#' @param FC Fold change threshold. Default is \code{1.5}.
#' @param parallelComputing If to use parallel computing. Default is \code{FALSE}.
#' @return The function outputs a \code{pdf} file for venn diagrams (total, up- and down-regulations). The function also exports overlapping gene or probe into a \code{csv} file.
#' @details When \code{"fdr"} set for DE, the p value threshold is set as \code{0.05}. When there is no significant genes or probes identified under \code{DE = "fdr"}, the threshold is set to \code{1}. If the arugments for \code{DE = "spikein"} are not complete, the function will automatically use \code{"fdr"}.
#' @import doParallel
#' @import foreach
#' @importFrom limma lmFit eBayes topTable contrasts.fit vennDiagram
#' @importFrom parallel detectCores makeCluster stopCluster
#' @examples
#' \dontrun{
#' rbioarray_venn_DE(plotName = "DE", cex = c(1, 2, 2), mar = rep(0.5,4), names = c("control", "stress1", "stress2"),
#'                   DEdata = fltdata_DE, geneName = TRUE, genesymbolVar = "GeneSymbol",
#'                   DE = "spikein", fltdata = fltdata, annot = annot, design = design, contra = contra, weights = fltdata$ArrayWeight,
#'                   parallelComputing = FALSE)
#' }
#' @export
rbioarray_venn_DE <- function(objTitle = "DE", plotName = "DE", plotWidth = 5, plotHeight = 5, ...,
                              annot = NULL,
                              DEdata = NULL, dataProbeVar = "ProbeName",
                              geneName = FALSE, annotProbeVar = "ProbeName", genesymbolVar = NULL,
                              DE = "fdr", fltdata = NULL, design = NULL, contra = NULL, weights = NULL, q.value = 0.05, FC = 1.5,
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
    if (DE == "fdr"){
      if (length(which(tmpdfm$adj.P.Val < q.value)) == 0){
        p_threshold <- 1
      } else {
        p_threshold <- max(tmpdfm[tmpdfm$adj.P.Val < q.value, ]$P.Value)
      }

    } else if (DE == "spikein") {
      # check arugments
      if (is.null(fltdata) | is.null(design) | is.null(contra) | is.null(weights)){
        warning(cat("Arguments not complete for \"spikein\" method. Proceed with \"fdr\" instead."))

        if (length(which(tmpdfm$adj.P.Val < q.value)) == 0){
          p_threshold <- 1
        } else {
          p_threshold <- max(tmpdfm[tmpdfm$adj.P.Val < q.value, ]$P.Value)
        }

      } else {
        # DE for extracting positive control
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
        cf <- names(m)
        PC <- fit[fit$genes$ControlType == 1, ]
        ifelse(min(PC$p.value[, cf[n]]) > q.value, p_threshold <- q.value, p_threshold <- min(PC$p.value[, cf[n]]))
      }

    } else {stop(cat("Please set p value thresholding method, \"fdr\" or \"spikein\"."))}
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
