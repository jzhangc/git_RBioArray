#' @title rbioseq_DE
#'
#' @description Legacy function. DE analysis function for RNA-seq data, with count filtering functionality.
#' @param objTitle Name for the output list. Default is \code{"seq_data"}.
#' @param dfm_count Dataframe contains the feature read counts, with rows as genomic featues (or genes) and column as samples. Default is \code{NULL}.
#' @param dfm_annot Dataframe contains the gene annotation information, with rows as genmic features and columns as annotation variables. The row lengths of this dataframe should be the same as \code{dfm_count}.
#' @param count_threshold Read count threshold. No filtering will be applied when set \code{"none"}. Otherwise, a numeric number can be set as the minimum read count for filtering. Default is \code{"none"}.
#' @param norm.method Normalization methods.Currently, the none-compositional methods are supported. Options are \code{"none", "TMM","RLE","upperquartile"}. Default is \code{"TMM"}. See \code{calcNormFactors} function from  \code{edgeR} pacakge for more.
#' @param qc_plot Wether or not to produce a QC plot upon filtering, normalization and weight calculation. Default is \code{TRUE}.
#' @param design Design matrix. Default is \code{NULL}.
#' @param contra Contrast matrix. Default is \code{NULL}.
#' @param ... arguments for \code{topTable()} from \code{limma} package.
#' @param plot If to generate volcano plots for the DE results. Defualt is \code{TRUE}. Plots are exported as \code{pdf} files.
#' @param geneName If to only plot probes with a gene name. Default is \code{FALSE}.
#' @param genesymbolVar The name of the variable for gene symbols from the \code{annot} object. Only set this argument when \code{geneName = TRUE}. Default is \code{NULL}.
#' @param topgeneLabel If to display the gene identification, i.e., probem name or gene name, on the plot. Default is \code{FALSE}.
#' @param nGeneSymbol When \code{topgeneLabel = TRUE}, to set how many genes to display. Default is \code{5}.
#' @param padding When \code{topgeneLabel = TRUE}, to set the distance between the dot and the gene symbol. Default is \code{0.5}.
#' @param FC Threshold for fold change (FC) for volcano plot. Default is \code{1.5}.
#' @param FDR Wether or not using FDR p value correction. Default is \code{TRUE}.
#' @param sig.p Threshold for the p value. Default is \code{0.05}.
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
#' @return The function outputs a \code{dataframe} object with filtered and normalized readcounts, a \code{list} object with DE results, a \code{dataframe} object for the F test results, merged with annotation. The function also exports DE reuslts to the working directory in \code{csv} format.
#' @details Sample weights are automatically caluclated during Voom normalization. When \code{count_threshold = 0}, the function will filter out all the genomic features with a total read count 0. When there is no significant genes or probes identified under \code{FDR = TRUE}, the function conducts DE analysis without p value correction and outputs an warning. Also note that both \code{geneName} and \code{genesymbolVar} need to be set to display gene sysmbols on the plot. Additionally, when set to display gene symbols, all the features without a gene symbol will be removed.
#' @import ggplot2
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @importFrom limma lmFit eBayes topTable contrasts.fit voomWithQualityWeights
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom grid grid.newpage grid.draw
#' @importFrom ggrepel geom_text_repel
#' @examples
#' \dontrun{
#' rbioseq_DE(objTitle = "test", dfm_count = Ann_Count_all[, 5:12], dfm_annot = Ann_Count_all[, 1:4], count_threshold = 50,
#'            design = design, contra = contra,
#'            FDR = TRUE, sig.p = 0.05,
#'            geneName = TRUE, genesymbolVar = "gene_name", topgeneLabel = TRUE, nGeneSymbol = 10)
#' }
#' @export
rbioseq_DE <- function(objTitle = "data_filtered", dfm_count = NULL, dfm_annot = NULL,
                       count_threshold = "none",
                       norm.method = "TMM",
                       qc_plot = TRUE,
                       design = NULL, contra = NULL,
                       ...,
                       plot = TRUE, geneName = FALSE, genesymbolVar = NULL, topgeneLabel = FALSE, nGeneSymbol = 5, padding = 0.5,
                       FC = 1.5, FDR = TRUE, sig.p = 0.05,
                       Title = NULL, xLabel = "log2(fold change)", yLabel = "-log10(p value)",
                       symbolSize = 2, sigColour = "red", nonsigColour = "gray",
                       xTxtSize = 10, yTxtSize =10,
                       plotWidth = 170, plotHeight = 150,
                       parallelComputing = FALSE, clusterType = "PSOCK"){
  ## check the key arguments
  if (!class(dfm_count) %in% c("data.frame", "matrix")){
    stop("dfm_count has to be either a data.frame or matrix object")
  }

  if (!class(dfm_annot) %in% c("data.frame", "matrix")){
    stop("dfm_annot has to be either a data.frame or matrix object")
  }

  if (nrow(dfm_count) != nrow(dfm_annot)){
    stop("Read count matrix doesn't have the same row number as the annotation matrix.")
  }

  if (is.null(design)){
    stop("Please set design matrix.")
  }

  if (is.null(contra)){
    stop("Please set contrast object.")
  }

  if (!norm.method %in% c("none", "TMM","RLE","upperquartile")){
    stop("Please set the norm.method with exactly one of the following: \"none\", \"TMM\", \"RLE\", \"upperquartile\".")
  }

  ## extract coefficients
  cf <- colnames(contra) # extract coefficient

  # set an empty matrix for exporting the threolding summery
  threshold_summary <- matrix(nrow = length(cf), ncol = 5)
  colnames(threshold_summary) <- c("coeffcient", "p.value.threshold", "fold.change.threshold", "True", "False")
  threshold_summary <- as.matrix(threshold_summary)

  ## temp func for plotting
  # i: outlist (object) listed below
  tmpfunc <- function(i, j){
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
    if (FDR){
      if (length(which(tmpdfm$adj.P.Val < sig.p)) == 0){
        warning("No significant results found using FDR correction. Please consider using another thresholding method. For now, sig.p is applied on raw p.values.")
        pcutoff <- sig.p
      } else {
        pcutoff <- max(tmpdfm[tmpdfm$adj.P.Val < sig.p, ]$P.Value)
      }
      cutoff <- as.factor(abs(tmpdfm$logFC) >= log2(FC) & tmpdfm$P.Value < pcutoff)
    } else  {
      pcutoff <- sig.p
      cutoff <- as.factor(abs(tmpdfm$logFC) >= log2(FC) & tmpdfm$P.Value < pcutoff)
    }

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
    pltgtb <- rightside_y(plt) # RBioplot::rightside_y() for displying rightside y-axis
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
  # create DGE object using edgeR
  cat("Data filtering and normalization...") # message
  dge <- DGEList(counts = dfm_count, genes = dfm_annot)

  if (count_threshold != "none"){ # set the count threshold for filtering
    count_s <- rowSums(dge$counts) # thresholdd
    isexpr <- count_s > count_threshold
    dge <- dge[isexpr, , keep.lib.size = FALSE] # filtering
  }

  # for data Voom normalization
  dgenormf <- calcNormFactors(dge, method = norm.method)
  vmwt <- voomWithQualityWeights(dgenormf, design = design, plot = qc_plot, normalization = "quantile") # Voom normalization with quality weights
  cat("DONE!\n") # message

  # fitting
  cat("Linear fitting...") # message
  fit <- lmFit(vmwt, design = design) # linear fitting
  fit <- contrasts.fit(fit, contrasts = contra)
  fit <- eBayes(fit)
  cat("DONE!\n") # message
  out <- fit

  ## output and plotting
  if(!parallelComputing){
    # compile resutls into a list
    outlist <- lapply(cf, function(i){
      tmp <- topTable(out, coef = i, number = Inf, ...)
      return(tmp)
    })

    names(outlist) <- cf
    # write DE results into files
    lapply(1:length(cf), function(j){
      write.csv(outlist[[j]], file = paste(cf[[j]], "_DE.csv", sep = ""), na = "NA", row.names = FALSE)
    })

    # volcano plot and output summary
    if (plot){
      threshold_summary[] <- t(sapply(1: length(cf), function(x)tmpfunc(i = outlist, j = x)))
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
          tmpfunc(i = outlist, j = j)
        }
      }
    } else { # macOS and Unix-like systems
      outlist <- mclapply(cf, FUN = function(i){
        tmp <- topTable(out, coef = i, number = Inf, ...)
        return(tmp)
      }, mc.cores = n_cores, mc.preschedule = FALSE)

      names(outlist) <- cf
      # write DE results into files
      mclapply(1:length(cf), FUN = function(j){
        write.csv(outlist[[j]], file = paste(cf[[j]], "_DE.csv", sep = ""), na = "NA", row.names = FALSE)
      }, mc.cores = n_cores, mc.preschedule = FALSE)

      # volcano plot and output summary
      if (plot){
        threshold_summary[] <- t(sapply(1: length(cf), function(x)tmpfunc(i = outlist, j = x)))
      }
    }
  }

  ## output the DE/fit objects to the environment, as well as the DE csv files into wd
  fitout <- topTable(out, number = Inf)
  assign(paste(objTitle, "_nrm", sep = ""), vmwt, envir = .GlobalEnv)
  assign(paste(objTitle, "_fit", sep = ""), fitout, envir = .GlobalEnv)
  write.csv(fitout, file = paste(objTitle, "_DE_Fstats.csv", sep = ""), row.names = FALSE)
  assign(paste(objTitle, "_DE", sep = ""), outlist, envir = .GlobalEnv)
  write.csv(threshold_summary, file = paste(objTitle, "_thresholding_summary.csv", sep = ""), row.names = FALSE)
  ## message
  if (geneName & !is.null(genesymbolVar)){
    print("Probes without a gene symbol are removed from the volcano plots")
  }
}


#' @title rbioarray_PreProc
#'
#' @description Legacy function. Data pre-processing function for the microarary data.
#' @param rawlist Input data, either a list, \code{EList} or \code{MAList} object.
#' @param logTrans If to perfom a log transformation on the data or not. Default is \code{FALSE}.
#' @param logTransMethod If \code{logTrans = TRUE}, set which method to use for the transformation, \code{"log2"} or \code{"log10"}. Default is \code{"log2"}.
#' @param logTransObjT If \code{logTrans = TRUE}, set the file name for the output \code{csv} file containing the log transformed data.
#' @param logTransParallelComputing If \code{logTrans = TRUE}, set if to use parallel computing for the transformation or not. Default is \code{FALSE}.
#' @param bgMethod Background correction method. Default is \code{"auto"}. See \code{backgroundCorrect()} function from \code{limma} package for details.
#' @param normMethod Normalization method. Default is \code{"quantile"}. See \code{normalizeBetweenArrays()} function from \code{limma} package for details.
#' @param ... arguments for \code{backgroundCorrect.matrix()} or \code{backgroundCorrect()} functions from \code{limma} package.
#' @details The function does not use design matrix for array weight calculation.
#'          Therefore, the DE analysis based on the output from this function will yeild slightly different resutls from the \code{\link{rbioarray_filter_combine}}.
#'
#' @return Depending on the input type, the function outputs a \code{list}, \code{Elist} or \code{MAList} object with corrected and normalized expression values.
#'         If \code{logTrans = TRUE}, the function also outputs a \code{csv} file containing the log transformed data.
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
#' @description Legacy function. Data filter function based on spike-in negative control.
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
      stop(cat("Since combineGeneDup = TRUE and normlst is not an Elist object, please set annotation dataframe for annot argument.
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
    neg <- apply(normlst$E, 2, function(x)quantile(x, p = percentile)) # 5% percentile of all the data
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


#' @title rbioarray_DE
#'
#' @description Legacy function. DE analysis function.
#' @param objTitle Name for the output list. Default is \code{"data_filtered"}.
#' @param output.mode Csv file output mode. Options are \code{"probe.all"}, \code{"probe.sig"}, \code{"geneName.all"} and \code{"geneName.sig"}. Default is \code{"probe.all"}. See details below.
#' @param fltlist filtered data, either a list, \code{EList} or \code{MAList} object. Default is \code{NULL}.
#' @param annot Annotation object, usually a \code{dataframe}. Make sure to name the probe ID variable \code{ProbeName}. Default is \code{NULL}.
#' @param design Design matrix. Default is \code{NULL}.
#' @param contra Contrast matrix. Default is \code{NULL}.
#' @param weights Array weights, determined by \code{arrayWeights()} function from \code{limma} package. Default is \code{NULL}.
#' @param ... arguments for \code{topTable()} from \code{limma} package.
#' @param plot If to generate volcano plots for the DE results. Defualt is \code{TRUE}. Plots are exported as \code{pdf} files.
#' @param geneName If to only plot probes with a gene name. Default is \code{FALSE}.
#' @param genesymbolVar The name of the variable for gene symbols from the \code{annot} object. Only set this argument when \code{geneName = TRUE}. Default is \code{NULL}.
#' @param topgeneLabel If to display the gene identification, i.e., probem name or gene name, on the plot. Default is \code{FALSE}.
#' @param nGeneSymbol When \code{topgeneLabel = TRUE}, to set how many genes to display. Default is \code{5}.
#' @param padding When \code{topgeneLabel = TRUE}, to set the distance between the dot and the gene symbol. Default is \code{0.5}.
#' @param FC Threshold for fold change (FC) for volcano plot. Default is \code{1.5}.
#' @param ctrlProbe Wether or not the data set has control type variable, with values \code{-1 (negative control)}, \code{0 (gene probes)} and \code{1 (positive control)}. Default is \code{TRUE}.
#' @param ctrlTypeVar Set only when \code{ctrlProbe = TRUE}, the control type variable. Default is the \code{Agilent} variable name \code{"ControlType"}.
#' @param sig.method DE methods set for p value thresholding. Values are \code{"fdr"}, \code{"spikein"} and \code{"none"}. Default is \code{"fdr"}. \code{"spikein"} can only be set when \code{ctrlProbe = TRUE} and \code{ctrlTypeVar} is properly set.
#' @param sig.p Only used when sig.method set as \code{"spikein"}, backup threshold for the p value if spikein p values is larger than \code{0.05}.
#' @param plotTitle Figure title. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
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
#' @details When \code{"fdr"} set for sig.method, the p value threshold is set as \code{0.05}.
#'
#'          When there is no significant genes or probes identified under \code{sig.method = "fdr"}, the threshold is set to \code{sig.p}. When set \code{sig.method = "none"}, the p cutoff will be \code{sig.p}.
#'          Also note that both \code{geneName} and \code{genesymbolVar} need to be set to display gene sysmbols on the plot. Otherwise, the labels will be probe names. Additionally, when set to display gene symbols,
#'          all the probes without a gene symbol will be removed.
#'
#'          The \code{output.mode} argument only affacts output csv files and also not for F test results. This means DE and DE F test resutls for all probes will always be stored to the R environment.
#'          And F test results for all probes will always be exported as csv file to the working directory.
#'
#' @import ggplot2
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @importFrom limma lmFit eBayes topTable contrasts.fit
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y
#' @importFrom ggrepel geom_text_repel
#' @examples
#' \dontrun{
#' rbioarray_DE(objTitle = "fltdata2", fltlist, annot = Anno, design, contra = contra,
#'              weights = fltdata$ArrayWeight, parallelComputing = TRUE,
#'              plot = TRUE, geneName = TRUE, genesymbolVar = "GeneSymbol",
#'              DE = "spikein")
#' }
#' @export
rbioarray_DE <- function(objTitle = "data_filtered", output.mode = "probe.all",
                         fltlist = NULL, annot = NULL,
                         design = NULL, contra = NULL, weights = NULL,
                         ...,
                         plot = TRUE,
                         geneName = FALSE, genesymbolVar = NULL, topgeneLabel = FALSE, nGeneSymbol = 5, padding = 0.5,
                         FC = 1.5,
                         ctrlProbe = TRUE, ctrlTypeVar = "ControlType", sig.method = "fdr", sig.p = 0.05,
                         plotTitle = NULL, xLabel = "log2(fold change)", yLabel = "-log10(p value)",
                         symbolSize = 2, sigColour = "red", nonsigColour = "gray",
                         xTxtSize = 10, yTxtSize =10,
                         plotWidth = 170, plotHeight = 150,
                         parallelComputing = FALSE, clusterType = "PSOCK"){
  ## check the key arguments
  if (is.null(fltlist)){
    stop("Please set input data object. Hint: either a list, EList or MAList object with pre-processed and flitered expression data. Function terminated.\n")
  }
  if (!is.null(annot)){ # check if the variable "ProbeName" is included in the annotation dataframe
    if (!"ProbeName" %in% names(annot)){
      stop(cat("For the annotation dataframe, make sure to name the variable containing probe name \"ProbeName\". Function terminated.\n"))
    }
  }
  if (is.null(design)){
    stop("Please set design matrix. Function terminated.\n")
  }
  if (is.null(contra)){
    stop("Please set contrast object. Function terminated.\n")
  }
  if (tolower(sig.method) == "spikein"){
    if (!ctrlProbe){
      stop(cat("\"spiken\" DE method can only be set when ctrlProbe = TRUE and ctrlTypeVar is properly set. Function terminated.\n"))
    }
  }
  if (ctrlProbe){
    if (!ctrlTypeVar %in% names(fltlist$genes)){
      stop(cat("ctrlTypeVar not found. Function terminated.\n"))
    }
  }

  if (!output.mode %in% c("probe.all", "probe.sig", "geneName.all", "geneName.sig")) {
    stop("output.mode should be one of \"probe.all\", \"probe.sig\", \"geneName.all\" or \"geneName.sig\".")
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
    # import DE data
    tmpdfm <- i[[j]]

    # set the cutoff
    if (tolower(sig.method) == "fdr"){
      if (length(which(tmpdfm$adj.P.Val < sig.p)) == 0){
        warning("No significant results found using FDR correction. Please consider using another thresholding method. For now, sig.p is applied on raw p.values.")
        pcutoff <- sig.p
      } else {
        pcutoff <- max(tmpdfm[tmpdfm$adj.P.Val < sig.p, ]$P.Value)
      }
    } else if (tolower(sig.method) == "spikein") {
      ifelse(min(PC$p.value[, cf[j]]) > sig.p, pcutoff <- sig.p, pcutoff <- min(PC$p.value[, cf[j]]))
    } else if (tolower(sig.method) == "none"){
      pcutoff <- sig.p
    } else {stop(cat("Please set p value thresholding method, \"fdr\" or \"spikein\"."))}

    # set the data frame
    if (geneName){
      if (!genesymbolVar %in% names(tmpdfm)){
        stop(cat("Invalid gene symbol variable"))
      }

      if (!is.null(genesymbolVar)){
        pltdfm <- tmpdfm[complete.cases(tmpdfm[, genesymbolVar]), ]
      } else {
        warning("No variable name for gene symbol set. Proceed with probe names with no probes removed.")
        pltdfm <- tmpdfm
      }
    } else {
      pltdfm <- tmpdfm
    }
    cutoff <- as.factor(abs(pltdfm$logFC) >= log2(FC) & pltdfm$P.Value < pcutoff)

    # plot
    loclEnv <- environment()
    plt <- ggplot(pltdfm, aes(x = logFC, y = -log10(P.Value)), environment = loclEnv) +
      geom_point(alpha = 0.4, size = symbolSize, aes(colour = cutoff)) +
      scale_color_manual(values = c(nonsigColour, sigColour)) +
      ggtitle(plotTitle) +
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
      tmpfltdfm <- pltdfm[abs(pltdfm$logFC) >= log2(FC) & pltdfm$P.Value < pcutoff, ]
      tmpfltdfm <- tmpfltdfm[order(tmpfltdfm$P.Value), ]
      plt <- plt + geom_text_repel(data = head(tmpfltdfm, n = nGeneSymbol),
                                   aes(x = logFC, y = -log10(P.Value), label = head(tmpfltdfm, n = nGeneSymbol)[, genesymbolVar]),
                                   point.padding = unit(padding, "lines"))
    }

    grid.newpage()
    pltgtb <- rightside_y(plt) # RBioplot::rightside_y() for displying rightside y-axis
    # export the file and draw a preview
    ggsave(filename = paste(objTitle, "_", cf[[j]],".volcano.pdf", sep = ""), plot = pltgtb,
           width = plotWidth, height = plotHeight, units = "mm",dpi = 600)
    grid.draw(pltgtb) # preview

    # save DE results files
    if (output.mode == "probe.all"){
      outdfm <- tmpdfm
    } else if (output.mode == "probe.sig"){
      # below: make sure to use all probes, i.e.tmpdfm
      outdfm <- tmpdfm[which(abs(tmpdfm$logFC) >= log2(FC) & tmpdfm$P.Value < pcutoff), ]
    } else if (output.mode == "geneName.all"){
      outdfm <- pltdfm
    } else {
      outdfm <- pltdfm[which(as.logical(cutoff)), ]  # it is ok to use cutoff straigtly here
    }
    write.csv(outdfm, file = paste(objTitle, "_", cf[[j]], "_DE.csv", sep = ""), na = "NA", row.names = FALSE)

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
  cat("Done!\n") # message

  if (ctrlProbe){
    out <- fit[fit$genes[, ctrlTypeVar] == 0, ] # remove control probes
  } else {
    out <- fit
  }

  if (tolower(sig.method) == "spikein"){ # extract PC stats for spikein method
    PCntl <- fit[fit$genes[, ctrlTypeVar] == 1, ]
  } else {
    PCntl <- NULL
  }

  ## output and plotting
  if(!parallelComputing){
    cat("Plotting and exporting files...") # message
    # compile resutls into a list
    outlist <- lapply(cf, function(i){
      tmp <- topTable(out, coef = i, number = Inf, ...)
      tmp$ProbeName <- rownames(tmp)
      if (!is.null(annot)){ # merge with annotation dataframe
        # tmp <- merge(tmp, annot, by = "ProbeName", all.x = TRUE)
        tmp <- merge(tmp, annot, all.x = TRUE)
      }
      return(tmp)
    })
    names(outlist) <- cf

    # volcano plot and output summary
    if (plot){
      threshold_summary[] <- t(sapply(1: length(cf), function(x)tmpfunc(i = outlist, j = x, PC = PCntl)))
    }
    cat("Done!\n") # message
  } else { ## parallel computing
    cat("Plotting and exporting files...") # message
    # check the cluster type
    if (clusterType != "PSOCK" & clusterType != "FORK"){
      stop("Please set the cluter type. Options are \"PSOCK\" (default) and \"FORK\".")
    }

    # set up cpu cores and cluster
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores, type = clusterType)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    outlist <- foreach(i = 1: length(cf), .packages = "limma") %dopar% {
      tmp <- limma::topTable(out, coef = cf[i], number = Inf, ...)
      tmp$ProbeName <- rownames(tmp)
      if (!is.null(annot)){ # merge with annotation dataframe
        tmp <- merge(tmp, annot, by = "ProbeName")
      }
      return(tmp)
    }
    names(outlist) <- cf

    # volcano plot and output summary
    if (plot){
      if (clusterType == "PSOCK"){
        threshold_summary[] <- foreach(j = 1: length(cf), .combine = "rbind",
                                       .packages = c("limma", "ggplot2", "gtable", "grid", "ggrepel", "RBioplot")) %dopar% {
                                         tmpfunc(i = outlist, j = j, PC = PCntl)
                                       }
      } else {
        threshold_summary[] <- foreach(j = 1: length(cf), .combine = "rbind",
                                       .packages = c("limma", "ggplot2", "gtable", "grid", "ggrepel", "RBioplot")) %do% {
                                         tmpfunc(i = outlist, j = j, PC = PCntl)
                                       }
      }
    }
    cat("Done!\n") # message
  }

  ## output the DE/fit objects to the environment, as well as the DE csv files into wd
  fitout <- topTable(out, number = Inf)
  if (!is.null(annot)){ # merge with annotation dataframe
    fitout <- merge(fitout, annot, by = "ProbeName")
  }
  assign(paste(objTitle, "_fit", sep = ""), fitout, envir = .GlobalEnv)
  write.csv(fitout, file = paste(objTitle, "_DE_Fstats.csv", sep = ""), row.names = FALSE)
  assign(paste(objTitle, "_DE", sep = ""), outlist, envir = .GlobalEnv)
  write.csv(threshold_summary, file = paste(objTitle, "_thresholding_summary.csv", sep = ""), row.names = FALSE)

  ## messages
  if (geneName & !is.null(genesymbolVar)){
    cat("Probes without a gene symbol are removed from the volcano plots.\n")
  }

  if (output.mode == "probe.all") cat("DE results for all probes exported into csv files. \n")
  if (output.mode == "probe.sig") cat("DE results for significant probes regardless of gene name saved to csv files. \n")
  if (output.mode == "geneName.all") {
    if (!geneName){
      cat("geneName = FALSE, DE results for all probes saved to csv files,
          even with output.mode = \"geneName.all\"\n")
    } else {
      cat("DE results for probes with gene name saved to csv files. \n")
    }
  }
  if (output.mode == "geneName.sig") {
    if (!geneName){
      cat("geneName = FALSE, DE results for significant probes regardless of gene name saved to csv files,
          even with output.mode = \"geneName.sig\"\n")
    } else {
      cat("DE results for significant probes with gene name saved to csv files. \n")
    }
  }
}


#' @title rbioarray_corcluster_super
#'
#' @description Legacy function. Wrapper for supervised (or unsupervised) Pearson correlation clustering analysis and heatmap visualization for both microarray and RNAseq, for gene co-expression analysis.
#' @param plotName File name for the export \code{pdf} plot file. Default is \code{"data"}.
#' @param fltlist Based on filtered data, a subset corresponding to the comparasion, either a list, \code{EList} or \code{MAList} object.
#' @param rmControl If to remove control probes (Agilent platform). Default is \code{TRUE}.
#' @param n_subgroup A vector of sample index (row number) for phenotype group. Default is \code{NULL}. The setting can be obtained from the corresponding condition summary object.
#' @param dfmDE A subset of the DE list, i.e. a \code{topTable} dataframe, corresponding to the comparasion (i.e., contrast).
#' @param FDR If to use FDR corrcted p value. Default is \code{TRUE}.
#' @param DE.sig.p P value cut off. Default is \code{0.05}. For unsupervised clustering, set \code{DE.sig.p = 1}.
#' @param FC Fold change (FC) filter for the heatmap. Default is \code{NULL}.
#' @param dataProbeVar \code{dfmDE} variable name for probe name. Default is \code{NULL}.
#' @param method The correlation method, options are "pearson", "spearman" and "pearson". Default is \code{"pearson"}.
#' @param sigPlot If to generate a significance heatmap. Default is \code{FALSE},
#' @param sigPlot.sig.FDR If to use FDR corrected correlation p value to plot sigPlot. Default is \code{FALSE}.
#' @param sigPlot.sig Only set when \code{sigPlot = TRUE}, the alpha value for correlation p value. Default is \code{0.05}
#' @param sigPlot.sigLabelColour Only set when \code{sigPlot = TRUE}, the colour for label for the significant pairs. Default is \code{"red"}.
#' @param sigPlot.sigLabelSize Only set when \code{sigPlot = TRUE}, the size for label for the significant pairs. Default is \code{3}.
#' @param sigPlot.labelColour Only set when \code{sigPlot = TRUE}, the significance heatmap axis label colour. Default is \code{"black"}.
#' @param sigPlot.labelSize Only set when \code{sigPlot = TRUE}, the significance heatmap axis label size. Default is \code{1}.
#' @param sigPlot.labelAngle Only set when \code{sigPlot = TRUE}, the significance heatmap axis label angle. Default is \code{90}.
#' @param sigPlot.keySize Only set when \code{sigPlot = TRUE}, the significance heatmap colour key size. Default is \code{1}.
#' @param axisLabel Whether to display label for both x- and y- axes. Default is \code{FALSE}.
#' @param annot The optional annotation matrix. Only needs to be set if \code{axisLabel = TRUE} AND if there is no genesymbolVar in the input data.
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
#'                            dataProbeVar = "gene_id", FDR = TRUE, DE.sig.p = 0.02,
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
                                       dfmDE = NULL, FDR = TURE, DE.sig.p = 0.05, FC = NULL,
                                       dataProbeVar = NULL,
                                       method = "pearson",
                                       sigPlot = FALSE, sigPlot.sig.FDR = FALSE, sigPlot.sig = 0.05,
                                       sigPlot.sigLabelColour = "red", sigPlot.sigLabelSize = 3,
                                       sigPlot.labelColour = "black", sigPlot.labelSize = 1, sigPlot.labelAngle = 90, sigPlot.keySize = 1,
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
      adj_corp <- matrix(p.adjust(corp, method = "fdr"), nrow = nrow(corp), byrow = T)  # fdr
      rownames(adj_corp) <- rownames(corp)
      colnames(adj_corp) <- colnames(corp)

      pdf(file = paste(plotName, "_corheatmap.pdf", sep = ""), width = plotWidth, height = plotHeight)
      heatmap.2(corcoef, symm = TRUE, trace = "none",
                col = brewer.pal(n_mapColour, mapColour), labRow = rownames(corcoef), labCol = colnames(corcoef), ...)
      garbage <- dev.off()

      if (sigPlot){
        if (sigPlot.sig.FDR){  # fdr correction or not
          sigplotmtx <- adj_corp
        } else {
          sigplotmtx <- corp
        }
        tryCatch(
          {
            pdf(file = paste(plotName, "_corheatmap.sigplot.pdf", sep = ""), width = plotWidth, height = plotHeight)
            corrplot(corr = corcoef, method = "color", type = "upper", p.mat = sigplotmtx, sig.level = sigPlot.sig,
                     insig = c("label_sig"), pch.col = sigPlot.sigLabelColour, pch.cex = sigPlot.sigLabelSize,
                     col = brewer.pal(n_mapColour, mapColour),
                     tl.col = sigPlot.labelColour, tl.cex = sigPlot.labelSize, tl.srt = sigPlot.labelAngle, cl.length = 3, cl.cex = sigPlot.keySize)
            garbage <- dev.off()
          },
          error = function(err){
            print("No significance found. Therefore no significance plot generated.")
            garbage <- dev.off()
          }
        )
      }
    } else {
      print("No gene symbol variable or annotation dataframe detected. Proceed without them.")

      mtx <- as.matrix(dfm[, s])
      rownames(mtx) <- dfm[, dataProbeVar]
      cormtx <- t(mtx)
      corcoef <- cor(cormtx[n_subgroup, ], method = method)
      corp <- foreach(i = corcoef, .combine = "cbind") %do% cor_pvalue(i, n = nrow(cormtx[n_subgroup, ])) # p value matrix
      diag(corp) <- NA
      adj_corp <- matrix(p.adjust(corp, method = "fdr"), nrow = nrow(corp), byrow = T)  # fdr
      rownames(adj_corp) <- rownames(corp)
      colnames(adj_corp) <- colnames(corp)

      pdf(file = paste(plotName, "_corheatmap.pdf", sep = ""), width = plotWidth, height = plotHeight)
      heatmap.2(corcoef, symm = TRUE, trace = "none",
                col = brewer.pal(n_mapColour, mapColour), labRow = rownames(corcoef), labCol = colnames(corcoef),...)
      garbage <- dev.off()

      if (sigPlot){
        if (sigPlot.sig.FDR){  # fdr correction or not
          sigplotmtx <- adj_corp
        } else {
          sigplotmtx <- corp
        }

        tryCatch(
          {
            pdf(file = paste(plotName, "_corheatmap.sigplot.pdf", sep = ""), width = plotWidth, height = plotHeight)
            corrplot(corr = corcoef, method = "color", type = "upper", p.mat = sigplotmtx, sig.level = sigPlot.sig,
                     insig = c("label_sig"), pch.col = sigPlot.sigLabelColour, pch.cex = sigPlot.sigLabelSize,
                     col = brewer.pal(n_mapColour, mapColour),
                     tl.col = sigPlot.labelColour, tl.cex = sigPlot.labelSize, tl.srt = sigPlot.labelAngle, cl.length = 3, cl.cex = sigPlot.keySize)
            garbage <- dev.off()
          },
          error = function(err){
            print("No significant correlation found. Therefore no significance plot generated.")
            garbage <- dev.off()
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
    adj_corp <- matrix(p.adjust(corp, method = "fdr"), nrow = nrow(corp), byrow = T)  # fdr
    rownames(adj_corp) <- rownames(corp)
    colnames(adj_corp) <- colnames(corp)

    pdf(file = paste(plotName, "_corheatmap.pdf", sep = ""), width = plotWidth, height = plotHeight)
    heatmap.2(corcoef, symm = TRUE, trace = "none",
              col = brewer.pal(n_mapColour, mapColour), labRow = FALSE, labCol = FALSE,...)
    garbage <- dev.off()

    if (sigPlot){
      if (sigPlot.sig.FDR){  # fdr correction or not
        sigplotmtx <- adj_corp
      } else {
        sigplotmtx <- corp
      }

      tryCatch(
        {
          pdf(file = paste(plotName, "_corheatmap.sigplot.pdf", sep = ""), width = plotWidth, height = plotHeight)
          corrplot(corr = corcoef, method = "color", type = "upper", p.mat = sigplotmtx, sig.level = sigPlot.sig,
                   insig = c("label_sig"), pch.col = sigPlot.sigLabelColour, pch.cex = sigPlot.sigLabelSize,
                   col = brewer.pal(n_mapColour, mapColour),
                   tl.col = sigPlot.labelColour, tl.cex = sigPlot.labelSize, tl.srt = sigPlot.labelAngle, cl.length = 3, cl.cex = sigPlot.keySize)
          garbage <- dev.off()
        },
        error = function(err){
          print("No significance found. Therefore no significance plot generated.")
          garbage <- dev.off()
        }
      )
    }
  }

  # export correlation matrix
  outdfm1 <- data.frame(
    group = paste(rownames(corcoef)[row(corcoef)[upper.tri(corcoef)]], colnames(corcoef)[col(corcoef)[upper.tri(corcoef)]], sep = "_"),
    row = rownames(corcoef)[row(corcoef)[upper.tri(corcoef)]],
    col = colnames(corcoef)[col(corcoef)[upper.tri(corcoef)]],
    coefficient = corcoef[upper.tri(corcoef)]
  )

  outdfm2 <- data.frame(
    group = paste(rownames(corp)[row(corp)[upper.tri(corp)]], colnames(corp)[col(corp)[upper.tri(corp)]], sep = "_"),
    row = rownames(corp)[row(corp)[upper.tri(corp)]],
    col = colnames(corp)[col(corp)[upper.tri(corp)]],
    p.value = corp[upper.tri(corp)]
  )

  outdfm3 <- data.frame(
    group = paste(rownames(adj_corp)[row(adj_corp)[upper.tri(adj_corp)]], colnames(adj_corp)[col(adj_corp)[upper.tri(adj_corp)]], sep = "_"),
    row = rownames(adj_corp)[row(adj_corp)[upper.tri(adj_corp)]],
    col = colnames(adj_corp)[col(adj_corp)[upper.tri(adj_corp)]],
    adj.p.value = adj_corp[upper.tri(adj_corp)]
  )

  out <- merge(outdfm1, outdfm2[, c("group", "p.value")], by = "group", x.all = TRUE)
  out <- merge(out, outdfm3[, c("group", "adj.p.value")], x.all = TRUE)

  write.csv(out, file = paste(plotName, ".cor.csv", sep = ""), row.names = FALSE)
}


#' @title rbioarray_hcluster
#'
#' @description Legacy function. Wrapper for hierarchical clustering analysis and heatmap visualization.
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
    stop(cat("Please provide sample index with argument fct. Function terminated.\nFunction terminated.\n"))
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
  garbage <- dev.off()
}


#' @title rbioseq_hcluster
#'
#' @description Legacy function. Wrapper for hierarchical clustering analysis and heatmap visualization for RNA seq data.
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
  garbage <- dev.off()
}


#' @title rbioarray_hcluster_super
#'
#' @description Legacy function. Wrapper for supervised hierarchical clustering analysis and heatmap visualization for both microarray and RNAseq.
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
      garbage <- dev.off()

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
      garbage <- dev.off()
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
    garbage <- dev.off()
  }
}

#' @title rbioarray_venn_DE
#'
#' @description (Legacy function) Venn diagrame for DE results. Uses \code{vennDiagram()} from \code{limma} package.
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
  garbage <- dev.off()

  pdf(file = paste(plotName, "_venn_up.pdf", sep = ""), width = plotWidth, height = plotHeight)
  vennDiagram(mtx, circle.col = 1:length(names(DEdata)), include = "up", ...)
  garbage <- dev.off()

  pdf(file = paste(plotName, "_venn_down.pdf", sep = ""), width = plotWidth, height = plotHeight)
  vennDiagram(mtx, circle.col = 1:length(names(DEdata)), include = "down", ...)
  garbage <- dev.off()

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

