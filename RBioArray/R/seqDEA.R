#' @title rbioseq_clr_ilr_transfo
#'
#' @description Log ratio tansformation function for read count data. Row: sample, column: features.
#' @param x Input read count data matrix.
#' @param offset Read count offset value added to avoid zero. Default is \code{0}.
#' @param mode Log ratio transformation method. Options are "clr" (centered log transformation) and "ilr" (isometric log transformation). Default is \code{"clr"}.
#' @param ilr.method.fast Useful only when \code{mode = "ilr"}. Default is \code{TRUE}.
#' @return A data matrix with log ratio transformed values.
#' @details This funciton is needed as part of the data pre-processing procedure to run LPS-DA analysis featured in \code{RBioFS} package. As per Quinn et al. (2018), NGS data can be considered as compositional data. As such, data must be transformed for usual statistical analysis visualization. Log ratio transfromation serves such purpose. Note that the number of features will be one less when using "ilr" method. It is not to be combined with the other normalization methods featured in the \code{\link{rbioseq_DE}}. Ref: Quinn TP, et al. 2018. Understanding sequencing data as compositions: an outlook and review. Bioinformatics. 2018: 1 - 9.
#' @examples
#' \dontrun{
#' tstX <- rbioseq_clr_ilr_transfo(tstdata, offset = 1, mode = "clr")
#' }
#' @export
rbioseq_clr_ilr_transfo <- function(x, offset = 0, mode = "clr", ilr.method.fast = TRUE){
  # data and arguments check
  if (!is.matrix(x))stop("x needs to b e a matrix")
  if (any(x == 0) & offset == 0)stop("zero detected in x. set offset to avoid it for ratio transformation")
  if (!tolower(mode) %in% c("clr", "ilr"))stop("choose the proper transformation mode: \"clr\" or \"ilr\"")

  # log ratio transformation
  if (tolower(mode) == "clr"){  # clr calculation
    if (dim(x)[2] == 1){
      cat()
      out <- list(x.clr = x, gm = rep(1, dim(x)[1]))
    } else {
      gm <- apply(x, 1, function(x) exp(mean(log(x + offset))))  # geometric mean = exp(mean(log(X)))
      clrX <- log((x + offset) / (gm))  # clr formula
      out <- clrX
    }
  } else if (tolower(mode) == "ilr"){  # ilr calculation, modified from ilr.transfo function from mixOmics package
    cat(paste0("ilr mode result in one less variable: ", ncol(x) - 1, " variables left for this dataset upon transformation.\n"))
    ilrX = matrix(NA, nrow = nrow(x), ncol = ncol(x) - 1)
    D = ncol(x)
    if (ilr.method.fast) {
      for (i in 1:ncol(ilrX)) {
        ilrX[, i] = sqrt((D - i)/(D - i + 1)) * log(((apply(as.matrix(x[,
                                                                        (i + 1):D, drop = FALSE]), 1, prod) + offset)^(1/(D -
                                                                                                                            i)))/(x[, i] + offset))
      }
    } else {
      for (i in 1:ncol(ilrX)) {
        ilrX[, i] = sqrt((D - i)/(D - i + 1)) * log(apply(as.matrix(x[,
                                                                      (i + 1):D]), 1, function(x) {
                                                                        exp(log(x))
                                                                      })/(x[, i] + offset) + offset)
      }
    }
    out <- as.matrix(ilrX)
  }

  return(out)
}


#' @title rbioseq_DE
#'
#' @description DE analysis function for RNA-seq data, with count filtering functionality.
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
#' @importFrom limma lmFit eBayes topTable contrasts.fit
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
