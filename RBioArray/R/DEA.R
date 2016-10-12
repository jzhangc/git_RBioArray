.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Written by Jing Zhang, Ph.D. Please direct questions to jzhangcad@gmail.com.")
  return(TRUE)
}

#' @title rbioSeq_DE
#'
#' @description The function below is an all-in-one solution to get DGE list and itrequires limma and edgeR packages for RNA-seq dataset.
#' @param AnnCountDfm Input data frame with annoatation merged.
#' @param DEGNum Number of targets to display. Default is \code{Inf}.
#' @param adjMethod P value correct methods. values: \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}.
#' @return Outputs a \code{dataframe} object DE resutls.
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma plotMDS voomWithQualityWeights lmFit eBayes plotSA plotMA topTable
#' @examples
#' \dontrun{
#' DE_dataframe <- rbioArray_DE(dataframe)
#' }
#' @export
rbioSeq_DE <- function(AnnCountDfm, DEGNum = Inf, adjMethod = "fdr"){
  # create the design matrix for Voom normalization
  Exp <- factor(SampleIndex$Conditions, levels= c("Pre", "Post"))
  design <- model.matrix( ~ Exp)
  # create DGE object using edgeR
  DGE_RNA <- DGEList(counts = AnnCountDfm[, 5:12], genes = AnnCountDfm[, 1:4])
  # for data Voom normalization
  DGE_RNA <- calcNormFactors(DGE_RNA)


  plotMDS(DGE_RNA) # MDS distribution before normalization
  VmWts_RNA <- voomWithQualityWeights(DGE_RNA, design, plot=TRUE,
                                      normalization = "quantile") # Voom normalization with quality weights

  Fit_RNA <- lmFit(VmWts_RNA, design) # linear fitting
  Bayes_RNA <- eBayes(Fit_RNA)

  plotSA(Bayes_RNA, xlab = "Average log-expression", ylab = "log2(sigma)",
         zero.weights=FALSE, pch=16, cex=0.2) # Simga vs Average plot
  limma::plotMA(Bayes_RNA) # log-ratios vs mean average plot
  DEGs_RNA <- topTable(Bayes_RNA, coef = "ExpPost", number = DEGNum, sort.by = "p",
                       adjust.method = adjMethod)
  # output
  return(DEGs_RNA)
}

#' @title rbioHtmap_unsuperv
#'
#' @description unsupvised heatmap: with or without voom normalization
#' @param AnnCountDfm Input data frame with annoatation merged.
#' @param voomNrm If the data needs to be voom normalized. Default is \code{TRUE}.
#' @param RNAtype RNAtype takes \code{"mRNA"}, \code{"miRNA"} and \code{"lncRNA"}. Default is \code{"mRNA"}.
#' @param nPlot sets the number of the features to be plotted in the heatmap, and takes \code{"all"} or a numeric value. Default is \code{"all"}.
#' @param cexRow Font size. Default is \code{0.6}.
#' @param cexCol Also font size. Default is \code{0.7}.
#' @return Outputs a heatmap based on the input dataframe.
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma voomWithQualityWeights
#' @importFrom gplots heatmap.2
#' @importFrom grid grid.newpage
#' @examples
#' \dontrun{
#' DE_dataframe <- rbioArray_DE(dataframe)
#' }
#' @export
####  ####
# make sure to load the library: gplots, grid, limma, edgeR
# cexRow and cexCol are the same functions in heatmap.2() that specify the font size
# rawDFM is the count data.frame for ALL the RNA types
#
# nPlot
rbioHtmap_unsuperv <- function(rawDFM, voomNrm = TRUE, RNAtype = "mRNA", nPlot = "all",
                               cexRow = 0.6, cexCol = 0.7){

  # normalization
  if (VoomNrm == TRUE){
    # create the design matrix for Voom normalization
    Exp <- factor(SampleIndex$Conditions, levels = c("Pre", "Post"))
    design <- model.matrix( ~ Exp)
    # create DGE object using edgeR
    DGE_RNA <- DGEList(counts = rawDFM[, 5:12], genes = rawDFM[, 1:4])
    # for data Voom normalization
    DGE_RNA <- calcNormFactors(DGE_RNA)
    VNrm <- voomWithQualityWeights(DGE_RNA, design, plot=TRUE,
                                   normalization = "none") # Voom normalization with quality weights

    if (RNAtype == "mRNA"){
      deDFM <- VNrm[VNrm$genes$gene_type == "protein_coding", ]$E # extract the mRNA portion
    } else if (RNAtype == "miRNA"){
      deDFM <- VNrm[VNrm$genes$gene_type == "miRNA", ]$E # extract the miRNA portion
    } else if (RNAtype == "lncRNA"){
      deDFM <- VNrm[VNrm$genes$gene_type == "lincRNA"
                    | VNrm$genes$gene_type == "antisense"
                    | VNrm$genes$gene_type == "processed_transcript"
                    | VNrm$genes$gene_type == "pseudogene", ]$E # extract the lncRNA portion
    } else {
      return("ERROR: Wrong RNA type")
    }

  } else {

    deDFM <- rawDFM # no normalization

  }

  if (nPlot == "all"){
    deDFM <- deDFM
  } else {
    deDFM <- deDFM[1:nPlot, ]
  }


  # the background dataframe. n = 10000 is to include all the entries. adjust accordingly should the number exceeds.
  htdfm <- apply(deDFM, c(1,2), function(x)(x + 1))

  grid.newpage()
  heatmap.2(as.matrix(htdfm), scale = "none", col = greenred(100),
            key = TRUE, symkey = FALSE, keysize = 1,
            density.info = "none", trace = "none",
            cexRow = cexRow, cexCol = cexCol)
}
