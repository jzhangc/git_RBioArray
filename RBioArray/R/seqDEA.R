.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Written by Jing Zhang, Ph.D. Please direct questions to jzhangcad@gmail.com.")
  suppressPackageStartupMessages(require(pathview))
  return(TRUE)
}

#' @title rbioseq_DE
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
rbioseq_DE <- function(AnnCountDfm, DEGNum = Inf, adjMethod = "fdr"){
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
