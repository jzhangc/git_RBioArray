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
