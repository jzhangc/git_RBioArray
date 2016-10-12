#' @title rbioArray_entrez2gene
#'
#' @description Add Human entrez ID to the DE dataframe
#' @param DEGdfm Input DE data frame.
#' @param cat Catergory, default is \code{"SYMBOL"}.
#' @param species Set the species, default is \code{"Hs"}.
#' @param pkg Name for the annotation package. The current package includes the human version.
#' @return Outputs a \code{dataframe} object with Entrez ID.
#' @importFrom pathview id2eg
#' @import org.Hs.eg.db
#' @examples
#' \dontrun{
#' DE_dataframe <- rbio_entrez2gene(dataframe)
#' }
#' @export
rbio_entrez2geneStats <- function(DEGdfm, cat = "SYMBOL", species = "Hs", pkg = "org.Hs.eg.db"){
  EnsemblVector <- DEGdfm$gene_name # create a vector conatining all the log-fold changes
  names(EnsemblVector) <- DEGdfm$gene_name # add names to each data point in the vector
  RNA_id <- id2eg(ids = names(EnsemblVector), category = cat, org = species, pkg.name = pkg)
  DEGdfm$SYMBOL <- DEGdfm$gene_name
  DEGdfm <- merge(DEGdfm, RNA_id, by = "SYMBOL")
  DEGdfm <- DEGdfm[complete.cases(DEGdfm), ]
  DEGdfm <- unique(DEGdfm)
  return(DEGdfm)
}



#' @title rbioArray_allGSA
#'
#' @description Add Human entrez ID to the DE dataframe
#' @param MyGS Pre-loaded gene set objects.
#' @param pVar Gene level p values. Could be, but not exclusive to, a variable of a dataframe.
#' @param logFCVar Gene level logFC (log fold change). Could be, but not exclusive to, a variable of a dataframe.
#' @param tVar Gene leve t values. Could be, but not exclusive to, a variable of a dataframe.#'
#' @return Outputs a \code{list} object with GSA results, both for p value based and t value based.
#' @details The function is based on piano package. It runs "fisher", "stouffer", "reporter", "tailStrength", "wilcoxon" for p value based GSA, and "page", "gsea", "maxmean" for t value based GSA.
#' @importFrom piano runGSA
#' @examples
#' \dontrun{
#' GOterm_bp <- loadGSC(file = "c5.bp.v5.0.entrez.gmt", type = "gmt") # load GO term biological process set
#'
#' pc_Kegg <- rbioArray_allGSA(GS_Kegg, pVar = pcGSdfm$p_value, logFCVar = pcGSdfm$logFC, tVar = pcGSdfm$t_value)
#'
#' }
#' @export
rbioArray_allGSA <- function(MyGS, pVar, logFCVar, tVar){

  gStats <- list(p_value = pVar,
                 logFC = logFCVar,
                 t_value = tVar)
  gStats <- lapply(gStats, function(x){names(x) <- get(paste(RNAtype, "_DE4GSA_flt", sep=""))$ENTREZID; x})

  GS_list_p <- list()
  GSigM_p <- c("fisher", "stouffer", "reporter", "tailStrength", "wilcoxon")
  for (i in 1:length(GSigM_p)){
    GS_S.i <- runGSA(gStats$p_value, gStats$logFC, geneSetStat = GSigM_p[i], signifMethod="geneSampling",
                     adjMethod= "fdr", gsc = MyGS, nPerm = 1000)
    GS_list_p[[i]] <- GS_S.i
  }
  names(GS_list_p) <- c("fisher", "stouffer", "reporter", "tailStrength", "wilcoxon")

  GS_list_t <- list()
  GSigM_t <- c("page", "gsea", "maxmean")
  for (j in 1:length(GSigM_t)){
    GS_S.j <- runGSA(gStats$t_value, geneSetStat = GSigM_t[j], signifMethod = "geneSampling",
                     adjMethod = "fdr", gsc = MyGS, nPerm = 1000)
    GS_list_t[[j]] <- GS_S.j
  }
  names(GS_list_t)<-c("page","gsea","maxmean")

  fullGS_list <- c(GS_list_p, GS_list_t)

  return(fullGS_list)
}
