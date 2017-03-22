#' @title rbioGS_sp2hsaEntrez
#'
#' @description Using up-to-date ensembl database, convert from mouse/rat ensemble transcript ID to and add Human entrez ID to the DE list from the DE functions, i.e. \code{\link{rbioarray_DE}} or \code{\link{rbioseq_DE}}, as human EntrezID is needed for GS analysis if using human gene sets.
#' @param DElst The list with DE reuslt, from functions \code{\link{rbioarray_DE}} or \code{\link{rbioseq_DE}}.
#' @param tgtSpecies The target species. Options are \code{"mmu"} and \code{"rno"}.
#' @param ensemblTransVar The name of the variable from DE list containing ensembl transcript ID.
#' @param parallelComputing If to use parallel computing. The cluster mode is \code{PSOCK} for now. Default is \code{FALSE}.
#' @details IMPORTANT: this function requires an internet connection as it retrieves information from ensembl website for human gene orthorlogs.
#' @return Outputs a DE \code{list} object with human Entrez ID for each dataframe. This list has the exact same format as the input DE list.
#' @import doParallel
#' @import foreach
#' @importFrom biomaRt useMart getBM
#' @examples
#' \dontrun{
#' rbioGS_sp2hsaEntrez(DElst = comparison_DE, tgtSpecies = "mmu", ensemblTransVar = "EnsemblID")
#' }
#' @export
rbioGS_sp2hsaEntrez <- function(DElst, tgtSpecies = "mmu", ensemblTransVar = NULL,
                                parallelComputing = FALSE){

  ## check the arguments
  if (is.null(ensemblTransVar)){
    stop("Please set the variable name for the ensembl transcript ID.")
  }

  ## prepare reference hsa entrezID
  # set the target species
  if (tgtSpecies == "mmu"){
    sp <- "mmusculus"
  } else if (tgtSpecies == "rno"){
    sp <- "rnorvegicus"
  }

  # extract hsa orthorlogy information
  sp_ensembl <- useMart("ensembl", dataset = paste0(sp, "_gene_ensembl"))
  attr <- c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene", "ensembl_transcript_id")
  sp_hsa_orth <- getBM(attr, filters = "with_hsapiens_homolog", values = TRUE,
                       mart = sp_ensembl)
  names(sp_hsa_orth)[names(sp_hsa_orth) == "ensembl_transcript_id"] <- paste0(tgtSpecies, "_ensembl_transcript_id") # generalized term for change column names

  # extract hsa entrezgene ID
  hsa_ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # establish the human set
  attr_hsa <- c("ensembl_gene_id", "entrezgene")
  hsa_entrez <- getBM(attr_hsa, filters = "", values = TRUE,
                      mart = hsa_ensembl)

  # merge the two dataframes
  sp_hsa_orth_entrez <- merge(sp_hsa_orth, hsa_entrez,
                              by.x = "hsapiens_homolog_ensembl_gene", by.y = "ensembl_gene_id",
                              all.x = TRUE)
  names(sp_hsa_orth_entrez)[names(sp_hsa_orth_entrez) == "entrezgene"] <- "hsa_entrezgene"
  sp_hsa_orth_entrez <- sp_hsa_orth_entrez[!duplicated(sp_hsa_orth_entrez$mmu_ensembl_transcript_id), ]

  ## add the hsa entrez ID to the non-hsa DElist
  # temp func for adding the variable, i is the DE dataframe
  tmpfunc <- function(i){
    j <- merge(i, sp_hsa_orth_entrez,
               by.x = ensemblTransVar, by.y = paste0(sp, "_ensembl_transcript_id"), all.x = TRUE)
    return(j)
  }

  # looping through the DElist
  # vectorize the output list
  out <- vector(mode = "list", length(names(DElst)))
  names(out) <- names(DElst)

  # looping
  if (!parallelComputing){

    out[] <- lapply(DElst, function(i)tmpfunc(i))

  } else { # parallel computing

    # set up clusters for PSOCK
    cl <- makeCluster(n_cores, type = cluster, outfile = "")
    registerDoParallel(cl) # part of doParallel package
    on.exit(stopCluster(cl)) # close connect when exiting the function

    # computing
    out[] <- foreach(i = DElist) %dopar% {
      tmpout <- tmpfunc(i)
    }

  }

  ## output
  assign(paste(deparse(substitute(DElist)), "_hsaEntrez",sep = ""), out, envir = .GlobalEnv)

  ## message
  message(cat("Human entrez ID has been added as variable \"hsa_entrezgene\". "))

}


#' @title rbioGS
#'
#' @description Add Human entrez ID to the DE dataframe
#' @param GS Pre-loaded gene set objects.
#' @param pVar Gene level p values. Could be, but not exclusive to, a variable of a dataframe. Must be the same length as \code{logFCVar}, \code{tVar} and \code{idVar}.
#' @param logFCVar Gene level logFC (log fold change). Could be, but not exclusive to, a variable of a dataframe. Must be the same length as \code{pVar}, \code{tVar} and \code{idVar}.
#' @param tVar Gene leve t values. Could be, but not exclusive to, a variable of a dataframe. Must be the same length as \code{pVar}, \code{logFCVar} and \code{idVar}.
#' @param idVar Gene IDs. Could be, but not exclusive to, a variable of a dataframe. Must be the same length as \code{pVar}, \code{logFCVar} and \code{tVar}. Currently only takes \code{Entrez ID}.
#' @param method_p Gene set ernichment methods that takes \code{p value} and \code{logFC}. Default is \code{c("fisher", "stouffer", "reporter", "tailStrength", "wilcoxon")}.
#' @param method_t Gene set ernichment methods that takes \code{t statistics}. Default is \code{c("page", "gsea", "maxmean")}.
#' @param ... Arguments to pass to \code{runGSA} function from \code{piano} pacakge. See the corresponding help page from of \code{piano} for details.
#' @param parallelComputing If to use parallel computing or not. Default is \code{FALSE}
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @details The function is based on piano package. It runs "fisher", "stouffer", "reporter", "tailStrength", "wilcoxon" for p value based GSA, and "page", "gsea", "maxmean" for t value based GSA.
#' @return Outputs a \code{list} object with GSA results, both for p value based and t value based.
#' @importFrom piano runGSA
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @examples
#' \dontrun{
#' gsoutput <- rbioGS(GS = kegg, pVar = dfm$p_value, logFCVar = dfm$logFC, tVar = dfm$t_value, idVar = dfm$EntrezID, nPerm = 1000, parallelComputing = TRUE, clusterType = "FORK")
#'
#'
#' }
#' @export
rbioGS <- function(GS, pVar, logFCVar, tVar, idVar,
                   method_p = c("fisher", "stouffer", "reporter", "tailStrength", "wilcoxon"),
                   method_t = c("page", "gsea", "maxmean"), ...,
                   parallelComputing = FALSE, clusterType = "PSOCK"){

  gStats <- list(p_value = pVar,
                 logFC = logFCVar,
                 t_value = tVar)
  gStats <- lapply(gStats, function(x){names(x) <- idVar; x})

  GSigM_p <- method_p
  GSigM_t <- method_t

  ## make empty output lists
  GS_list_p <- vector(mode = "list", length = length(GSigM_p))
  names(GS_list_p) <- GSigM_p

  GS_list_t <- vector(mode = "list", length = length(GSigM_t))
  names(GS_list_t) <- GSigM_t

  ## make tmp GS functions for parallel computing, as well as the
  tmpfunc_p <- function(i, GSmethod_p, ...){
    p <- gStats$p_value
    logfc <- gStats$logFC
    GS_S.i <- runGSA(p, logfc, geneSetStat = GSmethod_p[i], ...)
  }

  tmpfunc_t <- function(i, GSmethod_t, ...){
    t <- gStats$t_value
    GS_S.i <- runGSA(t, geneSetStat = GSmethod_t[i], ...)
  }

  if (!parallelComputing){

    for (m in 1:length(GSigM_p)){
      GS_list_p[[m]] <- tmpfunc_p(m, GSmethod_p = GSigM_p, gsc = GS, ...)
    }

    for (n in 1:length(GSigM_t)){
      GS_list_t[[n]] <- tmpfunc_t(n, GSmethod_t = GSigM_t, gsc = GS, ...)
    }

  } else {

    ## parallel computing
    # set up cpu core number
    n_cores <- detectCores() - 1

    # parallel computing
    if (clusterType == "FORK"){ # mac and linux only

      GS_list_p[] <- mclapply(1: length(GSigM_p), FUN = tmpfunc_p, GSmethod_p = GSigM_p, gsc = GS, ..., mc.cores = n_cores, mc.preschedule = FALSE)
      GS_list_t[] <- mclapply(1: length(GSigM_t), FUN = tmpfunc_t, GSmethod_t = GSigM_t, gsc = GS, ..., mc.cores = n_cores, mc.preschedule = FALSE)

    } else { # windows etc

      # set up cpu cluster for PSOCK
      cl <- makeCluster(n_cores, type = clusterType, outfile = "")
      registerDoParallel(cl) # part of doParallel package
      on.exit(stopCluster(cl)) # close connect when exiting the function

      # run functioins
      GS_list_p[] <- foreach(i = 1: length(GSigM_p), .packages = "piano") %dopar% {
        out <- tmpfunc_p(i, GSmethod_p = GSigM_p, gsc = GS, ...)
      }

      GS_list_t[] <- foreach(i = 1: length(GSigM_t), .packages = "piano") %dopar% {
        out <- tmpfunc_t(i, GSmethod_t =  GSigM_t, gsc = GS, ...)
      }

    }

  }

  fullGS_list <- c(GS_list_p, GS_list_t)

  return(fullGS_list)

}


#' @title rbioGS_boxplot
#'
#' @description Generate boxplot from piano GS rank object,
#' @param GSA_list piano GS results object.
#' @param fileName Output file name. Default is \code{"GS_list"}.
#' @param KEGG if \code{TRUE}, the function will remove the "KEGG_" string in the GS name variable. Default is \code{FALSE}.
#' @param pClass P class for the consensus score. Options are "distinct", "mixed" or "non". Default is \code{NULL}.
#' @param classDirection The directionality of the p class, if pClass is set to \code{"distinct"} or \code{"mixed"}. Options are \code{"up"} and \code{"down"}. Default is \code{NULL}.
#' @param ... Arguments passing to \code{consensusScores} function from \code{piano} package. See the corresponding \code{piano} package help page for details.
#' @param plotTitle Title of the plot. Default is \code{NULL}.
#' @param xLabel X-axis label. Default is \code{"rank"}.
#' @param yLabel Y-axis label. Default is \code{NULL}.
#' @param plotWidth Set the width of the plot. Default is \code{170}.
#' @param plotHeight Set the height of the plot. Default is \code{150}.
#' @details The function is a wrapper that takes resulted object from \code{runGSA} function from \code{piano} package.
#' @return Outputs a \code{pdf} boxplot figure file with allGSA results.
#' @importFrom reshape2 melt
#' @importFrom grid grid.newpage grid.draw
#' @import ggplot2
#' @importFrom piano consensusScores
#' @examples
#' \dontrun{
#'
#' GS_boxplot(GSA_list = experi_preVpost_GS, pClass = "non", n = 20, adjusted = TRUE, method = "median", KEGG = TRUE,fileName = "GS_analysis", plotWidth = 260, plotHeight = 240)
#'
#' }
#' @export
rbioGS_boxplot <- function(GSA_list, fileName = "GS_list", KEGG = FALSE, pClass = NULL, classDirection = NULL, ...,
                           plotTitle = NULL, xLabel = "rank", yLabel = NULL, plotWidth = 170, plotHeight = 150){

  # check the arguments
  if (is.null(pClass)) {
    stop(cat("Please set p class(es): \"distinct\", \"mixed\" or \"non\" "))
  }

  if (is.null(classDirection) & pClass != "non"){
    stop(cat("Please set direction: \"up\" or \"down\""))
  }

  # prepare consensus score list
  GSrank <- consensusScores(resList = GSA_list, plot = FALSE, class = pClass, direction = classDirection, ...)


  # prepare the dataframe for ggplot2
  DFM <- data.frame(GSrank$rankMat)
  DFM <- DFM[, c(-1, -2)]
  DFM <- data.frame(t(DFM)) # transpose the dataframe using t()
  DFM$enrichment <- rownames(DFM)
  DFM_mlt <- melt(DFM, id.vars = colnames(DFM)[length(colnames(DFM))])

  if (KEGG){ # remove "KEGG_" and the space between the terms
    DFM_mlt$variable <- sapply(DFM_mlt$variable, function(x){
      i <- substring(x, 6)
      gsub("_", " ", i)
    })
  } else {
    DFM_mlt$variable <- sapply(DFM_mlt$variable, function(x)gsub("_", " ", x))
  }

  DFM_mlt$variable <- factor(DFM_mlt$variable, levels = c(unique(DFM_mlt$variable)))

  # plot
  grid.newpage()
  plt <- ggplot(DFM_mlt, aes(x = variable, y = value)) +
    geom_boxplot() +
    guides(fill = FALSE) +
    scale_x_discrete(limits = with(DFM_mlt, rev(levels(variable)))) +
    ggtitle(plotTitle) +
    coord_flip() + # flip the axes
    xlab(yLabel) + # reverse the arguments because of the flipping
    ylab(xLabel) + # reverse the arguments because of the flipping
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          axis.title = element_text(face = "bold"),
          legend.position = "bottom")

  # export the file and draw a preview
  if (is.null(classDirection)){
    ggsave(filename = paste(fileName,"_boxplot_", pClass, ".pdf",sep = ""), plot = plt,
           width = plotWidth, height = plotHeight, units = "mm",dpi = 600)
  } else {
    ggsave(filename = paste(fileName,"_boxplot_", pClass, "_", classDirection , ".pdf",sep = ""), plot = plt,
           width = plotWidth, height = plotHeight, units = "mm",dpi = 600)
  }

  grid.draw(plt) # preview
}



#' @title rbioGS_scatter
#'
#' @description Generate scatter plot for piano GS rank heatmap obejct.
#' @param GSA_list GSA list generated from \code{\link{rbioArray_allGSA}}.
#' @param fileName Output file name.
#' @param ... Arguments passing to \code{consensusHeatmap} function from \code{piano} package. See the responding help page of \code{piano} for details.
#' @param rankCutoff Cutoff value for GS rank line.
#' @param pCutoff Cutoff value for GS p value line.
#' @param plotTitle Title of the plot. Default is \code{NULL}.
#' @param xLabel X-axis label. Default is \code{"median p value"}.
#' @param yLabel Y-axis label. Default is \code{"consensus score"}.
#' @param plotHeight Set the height of the plot. Default is \code{150}.
#' @return Outputs a \code{pdf} scatter figure file with allGSA results.
#' @details The function is based on piano package.
#' @importFrom reshape2 melt
#' @importFrom grid grid.newpage grid.draw
#' @importFrom piano consensusHeatmap
#' @import ggplot2
#' @examples
#' \dontrun{
#'
#' rbioGS_scatter(GSA_list = GS_Pos, fileName = "GS_pos", cutoff = 15, method = "median", adjusted = TRUE, rankCutoff = 50, pCutoff = 0.05)
#'
#' }
#' @export
rbioGS_scatter <- function(GSA_list, fileName = "GS_list",
                           ...,
                           plotTitle = NULL, xLabel = "median p value", yLabel = "consensus score",
                           rankCutoff = 20, pCutoff = 0.05,
                           plotWidth = 170, plotHeight = 150){

  HTmap <-consensusHeatmap(GSA_list, plot = FALSE, ...)

  ## data frame prep
  p_class_name_rank <- c("disdn_rank", "mixdn_rank", "nondir_rank", "mixup_rank",
                         "disup_rank")
  p_class_name_p<-c("disdn_p", "mixdn_p", "nondir_p", "mixup_p",
                    "disup_p")

  rank_tmp <- data.frame(HTmap$rankMat)
  p_value_tmp <- data.frame(HTmap$pMat)
  names(rank_tmp) <- p_class_name_rank
  names(p_value_tmp) <- p_class_name_p
  rank_tmp$GS <- rownames(rank_tmp)
  p_value_tmp$GS <- rownames(p_value_tmp)
  rp_tmp <- merge(rank_tmp, p_value_tmp, by="GS",sort=FALSE)

  # subset the total dfm4plot by p value class and merge all the sub dataframes into a new dataframe
  # with the new variables of "GS", "rank", "p_value", "p_value_class"
  p_class <- unlist(unique(lapply(strsplit(names(rp_tmp)[2:11], split="_", fixed=TRUE),
                                  function(x)x[1]))) # get the p class name from the character string from the column names

  for (n in 1:length(p_class)){
    assign(paste("dfm4plot_", p_class[n], sep = ""), data.frame(rp_tmp[, grep(paste(p_class[n], "+", sep=""),
                                                                              names(rp_tmp))],
                                                                GS = rp_tmp$GS, p_class = p_class[n]))
  } # create temp dataframes for all the p classes

  for (n in 1:length(p_class)){
    tmpDfm <- get(paste("dfm4plot_", p_class[n], sep = ""))
    names(tmpDfm)[1:2] <- c("rank", "p_value")
    assign(paste("dfm4plot_", p_class[n], sep = ""), tmpDfm)
    rm("tmpDfm")
  } # change variable names for all the temp data frames

  dfm4plot <- data.frame()

  for (n in 1:length(p_class)){
    dfm4plot <- rbind(dfm4plot, get(paste("dfm4plot_", p_class[n], sep = "")))
  } # merge all the temp data frames together


  ## plotting
  grid.newpage()
  ScatterP<-ggplot(dfm4plot, aes(x = p_value, y = rank)) +
    geom_point(aes(shape = factor(p_class)), size = 3) +
    ggtitle(plotTitle) +
    xlab(xLabel) +
    ylab(yLabel) +
    scale_x_continuous(breaks = c(0.1, 0.05, 0),
                       labels = c("0.1", "0.05", "0"),
                       expand = c(0, 0),
                       trans = "reverse", limits = c(0.1, NA)) +
    scale_y_continuous(breaks = c(100, 50, 0),
                       labels = c("100", "50", "0"),
                       expand = c(0, 0),
                       trans = "reverse", limits = c(100, NA)) +
    geom_vline(xintercept = pCutoff, linetype = "dashed") +
    geom_hline(yintercept = rankCutoff, linetype = "dashed") +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_rect(),
          legend.key = element_blank()) +
    scale_shape_manual(values = c(1:5))

  # export the files and draw a preview
  write.csv(dfm4plot, file = paste(fileName,".scatterplot.csv",sep = ""))
  ggsave(filename = paste(fileName,".scatterplot.pdf",sep = ""), plot = ScatterP,
         width = plotWidth, height = plotHeight, units = "mm",dpi = 600)
  grid.draw(ScatterP) # preview
}

#' @title rbioGS_kegg
#'
#' @description Download and generate DE results masked kegg pathway figures
#' @param dfm GS dataframe with \code{ENTREZID} and \code{logFC} variables.
#' @param entrezVar Name of the EntrezID variable in the \code{DElst} object.
#' @param keggID Make sure to have quotation marks around the ID number.
#' @param suffix Output file name suffix. Make sure to put it in quotation marks.
#' @param species Set the species. Default is \code{"hsa"}. Visit kegg website for details.
#' @details GSrank takes the resulted object generated from the piano function consensusScores().
#' @return Outputs a \code{list} kegg object, as well as masked figure files in \code{pdf} format.
#' @details Visit website \url{http://www.genome.jp/kegg/pathway.html} for kegg IDs.
#' @importFrom pathview pathview
#' @examples
#' \dontrun{
#'
#'  rbioGS_kegg(all_DE4GSA_flt, keggID = "04060") # Cytokine-cytokine receptor interaction, disdn + mixdn + nd
#'
#' }
#' @export
rbioGS_kegg <- function(dfm, entrezVar = NULL,
                       keggID, suffix, species = "hsa"){

  # check entrez ID variable name
  if (is.null(entrezVar)){
    stop("Please provide the name for the Entrez ID variable from the input dataframe")
  }

  # prepare objects
  logFC <- dfm$logFC
  names(logFC) <- dfm[, entrezVar]

  # visualize
  KEGG <- pathview(gene.data = logFC, pathway.id = keggID, species = species,
                   out.suffix = suffix, keys.align = "y", kegg.native = TRUE, match.data = FALSE,
                   key.pos = "topright")

  # set the .GlobalEnv to the envir argument so that the assign function will assign the value to a global object, aka outside the function
  return(assign(paste("kegg_", keggID, "_" , deparse(substitute(dfm)), sep = ""), KEGG, envir = .GlobalEnv))
}



#' @title rbioGS_all
#'
#' @description All-in-one wrapper for GSA and plotting.
#' @param objTitle Object title for the output GS analysis list from \code{piano} package.
#' @param DElst The input list with DE reuslt, from functions \code{\link{rbioarray_DE}} or \code{\link{rbioseq_DE}}.
#' @param entrezVar Name of the EntrezID variable in the \code{DElst} object.
#' @param GS Pre-loaded gene set objects.
#' @param ... Arguments to pass to \code{\link{rbioGS}}.
#' @param parallelComputing If to use parallel computing or not. Default is \code{FALSE}
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param boxplot If to plot boxplots. Default is \code{TRUE}.
#' @param boxplotKEGG When \code{boxplot = TRUE}, to set if the gene set is KEGG. Default is \code{FALSE}.
#' @param boxplotN When \code{boxplot = TRUE}, to set the \code{n} (rank cutoff) argument passed to \code{consensusScores} function from \code{piano} package. Default is \code{20}.
#' @param boxplotTitle When \code{boxplot = TRUE}, to set the title of the boxplots. Default is \code{NULL}.
#' @param boxplotXlabel When \code{boxplot = TRUE}, to set the boxplots x-axis label. Default is \code{"rank"}.
#' @param boxplotYlabel When \code{boxplot = TRUE}, to set the boxplots y-axis label. Default is \code{NULL}.
#' @param boxplotWidth When \code{boxplot = TRUE}, to set the boxplots width. Default is \code{170}.
#' @param boxplotHeight When \code{boxplot = TRUE}, to set the boxplots height. Default is \code{150}.
#' @param scatterplot When \code{scatterplot = TRUE}, to set if to plot a scatter plot. Default is \code{TRUE}.
#' @param scatterplotCutoff When \code{scatterplot = TRUE}, to set the rank cutoff of the scatter plot. Default is \code{20}.
#' @param scatterplotRankline When \code{scatterplot = TRUE}, to set the rank line on the scatter plot. Default is \code{20}.
#' @param scatterplotPline When \code{scatterplot = TRUE}, to set the p value line on the scatter plot. Default is \code{0.05}.
#' @param scatterTitle When \code{scatterplot = TRUE}, to set the title of the scatter plot. Default is \code{NULL}.
#' @param scatterXlabel When \code{scatterplot = TRUE}, to set the scatter plot x-axis label. Default is \code{"median p value"}.
#' @param scatterYlabel When \code{scatterplot = TRUE}, to set the scatter plot y-axis label. Default is \code{"consensus score"}.
#' @param scatterWidth When \code{scatterplot = TRUE}, to set the scatter plot width. Default is \code{170}.
#' @param scatterHeight When \code{scatterplot = TRUE}, to set the scatter plot height. Default is \code{150}.
#' @param plotMethod When \code{boxplot = TRUE} and/or \code{scatterplot = TRUE}, to set the p methods. Options are \code{"median"} and \code{"mean"}. Default is \code{"median"}.
#' @param plotPadjust When \code{boxplot = TRUE} and/or \code{scatterplot = TRUE}, to set if to use FDR adjusted p value or not. Default is \code{TRUE}.
#' @param plotGSname When \code{boxplot = TRUE} and/or \code{scatterplot = TRUE}, to set the GS name in the file name. Default is \code{"GS"}.
#' @details This is an all-in-one function for GS anlyasis based on piano package. It runs "fisher", "stouffer", "reporter", "tailStrength", "wilcoxon" for p value based GSA, and "page", "gsea", "maxmean" for t value based GSA (customizable). See arguments for \code{\link{rbioGS}} for details.
#' @return Outputs  \code{csv} files and \code{pdf} figure files with GSA results.
#' @import doParallel
#' @import foreach
#' @importFrom piano runGSA
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @examples
#' \dontrun{
#'
#' rbioGS_all(objTitle = "mydata", DElst = DElist, entrezVar = "EntrezID", method_p = c("stouffer", "fisher"), method_t = c("gsea", "maxmean"), nPerm = 1000, GS = kegg, parallelComputing = TRUE, clusterType = "FORK")
#'
#' }
#' @export
rbioGS_all <- function(objTitle = "DE", DElst, entrezVar = NULL,
                       GS = NULL, ...,
                       parallelComputing = FALSE, clusterType = "PSOCK",
                       boxplot = TRUE,
                       boxplotKEGG = FALSE, boxplotN = 20,
                       boxplotTitle = NULL, boxplotXlabel = "rank", boxplotYlabel = NULL,
                       boxplotWidth = 170, boxplotHeight = 150,
                       scatterplot = TRUE,
                       scatterplotCutoff = 20,
                       scatterRankline = 20, scatterPline = 0.05,
                       scatterTitle = NULL, scatterXlabel = "median p value", scatterYlabel = "consensus score",
                       scatterWidth = 170, scatterHeight = 150,
                       plotMethod = "median", plotPadjust = TRUE, plotGSname = "GS"){

  if (is.null(entrezVar)){
    stop("please set the name of the Entrez ID vaiable")
  }

  if (is.null(GS)){
    stop("please choose gene set(s)")
  }

  if (boxplot){
    # set up a temp function for boxplot with directinality
    tmpfunc <- function(x){

      pCl <- c("distinct", "mixed") # for boxplot
      classDirt <- c("up", "down") # for boxplot

      lapply(pCl, function(m)lapply(classDirt,
                                    function(n)RBioArray::rbioGS_boxplot(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""),
                                                                         KEGG = boxplotKEGG, pClass = m, classDirection = n, adjust = plotPadjust,
                                                                         n = boxplotN, xLabel = boxplotXlabel, yLabel = boxplotYlabel,
                                                                         plotTitle = boxplotTitle, plotWidth = boxplotWidth, plotHeight = boxplotHeight)))

    }
  }


  ## make an empty list to store the GS results
  GSlst <- vector(mode = "list", length(names(DElst)))
  names(GSlst) <- names(DElst)

  ## GSA
  if (!parallelComputing){
    # remove the rows with NA in the Entrez ID variable
    DElst <- lapply(DElst, function(x)x[complete.cases(x[, entrezVar]), ])

    # run GSA
    GSlst[] <- lapply(DElst, function(i)RBioArray::rbioGS(GS = GS, pVar = i$P.Value, logFCVar = i$logFC,
                                                          tVar = i$t, idVar = i[, entrezVar],
                                                          parallelComputing = parallelComputing, clusterType = clusterType, ...))


    if (boxplot){

      # boxplots
      lapply(1: length(GSlst), tmpfunc)

      lapply(1: length(GSlst), function(x)RBioArray::rbioGS_boxplot(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""),
                                                                    KEGG = boxplotKEGG, pClass = "non", adjust = plotPadjust,
                                                                    n = boxplotN, xLabel = boxplotXlabel, yLabel = boxplotYlabel,
                                                                    plotTitle = boxplotTitle, plotWidth = boxplotWidth, plotHeight = boxplotHeight))

    }

    if (scatterplot){

      # scatter plots
      lapply(1: length(GSlst), function(x)RBioArray::rbioGS_scatter(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""), cutoff = scatterplotCutoff,
                                                                    method = plotMethod, adjust = plotPadjust, rankCutoff = scatterRankline,
                                                                    pCutoff = scatterPline,
                                                                    plotTitle = scatterTitle, xLabel = scatterXlabel, yLabel = scatterYlabel,
                                                                    plotWidth = scatterWidth, plotHeight = scatterHeight))

    }

  } else { # parallel computing

    # set up cpu cores
    n_cores <- detectCores() - 1

    # parallel computing
    if (clusterType == "FORK"){ # mac and linux only
      # remove the rows with NA in the Entrez ID variable
      DElst <- mclapply(DElst, FUN = function(x)x[complete.cases(x[, entrezVar]), ], mc.cores = n_cores, mc.preschedule = FALSE)

      # run GSA
      GSlst[] <- mclapply(DElst, FUN = function(i)RBioArray::rbioGS(GS = GS, pVar = i$P.Value, logFCVar = i$logFC,
                                                                    tVar = i$t, idVar = i[, entrezVar],
                                                                    parallelComputing = FALSE, ...),
                          mc.cores = n_cores, mc.preschedule = FALSE)



      if (boxplot){

        # boxplots
        mclapply(1: length(GSlst), tmpfunc, mc.cores = n_cores, mc.preschedule = FALSE)

        mclapply(1: length(GSlst), FUN = function(x)RBioArray::rbioGS_boxplot(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""),
                                                                              KEGG = boxplotKEGG, pClass = "non", adjust = plotPadjust,
                                                                              n = boxplotN, xLabel = boxplotXlabel, yLabel = boxplotYlabel,
                                                                              plotTitle = boxplotTitle, plotWidth = boxplotWidth, plotHeight = boxplotHeight),
                 mc.cores = n_cores, mc.preschedule = FALSE)

      }

      if (scatterplot){

        mclapply(1: length(GSlst), FUN = function(x)RBioArray::rbioGS_scatter(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""),
                                                                              cutoff = scatterplotCutoff,
                                                                              method = plotMethod, adjust = plotPadjust, rankCutoff = scatterRankline,
                                                                              pCutoff = scatterPline,
                                                                              plotTitle = scatterTitle, xLabel = scatterXlabel, yLabel = scatterYlabel,
                                                                              plotWidth = scatterWidth, plotHeight = scatterHeight),
                 mc.cores = n_cores, mc.preschedule = FALSE)


      }

    } else { # windows etc

      # set up clusters for PSOCK
      cl <- makeCluster(n_cores, type = clusterType, outfile = "")
      registerDoParallel(cl) # part of doParallel package
      on.exit(stopCluster(cl)) # close connect when exiting the function

      # remove the rows with NA in the Entrez ID variable
      DElst <- foreach(x = DElst) %dopar% {
        out <- x[complete.cases(x[, entrezVar]), ]
      }

      GSlst[] <- foreach(i = DElst, .packages = c("RBioArray", "piano")) %dopar% {
        out <- rbioGS(GS = GS, pVar = i$P.Value, logFCVar = i$logFC,
                      tVar = i$t, idVar = i[, entrezVar], parallelComputing = FALSE, ...)
      }

      # plot
      if (boxplot){

        # boxplots
        foreach(x = 1: length(GSlst), .packages = c("RBioArray", "piano")) %dopar% {tmpfunc(x)}


        foreach(x = 1: length(GSlst), .packages = c("RBioArray", "piano")) %dopar% {
          RBioArray::rbioGS_boxplot(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""),
                                    KEGG = boxplotKEGG, pClass = "non", adjust = plotPadjust,
                                    n = boxplotN, xLabel = boxplotXlabel, yLabel = boxplotYlabel,
                                    plotTitle = boxplotTitle, plotWidth = boxplotWidth, plotHeight = boxplotHeight)
        }

      }

      if (scatterplot){
        foreach(x = 1: length(GSlst), .packages = c("RBioArray", "piano")) %dopar% {
          rbioGS_scatter(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""), cutoff = scatterplotCutoff,
                         method = plotMethod, adjust = plotPadjust, rankCutoff = scatterRankline,
                         pCutoff = scatterPline,
                         plotTitle = scatterTitle, xLabel = scatterXlabel, yLabel = scatterYlabel,
                         plotWidth = scatterWidth, plotHeight = scatterHeight)
        }
      }

    }

  }

  assign(paste(objTitle, "_GS_list_", deparse(substitute(GS)), sep = ""), GSlst, envir = .GlobalEnv)

}



#' @title rbioGS_all_noplot
#'
#' @description All-in-one wrapper for GSA (no plotting).
#' @param DElst The input list with DE reuslt, from functions \code{\link{rbioarray_DE}} or \code{\link{rbioseq_DE}}.
#' @param entrezVar Name of the EntrezID variable in the \code{DElst} object.
#' @param GS Pre-loaded gene set objects.
#' @param ... Arguments to pass to \code{\link{rbioGS}}.
#' @param parallelComputing If to use parallel computing or not. Default is \code{FALSE}
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @details This is an all-in-one function for GS anlyasis based on piano package. It runs "fisher", "stouffer", "reporter", "tailStrength", "wilcoxon" for p value based GSA, and "page", "gsea", "maxmean" for t value based GSA (customizable). See arguments for \code{\link{rbioGS}} for details.
#' @return Outputs  \code{csv} files and \code{pdf} figure files with GSA results.
#' @import doParallel
#' @import foreach
#' @importFrom piano runGSA
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @examples
#' \dontrun{
#'
#' rbioGS_all_noplot(DElst = DElist, entrezVar = "EntrezID", method_p = c("stouffer", "fisher"), method_t = c("gsea", "maxmean"), nPerm = 1000, GS = kegg, parallelComputing = TRUE, clusterType = "FORK")
#'
#' }
#' @export
rbioGS_all_noplot <- function(DElst, entrezVar = NULL,
                    GS = NULL, ...,
                    parallelComputing = FALSE, clusterType = "PSOCK"){

  if (is.null(entrezVar)){
    stop("please set the name of the Entrez ID vaiable")
  }

  if (is.null(GS)){
    stop("please choose gene set")
  }


  ## make an empty list to store the GS results
  GSlst <- vector(mode = "list", length(names(DElst)))
  names(GSlst) <- names(DElst)


  ## GSA
  if (!parallelComputing){
    # remove the rows with NA in the Entrez ID variable
    DElst <- lapply(DElst, function(x)x[complete.cases(x[, entrezVar]), ])

    # run GSA
    GSlst[] <- lapply(DElst, function(i)RBioArray::rbioGS(GS = GS, pVar = i$P.Value, logFCVar = i$logFC,
                                                          tVar = i$t, idVar = i[, entrezVar],
                                                          parallelComputing = parallelComputing, clusterType = clusterType, ...))


  } else { # parallel computing

    # set up cpu cores
    n_cores <- detectCores() - 1

    # parallel computing
    if (clusterType == "FORK"){ # mac and linux only
      # remove the rows with NA in the Entrez ID variable
      DElst <- mclapply(DElst, function(x)x[complete.cases(x[, entrezVar]), ], mc.cores = n_cores, mc.preschedule = FALSE)

      # run GSA
      GSlst[] <- mclapply(DElst, function(i)RBioArray::rbioGS(GS = GS, pVar = i$P.Value, logFCVar = i$logFC,
                                                              tVar = i$t, idVar = i[, entrezVar],
                                                              parallelComputing = FALSE, ...),
                          mc.cores = n_cores, mc.preschedule = FALSE)


    } else { # windows etc

      # set up clusters for PSOCK
      cl <- makeCluster(n_cores, type = clusterType, outfile = "")
      registerDoParallel(cl) # part of doParallel package
      on.exit(stopCluster(cl)) # close connect when exiting the function

      # remove the rows with NA in the Entrez ID variable
      DElst <- foreach(x = DElst) %dopar% {
        out <- x[complete.cases(x[, entrezVar]), ]
      }

      GSlst[] <- foreach(i = DElst, .packages = c("RBioArray", "piano")) %dopar% {
        out <- rbioGS(GS = GS, pVar = i$P.Value, logFCVar = i$logFC,
                      tVar = i$t, idVar = i[, entrezVar], parallelComputing = FALSE, ...)
      }

    }

  }

  return(GSlst)
}



#' @title rbioGS_plotting
#'
#' @description All-in-one wrapper for GSA plotting.
#' @param GSlst The input list with GS reuslt, i.e. a collection of dataframes from \code{topTable} function of \code{limma} package.
#' @param plotGSname When \code{boxplot = TRUE} and/or \code{scatterplot = TRUE}, to set the GS name in the file name. Default is \code{"GS"}.
#' @param boxplot If to plot boxplots. Default is \code{TRUE}.
#' @param boxplotKEGG When \code{boxplot = TRUE}, to set if the gene set is KEGG. Default is \code{FALSE}.
#' @param boxplotN When \code{boxplot = TRUE}, to set the \code{n} (rank cutoff) argument passed to \code{consensusScores} function from \code{piano} package. Default is \code{20}.
#' @param boxplotTitle When \code{boxplot = TRUE}, to set the title of the boxplots. Default is \code{NULL}.
#' @param boxplotXlabel When \code{boxplot = TRUE}, to set the boxplots x-axis label. Default is \code{"rank"}.
#' @param boxplotYlabel When \code{boxplot = TRUE}, to set the boxplots y-axis label. Default is \code{NULL}.
#' @param boxplotWidth When \code{boxplot = TRUE}, to set the boxplots width. Default is \code{170}.
#' @param boxplotHeight When \code{boxplot = TRUE}, to set the boxplots height. Default is \code{150}.
#' @param scatterplot When \code{scatterplot = TRUE}, to set if to plot a scatter plot. Default is \code{TRUE}.
#' @param scatterplotCutoff When \code{scatterplot = TRUE}, to set the rank cutoff of the scatter plot. Default is \code{20}.
#' @param scatterplotRankline When \code{scatterplot = TRUE}, to set the rank line on the scatter plot. Default is \code{20}.
#' @param scatterplotPline When \code{scatterplot = TRUE}, to set the p value line on the scatter plot. Default is \code{0.05}.
#' @param scatterTitle When \code{scatterplot = TRUE}, to set the title of the scatter plot. Default is \code{NULL}.
#' @param scatterXlabel When \code{scatterplot = TRUE}, to set the scatter plot x-axis label. Default is \code{"median p value"}.
#' @param scatterYlabel When \code{scatterplot = TRUE}, to set the scatter plot y-axis label. Default is \code{"consensus score"}.
#' @param scatterWidth When \code{scatterplot = TRUE}, to set the scatter plot width. Default is \code{170}.
#' @param scatterHeight When \code{scatterplot = TRUE}, to set the scatter plot height. Default is \code{150}.
#' @param plotMethod When \code{boxplot = TRUE} and/or \code{scatterplot = TRUE}, to set the p methods. Options are \code{"median"} and \code{"mean"}. Default is \code{"median"}.
#' @param plotPadjust When \code{boxplot = TRUE} and/or \code{scatterplot = TRUE}, to set if to use FDR adjusted p value or not. Default is \code{TRUE}.
#' @param parallelComputing If to use parallel computing or not. Default is \code{FALSE}
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @details The function takes the reuslted list from \code{\link{rbioGS_all_noplot}} function.
#' @return Outputs  \code{csv} files and \code{pdf} figure files, i.e. boxplots and scatter plot.
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @examples
#' \dontrun{
#'
#' rbioGS_plotting(GSlist = myGSlist, plotGSname = "kegg", parallelComputing = TRUE, clusterType = "FORK")
#'
#' }
#' @export
rbioGS_plotting <- function(GSlst, plotGSname = "GS",
                            boxplot = TRUE,
                            boxplotKEGG = FALSE, boxplotN = 20,
                            boxplotTitle = NULL, boxplotXlabel = "rank", boxplotYlabel = NULL,
                            boxplotWidth = 170, boxplotHeight = 150,
                            scatterplot = TRUE,
                            scatterplotCutoff = 20,
                            scatterRankline = 20, scatterPline = 0.05,
                            scatterTitle = NULL, scatterXlabel = "median p value", scatterYlabel = "consensus score",
                            scatterWidth = 170, scatterHeight = 150,
                            plotMethod = "median", plotPadjust = TRUE,
                            parallelComputing = FALSE, clusterType = "PSOCK"){


  # set up a temp function for boxplot with directinality
  tmpfunc <- function(x){

    pCl <- c("distinct", "mixed") # for boxplot
    classDirt <- c("up", "down") # for boxplot

    lapply(pCl, function(m)lapply(classDirt,
                                  function(n)RBioArray::rbioGS_boxplot(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""),
                                                                       KEGG = boxplotKEGG, pClass = m, classDirection = n, adjust = plotPadjust,
                                                                       n = boxplotN, xLabel = boxplotXlabel, yLabel = boxplotYlabel,
                                                                       plotTitle = boxplotTitle, plotWidth = boxplotWidth, plotHeight = boxplotHeight)))

  }

  if (!parallelComputing){

    if (boxplot){
      # boxplots
      lapply(1: length(GSlst), tmpfunc)

      lapply(1: length(GSlst), function(x)RBioArray::rbioGS_boxplot(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""),
                                                                    KEGG = boxplotKEGG, pClass = "non", adjust = plotPadjust,
                                                                    n = boxplotN, xLabel = boxplotXlabel, yLabel = boxplotYlabel,
                                                                    plotTitle = boxplotTitle, plotWidth = boxplotWidth, plotHeight = boxplotHeight))

    }

    if (scatterplot){

      # scatter plots
      lapply(1: length(GSlst), function(x)RBioArray::rbioGS_scatter(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""), cutoff = scatterplotCutoff,
                                                                    method = plotMethod, adjust = plotPadjust, rankCutoff = scatterRankline,
                                                                    pCutoff = scatterPline,
                                                                    plotTitle = scatterTitle, xLabel = scatterXlabel, yLabel = scatterYlabel,
                                                                    plotWidth = scatterWidth, plotHeight = scatterHeight))

    }

  } else { # parallel computing

    # set up cpu cores
    n_cores <- detectCores() - 1

    # parallel computing
    if (clusterType == "FORK"){ # mac and linux only

      if (boxplot){

        # boxplots
        mclapply(1: length(GSlst), tmpfunc, mc.cores = n_cores, mc.preschedule = FALSE)

        mclapply(1: length(GSlst), FUN = function(x)RBioArray::rbioGS_boxplot(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""),
                                                                              KEGG = boxplotKEGG, pClass = "non", adjust = plotPadjust,
                                                                              n = boxplotN, xLabel = boxplotXlabel, yLabel = boxplotYlabel,
                                                                              plotTitle = boxplotTitle, plotWidth = boxplotWidth, plotHeight = boxplotHeight),
                 mc.cores = n_cores, mc.preschedule = FALSE)

      }

      if (scatterplot){


        mclapply(1: length(GSlst), FUN = function(x)RBioArray::rbioGS_scatter(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""), cutoff = scatterplotCutoff,
                                                                              method = plotMethod, adjust = plotPadjust, rankCutoff = scatterRankline,
                                                                              pCutoff = scatterPline,
                                                                              plotTitle = scatterTitle, xLabel = scatterXlabel, yLabel = scatterYlabel,
                                                                              plotWidth = scatterWidth, plotHeight = scatterHeight),
                 mc.cores = n_cores, mc.preschedule = FALSE)


      }

    } else { # windows etc

      # set up clusters for PSOCK
      cl <- makeCluster(n_cores, type = clusterType, outfile = "")
      registerDoParallel(cl) # part of doParallel package
      on.exit(stopCluster(cl)) # close connect when exiting the function

      # plot
      if (boxplot){

        # boxplots
        foreach(x = 1: length(GSlst), .packages = c("RBioArray", "piano")) %dopar% {tmpfunc(x)}


        foreach(x = 1: length(GSlst), .packages = c("RBioArray", "piano")) %dopar% {
          RBioArray::rbioGS_boxplot(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""),
                                    KEGG = boxplotKEGG, pClass = "non", adjust = plotPadjust,
                                    n = boxplotN, xLabel = boxplotXlabel, yLabel = boxplotYlabel,
                                    plotTitle = boxplotTitle, plotWidth = boxplotWidth, plotHeight = boxplotHeight)
        }

      }

      if (scatterplot){
        foreach(x = 1: length(GSlst), .packages = c("RBioArray", "piano")) %dopar% {
          RBioArray::rbioGS_scatter(GSA_list = GSlst[[x]], fileName = paste(names(GSlst)[x], "_", plotGSname, sep = ""), cutoff = scatterplotCutoff,
                                    method = plotMethod, adjust = plotPadjust, rankCutoff = scatterRankline,
                                    pCutoff = scatterPline,
                                    plotTitle = scatterTitle, xLabel = scatterXlabel, yLabel = scatterYlabel,
                                    plotWidth = scatterWidth, plotHeight = scatterHeight)
        }
      }

    }

  }

}
