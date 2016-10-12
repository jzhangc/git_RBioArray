#' @title rbioGS_entrez2gene
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
#' DE_dataframe <- rbioGS_entrez2gene(dataframe)
#' }
#' @export
rbioGS_entrez2geneStats <- function(DEGdfm, cat = "SYMBOL", species = "Hs", pkg = "org.Hs.eg.db"){
  EnsemblVector <- DEGdfm$gene_name # create a vector conatining all the log-fold changes
  names(EnsemblVector) <- DEGdfm$gene_name # add names to each data point in the vector
  RNA_id <- id2eg(ids = names(EnsemblVector), category = cat, org = species, pkg.name = pkg)
  DEGdfm$SYMBOL <- DEGdfm$gene_name
  DEGdfm <- merge(DEGdfm, RNA_id, by = "SYMBOL")
  DEGdfm <- DEGdfm[complete.cases(DEGdfm), ]
  DEGdfm <- unique(DEGdfm)
  return(DEGdfm)
}



#' @title rbioGS_all
#'
#' @description Add Human entrez ID to the DE dataframe
#' @param MyGS Pre-loaded gene set objects.
#' @param pVar Gene level p values. Could be, but not exclusive to, a variable of a dataframe.
#' @param logFCVar Gene level logFC (log fold change). Could be, but not exclusive to, a variable of a dataframe.
#' @param tVar Gene leve t values. Could be, but not exclusive to, a variable of a dataframe.#'
#' @details The function is based on piano package. It runs "fisher", "stouffer", "reporter", "tailStrength", "wilcoxon" for p value based GSA, and "page", "gsea", "maxmean" for t value based GSA.
#' @return Outputs a \code{list} object with GSA results, both for p value based and t value based.
#' @importFrom piano runGSA
#' @examples
#' \dontrun{
#' GOterm_bp <- piano::loadGSC(file = "c5.bp.v5.0.entrez.gmt", type = "gmt") # load GO term biological process set
#'
#' pc_Kegg <- rbioGS_all(GOterm_bp, pVar = pcGSdfm$p_value, logFCVar = pcGSdfm$logFC, tVar = pcGSdfm$t_value)
#'
#' }
#' @export
rbioGS_all <- function(MyGS, pVar, logFCVar, tVar){

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


#' @title rbioGS_boxplot
#'
#' @description Generate boxplot from piano GS rank object,
#' @param GSrank piano GS rank object.
#' @param GS if GS = \code{"KEGG"}, the function will remove the "KEGG_" string in the GS name variable. Default is \code{"OTHER"}.
#' @param fileName Output file name.
#' @param plotWidth Set the width of the plot. Default is \code{170}.
#' @param plotHeight Set the height of the plot. Default is \code{150}.
#' @details GSrank takes the resulted object generated from the piano function consensusScores().
#' @return Outputs a \code{pdf} boxplot figure file with allGSA results.
#' @importFrom reshape2 melt
#' @importFrom grid grid.newpage grid.draw
#' @import ggplot2
#' @examples
#' \dontrun{
#' pcPos_rank_mxdn <- piano::consensusScores(pc_Pos, class = "mixed", direction="down",n = 15, adjusted = TRUE,
#'                                    method = "median", plot = TRUE, showLegend = F,
#'                                    rowNames = "names")
#'
#' rbioGS_boxplot(pcPos_rank_mxdn, fileName = "pcPos_rank_mxdn", plotWidth = 260, plotHeight = 240)
#'
#' }
#' @export
rbioGS_boxplot <- function(GSrank, GS = "OTHER", fileName,
                      plotWidth = 170, plotHeight = 150){
  # prepare the dataframe for ggplot2
  DFM <- data.frame(GSrank$rankMat)
  DFM <- DFM[, c(-1, -2)]
  DFM <- data.frame(t(DFM)) # transpose the dataframe using t()
  DFM$enrichment <- rownames(DFM)
  DFM_mlt <- melt(DFM, id.vars = colnames(DFM)[length(colnames(DFM))])

  if (GS == "KEGG"){ # remove "KEGG_" and the space between the terms
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
  plt <- ggplot(DFM_mlt, aes(x = variable, y = value)) + geom_boxplot() +
    guides(fill = FALSE) + scale_x_discrete(limits = with(DFM_mlt, rev(levels(variable)))) +
    ylab("rank") +
    xlab(NULL) +
    coord_flip() + # flip the axes
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          axis.title = element_text(face = "bold"),
          legend.position = "bottom")

  # export the file and draw a preview
  ggsave(filename = paste(fileName,".boxplot.pdf",sep = ""), plot = plt,
         width = plotWidth, height = plotHeight, units = "mm",dpi = 600)
  grid.draw(plt) # preview
}


#' @title rbioGS_scatter
#'
#' @description Generate scatter plot for piano GS rank heatmap obejct.
#' @param GSAList GSA list generated from \code{\link{rbioArray_allGSA}}.
#' @param GS if GS = \code{"KEGG"}, the function will remove the "KEGG_" string in the GS name variable. Default is \code{"OTHER"}.
#' @param rankCutoff Cutoff value for GS rank line.
#' @param pCutoff Cutoff value for GS p value line.
#' @param fileName Output file name.
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
#' rbioGS_scatter(GS_Pos, rankCutoff = 50, pCutoff = 0.05, fileName = "GS_pos")
#'
#' }
#' @export
rbioGS_scatter<-function(GSAList, rankCutoff, pCutoff,
                       fileName, plotWidth = 170, plotHeight = 150){
  HTmap <-consensusHeatmap(GSAList, cutoff = 15, method = "median",
                             colorkey = TRUE, colorgrad = c("blue","white"),
                             cellnote = "none")

  ## data frame prep
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


  ## ggplotting
  grid.newpage()
  ScatterP<-ggplot(dfm4plot, aes(x = p_value, y = rank))+
    geom_point(aes(shape = factor(p_class)), size=3)+
    scale_x_continuous(breaks = c(0.1, 0.05, 0),
                       labels = c("0.1", "0.05", "0"),
                       expand = c(0, 0),
                       trans = "reverse", limits = c(0.1, NA))+
    scale_y_continuous(breaks = c(100, 50, 0),
                       labels = c("100", "50", "0"),
                       expand = c(0, 0),
                       trans = "reverse", limits = c(100, NA))+
    geom_vline(xintercept = pCutoff)+
    geom_hline(yintercept = rankCutoff)+
    labs(x="median p value",y = "consensus score")+
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_rect(),
          legend.key = element_blank())+
    scale_shape_manual(values = c(1:5))

  # export the files and draw a preview
  write.csv(dfm4plot, file = paste(fileName,".scatterplot.csv",sep = ""))
  ggsave(filename = paste(fileName,".scatterplot.pdf",sep = ""), plot = ScatterP,
         width = plotWidth, height = plotHeight, units = "mm",dpi = 600)
  grid.draw(ScatterP) # preview

  # return the plot dataframe
  return(dfm4plot)
}

#' @title rbioGS_kegg
#'
#' @description Download and generate DE results masked kegg pathway figures
#' @param dfm GS dataframe with \code{ENTREZID} and \code{logFC} variables.
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
rbioGS_kegg<- function(dfm, keggID, suffix, species = "hsa"){
  logFC <- dfm$logFC
  names(logFC) <- dfm$ENTREZID

  KEGG <- pathview(gene.data = logFC, pathway.id = keggID, species = species,
                   out.suffix = suffix, keys.align = "y", kegg.native = TRUE, match.data = FALSE,
                   key.pos = "topright")

  # set the .GlobalEnv to the envir argument so that the assign function will assign the value to a global object, aka outside the function
  return(assign(paste("kegg_", keggID, "_" , deparse(substitute(dfm)), sep = "") ,KEGG, envir = .GlobalEnv))
}

