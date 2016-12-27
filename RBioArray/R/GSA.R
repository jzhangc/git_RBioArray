#' @title rbioGS_entrez2gene
#'
#' @description Add Human entrez ID to the DE dataframe
#' @param DEGdfm Input DE data frame.
#' @param cat Catergory, default is \code{"SYMBOL"}.
#' @param species Set the species, default is \code{"Hs"}.
#' @param pkg Name for the annotation package. The current package includes the human version.
#' @return Outputs a \code{dataframe} object with Entrez ID.
#' @importFrom pathview id2eg
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
#' @param multicore If to use parallel computing or not. Default is \code{FALSE}
#' @param clusterType Only set when \code{multicore = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @details The function is based on piano package. It runs "fisher", "stouffer", "reporter", "tailStrength", "wilcoxon" for p value based GSA, and "page", "gsea", "maxmean" for t value based GSA.
#' @return Outputs a \code{list} object with GSA results, both for p value based and t value based.
#' @importFrom piano runGSA
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @examples
#' \dontrun{
#' gsoutput <- rbioGS(GS = kegg, pVar = dfm$p_value, logFCVar = dfm$logFC, tVar = dfm$t_value, idVar = dfm$EntrezID, multicore = TRUE, clusterType = "FORK")
#'
#'
#' }
#' @export
rbioGS <- function(GS, pVar, logFCVar, tVar, idVar,
                   method_p = c("fisher", "stouffer", "reporter", "tailStrength", "wilcoxon"),
                   method_t = c("page", "gsea", "maxmean"),
                   multicore = FALSE, clusterType = "PSOCK", ...){

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

  if (!multicore){

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
      cl <- makeCluster(n_cores, type = clusterType)
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
#' @param GS_list piano GS results object.
#' @param ... Arguments passing to \code{consensusScores} function from \code{piano} package. See the corresponding \code{piano} package help page for details.
#' @param GS if GS = \code{"KEGG"}, the function will remove the "KEGG_" string in the GS name variable. Default is \code{"OTHER"}.
#' @param fileName Output file name. Default is \code{"GS_list"}.
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
#' rbioGS_boxplot(GS_object, fileName = "GS_analysis", class = "mixed", direction="down", n = 15, adjusted = TRUE, method = "median", plot = TRUE, rowNames = "names", plotWidth = 260, plotHeight = 240)
#'
#' }
#' @export
rbioGS_boxplot <- function(GS_list, ..., GS = "OTHER", fileName = "GS_list",
                           plotWidth = 170, plotHeight = 150){

  # prepare consensus score list
  GSrank <- consensusScores(resList = GS_list, plot = FALSE, ...)


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
#' @param GS_list GSA list generated from \code{\link{rbioArray_allGSA}}.
#' @param ... Arguments passing to \code{consensusHeatmap} function from \code{piano} package. See the responding help page of \code{piano} for details.
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
#' rbioGS_scatter(GS_Pos, cutoff = 15, method = "median", adjusted = TRUE, rankCutoff = 50, pCutoff = 0.05, fileName = "GS_pos")
#'
#' }
#' @export
rbioGS_scatter <- function(GSAList, ..., rankCutoff, pCutoff,
                       fileName, plotWidth = 170, plotHeight = 150){

  HTmap <-consensusHeatmap(GSAList, plot = FALSE, ...)

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

#' @title rbioGS_all
#'
#' @description Add Human entrez ID to the DE dataframe
#' @param fileName Output figure file name.
#' @param input Input list object for DE results.
#' @param entrezVar Name of the EntrezID variable in the \code{input} object.
#' @param GS Pre-loaded gene set objects.
#' @param multicore If to use parallel computing or not. Default is \code{FALSE}
#' @param clusterType Only set when \code{multicore = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @details This is an all-in-one function for GS anlyasis based on piano package. It runs "fisher", "stouffer", "reporter", "tailStrength", "wilcoxon" for p value based GSA, and "page", "gsea", "maxmean" for t value based GSA.
#' @return Outputs  \code{csv} files and \code{pdf} figure files with GSA results.
#' @importFrom piano runGSA
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @examples
#' \dontrun{
#' rbioGS_all(GS = kegg, pVar = dfm$p_value, logFCVar = dfm$logFC, tVar = dfm$t_value, idVar = dfm$EntrezID, multicore = TRUE, clusterType = "FORK")
#'
#'
#' }
#' @export
rbioGS_all <- function(fileName, input, entrezVar = NULL,
                   GS = NULL, plot = TRUE,
                   multicore = FALSE, clusterType = "PSOCK"){

  if (is.null(entrezVar)){
    stop("please tell the function the name of the Entrez ID vaiable")
  }

  if (is.null(GS)){
    stop("please choose a proper gene set")
  }

  ##
  GSlst <- vector(mode = "list", length(names(input)))
  names(GSlst) <- names(input)


  ## GSA
  if (!multicore){
    # remove the rows with NA in the Entrez ID variable
    DELst <- lapply(input, function(x)x[complete.cases(x[, entrezVar]), ])

    # run GSA
    GSlst[] <- lapply(DElst, function(i)RBioArray::rbioGS(GS = GS, pVar = i$P.Value, logFCVar = i$logFC,
                                                          tVar = i$t, idVar = i[, entrezVar]), multicore = multicore, clusterType = clusterType)

    if (plot){
      lapply(GSlst, RBioArray::rbioGS_boxplot)
      lapply(GSlst, RBioArray::rbioGS_scatter)
    }

  } else { # parallel computing

    # set up cpu cores
    n_cores <- detectCores() - 1

    # parallel computing
    if (clusterType == "FORK"){ # mac and linux only
      # remove the rows with NA in the Entrez ID variable
      DELst <- mclapply(input, function(x)x[complete.cases(x[, entrezVar]), ], mc.cores = n_cores, mc.preschedule = FALSE)

      # run GSA
      GSlst[] <- mclapply(DElst, function(i)RBioArray::rbioGS(GS = GS, pVar = i$P.Value, logFCVar = i$logFC,
                                                              tVar = i$t, idVar = i[, entrezVar]), mc.cores = n_cores, mc.preschedule = FALSE)

      if (plot){
        mclapply(GSlst, RBioArray::rbioGS_boxplot, mc.cores = n_cores, mc.preschedule = FALSE)
        mclapply(GSlst, RBioArray::rbioGS_scatter, mc.cores = n_cores, mc.preschedule = FALSE)
      }

    } else { # windows etc

      # set up clusters for PSOCK
      cl <- makeCluster(n_cores, type = clusterType)
      registerDoParallel(cl) # part of doParallel package
      on.exit(stopCluster(cl)) # close connect when exiting the function


      GS_list_p[] <- foreach(x = input) %dopar% {
        out <- function(x)x[complete.cases(x[, entrezVar]), ]
      }

      GS_list[] <- foreach(i = DElst, .packages = c("RBioArray", "piano")) %dopar% {
        out <- rbioGS(GS = GS, pVar = i$P.Value, logFCVar = i$logFC,
                      tVar = i$t, idVar = i[, entrezVar])
      }

      if (plot){
        foreach(x = GSlst, .packages = c("RBioArray", "piano")) %dopar% {
          out <- rbioGS_boxplot(x)
        }

        foreach(x = GSlst, .packages = c("RBioArray", "piano")) %dopar% {
          out <- rbioGS_scatter(x)
        }

      }

    }

  }

}
