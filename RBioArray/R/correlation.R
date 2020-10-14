#' @title cor_pvalue
#'
#' @description Calculate p value for the correlation analysis
#' @param r Correlation coefficient
#' @param n Sample size
#' @return P value for the correlation analysis
#' @examples
#' \dontrun{
#'
#' p <- cor_pvalue(r = 0.806687547, n = 6)
#'
#' }
#' @export
cor_pvalue <- function(r, n){
  t <- r/sqrt((1 - r^2)/(n - 2))
  p <- 2*pt(-abs(t), df = n - 2)
  return(p)
}


#' @title rbio_unsupervised_corcluster
#'
#' @description Generic unsupersived correlation clustering function.
#' @param object Input object in either \code{rbioseq_de} or \code{rbioarray_de} class.
#' @param ... Additional arguments for corresponding S3 class methods.
#' @details It behaves in the following ways:
#'          1. heatmap uses \code{heatmap.2} function from \code{gplots} package.
#'          2. significance plot uses \code{corrplot} package
#'
#'          Since the function depends on the \code{heatmap.2} function from \code{gplots} package for heatmap,
#'          arguments can be passed directly, as seen in the examples below.
#'
#' @return PDF files for unsupervised correlation heatmap and if set, significance plot. Also CSV files containing unsupervised correlation analysis results.
#'         The function also generates two list objects containing correlation results
#'
#'         It is generally not recommend to perform correlation analysis on the huge number of genes/probes/genomic features - extremely time consuming.
#' @examples
#'
#' \dontrun{
#' # rbioarray_de class input
#' rbio_unsupervised_corcluster(object = delist, n = 50, cor.method = "pearson",
#'                              rm.control = TRUE,
#'                              gene_symbol.only = TRUE,
#'                              map.colour = "PRGn", n.map.colour = 11,
#'                              heatmap.axis.label = FALSE,
#'                              cexCol = 0.6, cexRow = 0.6, margins = c(8, 8), key.title = NA,
#'                              heatmap.width = 7, heatmap.height = 7,
#'                              sigplot = TRUE,
#'                              sigplot.adj.p = FALSE, sigplot.alpha = 0.05,
#'                              sigplot.sigLabelColour = "red", sigplot.sigLabelSize = 1,
#'                              sigplot.labelColour = "black", sigplot.labelSize = 0.6, sigplot.labelAngle = 90,
#'                              sigplot.keySize = 1,
#'                              sigplot.mar = c(1.25, 1, 1, 0.5),
#'                              sigplot.width = 7, sigplot.height = 7)
#'
#' # rbioseq_de class input
#' rbio_unsupervised_corcluster(object = mrna_de, n = 50, cor.method = "pearson",
#'                              rm.control = TRUE,
#'                              gene_symbol.only = TRUE,
#'                              map.colour = "PRGn", n.map.colour = 11,
#'                              heatmap.axis.label = FALSE,
#'                              cexCol = 0.6, cexRow = 0.6, margins = c(8, 8), key.title = NA,
#'                              heatmap.width = 7, heatmap.height = 7,
#'                              sigplot = TRUE,
#'                              sigplot.adj.p = FALSE, sigplot.alpha = 0.05,
#'                              sigplot.sigLabelColour = "red", sigplot.sigLabelSize = 1,
#'                              sigplot.labelColour = "black", sigplot.labelSize = 0.6, sigplot.labelAngle = 90,
#'                              sigplot.keySize = 1,
#'                              sigplot.mar = c(1.25, 1, 1, 0.5),
#'                              sigplot.width = 7, sigplot.height = 7)
#' }
#'
#' @export
rbio_unsupervised_corcluster <- function(object, ...){
  ## check arguments
  if (!any(class(object) %in% c("rbioarray_de", "rbioseq_de"))) stop("The input object needs to be either \"rbioarray_de\" or \"rbioseq_de\" class object.")

  ## use methods
  UseMethod("rbio_unsupervised_corcluster", object)
}


#' @title rbio_unsupervised_corcluster.rbioarray_de
#'
#' @rdname rbio_unsupervised_corcluster
#' @method rbio_unsupervised_corcluster rbioarray_de
#' @param object Input object in \code{rbioarray_de} class.
#' @param ... Additional arguments for the default method.
#' @export
rbio_unsupervised_corcluster.rbioarray_de <- function(object, ...){
  ## variables
  export.name <- deparse(substitute(object))

  ## use methods
  rbio_unsupervised_corcluster.default(E = object$input_data$E, genes = object$input_data$genes,
                                       input.genes_annotation.control_type = object$genes_annotation.control_type,
                                       input.genes_annotation.gene_symbol.var_name = object$genes_annotation.gene_symbol.var_name,
                                       input.genes_annotation.gene_id.var_name = object$genes_annotation.gene_id.var_name,
                                       input.sample_groups = object$sample_groups,
                                       input.comparisons = object$comparisons,
                                       export.name = export.name, ...)
}


#' @title rbio_unsupervised_corcluster.rbioseq_de
#'
#' @rdname rbio_unsupervised_corcluster
#' @method rbio_unsupervised_corcluster rbioseq_de
#' @param object Input object in \code{rbioseq_de} class.
#' @param ... Additional arguments for the default method.
#' @details The function uses filtered count data, as opposed to normalized data.
#'          Due to the compositional nature of NGS data, the count data is transformed using CLR method prior to clustering.
#' @export
rbio_unsupervised_corcluster.rbioseq_de <- function(object, sample_id.var.name = NULL, ...){
  ## variables
  export.name <- deparse(substitute(object))
  # transform
  # cat("CLR transformation of filtered RNAseq count data...")
  # E_transfo <- rbioseq_clr_ilr_transfo(object$filter_results$filtered_counts$counts, offset = 1, mode = "clr")  # clr tranformation
  # cat("Done!\n")

  ## use methods
  # rbio_unsupervised_corcluster.default(E = E_transfo, genes = object$filter_results$filtered_counts$genes,
  #                 input.genes_annotation.control_type = NULL,
  #                 input.genes_annotation.gene_symbol.var_name = object$genes_annotation.gene_symbol.var_name,
  #                 input.genes_annotation.gene_id.var_name = object$genes_annotation.gene_id.var_name,
  #                 input.sample_groups = object$sample_groups,
  #                 input.comparisons = object$comparisons,
  #                 export.name = export.name, ...)
  rbio_unsupervised_corcluster.default(E = object$normalized_data$E, genes = object$filter_results$filtered_counts$genes,
                                       input.genes_annotation.control_type = NULL,
                                       input.genes_annotation.gene_symbol.var_name = object$genes_annotation.gene_symbol.var_name,
                                       input.genes_annotation.gene_id.var_name = object$genes_annotation.gene_id.var_name,
                                       input.sample_groups = object$sample_groups,
                                       input.comparisons = object$comparisons,
                                       export.name = export.name, ...)
}


#' @title rbio_unsupervised_corcluster.default
#'
#' @rdname rbio_unsupervised_corcluster
#' @method rbio_unsupervised_corcluster default
#' @param E Expression or count matrix, with rows for genes/probes/genomic features, columns for RNA samples.
#' @param genes Annotation data frame for genes/probes/genomic features.
#' @param input.sample_groups Input \code{factor} object for sample grouping labels.
#' @param n Number of genes/probes/genomic features to be clustered, numeric input or \code{"all"}. Default is \code{"all"}.
#' @param rm.control Whether to remove control probes (Agilent platform) or not. Default is \code{TRUE}.
#' @param input.genes_annotation.control_type Only set when \code{rm.control = TRUE}, input control type variable annotation list.
#' @param gene_symbol.only Whether or not to remove probes without gene symbol. Default is \code{FALSE}.
#' @param input.genes_annotation.gene_symbol.var_name Only set when \code{gene_symbol.only = TRUE}, variable name for gene symbol column in \code{genes} data frame.
#' @param input.genes_annotation.gene_id.var_name Only set when \code{gene_symbol.only = TRUE}, variable name for gene id column in \code{genes} data frame.
#' @param cor.method The correlation method, options are "spearman" and "pearson". Default is \code{"pearson"}.
#' @param export.name File name for the export \code{pdf} plot file.
#' @param hclust.heatmap Boolean. Whether or not to export a hcluster heatmap. Default is \code{TRUE}.
#' @param hclust.method String. Method for hclust function. Default is \code{"complete"}.
#' @param map.colour Heat map colour. Default is \code{"PRGn"}. See \code{RColorBrewer} package for more.
#' @param n.map.colour Number of colours displayed. Default is \code{11}. See \code{RColorBrewer} package for more.
#' @param heatmap.axis.label Whether to display label for both x- and y- axes. Default is \code{FALSE}.
#' @param ... Additional arguments for \code{heatmap.2} function from \code{gplots} package.
#' @param heatmap.width Width for correlation heatmap. Unit is \code{inch}. Default is \code{7}.
#' @param heatmap.height Height for correlation heatmap. Unit is \code{inch}. Default is \code{7}.
#' @param sigplot If to show significance plot. Default is \code{TRUE}.
#' @param sigplot.adj.p If to use FDR corrected correlation p value to plot sigPlot. Default is \code{FALSE}.
#' @param sigplot.alpha The alpha value for correlation p value. Default is \code{0.05}
#' @param sigplot.sigLabelColour The colour for label for the significant pairs. Default is \code{"red"}.
#' @param sigplot.sigLabelSize The size for label for the significant pairs. Default is \code{3}.
#' @param sigplot.labelColour The significance heatmap axis label colour. Default is \code{"black"}.
#' @param sigplot.labelSize The significance heatmap axis label size. Default is \code{1}.
#' @param sigplot.labelAngle The significance heatmap axis label angle. Default is \code{90}.
#' @param sigplot.keySize The significance heatmap colour key size. Default is \code{1}.
#' @param sigplot.mar The A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1.
#' @param sigplot.width Width for correlation significance plot. Unit is \code{inch}. Default is \code{7}.
#' @param sigplot.height Height for correlation significance plot. Unit is \code{inch}. Default is \code{7}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return A heatmap based on hierarchical clustering analysis in \code{pdf} format.
#' @import corrplot
#' @import foreach
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @export
rbio_unsupervised_corcluster.default <- function(E, genes, input.sample_groups, n = "all",
                                                 input.comparisons = NULL,
                                                 cor.method = c("pearson", "spearman"),
                                                 rm.control = FALSE, input.genes_annotation.control_type,
                                                 gene_symbol.only = FALSE,
                                                 input.genes_annotation.gene_symbol.var_name = NULL,
                                                 input.genes_annotation.gene_id.var_name = NULL,
                                                 export.name = NULL,
                                                 hclust.heatmap = TRUE,
                                                 hclust.method = c("complete", "ward.D", "ward.D2", "single",  "average", "mcquitty", "median", "centroid"),
                                                 map.colour = "PRGn", n.map.colour = 11, ...,
                                                 heatmap.axis.label = FALSE, heatmap.width = 7, heatmap.height = 7,
                                                 sigplot = TRUE,
                                                 sigplot.adj.p = FALSE, sigplot.alpha = 0.05,
                                                 sigplot.sigLabelColour = "red", sigplot.sigLabelSize = 3,
                                                 sigplot.labelColour = "black", sigplot.labelSize = 1, sigplot.labelAngle = 90, sigplot.keySize = 1,
                                                 sigplot.mar = c(5, 4, 4, 2) + 0.1,
                                                 sigplot.width = 7, sigplot.height = 7, verbose = TRUE){
  ## check arguments
  cor.method <- match.arg(tolower(cor.method), c("pearson", "spearman"))

  if (n != "all" && n %% 1 != 0) stop("Argument n needs to be either \"all\" or an integer number.")
  if (n.map.colour %% 1 != 0) stop("Argument n.map.colour needs to be an integer number.")
  if (rm.control && is.null(input.genes_annotation.control_type)) {
    cat("Argument input.genes_annotation.control_type is NULL when rm.control = TRUE, automatically set rm.control = FALSE.\n\n")
    rm.control <- FALSE
  }
  if (gene_symbol.only && ! input.genes_annotation.gene_symbol.var_name %in% names(genes)) {
    if (verbose) cat("Argument input.genes_annotation.gene_symbol.var_name not found in genes data frame when gene_symbol.only = TRUE, automatically set gene_symbol.only = FALSE.\n\n")
    gene_symbol.only <- FALSE
  }
  # if (!cor.method %in% c("pearson", "spearman")) stop("Argument cor.method needs to be \"pearson\" or \"spearman\".")
  if (missing(export.name) || is.null(export.name)) stop("Please set value for argument export.name.")

  # check contrast levels against sample groups
  contra_levels_all <- unique(foreach(i = 1:length(input.comparisons$comparison_levels), .combine = "c") %do% {
    input.comparisons$comparison_levels[[i]]
  })
  if (!all(contra_levels_all %in% unique(levels(input.sample_groups)))) stop("Contrast levels not matching sample groups. Please check the input.")

  ## variables
  row.lab.var_name <- input.genes_annotation.gene_id.var_name

  # hclust function
  clustfunc <- function(x)hclust(x, method = hclust.method)

  ## prepare dfm for clustering
  dfm <- data.frame(genes, E, check.names = FALSE)

  if (rm.control){ # remove control
    dfm <- dfm[dfm[, input.genes_annotation.control_type$control_type.var_name] == input.genes_annotation.control_type$exp_type.value, ]
  }
  if (gene_symbol.only) {
    dfm <- dfm[complete.cases(dfm[, input.genes_annotation.gene_symbol.var_name]),]
    row.lab.var_name <- input.genes_annotation.gene_symbol.var_name
  }
  if (n != "all"){ # subset
    dfm <- dfm[1:n, ]
  }

  mtx <- as.matrix(dfm[, !names(dfm) %in% names(genes)])
  row.lab <- dfm[, row.lab.var_name]
  if (heatmap.axis.label){
    heatmap.labRow <- row.lab
    heatmap.labCol <- row.lab
  } else {
    heatmap.labRow <- NA
    heatmap.labCol <- NA
  }

  ## cluster
  if (verbose) cat(paste0("Unsupervised correlation analysis method: ", cor.method, "\n"))
  if (verbose) cat(paste0("Genes/Probes/Genomic features assessed: ", n, "\n"))
  if (verbose) cat("-------------------\n")
  corcoef_list <- vector(mode = "list", length = length(contra_levels_all))
  corp_list <- vector(mode = "list", length = length(contra_levels_all))
  adj_corp_list <- vector(mode = "list", length = length(contra_levels_all))
  cor_names_vector <- vector(length = length(contra_levels_all))
  for (i in seq(length(contra_levels_all))) {
    cor_sample_level <- contra_levels_all[i]
    cor_mtx <- mtx[, which(input.sample_groups %in% cor_sample_level)]
    cor_mtx <- t(cor_mtx)
    corcoef <- cor(cor_mtx, method = cor.method)
    rownames(corcoef) <- row.lab
    colnames(corcoef) <- row.lab
    corp <- foreach(m = corcoef, .combine = "cbind") %do% cor_pvalue(m, n = nrow(cor_mtx)) # p value matrix
    diag(corp) <- NA
    adj_corp <- matrix(p.adjust(corp, method = "fdr"), nrow = nrow(corp), byrow = T)  # fdr
    rownames(adj_corp) <- rownames(corp)
    colnames(adj_corp) <- colnames(corp)
    corname <- cor_sample_level

    if (hclust.heatmap){
      if (verbose) cat(paste0("Correlation heatmap saved to file: ", cor_sample_level, "_cor.unsuper.heatmap.pdf..."))
      pdf(file = paste0(cor_sample_level, "_cor.unsuper.heatmap.pdf"), width = heatmap.width, height = heatmap.width)
      heatmap.2(corcoef, symm = TRUE, trace = "none",
                distfun = function(x)as.dist(1-corcoef),
                hclustfun = clustfunc,
                col = brewer.pal(n.map.colour, map.colour), labRow = heatmap.labRow, labCol = heatmap.labCol,
                ...)
      dev.off()
      if (verbose) cat("Done!\n")
    }
    corcoef_list[[i]] <- corcoef
    corp_list[[i]] <- corp
    adj_corp_list[[i]] <- adj_corp
    cor_names_vector[i] <- corname
  }
  names(corcoef_list) <- cor_names_vector
  names(corp_list) <- cor_names_vector
  names(adj_corp_list) <- cor_names_vector

  # sig plot
  if (sigplot) {
    if (sigplot.adj.p) {
      sigplot.corp_list <- adj_corp_list
    } else {
      sigplot.corp_list <- corp_list
    }

    if (verbose) cat("\n")
    for (i in seq(length(corcoef_list))) {
      tryCatch(
        {
          pdf(file = paste(names(corcoef_list)[i], "_cor.unsuper.sigplot.pdf", sep = ""), width = sigplot.width, height = sigplot.height)
          corrplot(corr = corcoef_list[[i]], method = "color", type = "upper", p.mat = sigplot.corp_list[[i]], sig.level = sigplot.alpha,
                   insig = c("label_sig"), pch.col = sigplot.sigLabelColour, pch.cex = sigplot.sigLabelSize,
                   col = brewer.pal(n.map.colour, map.colour), mar = sigplot.mar,
                   tl.col = sigplot.labelColour, tl.cex = sigplot.labelSize, tl.srt = sigplot.labelAngle, cl.length = 3, cl.cex = sigplot.keySize)
          if (verbose) cat(paste0("Correlation significance plot saved to file: ", names(corcoef_list)[i], "_cor.unsuper.sigplot.pdf..."))
          dev.off()
          if (verbose) cat("Done!\n")
        },
        error = function(err){
          if (verbose) cat(paste0("No significant correlation found for ", names(corcoef_list)[i], ". ", "Therefore no significance plot generated.\n"))
          dev.off()
        }
      )
    }
  }

  # csv export
  if (verbose) cat("\n")
  out_env <- vector(mode = "list", length = length(corcoef_list))
  for (i in seq(length(corcoef_list))) {
    out_corcoef <- corcoef_list[[i]]
    out_corp <- corp_list[[i]]
    out_adj_corp <- adj_corp_list[[i]]
    outdfm1 <- data.frame(
      group = paste(rownames(out_corcoef)[row(out_corcoef)[upper.tri(out_corcoef)]], colnames(out_corcoef)[col(out_corcoef)[upper.tri(out_corcoef)]], sep = "_"),
      row = rownames(out_corcoef)[row(out_corcoef)[upper.tri(out_corcoef)]],
      col = colnames(out_corcoef)[col(out_corcoef)[upper.tri(out_corcoef)]],
      coefficient = out_corcoef[upper.tri(out_corcoef)]
    )

    outdfm2 <- data.frame(
      group = paste(rownames(out_corp)[row(out_corp)[upper.tri(out_corp)]], colnames(out_corp)[col(out_corp)[upper.tri(out_corp)]], sep = "_"),
      row = rownames(out_corp)[row(out_corp)[upper.tri(out_corp)]],
      col = colnames(out_corp)[col(out_corp)[upper.tri(out_corp)]],
      p.value = out_corp[upper.tri(out_corp)]
    )

    outdfm3 <- data.frame(
      group = paste(rownames(out_adj_corp)[row(out_adj_corp)[upper.tri(out_adj_corp)]], colnames(out_adj_corp)[col(out_adj_corp)[upper.tri(out_adj_corp)]], sep = "_"),
      row = rownames(out_adj_corp)[row(out_adj_corp)[upper.tri(out_adj_corp)]],
      col = colnames(out_adj_corp)[col(out_adj_corp)[upper.tri(out_adj_corp)]],
      adj.p.value = out_adj_corp[upper.tri(out_adj_corp)]
    )

    out <- merge(outdfm1, outdfm2[, c("group", "p.value")], by = "group", x.all = TRUE)
    out <- merge(out, outdfm3[, c("group", "adj.p.value")], x.all = TRUE)
    out_env[[i]] <- out

    # if (verbose) cat(paste0("Correlation analysis results saved to file: ", names(corcoef_list)[i], ".cor.unsuper.csv..."))
    # write.csv(out, file = paste0(names(corcoef_list)[i], ".cor.unsuper.csv"), row.names = FALSE)
    # if (verbose) cat("Done!\n")
  }
  names(out_env) <- names(corcoef_list)
  assign(paste0(export.name, "_cor_unsuper"), out_env, envir = .GlobalEnv)
  assign(paste0(export.name, "_cor_adjacency_unsuper"), corcoef_list, envir = .GlobalEnv)
}


#' @title rbio_supervised_corcluster
#'
#' @description Supervised correlation clustering analysis and heatmap visualization for both microarray and RNAseq, for gene co-expression analysis.
#' @param object Input \code{sig} class object.
#' @param gene_symbol.only If to only use probes/genes/genomic features with a gene sybmol id. Default is \code{FALSE}.
#' @param cor.method The correlation method, options are "pearson", "spearman" and "pearson". Default is \code{"pearson"}.
#' @param hclust.heatmap Boolean. Whether or not to export a hcluster heatmap. Default is \code{TRUE}.
#' @param hclust.method String. Method for hclust function. Default is \code{"complete"}.
#' @param map.colour Heat map colour. Default is \code{"PRGn"}. See \code{RColorBrewer} package for more.
#' @param n.map.colour Number of colours displayed. Default is \code{11}. See \code{RColorBrewer} package for more.
#' @param ... Additional arguments for \code{heatmap.2} function from \code{gplots} package.
#' @param heatmap.width Width for correlation heatmap. Unit is \code{inch}. Default is \code{7}.
#' @param heatmap.height Height for correlation heatmap. Unit is \code{inch}. Default is \code{7}.
#' @param sigplot If to show significance plot. Default is \code{TRUE}.
#' @param sigplot.adj.p If to use FDR corrected correlation p value to plot sigPlot. Default is \code{FALSE}.
#' @param sigplot.alpha The alpha value for correlation p value. Default is \code{0.05}
#' @param sigplot.sigLabelColour The colour for label for the significant pairs. Default is \code{"red"}.
#' @param sigplot.sigLabelSize The size for label for the significant pairs. Default is \code{3}.
#' @param sigplot.labelColour The significance heatmap axis label colour. Default is \code{"black"}.
#' @param sigplot.labelSize The significance heatmap axis label size. Default is \code{1}.
#' @param sigplot.labelAngle The significance heatmap axis label angle. Default is \code{90}.
#' @param sigplot.keySize The significance heatmap colour key size. Default is \code{1}.
#' @param sigplot.mar The A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1.
#' @param sigplot.width Width for correlation significance plot. Unit is \code{inch}. Default is \code{7}.
#' @param sigplot.height Height for correlation significance plot. Unit is \code{inch}. Default is \code{7}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @details The function takes \code{sig} class object. It behaves in the following ways:
#'          1. heatmap uses \code{heatmap.2()} function from \code{gplots} package.
#'          2. significance plot uses \code{corrplot} package
#'          3. the function will use normalized and filtered expression data for both RNAseq and microaorray
#'          4. the function will automatically subset data using \code{thresholding_summary} from \code{sig} class input
#'          5. the function will automatically subset and plot for each group under each comparison
#'
#' @return A supervised heatmap based on hierarchical clustering analysis in \code{pdf} format, along with correlation coefficient and p value matrices.
#'         The function also outputs a significance plot.
#'         The function also generates two list objects containing correlation results
#' @import corrplot
#' @import foreach
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' rbio_supervised_corcluster(object = sig_input, gene_symbol.only = TRUE, cor.method = "pearson",
#'                            map.colour = "PRGn", n.map.colour = 11,
#'                            heatmap.width = 7, heatmap.height = 7,
#'                            sigplot.adj.p = FALSE, sigplot.alpha = 0.05, sigplot.sigLabelColour = "red", sigplot.sigLabelSize = 3,
#'                            sigplot.labelColour = "black", sigplot.labelSize = 1, sigplot.labelAngle = 90, sigplot.keySize = 1,
#'                            sigplot.width = 7, sigplot.height = 7)
#' }
#' @export
rbio_supervised_corcluster <- function(object,
                                       gene_symbol.only = FALSE,
                                       cor.method = c("pearson", "spearman"),
                                       hclust.heatmap = TRUE,
                                       hclust.method = c("complete", "ward.D", "ward.D2", "single",  "average", "mcquitty", "median", "centroid"),
                                       map.colour = "PRGn", n.map.colour = 11, ...,
                                       heatmap.width = 7, heatmap.height = 7,
                                       sigplot = TRUE,
                                       sigplot.adj.p = FALSE, sigplot.alpha = 0.05,
                                       sigplot.sigLabelColour = "red", sigplot.sigLabelSize = 3,
                                       sigplot.labelColour = "black", sigplot.labelSize = 1, sigplot.labelAngle = 90, sigplot.keySize = 1,
                                       sigplot.mar = c(5, 4, 4, 2) + 0.1,
                                       sigplot.width = 7, sigplot.height = 7, verbose = TRUE){
  ## argument check
  cor.method <- match.arg(tolower(cor.method), c("pearson", "spearman"))

  if (any(class(object) != "sig")) stop("The input object has to be a \"sig\" class.")
  if (n.map.colour %% 1 != 0) stop("Argument n.map.colour needs to be an integer number.")
  # if (!cor.method %in% c("pearson", "spearman")) stop("Argument cor.method needs to be \"pearson\" or \"spearman\".")

  ## variables
  E <- object$input_data$norm_E
  genes <- object$input_data$genes
  input.genes_annotation.gene_symbol.var_name = object$input_data$input.genes_annotation.gene_symbol.var_name
  input.genes_annotation.gene_id.var_name = object$input_data$input.genes_annotation.gene_id.var_name
  export.name <- deparse(substitute(object))
  input.sample_groups <- object$input_data$sample_groups
  input.genes_annotation.control_type <- object$input_data$input.genes_annotation.control_type
  comparisons <- object$input_data$comparisons$comparisons
  comparison_levels <- object$input_data$comparisons$comparison_levels
  comp_to_remove <- which(as.numeric(object$significant_change_summary[, "True"]) < 2)
  thresholding_summary <- object$thresholding_summary
  row.lab.var_name <- input.genes_annotation.gene_id.var_name

  if (gene_symbol.only && ! input.genes_annotation.gene_symbol.var_name %in% names(genes)) {
    cat("Argument input.genes_annotation.gene_symbol.var_name not found in genes data frame when gene_symbol.only = TRUE, automatically set gene_symbol.only = FALSE.\n\n")
    gene_symbol.only <- FALSE
  }
  if (is.null(input.genes_annotation.control_type)) {
    cat("Argument input.genes_annotation.control_type is NULL, no control probes are removed.\n\n")
    rm.control <- FALSE
  } else {
    rm.control <- TRUE
  }

  if (length(comp_to_remove) > 0) {
    cat("Comparisons with less than two significant changes were removed: ", comparisons[comp_to_remove], "\n")
    comparisons <- comparisons[-comp_to_remove]
    comparison_levels <- comparison_levels[-comp_to_remove]
    thresholding_summary <- thresholding_summary[-comp_to_remove]
  }

  # hclust function
  clustfunc <- function(x)hclust(x, method = hclust.method)

  # prepare dfm for clustering
  dfm <- data.frame(genes, E, check.names = FALSE)

  if (rm.control){ # remove control
    dfm <- dfm[dfm[, input.genes_annotation.control_type$control_type.var_name] == input.genes_annotation.control_type$exp_type.value, ]
  }

  if (gene_symbol.only) {
    dfm <- dfm[complete.cases(dfm[, input.genes_annotation.gene_symbol.var_name]),]
    row.lab.var_name <- input.genes_annotation.gene_symbol.var_name
  }

  ## cluster
  if (verbose) cat(paste0("Supervised correlation analysis method: ", cor.method, "\n"))
  if (verbose) cat("-------------------\n")
  cor_length <- sum(foreach(i = comparison_levels, .combine = "c") %do% length(i)) # total length for correlation analysis
  corcoef_list <- vector(mode = "list", length = cor_length)
  corp_list <- vector(mode = "list", length = cor_length)
  adj_corp_list <- vector(mode = "list", length = cor_length)
  cor_names_vector <- vector(length = cor_length)

  # cluster and heatmap
  n <- 0
  for (i in seq(length(comparisons))) {
    # set up plotting matrix for comparison i
    plt_dfm <- dfm[as.logical(thresholding_summary[[i]]), ]
    plt_mtx <- as.matrix(plt_dfm[, !names(plt_dfm) %in% names(genes)])
    # plt_mtx <- plt_mtx[, which(input.sample_groups %in% comparison_levels[[i]])]  # subsetting samples for the comparison levels
    row.lab <- plt_dfm[, row.lab.var_name]

    # correlation clustering for each group j under comparison i
    for (j in seq(length(comparison_levels[[i]]))) {
      n <- n + 1
      cor_sample_level <- comparison_levels[[i]][j]
      cor_mtx <- plt_mtx[, which(input.sample_groups %in% cor_sample_level)]
      cor_mtx <- t(cor_mtx)
      corcoef <- cor(cor_mtx, method = cor.method)
      rownames(corcoef) <- row.lab
      colnames(corcoef) <- row.lab
      corp <- foreach(m = corcoef, .combine = "cbind") %do% cor_pvalue(m, n = nrow(cor_mtx)) # p value matrix
      diag(corp) <- NA
      adj_corp <- matrix(p.adjust(corp, method = "fdr"), nrow = nrow(corp), byrow = T)  # fdr
      rownames(adj_corp) <- rownames(corp)
      colnames(adj_corp) <- colnames(corp)
      corname <- paste0(comparisons[[i]], "_", comparison_levels[[i]][j])

      if (hclust.heatmap) {
        if (verbose) cat(paste0("Correlation heatmap saved to file: ", comparisons[[i]], "_", comparison_levels[[i]][j], "_cor.heatmap.pdf..."))
        pdf(file = paste0(comparisons[[i]], "_", comparison_levels[[i]][j], "_cor.heatmap.pdf"), width = heatmap.width, height = heatmap.width)
        heatmap.2(corcoef, symm = TRUE, trace = "none",
                  col = brewer.pal(n.map.colour, map.colour),
                  distfun = function(x)as.dist(1-corcoef),
                  hclustfun = clustfunc,
                  labRow = rownames(corcoef), labCol = colnames(corcoef),
                  ...)
        dev.off()
        if (verbose) cat("Done!\n")
      }

      corcoef_list[[n]] <- corcoef
      corp_list[[n]] <- corp
      adj_corp_list[[n]] <- adj_corp
      cor_names_vector[n] <- corname
    }
  }
  names(corcoef_list) <- cor_names_vector
  names(corp_list) <- cor_names_vector
  names(adj_corp_list) <- cor_names_vector

  # sig plot
  if (sigplot) {
    if (sigplot.adj.p) {
      sigplot.corp_list <- adj_corp_list
    } else {
      sigplot.corp_list <- corp_list
    }

    if (verbose) cat("\n")
    for (i in seq(length(corcoef_list))) {
      tryCatch(
        {
          pdf(file = paste(names(corcoef_list)[i], "_cor.sigplot.pdf", sep = ""), width = sigplot.width, height = sigplot.height)
          corrplot(corr = corcoef_list[[i]], method = "color", type = "upper", p.mat = sigplot.corp_list[[i]], sig.level = sigplot.alpha,
                   insig = c("label_sig"), pch.col = sigplot.sigLabelColour, pch.cex = sigplot.sigLabelSize,
                   col = brewer.pal(n.map.colour, map.colour), mar = sigplot.mar,
                   tl.col = sigplot.labelColour, tl.cex = sigplot.labelSize, tl.srt = sigplot.labelAngle, cl.length = 3, cl.cex = sigplot.keySize)
          if (verbose) cat(paste0("Correlation significance plot saved to file: ", names(corcoef_list)[i], "_cor.sigplot.pdf..."))
          dev.off()
          if (verbose) cat("Done!\n")
        },
        error = function(err){
          if (verbose) cat(paste0("No significant correlation found for ", names(corcoef_list)[i], ". ", "Therefore no significance plot generated.\n"))
          dev.off()
        }
      )
    }
  }

  # csv export
  if (verbose) cat("\n")
  out_env <- vector(mode = "list", length = length(corcoef_list))
  for (i in seq(length(corcoef_list))) {
    out_corcoef <- corcoef_list[[i]]
    out_corp <- corp_list[[i]]
    out_adj_corp <- adj_corp_list[[i]]
    outdfm1 <- data.frame(
      group = paste(rownames(out_corcoef)[row(out_corcoef)[upper.tri(out_corcoef)]], colnames(out_corcoef)[col(out_corcoef)[upper.tri(out_corcoef)]], sep = "_"),
      row = rownames(out_corcoef)[row(out_corcoef)[upper.tri(out_corcoef)]],
      col = colnames(out_corcoef)[col(out_corcoef)[upper.tri(out_corcoef)]],
      coefficient = out_corcoef[upper.tri(out_corcoef)]
    )

    outdfm2 <- data.frame(
      group = paste(rownames(out_corp)[row(out_corp)[upper.tri(out_corp)]], colnames(out_corp)[col(out_corp)[upper.tri(out_corp)]], sep = "_"),
      row = rownames(out_corp)[row(out_corp)[upper.tri(out_corp)]],
      col = colnames(out_corp)[col(out_corp)[upper.tri(out_corp)]],
      p.value = out_corp[upper.tri(out_corp)]
    )

    outdfm3 <- data.frame(
      group = paste(rownames(out_adj_corp)[row(out_adj_corp)[upper.tri(out_adj_corp)]], colnames(out_adj_corp)[col(out_adj_corp)[upper.tri(out_adj_corp)]], sep = "_"),
      row = rownames(out_adj_corp)[row(out_adj_corp)[upper.tri(out_adj_corp)]],
      col = colnames(out_adj_corp)[col(out_adj_corp)[upper.tri(out_adj_corp)]],
      adj.p.value = out_adj_corp[upper.tri(out_adj_corp)]
    )

    out <- merge(outdfm1, outdfm2[, c("group", "p.value")], by = "group", x.all = TRUE)
    out <- merge(out, outdfm3[, c("group", "adj.p.value")], x.all = TRUE)
    out_env[[i]] <- out

    if (verbose) cat(paste0("Correlation analysis results saved to file: ", names(corcoef_list)[i], ".cor.csv..."))
    write.csv(out, file = paste0(names(corcoef_list)[i], ".cor.csv"), row.names = FALSE)
    if (verbose) cat("Done!\n")
  }
  names(out_env) <- names(corcoef_list)
  assign(paste0(deparse(substitute(object)), "_cor_sig"), out_env, envir = .GlobalEnv)
  assign(paste0(deparse(substitute(object)), "_cor_adjacency_sig"), corcoef_list, envir = .GlobalEnv)
}
