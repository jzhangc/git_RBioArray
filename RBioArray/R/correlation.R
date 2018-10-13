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


#' @title rbio_supervised_corcluster
#'
#' @description Supervised correlation clustering analysis and heatmap visualization for both microarray and RNAseq, for gene co-expression analysis.
#' @param object Input \code{sig} class object.
#' @param gene_symbol.only If to only use probes/genes/genomic features with a gene sybmol id. Default is \code{FALSE}.
#' @param cor.method The correlation method, options are "pearson", "spearman" and "pearson". Default is \code{"pearson"}.
#' @param map.colour Heat map colour. Default is \code{"PRGn"}. See \code{RColorBrewer} package for more.
#' @param n.map.colour Number of colours displayed. Default is \code{11}. See \code{RColorBrewer} package for more.
#' @param heatmap.axis.label Whether to display label for both x- and y- axes. Default is \code{FALSE}.
#' @param ... Additional arguments for \code{heatmap.2} function from \code{gplots} package.
#' @param heatmap.width Width for correlation heatmap. Unit is \code{inch}. Default is \code{7}.
#' @param heatmap.height Height for correlation heatmap. Unit is \code{inch}. Default is \code{7}.
#' @param sigplot.adj.p If to use FDR corrected correlation p value to plot sigPlot. Default is \code{FALSE}.
#' @param sigplot.alpha The alpha value for correlation p value. Default is \code{0.05}
#' @param sigplot.sigLabelColour The colour for label for the significant pairs. Default is \code{"red"}.
#' @param sigplot.sigLabelSize The size for label for the significant pairs. Default is \code{3}.
#' @param sigplot.labelColour The significance heatmap axis label colour. Default is \code{"black"}.
#' @param sigplot.labelSize The significance heatmap axis label size. Default is \code{1}.
#' @param sigplot.labelAngle The significance heatmap axis label angle. Default is \code{90}.
#' @param sigplot.keySize The significance heatmap colour key size. Default is \code{1}.
#' @param sigplot.width Width for correlation significance plot. Unit is \code{inch}. Default is \code{7}.
#' @param sigplot.height Height for correlation significance plot. Unit is \code{inch}. Default is \code{7}.
#' @details The function takes \code{sig} class object. It behaves in the following ways:
#'          1. heatmap uses \code{heatmap.2} function from \code{gplots} package.
#'          2. significance plot uses \code{corrplot} package
#'          3. the function will use normalized and filtered expression data for both RNAseq and micaorray
#'          4. the function will automatically subset data using \code{thresholding_summary} from \code{sig} class input
#'          5. the function will automatically subset and plot for each group under each comparison
#' @return A supervised heatmap based on hierarchical clustering analysis in \code{pdf} format, along with correaltion coefficient and p value matrices. The function also outputs a significance plot.
#' @import corrplot
#' @import foreach
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' rbio_supervised_corcluster(object = sig_input, gene_symbol.only = TRUE, cor.method = "pearson",
#'                            map.colour = "PRGn", n.map.colour = 11,
#'                            heatmap.axis.label = TRUE, heatmap.width = 7, heatmap.height = 7,
#'                            sigplot.adj.p = FALSE, sigplot.alpha = 0.05, sigplot.sigLabelColour = "red", sigplot.sigLabelSize = 3,
#'                            sigplot.labelColour = "black", sigplot.labelSize = 1, sigplot.labelAngle = 90, sigplot.keySize = 1,
#'                            sigplot.width = 7, sigplot.height = 7)
#' }
#' @export
rbio_supervised_corcluster <- function(object,
                                       gene_symbol.only = FALSE,
                                       cor.method = "pearson",
                                       map.colour = "PRGn", n.map.colour = 11, ...,
                                       heatmap.axis.label = FALSE, heatmap.width = 7, heatmap.height = 7,
                                       sigplot.adj.p = FALSE, sigplot.alpha = 0.05,
                                       sigplot.sigLabelColour = "red", sigplot.sigLabelSize = 3,
                                       sigplot.labelColour = "black", sigplot.labelSize = 1, sigplot.labelAngle = 90, sigplot.keySize = 1,
                                       sigplot.width = 7, sigplot.height = 7){
  ## argument check
  if (class(object) != "sig") stop("The input object has to be a \"sig\" class.")
  if (n.map.colour %% 1 != 0) stop("Argument n.map.colour needs to be an integer number.")
  if (!cor.method %in% c("pearson", "spearman")) stop("Argument cor.method needs to be \"pearson\" or \"spearman\".")

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
  cat(paste0("Supervised correlation analysis method: ", cor.method, "\n"))
  cat("-------------------\n")
  cor_length = sum(foreach(i = comparison_levels, .combine = "c") %do% length(i)) # total length for correlation analysis
  corcoef_list <- vector(mode = "list", length = cor_length)
  corp_list <- vector(mode = "list", length = cor_length)
  adj_corp_list <- vector(mode = "list", length = cor_length)
  cor_names_vector <- vector(length = cor_length)
  n = 0
  # cluster and heatmap
  for (i in seq(length(comparisons))) {
    # set up plotting matrix for comparison i
    plt_dfm <- dfm[as.logical(thresholding_summary[[i]]), ]
    plt_mtx <- as.matrix(plt_dfm[, !names(plt_dfm) %in% names(genes)])
    # plt_mtx <- plt_mtx[, which(input.sample_groups %in% comparison_levels[[i]])]  # subsetting samples for the comparison levels
    row.lab <- plt_dfm[, row.lab.var_name]

    # correlation clustering for each group j under comparison i
    for (j in seq(length(comparison_levels[[i]]))) {
      n = n + 1
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

      cat(paste0("Correlation heatmap saved to file: ", comparisons[[i]], "_", comparison_levels[[i]][j], "_cor.heatmap.pdf..."))
      pdf(file = paste0(comparisons[[i]], "_", comparison_levels[[i]][j], "_cor.heatmap.pdf"), width = heatmap.width, height = heatmap.width)
      heatmap.2(corcoef, symm = TRUE, trace = "none",
                col = brewer.pal(n.map.colour, map.colour), labRow = rownames(corcoef), labCol = colnames(corcoef),
                ...)
      dev.off()
      cat("Done!\n")

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
  if (sigplot.adj.p) {
    sigplot.corp_list <- adj_corp_list
  } else {
    sigplot.corp_list <- corp_list
  }

  cat("\n")
  for (i in seq(length(corcoef_list))) {
    tryCatch(
      {
        pdf(file = paste(names(corcoef_list)[i], "_cor.sigplot.pdf", sep = ""), width = sigplot.width, height = sigplot.height)
        corrplot(corr = corcoef_list[[i]], method = "color", type = "upper", p.mat = sigplot.corp_list[[i]], sig.level = sigplot.alpha,
                 insig = c("label_sig"), pch.col = sigplot.sigLabelColour, pch.cex = sigplot.sigLabelSize,
                 col = brewer.pal(n.map.colour, map.colour),
                 tl.col = sigplot.labelColour, tl.cex = sigplot.labelSize, tl.srt = sigplot.labelAngle, cl.length = 3, cl.cex = sigplot.keySize)
        cat(paste0("Correlation significance plot saved to file: ", names(corcoef_list)[i], "_cor.sig.pdf..."))
        dev.off()
        cat("Done!\n")
      },
      error = function(err){
        cat(paste0("No significant correlation found for ", names(corcoef_list)[i], ". ", "Therefore no significance plot generated.\n"))
        dev.off()
      }
    )
  }

  # csv export
  cat("\n")
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

    cat(paste0("Correlation analysis results saved to file: ", names(corcoef_list)[i], ".cor.csv..."))
    write.csv(out, file = paste0(names(corcoef_list)[i], ".cor.csv"), row.names = FALSE)
    cat("Done!\n")
  }
}
