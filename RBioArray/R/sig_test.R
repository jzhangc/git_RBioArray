#' Title sig
#'
#' @description Generic significance test function
#' @param object Object containing DE inforamtion. Currently the function supports \code{rbioseq_de} and \code{rbioarray_de} objects.
#' @param ... Additional arguments for corresponding S3 class methods.
#' @return A differential expression result object.
#' @examples
#'
#' \dontrun{
#' sig(object = mrnade)
#' }
#'
#' @export
sig <- function(object, ...){
  ## check object
  if (!class(object) %in% c("rbioseq_de", "rbioarray_de")) stop("object needs to be either a \"rbioseq_de\" or \"rbioarray_de\" object")

  ## use methods
  UseMethod("sig", object)
}


#' Title sig.rbioseq_de
#'
#' @description The \code{sig} function for \code{rbioseq_de} object from \code{\link{rnaseq_de}} function.
#' @param object A \code{rbioseq_de} object \code{\link{rnaseq_de}} function.
#' @param export.name Optional name used for output objects to the environment and directory. Default is \code{NULL}.
#' @param ... Additional arguments for \code{\link{sig.defuault}}.
#'
#' @export
sig.rbioseq_de <- function(object, export.name = NULL, ...){
  if (is.null(export.name)){
    export.name = deparse(substitute(object))
  }

  sig.default(input.de.list = object$DE_results, input.gene_symbol.var.name = object$gene_symbol_var_name,
              input.Fstats.matrix = object$F_stats, experiment = "RNAseq", export.name = export.name, ...)
}


#' Title sig.default
#'
#' @description The default \code{sig} function.
#' @param input.de.list  Input list cantaining DE dataframes for each comparison.
#' @param input.gene_symbol.var.name Input gene sysmbol variable name from teh DE dataframes.
#' @param input.Fstats.matrix Input dataframe containing F stats
#' @param experiment Character string describing the experiment used to generate data, e.g. "microarray", "RNAseq".
#' @param FC Threshold for fold change. Default is \code{1.5}.
#' @param alpha Threshold for p values. Default is \code{0.05}.
#' @param FDR If to use FDR correction for p values. Default is \code{TRUE}.
#' @param gene_symbol If to apply gene symbols in the plot and exported results. Default is \code{TRUE}.
#' @param export.name Name used for output objects to the environment and directory. Not optional. Default is \code{NULL}.
#' @param export.mode Mode used to export results to the directory. Options are \code{"all"}, \code{"all.gene_symbol"}, \code{"sig"}. Default is \code{"all"}. See details.
#' @param plot If to plot volcano plot. Default is \code{TRUE}.
#' @param plot.top.gene If to display the top gene identification, i.e., probem name or gene name, on the plot. Default is \code{FALSE}.
#' @param plot.top.gene.n When \code{plot.top.gene = TRUE}, to set how many genes to display. Default is \code{5}.
#' @param plot.top.gene.padding When \code{plot.top.gene = TRUE}, to set the distance between the dot and the gene symbol. Default is \code{0.5}.
#' @param plot.Title Figure title. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param plot.xLabel X-axis label. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param plot.yLabel Y-axis label. Make sure to use quotatio marks. Use \code{NULL} to hide. Default is \code{"Mean Decrease in Accurac"}
#' @param plot.symbolSize Size of the symbol. Default is \code{2}.
#' @param plot.sigColour Colour of the significant genes or probes. Default is \code{"red"}.
#' @param plot.nonsigColour Colour of the non-significant genes or probes. Default is \code{"gray"}.
#' @param plot.xTxtSize Font size for the x-axis text. Default is \code{10}.
#' @param plot.yTxtSize Font size for the y-axis text. Default is \code{10}.
#' @param plot.Width The width of the figure for the final output figure file. Default is \code{170}.
#' @param plot.Height The height of the figure for the final output figure file. Default is \code{150}.
#' @param genesymbolVar The name of the variable for gene symbols from the \code{annot} object. Only set this argument when \code{geneName = TRUE}. Default is \code{NULL}.
#' @return Signifiance test results as \code{csv} files to the working directory, as well as a \code{sig} object to the environment.
#'         The \code{sig} object contains the following items;
#'
#'         \code{significant_change_summary}
#'
#'         \code{export.mode}
#'
#'         \code{signifianct_change_results}: this only contains significantly changed genes.
#'
#'         \code{experiment}
#'
#' @details Explanation for \code{export.mode} options:
#'
#'         \code{all}: export all results for all probes/features with or without name (e.g. gene symbol, gene name, etc.) annotations.
#'
#'         \code{all.gene_symbol}: export all results for probes/features only with name (e.g. gene symbol, gene name, etc.) annotations.
#'
#'         \code{sig}: export only the signifiant changes. Gene symbol settings depends on argument \code{gene_symbol}.
#'
#' @import ggplot2
#' @importFrom RBioplot rightside_y
#' @importFrom grid grid.newpage grid.draw
#' @importFrom ggrepel geom_text_repel
#'
#' @export
sig.default <- function(input.de.list, input.gene_symbol.var.name, input.Fstats.matrix,
                        experiment = NULL,
                        FC = 1.5, alpha = 0.05, FDR = TRUE, gene_symbol = TRUE,
                        export.name = NULL, export.mode = "all",
                        plot = TRUE,
                        plot.top.gene = FALSE, plot.top.gene.n = 5,  plot.top.gene.padding = 0.5,
                        plot.Title = NULL,  plot.xLabel = "log2(fold change)", plot.yLabel = "-log10(p value)",
                        plot.symbolSize = 2, plot.sigColour = "red", plot.nonsigColour = "gray",
                        plot.xTxtSize = 10, plot.yTxtSize = 10, plot.Width = 170, plot.Height = 150){
  ## check arguments
  if (is.null(export.name)) stop("export.name is needed.")
  if (!export.mode %in% c("all", "all.gene_symbol", "sig", "sig.gene_sybmol")) {
    stop("export.mode should be one of \"all\", \"all.gene_symbol\", \"sig\", and \"sig.gene_symbol\"")
  }
  if (is.null(names(input.de.list))) names(input.de.list) <- seq(length(input.de.list))

  ## sig test
  # set up data.frames
  de_list <- vector(mode = "list", length = length(input.de.list))
  de_list[] <- foreach(i = seq(length(input.de.list))) %do% {
    # set up DE dfm
    de_dfm <- input.de.list[[i]]
    if (gene_symbol){  # if to display gene name, then subset
      de_dfm <- de_dfm[complete.cases(de_dfm[, input.gene_symbol.var.name]), ]
    }
    de_dfm
  }
  names(de_list) <- names(input.de.list)

  # set the cutoff and summary matrix
  sig_summary_mtx <- matrix(nrow = length(input.de.list), ncol = 7)  # initiate summary matrix
  colnames(sig_summary_mtx) <- c("comparisons", "raw.p.value.threshold", "fold.change.threshold", "alpha", "FDR", "True", "False")

  pcutoff_vector <- vector(length = length(input.de.list))

  cutoff_list <- vector(mode = "list", length = length(input.de.list))
  names(cutoff_list) <- names(input.de.list)

  for (i in seq(length(input.de.list))) {
    sig_dfm <- de_list[[i]]
    # cut off
    if (FDR){
      if (length(which(sig_dfm$adj.P.Val < alpha)) == 0){
        cat(paste0("No FDR corrected p-values found less than alpha for the comparison: ", names(input.de.list)[i],
                   ". \nPlease consider using another thresholding method. For now, alpha is applied to the raw p.values.\n"))
        pcutoff <- alpha
        fdr.stats <- FALSE
      } else {
        pcutoff <- max(sig_dfm[sig_dfm$adj.P.Val < alpha, ]$P.Value)
        fdr.stats <- TRUE
      }
      cutoff <- as.factor(abs(sig_dfm$logFC) >= log2(FC) & sig_dfm$P.Value < pcutoff)
    } else  {
      pcutoff <- alpha
      fdr.stats <- FALSE
      cutoff <- as.factor(abs(sig_dfm$logFC) >= log2(FC) & sig_dfm$P.Value < pcutoff)
    }
    cutoff_list[[i]] <- cutoff

    # store pcutoffs
    pcutoff_vector[i] <- pcutoff

    # summary matrix
    if (length(levels(cutoff)) == 1){
      if (levels(cutoff) == "TRUE"){
        sig.summary <- c(names(input.de.list)[i], signif(pcutoff, digits = 4), FC, alpha, fdr.stats, summary(cutoff)[["TRUE"]], 0)
      } else {
        sig.summary <- c(names(input.de.list)[i], signif(pcutoff, digits = 4), FC, alpha, fdr.stats, 0, summary(cutoff)[["FALSE"]])
      }
    } else {
      sig.summary <- c(names(input.de.list)[i], signif(pcutoff, digits = 4), FC, alpha, fdr.stats, summary(cutoff)[["TRUE"]], summary(cutoff)[["FALSE"]])
    }
    sig_summary_mtx[i, ] <- sig.summary
  }

  ## plot
  if (plot){
    plt_list <- vector(mode = "list", length = length(input.de.list))
    plt_list[] <- foreach(i = seq(length(input.de.list))) %do% {
      plt_dfm <- de_list[[i]]
      plt_cutoff <- cutoff_list[[i]]
      plt_pcutoff <- pcutoff_vector[i]

      loclEnv <- environment()
      plt <- ggplot(plt_dfm, aes(x = logFC, y = -log10(P.Value)), environment = loclEnv) +
        geom_point(alpha = 0.4, size = plot.symbolSize , aes(colour = plt_cutoff)) +
        scale_color_manual(values = c(plot.nonsigColour, plot.sigColour)) +
        ggtitle(plot.Title) +
        scale_y_continuous(expand = c(0.02, 0)) +
        xlab(plot.xLabel) +
        ylab(plot.yLabel) +
        geom_vline(xintercept = log2(FC), linetype = "dashed") +
        geom_vline(xintercept = - log2(FC), linetype = "dashed") +
        geom_hline(yintercept = - log10(plt_pcutoff), linetype = "dashed") +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(hjust = 0.5),
              legend.position = "none",
              legend.title = element_blank(),
              axis.text.x = element_text(size = plot.xTxtSize),
              axis.text.y = element_text(size = plot.yTxtSize, hjust = 0.5))
      if (plot.top.gene){
        plt_fltdfm <- plt_dfm[abs(plt_dfm$logFC) >= log2(FC) & plt_dfm$P.Value < plt_pcutoff, ]
        plt_fltdfm <- plt_fltdfm[order(plt_fltdfm$P.Value), ]
        plt <- plt + geom_text_repel(data = head(plt_fltdfm, n = plot.top.gene.n),
                                     aes(x = logFC, y = -log10(P.Value), label = head(plt_fltdfm, n = plot.top.gene.n)[, input.gene_symbol.var.name]),
                                     point.padding = unit(plot.top.gene.padding, "lines"))
      }
      plt <- RBioplot::rightside_y(plt)
      plt
    }
    names(plt_list) <- names(input.de.list)

    cat("\n")
    for (i in seq(length(input.de.list))) {
      cat(paste0("Saving volcano plots to file: ", names(input.de.list)[i], ".volcano.pdf..."))
      grid.newpage()
      ggsave(filename = paste(names(input.de.list)[i],".volcano.pdf", sep = ""), plot = plt_list[[i]],
             width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
      grid.draw(plt_list[[i]])
      cat("Done!\n")
    }
  }

  ## output
  sig_out_list <- vector(mode = "list", length = length(input.de.list))
  sig_out_list[] <- foreach(i = seq(length(input.de.list))) %do% {
    out_dfm <- de_list[[i]][as.logical(cutoff_list[[i]]), ]
    out_dfm
  }
  names(sig_out_list) <- names(input.de.list)

  out <- list(significant_change_summary = sig_summary_mtx,
              export.mode = export.mode,
              significant_change_results = sig_out_list,
              experiment = experiment)
  class(out) <- "sig"
  assign(paste0(export.name, "_sig"), out, envir = .GlobalEnv)

  ## export
  cat("\n")
  cat(paste0("Exporting DE F-stats to file: ", export.name, "_Fstats.csv..."))
  write.csv(file = paste0(export.name, "_Fstats.csv"), input.Fstats.matrix, row.names = FALSE)
  cat("Done!\n")
  cat(paste0("Saving signifiance test summary to file: ", export.name, "_sig_summary.csv..."))
  write.csv(file = paste0(export.name, "_sig_summary.csv"), sig_summary_mtx, row.names = FALSE)
  cat("Done!\n")
  if (export.mode == "all") {
    for (i in seq(length(input.de.list))) {  # use the input DE dataframes
      cat(paste0("Exporting all DE results for all genes/probes/features to file: ", export.name, "_", names(input.de.list)[i], "_de_all.csv..."))
      write.csv(file = paste0(export.name, "_", names(input.de.list)[i], "_de_all.csv"), input.de.list[[i]], row.names = FALSE)
      cat("Done!\n")
    }
  } else if (export.mode == "all.gene_symbol") {  # use the subsetted DE dataframes
    for (i in seq(length(input.de.list))) {
      cat(paste0("Exporting all results for genes/probes/features with a symbol to file: ", export.name, "_", names(de_list)[i], "_de_gene_symbol.csv..."))
      write.csv(file = paste0(export.name, "_", names(de_list)[i], "_de_gene_symbol.csv"), de_list[[i]], row.names = FALSE)
      cat("Done!\n")
    }
  } else if (export.mode == "sig") {
    for (i in seq(length(input.de.list))) {
      cat(paste0("Exporting only the significantly changes to file: ", export.name, "_", names(de_list)[i], ifelse(gene_symbol, "_sig_gene_symbol.csv...", "_sig.csv...")))
      write.csv(file = paste0(export.name, "_", names(de_list)[i], ifelse(gene_symbol, "_sig_gene_symbol.csv", "_sig.csv")), sig_out_list[[i]], row.names = FALSE)
      cat("Done!\n")
    }
  }
}


#' @export
print.sig <- function(x, ...){
  cat("\n")
  cat(paste0("--- ", x$experiment, " significant test results ---\n"))
  print(as.data.frame(x$significant_change_summary, row.names = NULL))
  cat("\n")
}
