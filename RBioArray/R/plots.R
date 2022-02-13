#' @title rbio_venn_de
#'
#' @description Venn diagram for DE results for the \code{sig} class input object.
#' @param object Input \code{sig} class object.
#' @param gene_symbol.only Whether or not to remove probes without gene symbol. Default is \code{FALSE}.
#' @param export.name Optional name used for output objects to the environment and directory. Default is \code{NULL}.
#' @param ... arguments for \code{vennDiagram()} from \code{limma} package.
#' @param plot.Width The width of the figure for the final output figure file. Default is \code{170}.
#' @param plot.Height The height of the figure for the final output figure file. Default is \code{150}.
#' @param verbose If to turn on messaging. Default is \code{TRUE}.
#' @return The function outputs a \code{pdf} file for Venn diagrams (total, up- and down-regulations). The function also exports overlapping gene or probe into a \code{csv} file.
#' @details If exists, the F stats information will be removed from the \code{object} prior to the Venn analysis.
#' @import foreach
#' @importFrom limma vennDiagram
#' @examples
#' \dontrun{
#' rbio_venn_de(object = siglist, gene_symbol.only = TRUE)
#' }
#' @export
rbio_venn_de <- function(object, gene_symbol.only = FALSE,
                         export.name = NULL,
                         ...,
                         plot.Width = 5, plot.Height = 5,
                         verbose = TRUE){
  ## argument check
  if (any(class(object) != "sig")) stop("The input object has to be a \"sig\" class.")
  if (gene_symbol.only && ! input.genes_annotation.gene_symbol.var_name %in% names(object$input_data$genes)) {
    cat("Argument input.genes_annotation.gene_symbol.var_name not found in genes data frame when gene_symbol.only = TRUE, automatically set gene_symbol.only = FALSE.\n\n")
    gene_symbol.only <- FALSE
  }

  # check contrast levels against sample groups
  contra_levels_all <- unique(foreach(i = 1:length(object$input_data$comparisons$comparison_levels), .combine = "c") %do% {
    object$input_data$comparisons$comparison_levels[[i]]
  })
  if (!all(contra_levels_all %in% unique(levels(object$input_data$sample_groups)))) stop("Contrast levels not matching sample groups. Please check the input.")

  ## variables
  if (is.null(export.name)){
    export.name = deparse(substitute(object))
  }

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

  venn_de_list <- object$input_data$full_de_results
  stats_summary <- object$significant_change_summary
  f_stats_indices <- which(stats_summary == "F stats", arr.ind=TRUE)
  if (length(f_stats_indices) > 0) {
    if (verbose) cat("F stats inforamtion removed prior to Venn analysis.\n")
    stats_summary <- stats_summary[-f_stats_indices[1],]
    venn_de_list <- venn_de_list[!names(venn_de_list) %in% "F_stats"]
  }
  pcutoff.vector <- as.numeric(stats_summary[, "raw.p.value.threshold"])
  fc.vector <- as.numeric(complete.cases(stats_summary[, "fold.change.threshold"]))
  venn_display_var_name <- input.genes_annotation.gene_id.var_name
  if (gene_symbol.only) {
    for (i in 1:length(venn_de_list)) {
      venn_de_list[[i]] <- venn_de_list[[i]][complete.cases(venn_de_list[[i]][, input.genes_annotation.gene_symbol.var_name]), ]
    }
    if (verbose) cat("Probes without a gene symbol removed.\n\n")
    venn_display_var_name <- input.genes_annotation.gene_symbol.var_name
  }
  lfc <- matrix(nrow = nrow(venn_de_list[[1]]),
                ncol = length(venn_de_list))  # same size as the DE dataframes
  colnames(lfc) <- names(venn_de_list)
  rownames(lfc) <- venn_de_list[[1]][, venn_display_var_name]
  p <- array(NA, dim(lfc), dimnames = dimnames(lfc)) # this is another way to create a matrix. Give the value of NA.
  venn_mtx <- array(NA, dim(lfc), dimnames = dimnames(lfc))

  ## processing
  # populate the matrices
  for (i in 1:length(names(venn_de_list))){
    lfc[, i] <- venn_de_list[[i]]$logFC # extract log fold change
    p[, i] <- venn_de_list[[i]]$P.Value # extract p value (p) to a matrix
    pcutoff <- pcutoff.vector[i]
    fc <- fc.vector[i]
    # note we are using factors here e.g. "-1L, 0L, 1L". the starting value is 0L
    venn_mtx[, i] <- ifelse(p[, i] >= pcutoff | abs(lfc[, i]) < log2(fc), 0L, ifelse(lfc[, i] > 0, 1L, -1L))
  }

  ## output
  # diagrams
  if (verbose) cat("Saving Venn diagram files...")
  pdf(file = paste(export.name, "_venn_total.pdf", sep = ""), width = plot.Width, height = plot.Height)
  vennDiagram(venn_mtx, circle.col = 1:length(names(venn_de_list)), ...)
  dev.off()

  pdf(file = paste(export.name, "_venn_up.pdf", sep = ""), width = plot.Width, height = plot.Height)
  vennDiagram(venn_mtx, circle.col = 1:length(names(venn_de_list)), include = "up", ...)
  dev.off()

  pdf(file = paste(export.name, "_venn_down.pdf", sep = ""), width = plot.Width, height = plot.Height)
  vennDiagram(venn_mtx, circle.col = 1:length(names(venn_de_list)), include = "down", ...)
  dev.off()
  if (verbose) cat("Done!\n")

  # file
  if (verbose) cat("Saving Venn table...")
  if (is.null(gene_symbol.only)){
    outdfm <- data.frame(GeneSymbol = rownames(venn_mtx), venn_mtx)
    write.csv(outdfm, file = paste(export.name, "_venn_table.csv", sep = ""), row.names = FALSE)
  } else {
    outdfm <- data.frame(ProbeName = rownames(venn_mtx), venn_mtx)
    write.csv(outdfm, file = paste(export.name, "_venn_table.csv", sep = ""), row.names = FALSE)
  }
  if (verbose) cat("Done!\n")
}

