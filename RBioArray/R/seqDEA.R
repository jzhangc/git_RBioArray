#' @title rbioseq_clr_ilr_transfo
#'
#' @description Log ratio tansformation function for read count data. Row: sample, column: features.
#' @param x Input read count data matrix.
#' @param offset Read count offset value added to avoid zero. Default is \code{1}.
#' @param mode Log ratio transformation method. Options are "clr" (centered log transformation) and "ilr" (isometric log transformation). Default is \code{"clr"}.
#' @param ilr.method.fast Useful only when \code{mode = "ilr"}. Default is \code{TRUE}.
#' @return A data matrix with log ratio transformed values.
#' @details This funciton is needed as part of the data pre-processing procedure to run multivariate and machine learning analysis featured in \code{RBioFS} package.
#'
#'          As per Quinn et al. (2018), NGS data can be considered as compositional data. As such, data must be transformed for usual statistical analysis and visualization.
#'          Log ratio transfromation serves such purpose. Note that the number of features will be one less when using "ilr" method.
#'
#'          It is not to be combined with the other normalization methods featured in the \code{\link{rbioseq_DE}}.
#'
#'          Ref: Quinn TP, et al. 2018. Understanding sequencing data as compositions: an outlook and review. Bioinformatics. 2018: 1 - 9.
#' @examples
#' \dontrun{
#' tstX <- rbioseq_clr_ilr_transfo(tstdata, offset = 1, mode = "clr")
#' }
#' @export
rbioseq_clr_ilr_transfo <- function(x, offset = 1, mode = "clr", ilr.method.fast = TRUE){
  # data and arguments check
  if (!is.matrix(x))stop("x needs to b e a matrix")
  if (any(x == 0) & offset == 0)stop("zero detected in x. set offset to avoid it for ratio transformation")
  if (!tolower(mode) %in% c("clr", "ilr"))stop("choose the proper transformation mode: \"clr\" or \"ilr\"")

  # log ratio transformation
  if (tolower(mode) == "clr"){  # clr calculation
    if (dim(x)[2] == 1){
      cat()
      out <- list(x.clr = x, gm = rep(1, dim(x)[1]))
    } else {
      gm <- apply(x, 1, function(x) exp(mean(log(x + offset))))  # geometric mean = exp(mean(log(X)))
      clrX <- log((x + offset) / (gm))  # clr formula
      out <- clrX
    }
  } else if (tolower(mode) == "ilr"){  # ilr calculation, modified from ilr.transfo function from mixOmics package
    cat(paste0("ilr mode result in one less variable: ", ncol(x) - 1, " variables left for this dataset upon transformation.\n"))
    ilrX = matrix(NA, nrow = nrow(x), ncol = ncol(x) - 1)
    D = ncol(x)
    if (ilr.method.fast) {
      for (i in 1:ncol(ilrX)) {
        ilrX[, i] = sqrt((D - i)/(D - i + 1)) * log(((apply(as.matrix(x[,
                                                                        (i + 1):D, drop = FALSE]), 1, prod) + offset)^(1/(D -
                                                                                                                            i)))/(x[, i] + offset))
      }
    } else {
      for (i in 1:ncol(ilrX)) {
        ilrX[, i] = sqrt((D - i)/(D - i + 1)) * log(apply(as.matrix(x[,
                                                                      (i + 1):D]), 1, function(x) {
                                                                        exp(log(x))
                                                                      })/(x[, i] + offset) + offset)
      }
    }
    out <- as.matrix(ilrX)
  }

  return(out)
}


