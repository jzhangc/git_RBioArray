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
