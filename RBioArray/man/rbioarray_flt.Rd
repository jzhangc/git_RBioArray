% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arrayDEA.R
\name{rbioarray_flt}
\alias{rbioarray_flt}
\title{rbioarray_flt}
\usage{
rbioarray_flt(normlst, percentile = 0.95)
}
\arguments{
\item{normlst}{Normalized data, either a list, \code{EList} or \code{MAList} object.}

\item{percentile}{The percentile cutoff, only used when muliptle negative control probes are detected. Default is \code{0.95}.}

\item{...}{arguments for \code{backgroundCorrect.matrix()} or \code{backgroundCorrect()} functions from \code{limma} package.}
}
\value{
Depending on the input type, the function outputs a \code{list}, \code{Elist} or \code{MAList} object with filtered expression values.
}
\description{
data filter function based on spike-in negative control.
}
\examples{
\dontrun{
fltdata <- rbioarray_flt(normdata)
}
}
