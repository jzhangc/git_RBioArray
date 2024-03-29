% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSA.R
\name{rbioGS_all}
\alias{rbioGS_all}
\title{rbioGS_all}
\usage{
rbioGS_all(
  objTitle = "DE",
  DElst,
  entrezVar = NULL,
  GS = NULL,
  GSfile = NULL,
  ...,
  export.name = NULL,
  parallelComputing = FALSE,
  clusterType = "PSOCK",
  boxplot = TRUE,
  boxplotKEGG = FALSE,
  boxplotN = 20,
  boxplotTitle = NULL,
  boxplotXlabel = "rank",
  boxplotYlabel = NULL,
  boxplotYlabelsize = 7,
  boxplotWidth = 170,
  boxplotHeight = 150,
  scatterplot = TRUE,
  scatterplotCutoff = 20,
  scatterRankline = 20,
  scatterPline = 0.05,
  scatterTitle = NULL,
  scatterXlabel = "median p value",
  scatterYlabel = "consensus score",
  scatterWidth = 170,
  scatterHeight = 150,
  plotMethod = "median",
  plotPadjust = TRUE,
  plotGSname = "GS"
)
}
\arguments{
\item{objTitle}{Object title for the output GS analysis list from \code{piano} package.}

\item{DElst}{The input list with DE result, from legacy functions \code{\link{rbioarray_DE}} or \code{\link{rbioseq_DE}}.}

\item{entrezVar}{Name of the EntrezID variable in the \code{DElst} object.}

\item{GS}{Pre-loaded gene set objects. Set only if \code{GSfile} argument is \code{NULL}. Default is \code{NULL}.}

\item{GSfile}{GS database file. Set only if \code{GS} argument is \code{NULL}. File format should be \code{gmt}. If the working directory isn't set, make sure to include the full path. Default is \code{NULL}.}

\item{...}{Arguments to pass to \code{\link{rbioGS}}.}

\item{export.name}{Test string for export filer prefix. Default is \code{NULL}.}

\item{parallelComputing}{If to use parallel computing or not. Default is \code{FALSE}}

\item{clusterType}{Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.}

\item{boxplot}{If to plot boxplots. Default is \code{TRUE}.}

\item{boxplotKEGG}{When \code{boxplot = TRUE}, to set if the gene set is KEGG. Default is \code{FALSE}.}

\item{boxplotN}{When \code{boxplot = TRUE}, to set the \code{n} (rank cutoff) argument passed to \code{consensusScores} function from \code{piano} package. Default is \code{20}.}

\item{boxplotTitle}{When \code{boxplot = TRUE}, to set the title of the boxplots. Default is \code{NULL}.}

\item{boxplotXlabel}{When \code{boxplot = TRUE}, to set the boxplots x-axis label. Default is \code{"rank"}.}

\item{boxplotYlabel}{When \code{boxplot = TRUE}, to set the boxplots y-axis label. Default is \code{NULL}.}

\item{boxplotYlabelsize}{When \code{boxplot = TRUE}, to set the boxplots y-axis label size. Default is \code{7}.}

\item{boxplotWidth}{When \code{boxplot = TRUE}, to set the boxplots width. Default is \code{170}.}

\item{boxplotHeight}{When \code{boxplot = TRUE}, to set the boxplots height. Default is \code{150}.}

\item{scatterplot}{When \code{scatterplot = TRUE}, to set if to plot a scatter plot. Default is \code{TRUE}.}

\item{scatterplotCutoff}{When \code{scatterplot = TRUE}, to set the rank cutoff of the scatter plot. Default is \code{20}.}

\item{scatterTitle}{When \code{scatterplot = TRUE}, to set the title of the scatter plot. Default is \code{NULL}.}

\item{scatterXlabel}{When \code{scatterplot = TRUE}, to set the scatter plot x-axis label. Default is \code{"median p value"}.}

\item{scatterYlabel}{When \code{scatterplot = TRUE}, to set the scatter plot y-axis label. Default is \code{"consensus score"}.}

\item{scatterWidth}{When \code{scatterplot = TRUE}, to set the scatter plot width. Default is \code{170}.}

\item{scatterHeight}{When \code{scatterplot = TRUE}, to set the scatter plot height. Default is \code{150}.}

\item{plotMethod}{When \code{boxplot = TRUE} and/or \code{scatterplot = TRUE}, to set the p methods. Options are \code{"median"} and \code{"mean"}. Default is \code{"median"}.}

\item{plotPadjust}{When \code{boxplot = TRUE} and/or \code{scatterplot = TRUE}, to set if to use FDR adjusted p value or not. Default is \code{TRUE}.}

\item{plotGSname}{When \code{boxplot = TRUE} and/or \code{scatterplot = TRUE}, to set the GS name in the file name. Default is \code{"GS"}.}

\item{scatterplotRankline}{When \code{scatterplot = TRUE}, to set the rank line on the scatter plot. Default is \code{20}.}

\item{scatterplotPline}{When \code{scatterplot = TRUE}, to set the p value line on the scatter plot. Default is \code{0.05}.}
}
\value{
Outputs  \code{csv} files and \code{pdf} figure files with GSA results.
}
\description{
All-in-one wrapper for GSA and plotting.
}
\details{
This is an all-in-one function for GS analysis based on piano package.
         It runs "fisher", "stouffer", "reporter", "tailStrength", "wilcoxon" for p value based GSA, and "page", "gsea", "maxmean" for t value based GSA (customizable).

         See arguments for \code{\link{rbioGS}} for details.

         For the new DE classes from the current version of package, such as \code{rbioseq_de} class, the DElist is the \code{DE_results} item.
}
\examples{
\dontrun{

rbioGS_all(objTitle = "mydata", DElst = DElist, entrezVar = "EntrezID", method_p = c("stouffer", "fisher"), method_t = c("gsea", "maxmean"), nPerm = 1000, GS = kegg, parallelComputing = TRUE, clusterType = "FORK")

}
}
