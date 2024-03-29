% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{relativizeGeneExpresion}
\alias{relativizeGeneExpresion}
\title{relativizeGeneExpresion}
\usage{
relativizeGeneExpresion(geneDataFrame, houseKeepingNames)
}
\arguments{
\item{geneDataFrame}{Dataframe produced by readGenes with averaged duplicates.}

\item{houseKeepingNames}{Vector with names of housekeepings that will be used for relativization.}
}
\value{
A dataframe with relativized gene expression levels
}
\description{
This function relativizes averaged gene expression to the geometric mean of a set of housekeepings.
}
