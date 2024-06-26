% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{genes}
\alias{genes}
\title{Information on human genes within the Ensembl database.}
\format{
A \link{data.table} with the following columns:
\describe{
\item{id}{Ensembl gene ID}
\item{name}{The gene's HGNC name (if available)}
\item{chromosome}{The human chromosome the gene is located on}
}
}
\usage{
genes
}
\description{
This includes only genes on the primary suggested assembly of the human
nuclear DNA.
}
\keyword{datasets}
