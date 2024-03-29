% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_positions}
\alias{plot_positions}
\title{Plot gene positions.}
\usage{
plot_positions(species_ids, gene_sets, reference_gene_ids = NULL)
}
\arguments{
\item{species_ids}{IDs of species to show in the plot.}

\item{gene_sets}{A list of gene sets (containing vectors of gene IDs) that
will be highlighted in the plot. The names will be used as labels.}

\item{reference_gene_ids}{Optionally, a set of reference genes that will be
used to reorder the species.}
}
\description{
This function requires the package \code{plotly}.
}
