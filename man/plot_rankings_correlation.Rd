% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_rankings_correlation}
\alias{plot_rankings_correlation}
\title{Plot a scatter plot to compare two rankings.}
\usage{
plot_rankings_correlation(
  ranking_x,
  ranking_y,
  name_x,
  name_y,
  gene_sets = NULL,
  use_ranks = TRUE,
  use_sample = TRUE
)
}
\arguments{
\item{ranking_x}{The ranking to be shown on the X axis.}

\item{ranking_y}{The ranking to be shown on the Y axis.}

\item{name_x}{Title of the X axis.}

\item{name_y}{Title of the Y axis.}

\item{gene_sets}{A named list of vectors of gene IDs to highlight. The names
will be used to distinguish the sets and in the legend.}

\item{use_ranks}{Show ranks instead of scores.}

\item{use_sample}{Limit genes outside of the gene sets to a random sample.}
}
\description{
This function requires the package \code{plotly}.
}
