% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{densest}
\alias{densest}
\title{Find the densest value in the data.}
\usage{
densest(data)
}
\arguments{
\item{data}{The input data.}
}
\value{
The densest value of data.
}
\description{
This function assumes that data represents a continuous variable and finds
a single value with the highest estimated density. This can be used to
estimate the mode of the data. If there is only one value that value is
returned. If multiple density maxima with the same density exist, their mean
is returned.
}
