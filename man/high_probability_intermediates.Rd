% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mainlib.R
\name{high_probability_intermediates}
\alias{high_probability_intermediates}
\title{Obtain high probability intermediates}
\usage{
high_probability_intermediates(x, intermediates, percentile)
}
\arguments{
\item{x}{A data frame of gene expression values. Has to contain a column called Gene}

\item{intermediates}{A one-column data frame of intermediate signaling molecules.}

\item{percentile}{The percentile above which an intermediate molecule is considered significant (range: 0-100)}
}
\value{
A vector of high probability intermediates
}
\description{
The function computes high probability intermediates
}
