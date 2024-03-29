% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mainlib.R
\name{SigHotSpotter_pipeline}
\alias{SigHotSpotter_pipeline}
\title{General Pipeline Call to SigHotSpotter}
\usage{
SigHotSpotter_pipeline(
  species,
  idata,
  cutoff,
  DE_Genes,
  percentile,
  invert_DE = FALSE,
  showprogress = TRUE,
  ncores = 4
)
}
\arguments{
\item{species}{Currently supported species: "HUMAN", "MOUSE"}

\item{idata}{Data frame of input data. Rows correspond to genes, columns to cells}

\item{cutoff}{Maximum number of zero-value genes, above this cutoff the genes are excluded}

\item{DE_Genes}{(Differential) expression dataset as a one column data frame (1 for up-regulated, -1 for down-regulated genes)}

\item{percentile}{Predicted intermediates are taken into account above this threshold}

\item{invert_DE}{If the differential expression should be inverted (default: FALSE)}

\item{showprogress}{shows progress bar in shiny app if set to TRUE, set it to FALSE in batch mode without GUI (default: TRUE)}

\item{ncores}{Number of cores to use (default: 4)}
}
\value{
Compatibility scores
}
\description{
The function computes compatibility scores for signaling intermediates
}
