% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{InterCom}
\alias{InterCom}
\title{Reconstructs functional cell-cell communication networks}
\usage{
InterCom(
  data,
  anno.tbl,
  species,
  sighot.cutoff = 0.1,
  sighot.percentile = 70,
  consv.thrs = 0.05,
  ncores = 4,
  sig.cutoff = 0.9,
  z.score.cutoff = 2,
  min.cells = 10,
  tissue.name,
  temp.folder.name = "temp",
  out.path
)
}
\arguments{
\item{data}{Matrix or data frame of expression values; rows have to correspond to genes and columns
to cells. Row names \strong{must} be gene symbols while column names can be any identifier for a cell.
Raw counts or normalized data is permitted.}

\item{anno.tbl}{Data frame annotating each cell with a cluster/cell type. Cell names have to be provided in
the first and cluster/cell type information in the second column.}

\item{species}{The organism from which the data was obtained. Currently, only "MOUSE" and "HUMAN" are supported}

\item{sighot.cutoff}{Cutoff parameter for SigHotSpotter. Can be between 0 and 1.}

\item{sighot.percentile}{Percentile parameter for SigHotSpotter. Can be between 0 and 100.}

\item{consv.thrs}{Fraction of cells per cluster/cell type that must express a ligand, receptor or TF to be considered.}

\item{ncores}{Number of cores to use. Depending on the machine, can be any integer value greater or equal to 1.}

\item{sig.cutoff}{Significance cutoff between 0 (weakest) and 1 (strictest). Default: 0.9}

\item{z.score.cutoff}{Cutoff parameter to determine significant associations between receptors and interface TFs.
Default: 2.}

\item{min.cells}{The minimum number of cells/samples per cell type. Default: 10}

\item{tissue.name}{A name of the dataset}

\item{temp.folder.name}{Name of the temporary folder to be created within the output path. Default: "temp"}

\item{out.path}{Path to a folder where the output should be stored.}
}
\value{
List object containing the final interactome and auxiliary information
}
\description{
Reconstructs functional cell-cell communication networks for a given dataset.
}
