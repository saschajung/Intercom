% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mainlib.R
\name{get.gene.expr}
\alias{get.gene.expr}
\title{Get Gene Expression Subset}
\usage{
get.gene.expr(exp.tbl, genes, cell.type = NULL)
}
\arguments{
\item{exp.tbl}{A data frame of gene expression values. Row names must be genes, column names must be cell types/clusters}

\item{genes}{A character vector of genes}

\item{cell.type}{A character vector of cell types/clusters}
}
\value{
Data frame of gene expression values
}
\description{
The function derives a subset of a gene expression data frame for a given list of genes and optionally
cell types/clusters
}
