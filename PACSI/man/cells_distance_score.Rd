% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cells_distance_score.R
\name{cells_distance_score}
\alias{cells_distance_score}
\title{Title}
\usage{
cells_distance_score(
  graph,
  distance_matrix,
  cell_signature.list,
  phenotye_sample_signature.list
)
}
\arguments{
\item{graph}{PPI data}

\item{distance_matrix}{matrix calculated by igraph}

\item{cell_signature.list}{all cell signatures}

\item{phenotye_sample_signature.list}{signatures of Samples belonging to the same phenotype}
}
\value{
a list of proximity between all cells and the phenotype
}
\description{
Title
}
