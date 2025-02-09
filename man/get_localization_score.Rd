% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvo.R
\name{get_localization_score}
\alias{get_localization_score}
\title{Calculate Localization Score for Genes Based on ST data}
\usage{
get_localization_score(
  st_data,
  location,
  th_vol,
  data_struct = "deldir",
  deldir.nb = 3,
  dist.th = 2
)
}
\arguments{
\item{st_data}{A matrix or data frame containing spatial transcriptomics data
with genes as rows and spatial spots as columns.}

\item{location}{A matrix or data frame of coordinates, with each row representing
a point and two columns for x and y coordinates.}

\item{th_vol}{A threshold for the volume ratio to filter genes based on expression
level. Default is 0.025.}

\item{data_struct}{A character string specifying the structure to use.
Options are "deldir" (default) for Delaunay triangulation
or "dist" for distance threshold.}

\item{deldir.nb}{An integer specifying the threshold for the number of neighbors to
consider in the Gabriel graph construction. Defaults to 3.}

\item{dist.th}{A numeric threshold for the distance used when data_struct is set
to "dist". Defaults to 2.}
}
\value{
A named numeric vector where each element is the calculated localization
score for a gene.
}
\description{
Calculate Localization Score for Genes Based on ST data
}
