% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interaction.R
\name{get_receptor_sig}
\alias{get_receptor_sig}
\title{Retrieve direct signaling targets of receptors using Omnipath}
\usage{
get_receptor_sig(lr_sp_pairs_all, sp_genes, database)
}
\arguments{
\item{lr_sp_pairs_all}{A data frame containing ligand-receptor pairs}

\item{sp_genes}{A character vector of spatially expressed genes}

\item{database}{Prior knowledge of signaling pathways (Omnipath_interactions),
which is a tibble with columns named "source", "target",
"source_genesymbol", "target_genesymbol", "is_directed",
"is_stimulation"}
}
\value{
A named list where each receptor of interest is a name, and the values
are character vectors of signaling target genes.
}
\description{
Retrieve direct signaling targets of receptors using Omnipath
}
