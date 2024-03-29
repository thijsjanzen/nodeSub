% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_rel_contr_node_time.R
\name{calc_expected_hidden_nodes}
\alias{calc_expected_hidden_nodes}
\title{Calculate the number of expected hidden nodes in a phylogenetic tree}
\usage{
calc_expected_hidden_nodes(phy, lambda = NULL, mu = NULL)
}
\arguments{
\item{phy}{phylogenetic tree}

\item{lambda}{birth rate}

\item{mu}{death rate}
}
\value{
expected number of hidden nodes
}
\description{
Calculate the number of expected hidden nodes
using equation 1 in Manceau et al. 2020
}
\references{
Manceau, M., Marin, J., Morlon, H., & Lambert, A. (2020).
Model-based inference of punctuated molecular evolution.
Molecular Biology and Evolution, 37(11), 3308-3323.
}
