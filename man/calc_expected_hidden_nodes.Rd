% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_rel_contr_node_time.R
\name{calc_expected_hidden_nodes}
\alias{calc_expected_hidden_nodes}
\title{calculate the number of expected hidden nodes in a phylogenetic tree}
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
calculate the number of expected hidden nodes
using equation 1 in Maceau et al. 2019
}
