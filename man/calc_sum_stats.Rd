% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_summary_statistics.R
\name{calc_sum_stats}
\alias{calc_sum_stats}
\title{calculate summary statistics of a phylogenetic tree,
compared with a reference tree. The following statistics are calculated:
the beta statistic, gamma statistic, crown age, mean branch length,
number of tips, the nLTT statistic and the laplacian difference, given by
RPANDA's JSDtree. Because JSDtree can sometimes cause issues, some additional
checks are performed to ensure that is possible to run this function.}
\usage{
calc_sum_stats(trees, true_tree, verbose = FALSE)
}
\arguments{
\item{trees}{a phyloList object containing multiple trees}

\item{true_tree}{a phylo object containing the reference tree, preferably
without extinct lineages. If extinct lineages are found,
these are dropped.}

\item{verbose}{verbose output if true (e.g. progressbars)}
}
\value{
list with two tibbles
1) containing the summary statistics of all trees and
2) containing the difference with the true tree
}
\description{
calculate summary statistics of a phylogenetic tree,
compared with a reference tree. The following statistics are calculated:
the beta statistic, gamma statistic, crown age, mean branch length,
number of tips, the nLTT statistic and the laplacian difference, given by
RPANDA's JSDtree. Because JSDtree can sometimes cause issues, some additional
checks are performed to ensure that is possible to run this function.
}
