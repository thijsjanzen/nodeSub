#' @keywords internal
calc_beta <- function(focal_tree) {
  beta <- apTreeshape::maxlik.betasplit(focal_tree)$max_lik
  return(beta)
}

#' @keywords internal
calc_gamma <- function(focal_tree) {
  gamma <- ape::gammaStat(focal_tree)
  return(gamma)
}

#' @keywords internal
calc_tree_height <- function(focal_tree) {
  return(max(ape::branching.times(focal_tree)))
}

#' @keywords internal
calc_mean_branch_length <- function(focal_tree) {
  return(mean(focal_tree$edge.length, na.rm = T))
}

#' @keywords internal
calc_all_stats <- function(focal_tree) {
  output <- c( calc_beta(focal_tree),
               calc_gamma(focal_tree),
               calc_tree_height(focal_tree),
               calc_mean_branch_length(focal_tree))
  return(output)
}

#' calculate summary statistics compared with a reference tree
#' @param trees a phyloList object containing multiple trees
#' @param true_tree a phylo object containing the reference tree
#' @param verbose verbose output if true (e.g. progressbars)
#' @return list with two tibbles 1) containing the summary statistics of all trees and 2) containing the difference with the true tree
#' @export
calc_sum_stats <- function(trees,
                           true_tree,
                           verbose = FALSE) {

  sum_stats_true_tree <- calc_all_stats(true_tree)
  sum_stats_trees <- list()
  if(!verbose) sum_stats_trees <- lapply(trees,calc_all_stats)
  if(verbose) sum_stats_trees <-  pbapply::pblapply(trees,calc_all_stats)

  all_sum_stats <- matrix(NA, nrow = length(trees), ncol = 4)
  all_differences <- matrix(NA, nrow = length(trees), ncol = 5)
  if(verbose) pb <- utils::txtProgressBar(max = length(trees), style = 3)
  for(i in 1:length(sum_stats_trees)) {
    # this for loop could be optimized later.
    to_add <- sum_stats_trees[[i]]
    local_diff <- abs(to_add - sum_stats_true_tree)
    local_nltt <- nLTT::nltt_diff_exact(true_tree, trees[[i]])
    #local_rf <- phangorn::RF.dist(true_tree, trees[[i]])
   # local_diff <- c(local_diff, local_nltt, local_rf)
    local_diff <- c(local_diff, local_nltt)

    all_differences[i, ] <- local_diff
    all_sum_stats[i, ] <- to_add
    if(verbose) utils::setTxtProgressBar(pb, i)
  }
  colnames(all_sum_stats) <- c("beta","gamma","tree_height","mean_branch_length")
  colnames(all_differences) <- c("beta","gamma","tree_height","mean_branch_length", "nLTT")

  all_sum_stats <- tibble::as_tibble(all_sum_stats)
  all_differences <- tibble::as_tibble(all_differences)

  return(list("stats" = all_sum_stats,
              "differences" = all_differences))
}