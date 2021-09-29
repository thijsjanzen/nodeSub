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
calc_num_tips <- function(focal_tree) {
  return(length(focal_tree$tip.label))
}

#' @keywords internal
calc_tree_height <- function(focal_tree) {
  return(beautier::get_crown_age(focal_tree))
}

#' @keywords internal
calc_mean_branch_length <- function(focal_tree) {
  return(mean(focal_tree$edge.length, na.rm = TRUE))
}

#' @keywords internal
calc_all_stats <- function(focal_tree) {
  focal_tree <- ape::multi2di(focal_tree)
  output <- c(calc_beta(focal_tree),
              calc_gamma(focal_tree),
              calc_tree_height(focal_tree),
              calc_mean_branch_length(focal_tree),
              calc_num_tips(focal_tree))
  return(output)
}

#' calculate summary statistics of a phylogenetic tree,
#' compared with a reference tree
#' @param trees a phyloList object containing multiple trees
#' @param true_tree a phylo object containing the reference tree, preferably
#'                  without extinct lineages. If extinct lineages are found,
#'                  these are dropped.
#' @param verbose verbose output if true (e.g. progressbars)
#' @return list with two tibbles
#' 1) containing the summary statistics of all trees and
#' 2) containing the difference with the true tree
#' @export
calc_sum_stats <- function(trees,
                           true_tree,
                           verbose = FALSE) {

  if (class(trees) != "multiPhylo") {
    if (class(trees) == "phylo") {
      trees <- list(trees)
      class(trees) <- "multiPhylo"
    } else {
      stop("input needs to be correct phylo object or multiPhylo")
    }
  }

  true_tree <- ape::multi2di(true_tree)

  if (length(geiger::is.extinct(true_tree) > 0)) {
    warning("Found extinct lineages, removed these from tree\n")
    true_tree <- geiger::drop.extinct(true_tree)
  }
  sum_stats_true_tree <- calc_all_stats(true_tree)
  sum_stats_trees <- list()
  if (!verbose) sum_stats_trees <- lapply(trees, calc_all_stats)
  if (verbose) sum_stats_trees <-  pbapply::pblapply(trees, calc_all_stats)

  all_sum_stats <- matrix(NA,
                          nrow = length(trees),
                          ncol = length(sum_stats_true_tree))
  all_differences <- matrix(NA,
                            nrow = length(trees),
                            ncol = length(sum_stats_true_tree) + 2)

  if (verbose) pb <- utils::txtProgressBar(max = length(trees), style = 3)
  for (i in seq_along(sum_stats_trees)) {
    # this for loop could be optimized later.
    to_add <- sum_stats_trees[[i]]
    local_diff <- abs(to_add - sum_stats_true_tree)
    local_nltt <- nLTT::nltt_diff_exact(true_tree, trees[[i]])

    # this is rather inefficient off course,
    # RPANDA can do all pairwise in one go.
    local_jsd <- tryCatch(RPANDA::JSDtree(list(true_tree, trees[[i]]),
                                          meth = "standard")[1, 2],
                          error = NA)


    local_diff <- c(local_diff, local_nltt, local_jsd)

    all_differences[i, ] <- local_diff
    all_sum_stats[i, ] <- to_add
    if (verbose) utils::setTxtProgressBar(pb, i)
  }
  colnames(all_sum_stats) <- c("beta", "gamma", "crown_age",
                               "mean_branch_length", "num_tips")
  colnames(all_differences) <- c(colnames(all_sum_stats), "nLTT",
                                 "jsd")

  all_sum_stats <- tibble::as_tibble(all_sum_stats)
  all_differences <- tibble::as_tibble(all_differences)

  return(list("stats" = all_sum_stats,
              "differences" = all_differences))
}
