#' Calculate the expected fraction of substitutions at the nodes,
#' relative to the fraction at the branches
#' @description calculates the relative contribution of substitutions at
#' the nodes
#' @param phy phylogenetic tree (optional)
#' @param node_time time spent at the node
#' @param model node substitution model
#' @return expected fraction
#' @export
calc_fraction <- function(phy = NULL,
                          node_time = 0,
                          model = "unlinked") {

  if (is.null(phy) || class(phy) != "phylo") {
    stop("phy needs to be a valid phylo object")
  }

  num_hidden_nodes <- nodeSub::count_hidden(phy)
  phy_no_extinct <- geiger::drop.extinct(phy)
  total_bl <- sum(phy_no_extinct$edge.length)
  num_nodes <- phy_no_extinct$Nnode

  num_node_subs <- 0
  if (model == "linked")   num_node_subs <-
                                    (num_nodes + num_hidden_nodes) * node_time
  if (model == "unlinked") num_node_subs <-
                                (2 * num_nodes + num_hidden_nodes) * node_time

  s <- num_node_subs / (num_node_subs + total_bl)
  return(s)
}


#' Calculate the required node time to obtain a desired fraction of
#' substitutions at the node
#' @description calculates the required node time to obtain a desired fraction
#' of substitutions at the node
#' @param phy phylogenetic tree
#' @param s desired fraction
#' @param model node substitution model, either "linked" or "unlinked".
#' @return expected fraction
#' @export
calc_required_node_time <- function(phy = NULL,
                                    s = 0.5,
                                    model = "unlinked") {

  if (is.null(phy) || class(phy) != "phylo") {
    stop("phy needs to be a valid phylo object")
  }

  if (!(model %in% c("linked", "unlinked"))) {
    stop("model needs to be either 'linked' or 'unlinked'")
  }

  num_hidden_nodes <- nodeSub::count_hidden(phy)
  phy_no_extinct <- geiger::drop.extinct(phy)
  total_bl <- sum(phy_no_extinct$edge.length)
  num_nodes <- phy_no_extinct$Nnode

  node_time <- -1
  if (model == "linked")  node_time <- (s * total_bl) /
                               ((1 - s) * (num_nodes + num_hidden_nodes))
  if (model == "unlinked") node_time <- (s * total_bl) /
                               ((1 - s) * (2 * num_nodes + num_hidden_nodes))
  return(node_time)
}

#' Calculate the number of expected hidden nodes in a phylogenetic tree
#' @description Calculate the number of expected hidden nodes
#' using equation 1 in Manceau et al. 2020
#' @param phy phylogenetic tree
#' @param lambda birth rate
#' @param mu death rate
#' @return expected number of hidden nodes
#' @references Manceau, M., Marin, J., Morlon, H., & Lambert, A. (2020).
#' Model-based inference of punctuated molecular evolution.
#' Molecular Biology and Evolution, 37(11), 3308-3323.
#' @export
calc_expected_hidden_nodes <- function(phy,
                                       lambda = NULL,
                                       mu = NULL) {

  if (class(phy) != "phylo") {
    stop("requires a valid phylogeny as input")
  }
  if (is.null(lambda) || is.null(mu)) {
    stop("requires valid input for lambda and mu")
  }

  if (length(geiger::is.extinct(phy)) > 0) {
    warning("can not calculate number of hidden nodes for a tree
          with extinct branches, removing extinct branches!")
    phy <- geiger::drop.extinct(phy)
  }

  if (mu == 0) return(0)

  bt <- ape::branching.times(phy)
  branches <- c()
  for (i in seq_along(bt)) {
    a <- as.numeric(names(bt)[[i]])
    desc <- which(phy$edge[, 1] == a)
    for (j in seq_along(desc)) {
      start_node <- phy$edge[desc[j], 1]
      end_node <- phy$edge[desc[j], 2]
      t0 <- bt[names(bt) == start_node]
      t1 <- 0
      if (end_node %in% names(bt)) {
        t1 <- bt[names(bt) == end_node]
      }
      branches <- rbind(branches, c(t0, t1))
    }
  }

  #now we integrate over each branch
  calc_exp_hidden_nodes_per_dt <- function(t, lambda, mu) {
    t0 <- t[1]
    t1 <- t[2]
    numerator <- 1 - lambda / (mu) * exp((lambda - mu) * t0)
    denominator <- 1 - lambda / (mu) * exp((lambda - mu) * t1)
    return(2 * lambda * (t0 - t1) - 2 * log(numerator / denominator))
  }
  exp_n <- apply(branches, 1, calc_exp_hidden_nodes_per_dt, lambda, mu)
  return(sum(exp_n))
}
