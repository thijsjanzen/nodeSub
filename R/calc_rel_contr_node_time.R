#' calculate the fraction of substitutions at the nodes, relative to the fraction
#' at the branches
#' @description calculates the relative contribution of substitutions at the nodes
#' @param phy phylogenetic tree
#' @param node_time time spent at the node
#' @param lambda the birth rate to generate the phylogenetic tree
#' @param mu the death rate to generate the phylogenetic tree
#' @return required fraction
#' @export
calc_rel_contr_node_time <- function(phy,
                                     node_time,
                                     lambda = NULL,
                                     mu = NULL) {

  if(is.null(lambda) || is.null(mu)) {
    return(NA)
  }

  calc_u <- function(t) {
    return( 1 - (lambda / mu) * exp((lambda - mu) * t))
  }

  n <- phy$Nnode

  expected_branch_length <- (mu + (lambda - mu) * log(1 - mu / lambda)) /
                            (mu * mu)

  expected_hidden_per_bl <- 2 * lambda * expected_branch_length -
                            2 * log(  calc_u(expected_branch_length) / calc_u(0))

  total_normal <- n * expected_branch_length
  total_node   <- (n + expected_hidden_per_bl) * node_time

  fraction <- total_node / (total_node + total_normal)

  return(fraction)
}

#' calculate the fraction of substitutions at the nodes, relative to the fraction
#' at the branches
#' @description calculates the relative contribution of substitutions at the nodes
#' @param phy phylogenetic tree
#' @param node_time time spent at the node
#' @param lambda the birth rate to generate the phylogenetic tree
#' @param mu the death rate to generate the phylogenetic tree
#' @return list with things
#' @export
estim_rel_contr_node_time <- function(phy,
                                     node_time,
                                     lambda = NULL,
                                     mu = NULL) {
  if(is.null(lambda) || is.null(mu)) {
    return(NA)
  }

  calc_u <- function(t) {
    return( 1 - (lambda / mu) * exp((lambda - mu) * t))
  }

  purged_phy <- geiger::drop.extinct(phy)

  total_num_hidden_nodes <- 0
  for(tt in purged_phy$edge.length) {
    total_num_hidden_nodes <- total_num_hidden_nodes +
                              2 * lambda * tt -
                                  2 * log(  calc_u(tt) / calc_u(0))
  }


  total_normal <- sum(purged_phy$edge.length)
  total_node   <- (purged_phy$Nnode + total_num_hidden_nodes) * node_time

  exp_fraction <- total_node / (total_node + total_normal)


  obs_hidden_nodes <- nodeSub::count_hidden(phy)
  obs_total_node   <- (purged_phy$Nnode + obs_hidden_nodes) * node_time

  obs_fraction <- obs_total_node / (obs_total_node + total_normal)

  output <- list("total_bl" = total_normal,
                 "exp_num_hidden_nodes" = total_num_hidden_nodes,
                 "exp_fraction" = exp_fraction,
                 "obs_num_hidden_nodes" = obs_hidden_nodes,
                 "obs_fraction" = obs_fraction)

  return(output)
}
