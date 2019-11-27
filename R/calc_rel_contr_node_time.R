#' calculate the expected fraction of substitutions at the nodes, relative to the fraction
#' at the branches
#' @description calculates the relative contribution of substitutions at the nodes
#' @param phy phylogenetic tree (optional)
#' @param num_tips number of tips (optional)
#' @param node_time time spent at the node
#' @param model node substitution model
#' @param lambda the birth rate to generate the phylogenetic tree
#' @param mu the death rate to generate the phylogenetic tree
#' @return expected fraction
#' @export
calc_expected_fraction <- function(phy = NULL,
                                   num_tips = NULL,
                                   node_time,
                                   model = "unlinked",
                                   lambda = NULL,
                                   mu = NULL) {

  if(is.null(mu) || mu == 0) {
    # pure birth tree
    s <- 0
    if(model == "unlinked") s <- (2 * lambda * node_time) / (2 * lambda * node_time + 1)
    if(model == "linked")   s <- (lambda * node_time) / (lambda * node_time + 1)
    return(s)
  } else {
    total_bl <- 0
    num_nodes <- 0
    if(!is.null(phy)) {
      total_bl <- sum(phy$edge.length)
      num_nodes <- phy$Nnode
    } else {
      exp_bl <- (mu + (lambda - mu) * log(1 - mu/lambda) ) / (mu * mu)
      total_bl <- (num_tips - 1) * exp_bl
      num_nodes <- num_tips - 1
    }

    func_u <- function(x) {
      return(1 - (lambda/mu) * exp((lambda - mu) * x))
    }

    num_hidden_nodes <- 2 * lambda * total_bl - 2 *  log(func_u(total_bl) / func_u(0))

    num_node_subs <- 0
    if(model == "linked")   num_node_subs <-     (num_nodes + num_hidden_nodes) * node_time
    if(model == "unlinked") num_node_subs <- (2 * num_nodes + num_hidden_nodes) * node_time

    s <- num_node_subs / (num_node_subs + total_bl)
    return(s)
  }
  return(-1)
}

#' calculate the required node time to obtain a desired fraction of substitions at the node
#' @description calculates the required node time to obtain a desired fraction of susbstitutions at the node
#' @param phy phylogenetic tree (optional)
#' @param num_tips number of tips (optional)
#' @param s desired fraction
#' @param model node substitution model
#' @param lambda the birth rate to generate the phylogenetic tree
#' @param mu the death rate to generate the phylogenetic tree
#' @return expected fraction
#' @export
calc_required_node_time <- function(phy = NULL,
                                    num_tips = NULL,
                                    s,
                                    model = "unlinked",
                                    lambda = NULL,
                                    mu = NULL) {
  if(is.null(mu) || mu == 0) {
    # pure birth tree
    node_time <- -1
    if(model == "unlinked") node_time <- s / (2 * lambda - 2 * lambda * s)
    if(model == "linked")   node_time <- s / (lambda - lambda * s)
    return(node_time)
  } else {
    total_bl <- 0
    num_nodes <- 0
    if(!is.null(phy)) {
      total_bl <- sum(phy$edge.length)
      num_nodes <- phy$Nnode
    } else {
      num_nodes <- num_tips - 1
      exp_bl <- (mu + (lambda - mu) * log(1 - mu/lambda) ) / (mu * mu)
      total_bl <- 2 * num_nodes * exp_bl
    }

    func_u <- function(x) {
      return(1 - (lambda/mu) * exp((lambda - mu) * x))
    }

    num_hidden_nodes <- 2 * lambda * total_bl - 2 *  log(func_u(total_bl) / func_u(0))

    node_time <- -1
    if(model == "linked")   node_time <- (s * total_bl) / ((1-s) * (    num_nodes + num_hidden_nodes))
    if(model == "unlinked") node_time <- (s * total_bl) / ((1-s) * (2 * num_nodes + num_hidden_nodes))
    return(node_time)
  }
  return(-1)
}
