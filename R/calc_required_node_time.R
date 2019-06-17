#' calculate the required node time to obtain a % of all substitutions
#' in the nodes
#' @description calculates the required node time
#' @param phy phylogenetic tree
#' @param fraction the targeted fraction
#' @param birth_rate the birth rate to generate the phylogenetic tree
#' @param death_Rate the death rate to generate the phylogenetic tree
#' @param is_birth_death if false, a yule tree is assumed
#' @param model the used node substitution model
#' @details if no birth or death rates are provided, birth and death rates are
#' estimated by fitting a birth-death model to the tree (or a yule tree)
#' @return required fraction
#' @export
calc_required_node_time <- function(phy,
                                    fraction,
                                    lambda = NULL,
                                    mu = NULL,
                                    is_birth_death = FALSE,
                                    model = "parent") {

  if(is.null(lambda) && !is_birth_death) {
    birth_rate <- ape::yule(phy)
  }
  if(is.null(lambda) && is.null(mu)) {
    birth_rate <- ape::birthdeath(phy)
  }

  node_time <- -1
  if(model == "parent") {
    if(!is_birth_death) {
      node_time <- fraction / (lambda - lambda * fraction)
    }
    if(is_birth_death) {
      t <- max(ape::branching.times(phy))
      tips <- length(phy$tip.label)
      nodes <- phy$Nnode
      a <- lambda / (lambda - mu) - mu / ((lambda - mu) * exp((lambda - mu) * t))
      node_time <- fraction/(1 - fraction) * tips/(mu * nodes) * log(a)
    }
  }


  if(model == "independent") {
    if(!is_birth_death) {
      node_time <- fraction / (lambda - lambda * fraction)
    }
    if(is_birth_death) {
      t <- max(ape::branching.times(phy))
      tips <- length(phy$tip.label)
      nodes <- phy$Nnode
      a <- lambda / (lambda - mu) - mu / ((lambda - mu) * exp((lambda - mu) * t))
      node_time <- fraction/(1 - fraction) * tips/(mu * nodes) * log(a)
    }
  }
  return(node_time)
}
