#' calculate the expected time spent on the node
#' in the nodes
#' @description calculates the required node time
#' @param phy phylogenetic tree
#' @param lambda the birth rate to generate the phylogenetic tree
#' @param mu the death rate to generate the phylogenetic tree
#' @param is_birth_death if false, a yule tree is
#' @param model the used node substitution model
#' @details if no birth or death rates are provided, birth and death rates are
#' estimated by fitting a birth-death model to the tree (or a yule tree)
#' @return expected fraction of time spent on the nodes
#' @export
calc_time_spent_at_node <- function(phy,
                                    node_time = NULL,
                                    lambda = NULL,
                                    mu = NULL,
                                    is_birth_death = FALSE,
                                    model = "parent") {

  if(is.null(node_time)) {
    stop("Have to provide a node time to calculate the expected fraction of \n
         time spent on the nodes.")
  }

  if(is.null(lambda) && !is_birth_death) {
    birth_rate <- ape::yule(phy)
  }
  if(is.null(lambda) && is.null(mu)) {
    birth_rate <- ape::birthdeath(phy)
  }

  rel_node_spent <- -1


  if(model == "parent") {

    if(!is_birth_death) {
      rel_node_spent <- lambda * node_time / (lambda * node_time + 1)
    }
    if(is_birth_death) {
      t <- max(ape::branching.times(phy))
      nodes <- phy$Nnode
      a <- lambda / (lambda - mu) - mu /
                                  ((lambda - mu) * exp((lambda - mu) * t))
      rel_node_spent <- mu * nodes * node_time /
                       (mu * nodes * node_time + (1 + nodes) * log(a))
    }
  }

  if(model == "independent") {
    if(!is_birth_death) {
      rel_node_spent <- 2 * lambda * node_time / (2 * lambda * node_time + 1)
    }
    if(is_birth_death) {
      t <- max(ape::branching.times(phy))
      nodes <- phy$Nnode
      a <- lambda / (lambda - mu) - mu /
                                   ((lambda - mu) * exp((lambda - mu) * t))
      rel_node_spent <- 2 * mu * nodes * node_time /
        (2 * mu * nodes * node_time + (1 + nodes) * log(a))
    }
  }

  if(rel_node_spent == -1) {
    stop("An error occurred, relative time spent at the node was not calculated")
  }

  return(rel_node_spent)
}
