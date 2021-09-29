#' Function to remove speciation events occuring after an extinction event.
#' Extinct species are pruned randomly, such that only a single extinct species
#' per branching event (if any extinct species) remains.
#' @param tree phylo object
#' @return pruned tree
#' @export
reduce_tree <- function(tree) {
  extant_species <- geiger::drop.extinct(tree)$tip.label
  num_extant_species <- length(extant_species)

  if (num_extant_species == length(tree$tip.label)) {
    return(tree)
  }

  extinct_species <- tree$tip.label[!(tree$tip.label %in% extant_species)]

  root_node <- min(tree$edge[, 1])
  all_nodes <- unique(tree$edge[, 1])

  for (n in all_nodes) {
    desc <- names(phylobase::descendants(phylobase::phylo4(tree),
                                         n,
                                         "children"))

    b <- desc[desc %in% extinct_species]
    if (n == root_node && length(b) == 1)  {
      tree <- ape::drop.tip(tree, b[[1]])
      return(reduce_tree(tree))
    } else {
      if (length(b) >= 2) {
        for (i in 2:length(b)) {
          tree <- ape::drop.tip(tree, b[[i]])
        }
        return(reduce_tree(tree))
      }
    }
  }
  return(tree)
}

#' Function to calculate the number of hidden speciation events, e.g. speciation
#' events that have lead to an extinct species. Thus, these hidden speciation
#' events can only be detected in complete trees (as opposed to reconstructed
#' trees).
#' @param tree phylo object
#' @return number of hidden speciation events
#' @export
count_hidden <- function(tree) {
  reduced_tree <- reduce_tree(tree)
  num_extinct_taxa <- length(geiger::is.extinct(reduced_tree))
  return(num_extinct_taxa)
}
