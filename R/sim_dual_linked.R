#' @keywords internal
make_transition_matrix <- function(mu) {
  output <- matrix(NA, nrow = 4, ncol = 10)

  one_sub <- mu * (1-mu)
  two_sub <- (mu^2)
  no_sub <- (1 - mu)^2

  order_events <- c("AA", "TT","CC","GG","AC","AT","AG","TC","TG","CG")

  #for(focal in c("A","C","T","G") {
  bases <- c("A","T","C","G")
  for(j in 1:4) {
    focal <- bases[j]
    to_add <- c()
    for(i in 1:length(order_events)) {
      focal_target <- order_events[i]
      matching_bases <- stringr::str_count(focal_target, focal)
      if(matching_bases == 2) {
        to_add[i] <- no_sub
      }
      if(matching_bases == 1) {
        to_add[i] <- one_sub
      }
      if(matching_bases == 0) {
        a <- substr(focal_target, 1, 1)
        b <- substr(focal_target, 2, 2)
        if(a != b) {
          to_add[i] <- two_sub
        } else {
          to_add[i] <- 0
        }
      }
    }
    to_add <- to_add / sum(to_add)
    output[j, ] <- to_add
  }

  return(output)
}

#' @keywords internal
get_index <- function(local_matrix, parent, offspring) {
  candidate_indices <- which(local_matrix[,1] == parent)
  for(i in candidate_indices) {
    if(local_matrix[i,2] == offspring) return(i)
  }
  return(-1)
}

#' @keywords internal
draw_bases <- function(focal_base, trans_matrix) {
  bases <- c("a","t","c","g")
  output_table <- matrix(nrow = 10, ncol = 2)
  output_table[1,] <- c("a","a")
  output_table[2,] <- c("t","t")
  output_table[3,] <- c("c","c")
  output_table[4,] <- c("g","g")
  output_table[5,] <- c("a","c")
  output_table[6,] <- c("a","t")
  output_table[7,] <- c("a","g")
  output_table[8,] <- c("t","c")
  output_table[9,] <- c("t","g")
  output_table[10, ] <- c("c","g")

  focus <- which(focal_base == bases)
  output_bases <- sample(1:10, 1, prob = trans_matrix[focus, ], replace = T)
  return(output_table[output_bases,])
}


#' @keywords internal
get_mutated_sequences <- function(parent_seq, trans_matrix) {

  child1_seq <- parent_seq
  child2_seq <- parent_seq

  for(i in 1:length(parent_seq)) {
    bases <- draw_bases(parent_seq[i], trans_matrix)
    child1_seq[i] <- bases[[1]]
    child2_seq[i] <- bases[[2]]
  }
  return(list(child1_seq, child2_seq))
}

#' simulate a sequence assuming conditional substitutions on the node
#' @param phy tree for which to simulate sequences
#' @param Q substitution matrix along the branches, default = JC
#' @param rate mutation rate , default = 1
#' @param mu mutation rate on the node, default = 1e-9
#' @param l number of base pairs to simulate
#' @param bf base frequencies, default = c(0.25, 0.25, 0.25, 0.25)
#' @param rootseq sequence at the root, simulated by default
#' @return phyDat object
#' @export
sim_dual_subs_linked <- function(phy,
                                 Q = NULL,
                                 rate = 1,
                                 mu = 1e-9,
                                 l = 1000,
                                 bf = NULL,
                                 rootseq = NULL) {
  levels <- c("a", "c", "g", "t")
  lbf <- length(levels)
  if (is.null(bf)) bf <- rep(1 / lbf, lbf) # default is c(0.25, 0.25, 0.25, 0.25)
  if (is.null(Q)) Q <- rep(1, lbf * (lbf - 1) / 2) # default is JC69

  # only extract the 6 important rates.
  if (is.matrix(Q)) Q <- Q[lower.tri(Q)]

  eigQ <- phangorn::edQt(Q, bf) # eigen values

  m <- length(levels) # always 4 (bases)

  if (is.null(rootseq)) rootseq <- sample(levels, l, replace = TRUE, prob = bf)

  phy <- stats::reorder(phy)
  edge <- phy$edge
  nNodes <- max(edge)

  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parent[!match(parent, child, 0)][1])

  res <- matrix(NA, l, nNodes)
  res[, root] <- rootseq

  parents <- sort(unique(as.integer(edge[,1])))
  testit::assert(parents[1] == root) # the first parent should be the root, otherwise the algorithm doesn't work
  testit::assert(mu >= 0) # if mu < 0, the model is undefined
  node_transition_matrix <-  make_transition_matrix(mu)

  for(focal_parent in parents) {
    # given parent alignment
    # generate two children aligments
    offspring = edge[which(parent == focal_parent), 2]
    # first we do substitutions due to the node model:
    result <- get_mutated_sequences(res[,focal_parent], node_transition_matrix)

    indices <- which(parent == focal_parent)
    for(i in 1:2) {
      branch_length <- phy$edge.length[indices[i]]
      P <- getP(branch_length, eigQ, rate)[[1]]

      # avoid numerical problems for larger P and small t
      if (any(P < 0)) P[P < 0] <- 0
      before_mut_seq <- result[[i]]
      after_mut_seq <- c()
      for (j in 1:m) {
        ind <- before_mut_seq == levels[j]
        after_mut_seq[ind] <- sample(levels, sum(ind), replace = TRUE, prob = P[, j])
      }
      res[, offspring[i] ] <- after_mut_seq
    }
  }

  k <- length(phy$tip.label)
  label <- c(phy$tip.label, as.character( (k + 1):nNodes))
  colnames(res) <- label
  res <- res[, phy$tip.label, drop = FALSE]
  alignment_phydat <- phyDat.DNA(as.data.frame(res, stringsAsFactors = FALSE), return.index = TRUE)
  return(alignment_phydat)
}
