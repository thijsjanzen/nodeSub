#' @keywords internal
get_index <- function(local_matrix, parent, offspring) {
  candidate_indices <- which(local_matrix[, 1] == parent)
  for (i in candidate_indices) {
    if (local_matrix[i, 2] == offspring) return(i)
  }
  return(-1)
}

#' @keywords internal
draw_bases <- function(focal_base, trans_matrix) {
  bases <- c("a", "t", "c", "g")

  focus <- which(focal_base == bases)
  picked_index <- sample(1:10, 1, prob = trans_matrix[focus, ], replace = T)

  output <- switch(picked_index,
         c("a","a"),
         c("t","t"),
         c("c","c"),
         c("g","g"),
         c("a", "c"),
         c("a", "t"),
         c("a", "g"),
         c("t", "c"),
         c("t", "g"),
         c("c", "g"))

  if(runif(1, 0, 1) < 0.5) return(rev(output))

  return(output)
}

#' @keywords internal
make_transition_matrix <- function(mut_double) {

  output <- matrix(NA, nrow = 4, ncol = 10)

  one_sub <- 1
  two_sub <- mut_double
  no_sub <- NA

  order_events <- c("aa", "tt", "cc", "gg",
                    "ac", "at", "ag",
                    "tc", "tg",
                    "cg")

  bases <- c("a", "t", "c", "g")
  for (j in 1:4) {
    focal <- bases[j]
    to_add <- c()
    for (i in seq_along(order_events)) {
      focal_target <- order_events[i]
      matching_bases <- stringr::str_count(focal_target, focal)
      if (matching_bases == 2) {
        to_add[i] <- no_sub
      }
      if (matching_bases == 1) {
        to_add[i] <- one_sub
      }
      if (matching_bases == 0) {
        a <- substr(focal_target, 1, 1)
        b <- substr(focal_target, 2, 2)
        if (a != b) {
          to_add[i] <- two_sub
        } else {
          to_add[i] <- 0
        }
      }
    }
    output[j, ] <- to_add
  }

  for(i in 1:4) {
    output[i, ] <- output[i, ] / sum(output[i, ], na.rm = T)
  }

  for (i in 1:4) {
    output[i, i] <-  - sum(output[i, ], na.rm = T)
  }
  for (j in 1:6) {
    output <- rbind(output, rep(0, length(output[1, ])))
  }

  return(output)
}

#' @keywords internal
get_mutated_sequences <- function(parent_seq, trans_matrix) {

  vx <- sapply(parent_seq, draw_bases, trans_matrix)

  child1_seq <- vx[1,]
  child2_seq <- vx[2,]

  return(list(child1_seq, child2_seq))
}


#' simulate a sequence assuming conditional substitutions on the node
#' @param phy tree for which to simulate sequences
#' @param Q substitution matrix along the branches, default = JC
#' @param rate mutation rate , default = 1
#' @param node_mut_rate_double mutation rate on the node, default = 1e-9
#' @param l number of base pairs to simulate
#' @param bf base frequencies, default = c(0.25, 0.25, 0.25, 0.25)
#' @param rootseq sequence at the root, simulated by default
#' @param node_time time spent at the node
#' @return phyDat object
#' @export
sim_dual_linked <- function(phy,
                            Q = NULL,
                            rate = 0.1,
                            node_mut_rate_double = 1e-9,
                            l = 1000,
                            bf = NULL,
                            rootseq = NULL,
                            node_time = 0.01) {
  levels <- c("a", "c", "g", "t")
  lbf <- length(levels)

  # default is c(0.25, 0.25, 0.25, 0.25)
  if (is.null(bf)) bf <- rep(1 / lbf, lbf)
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

  parents <- sort(unique(as.integer(edge[, 1])))

  # the first parent should be the root, otherwise the algorithm doesn't work
  testit::assert(parents[1] == root)

  # testit::assert(node_mut_rate_single >= 0) # if mu < 0, the model is undefined
  testit::assert(node_mut_rate_double >= 0) # if mu < 0, the model is undefined

  node_transition_matrix <- make_transition_matrix(node_mut_rate_double)
  eigen_obj <- eigen(node_transition_matrix, FALSE)
  eigen_obj$inv <- solve.default(eigen_obj$vec)

  total_node_subs <- 0
  total_branch_subs <- 0

  for (focal_parent in parents) {
    # given parent alignment
    # generate two children aligments
    offspring <- edge[which(parent == focal_parent), 2]
    # first we do substitutions due to the node model:
    p_matrix <- get_p_matrix(node_time,
                             eig = eigen_obj,
                             rate = rate)

    result <- get_mutated_sequences(res[, focal_parent],
                                      p_matrix)

    node_subs_1 <- sum(res[,focal_parent] != result[[1]])
    node_subs_2 <- sum(res[,focal_parent] != result[[2]])
    total_node_subs <- total_node_subs + node_subs_1 + node_subs_2

    indices <- which(parent == focal_parent)
    for (i in 1:2) {
      branch_length <- phy$edge.length[indices[i]]
      P <- get_p_matrix(branch_length, eigQ, rate)

      # avoid numerical problems for larger P and small t
      if (any(P < 0)) P[P < 0] <- 0
      before_mut_seq <- result[[i]]
      after_mut_seq <- c()
      for (j in 1:m) {
        ind <- before_mut_seq == levels[j]
        after_mut_seq[ind] <- sample(levels, sum(ind), replace = TRUE,
                                     prob = P[, j])
      }
      res[, offspring[i]] <- after_mut_seq
      branch_subs <- sum(after_mut_seq != before_mut_seq)
      total_branch_subs <- total_branch_subs + branch_subs
    }
  }

  phy_no_extinct <- geiger::drop.extinct(phy)

  k <- length(phy$tip.label)
  label <- c(phy$tip.label, as.character( (k + 1):nNodes))
  colnames(res) <- label
  res <- res[, phy_no_extinct$tip.label, drop = FALSE]
  alignment_phydat <- phyDat.DNA(as.data.frame(res, stringsAsFactors = FALSE))
  output <- list("alignment" = alignment_phydat,
                "root_seq" = rootseq,
                "total_branch_substitutions" = total_branch_subs,
                "total_node_substitutions" = total_node_subs)
  return(output)
}
