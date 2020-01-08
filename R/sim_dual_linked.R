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
  picked_index <- sample(1:10, 1, prob = trans_matrix[focus, ], replace = TRUE)

  output <- switch(picked_index,
         c("a", "a"),
         c("t", "t"),
         c("c", "c"),
         c("g", "g"),
         c("a", "c"),
         c("a", "t"),
         c("a", "g"),
         c("t", "c"),
         c("t", "g"),
         c("c", "g"))

  if (stats::runif(1, 0, 1) < 0.5) return(rev(output))

  return(output)
}

#' @keywords internal
pick_rate <- function(matching_bases, double_rate, focal_target) {
  output <- NA

  if (matching_bases == 2) {
    output <- NA  # e.g. A -> AA, to be filled in later
  }
  if (matching_bases == 1) {
    output <- 1 # e.g. A -> AT, one mutation
  }
  if (matching_bases == 0) { # double mutation!
    a <- substr(focal_target, 1, 1)
    b <- substr(focal_target, 2, 2)
    if (a != b) {
      output <- double_rate  # e.g. A -> TG
    } else {
      output <- 0            # e.g. A -> TT, impossible by definition
    }
    # if a == b, use default output = 0
  }
  return(output)
}

#' @keywords internal
make_transition_matrix <- function(mut_double) {

  output <- matrix(NA, nrow = 4, ncol = 10)



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

      to_add[i] <- pick_rate(matching_bases, mut_double, focal_target)
    }
    output[j, ] <- to_add
  }

  for (i in 1:4) {
    output[i, ] <- output[i, ] / sum(output[i, ], na.rm = TRUE)
  }

  for (i in 1:4) {
    output[i, i] <-  - sum(output[i, ], na.rm = TRUE)
  }
  for (j in 1:6) {
    output <- rbind(output, rep(0, length(output[1, ])))
  }

  return(output)
}

#' @keywords internal
get_mutated_sequences <- function(parent_seq, trans_matrix) {

  vx <- sapply(parent_seq, draw_bases, trans_matrix)

  child1_seq <- vx[1, ]
  child2_seq <- vx[2, ]

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
                            Q = rep(1, 6),  # nolint
                            rate = 0.1,
                            node_mut_rate_double = 1e-9,
                            l = 1000,
                            bf = rep(0.25, 4),
                            rootseq = NULL,
                            node_time = 0.01) {
  levels <- c("a", "c", "g", "t")
  if (is.null(rootseq)) {
    rootseq <- sample(levels, l, replace = TRUE, prob = bf)
  }
  if (length(rootseq) != l) {
    stop(
      "'rootseq' must have the same length as 'l'. \n",
      "length 'rootseq': ", length(rootseq), " \n",
      "value of 'l': ", l, " \n"
    )
  }



  # only extract the 6 important rates.
  if (is.matrix(Q)) Q <- Q[lower.tri(Q)] # nolint
  # capital Q is retained to conform to mathematical notation on wikipedia
  # and in the literature

  eig_q <- phangorn::edQt(Q, bf) # eigen values

  m <- length(levels) # always 4 (bases)

  phy <- stats::reorder(phy)
  edge <- phy$edge
  num_nodes <- max(edge)

  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parent[!match(parent, child, 0)][1])

  res <- matrix(NA, l, num_nodes)
  res[, root] <- rootseq

  parents <- sort(unique(as.integer(edge[, 1])))

  # the first parent should be the root, otherwise the algorithm doesn't work

  assertthat::assert_that(parents[1] == root)

  # if mu < 0, the model is undefined
  assertthat::assert_that(node_mut_rate_double >= 0)

  node_transition_matrix <- make_transition_matrix(node_mut_rate_double)
  eigen_obj <- base::eigen(node_transition_matrix, FALSE)
  eigen_obj$inv <- base::solve.default(eigen_obj$vec)

  total_node_subs   <- 0
  total_branch_subs <- 0

  phy_no_extinct <- geiger::drop.extinct(phy)

  check_to_extinct_tip <- function(number) {
    if (number > length(phy$tip.label)) return(TRUE)
    if (number <= length(phy_no_extinct$tip.label)) return(TRUE)
    return(FALSE)
  }

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

    node_subs_1 <- sum(res[, focal_parent] != result[[1]])
    node_subs_2 <- sum(res[, focal_parent] != result[[2]])
    node_subs <- c(node_subs_1, node_subs_2)

    indices <- which(parent == focal_parent)
    for (i in 1:2) {
      branch_length <- phy$edge.length[indices[i]]
      P <- get_p_matrix(branch_length, eig_q, rate)   # nolint
      # capital P is retained to conform to mathematical notation on wikipedia
      # and in the literature

      before_mut_seq <- result[[i]]
      after_mut_seq <- c()
      for (j in 1:m) {
        ind <- before_mut_seq == levels[j]
        after_mut_seq[ind] <- sample(levels, sum(ind), replace = TRUE,
                                     prob = P[, j])
      }
      res[, offspring[i]] <- after_mut_seq
      branch_subs <- sum(after_mut_seq != before_mut_seq)

      if (check_to_extinct_tip(offspring[i])) {
            total_branch_subs <- total_branch_subs + branch_subs
            total_node_subs   <- total_node_subs + node_subs[i]
      }
    }
  }

  phy_no_extinct <- geiger::drop.extinct(phy)

  k <- length(phy$tip.label)
  label <- c(phy$tip.label, as.character((k + 1):num_nodes))
  colnames(res) <- label
  res <- res[, phy_no_extinct$tip.label, drop = FALSE]
  alignment_phydat <- phyDat.DNA(as.data.frame(res, stringsAsFactors = FALSE))
  output <- list("alignment" = alignment_phydat,
                "root_seq" = rootseq,
                "total_branch_substitutions" = total_branch_subs,
                "total_node_substitutions" = total_node_subs)
  return(output)
}
