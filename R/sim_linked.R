#' simulate a sequence assuming conditional substitutions on the node.
#' @param phy tree for which to simulate sequences
#' @param Q substitution matrix along the branches, default = JC
#' @param rate mutation rate , default = 1
#' @param node_mut_rate_double mutation rate on the node, default = 1e-9
#' @param l number of base pairs to simulate
#' @param bf base frequencies, default = c(0.25, 0.25, 0.25, 0.25)
#' @param rootseq sequence at the root, simulated by default
#' @param node_time time spent at the node
#' @return list with four items \enumerate{
#' \item{alignment} Phydat object with the resulting alignment
#' \item{rootseq} the rootsequence used
#' \item{total_branch_substitutions} total number of substitutions accumulated
#' on the branches
#' \item{total_node_substitutions} total number of substitutions accumulated at
#' the nodes}
#' @export
sim_linked <- function(phy,
                       Q = rep(1, 6), # nolint
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
  if (is.matrix(Q)) Q <- Q[lower.tri(Q)]  # nolint
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

  eigen_obj <- eigen(node_transition_matrix, FALSE)
  eigen_obj$inv <- solve.default(eigen_obj$vec)


  branch_subs_all <- cbind(0, 0, rep(0, max(parents)))
  node_subs_all   <- cbind(0, 0, rep(0, max(parents)))

  for (focal_parent in parents) {
    # given parent alignment
    # generate two children alignments
    offspring <- edge[which(parent == focal_parent), 2]
    # first we do substitutions due to the node model:
    p_matrix <- get_p_matrix(node_time,
                             eig = eigen_obj,
                             rate = rate)

    result <- get_mutated_sequences(res[, focal_parent],
                                    p_matrix)

    node_subs_1 <- sum(res[, focal_parent] != result[[1]])
    node_subs_2 <- sum(res[, focal_parent] != result[[2]])
    all_node_subs <- c(node_subs_1, node_subs_2)

    indices <- which(parent == focal_parent)
    for (i in 1:2) {
      branch_length <- phy$edge.length[indices[i]]
      P <- get_p_matrix(branch_length, eig_q, rate)  # nolint
      # capital P is retained to conform to mathematical notation on wikipedia
      # and in the literature

      before_mut_seq <- result[[i]] # sequence 1 after node substitutions
      after_mut_seq <- c()
      for (j in 1:m) {
        ind <- before_mut_seq == levels[j]
        after_mut_seq[ind] <- sample(levels, sum(ind), replace = TRUE,
                                     prob = P[, j])
      }
      res[, offspring[i]] <- after_mut_seq
      branch_subs <- sum(after_mut_seq != before_mut_seq)

      from <- focal_parent
      to   <- offspring[i]
      a <- which(edge[, 1] == from & edge[, 2] == to)

      branch_subs_all[a, ] <- c(focal_parent, offspring[i], branch_subs)

      node_subs_all[a, ] <- c(focal_parent, offspring[i], all_node_subs[i])
    }
  }

  branch_subs_all <- branch_subs_all[-which(branch_subs_all[, 1] == 0), ]
  node_subs_all   <- node_subs_all[-which(node_subs_all[, 1] == 0), ]

  updated_subs <- calc_accumulated_substitutions(phy,
                                                 branch_subs_all[, 3],
                                                 node_subs_all[, 3])
  phy_no_extinct <- geiger::drop.extinct(phy)

  k <- length(phy$tip.label)
  label <- c(phy$tip.label, as.character((k + 1):num_nodes))
  colnames(res) <- label
  res <- res[, phy_no_extinct$tip.label, drop = FALSE]
  alignment_phydat <- phyDat.DNA(as.data.frame(res, stringsAsFactors = FALSE))

  output <- list("alignment" = alignment_phydat,
                 "root_seq" = rootseq,
                 "total_branch_substitutions" = updated_subs$total_branch_subs,
                 "total_node_substitutions" = updated_subs$total_node_subs,
                 "total_accumulated_substitutions" =
                   updated_subs$total_accumulated_substitutions)
  return(output)
}
