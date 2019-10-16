#' simulate a sequence assuming node substitutions are shared among the offspring
#' @param phy tree for which to simulate sequences
#' @param Q1 substitution matrix along the branches, default = JC
#' @param Q2 substitution matrix on the nodes, default = JC
#' @param rate1 mutation rate along the branch, default = 1
#' @param rate2 mutation rate on the node, default = 1
#' @param l number of base pairs to simulate
#' @param bf base frequencies, default = c(0.25, 0.25, 0.25, 0.25)
#' @param rootseq sequence at the root, simulated by default
#' @param node_time amount of time spent at the nodes
#' @return phyDat object
#' @export
sim_unlinked <- function(phy,
                         Q1 = NULL,
                         Q2 = NULL,
                         rate1 = 1,
                         rate2 = 1,
                         l = 1000,
                         bf = NULL,
                         rootseq = NULL,
                         node_time = 1e-3) {

  levels <- c("a", "c", "g", "t")
  lbf <- length(levels)

  # default is c(0.25, 0.25, 0.25, 0.25)
  if (is.null(bf)) bf <- rep(1 / lbf, lbf)

  if (is.null(Q1)) Q1 <- rep(1, lbf * (lbf - 1) / 2) # default is JC69
  if (is.null(Q2)) Q2 <- rep(1, lbf * (lbf - 1) / 2) # default is JC69

  # only extract the 6 important rates.
  if (is.matrix(Q1)) Q1 <- Q1[lower.tri(Q1)]
  if (is.matrix(Q2)) Q2 <- Q2[lower.tri(Q2)]

  eig_q1 <- phangorn::edQt(Q1, bf) # eigen values
  eig_q2 <- phangorn::edQt(Q2, bf) # eigen values

  m <- length(levels) # always 4 (bases)

  if (is.null(rootseq)) rootseq <- sample(levels, l, replace = TRUE, prob = bf)

  phy <- stats::reorder(phy)
  edge <- phy$edge
  num_nodes <- max(edge)

  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parent[!match(parent, child, 0)][1])

  res <- matrix(NA, l, num_nodes)
  res[, root] <- rootseq
  tl <- phy$edge.length

  total_node_subs <- 0
  total_branch_subs <- 0

  for (i in seq_along(tl)) {
    from <- parent[i]
    to <- child[i]

    # first we do substitutions due to the node model:
    P <- get_p_matrix(node_time, eig_q2, rate2)
    # avoid numerical problems for larger P and small t
    if (any(P < 0)) P[P < 0] <- 0
    for (j in 1:m) {
      ind <- res[, from] == levels[j]
      res[ind, to] <- sample(levels, sum(ind), replace = TRUE, prob = P[, j])
    }

    node_subs <- sum(res[, to] != res[, from])
    total_node_subs <- total_node_subs + node_subs

    # and then we add extra substitutions
    from <- to # the parent is now the individual again
    P <- get_p_matrix(tl[i], eig_q1, rate1)
    # avoid numerical problems for larger P and small t
    if (any(P < 0)) P[P < 0] <- 0
    before_mut_seq <- res[, from]
    after_mut_seq <- before_mut_seq
    for (j in 1:m) {
      ind <- before_mut_seq == levels[j]
      after_mut_seq[ind] <- sample(levels, sum(ind), replace = TRUE,
                                   prob = P[, j])
    }

    branch_subs <- sum(after_mut_seq != before_mut_seq)
    total_branch_subs <- total_branch_subs + branch_subs

    res[, to] <- after_mut_seq
  }

  phy_no_extinct <- geiger::drop.extinct(phy)

  k <- length(phy$tip.label)
  label <- c(phy$tip.label, as.character( (k + 1):num_nodes))
  colnames(res) <- label
  res <- res[, phy_no_extinct$tip.label, drop = FALSE]
  alignment_phydat <- phyDat.DNA( as.data.frame(res, stringsAsFactors = FALSE))

  return(alignment_phydat)
}
