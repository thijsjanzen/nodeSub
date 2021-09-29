#' Simulate a sequence assuming node substitutions are not shared amongst
#' offspring, using the explicit simulation method (e.g. reverse substitutions
#' are modeled explicitly)
#' @param phy phylogenetic tree for which to simulate sequences
#' @param Q1 substitution matrix along the branches, default = JC
#' @param Q2 substitution matrix on the nodes, default = JC
#' @param rate1 mutation rate along the branch, default = 0.1
#' @param rate2 mutation rate on the node, default = 0.1
#' @param l number of base pairs to simulate
#' @param bf base frequencies, default = c(0.25, 0.25, 0.25, 0.25)
#' @param rootseq sequence at the root, simulated by default
#' @param node_time amount of time spent at the nodes
#' @return list with four items \enumerate{
#' \item{alignment} Phydat object with the resulting alignment
#' \item{rootseq} the rootsequence used
#' \item{total_branch_substitutions} total number of substitutions accumulated
#' on the branches
#' \item{total_node_substitutions} total number of substitutions accumulated at
#' the nodes}
#' @export
sim_unlinked_explicit <- function(phy,
                         Q1 = rep(1, 6),   # nolint
                         Q2 = rep(1, 6),   # nolint
                         rate1 = 0.1,
                         rate2 = 0.1,
                         l = 1000,
                         bf = rep(0.25, 4),
                         rootseq = NULL,
                         node_time = 1e-3) {

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
  if (is.matrix(Q1)) Q1 <- Q1[lower.tri(Q1)]  # nolint
  if (is.matrix(Q2)) Q2 <- Q2[lower.tri(Q2)]  # nolint
  # capital Q is retained to conform to mathematical notation on wikipedia
  # and in the literature

  phy <- stats::reorder(phy)
  edge <- phy$edge
  num_nodes <- max(edge)

  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parent[!match(parent, child, 0)][1])

  res <- matrix(NA, l, num_nodes)
  res[, root] <- rootseq
  tl <- phy$edge.length

  branch_subs_all <- rep(0, length(parent))
  node_subs_all   <- rep(0, length(parent))

  for (i in seq_along(tl)) {
    from <- parent[i]
    to <- child[i]

    # first we do substitutions due to the node model:
    P <-  p_n(100, rate2, node_time)  # get_p_matrix(tl[i], eig, rate)  # nolint

    # capital P is retained to conform to mathematical notation on wikipedia
    # and in the literature

    seq_before_mut <- res[, from]
    seq_after_mut <- mutate_seq_explicit(seq_before_mut, P)
    node_subs <- seq_after_mut$num_mut
    res[, to] <- seq_after_mut$seq

    # and then we add extra substitutions
    from <- to # the parent is now the individual again
    P <-  p_n(100, rate1, tl[i]) # nolint
    # capital P is retained to conform to mathematical notation on wikipedia
    # and in the literature

    seq_before_mut <- res[, from]
    seq_after_mut <- mutate_seq_explicit(seq_before_mut, P)
    branch_subs <- seq_after_mut$num_mut
    res[, to] <- seq_after_mut$seq

    branch_subs_all[i] <- branch_subs_all[i] + branch_subs
    node_subs_all[i]   <- node_subs_all[i] + node_subs
  }

  updated_subs <- calc_accumulated_substitutions(phy, branch_subs_all,
                                                 node_subs_all)

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
