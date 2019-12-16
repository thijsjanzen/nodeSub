#' calculate the expected number of hidden nodes, and draw from Poisson
#' distribution to obtain the realized number of hidden nodes
#' @param phy phylogeny
#' @param lambda birth rate
#' @param mu death rate
add_hidden_nodes <- function(phy,
                             lambda,
                             mu) {

  if(class(phy) != "phylo") {
    stop("need a phylo object as input")
  }

  if (length(geiger::is.extinct(phy)) > 0) {
    warning("this function should only be used with reconstructed trees, \n
            removing extinct branches")
    phy <- geiger::drop.extinct(phy)
  }

  bt <- ape::branching.times(phy)
  branches <- c()
  for (i in seq_along(bt)) {
    a <- as.numeric(names(bt)[[i]])
    desc <- which(phy$edge[, 1] == a)
    for (j in seq_along(desc)) {
      start_node <- phy$edge[desc[j], 1]
      end_node <- phy$edge[desc[j], 2]
      t0 <- bt[names(bt) == start_node]
      t1 <- 0
      if (end_node %in% names(bt)) {
        t1 <- bt[names(bt) == end_node]
      }

      to_add <- c(start_node, end_node, t0, t1, t0-t1)

      branches <- rbind(branches, to_add)
    }
  }

  calc_expected_hidden_nodes_per_dt <- function(t, lambda, mu) {
    t0 <- t[1]
    t1 <- t[2]
    numerator <- 1 - lambda / (mu) * exp((lambda - mu) * t0)
    denominator <- 1 - lambda / (mu) * exp((lambda - mu) * t1)
    mean_val <- 2 * lambda * (t0 - t1) - 2 * log(numerator / denominator)
    return(mean_val)
  }

  exp_hidden <- as.numeric(
                  apply(branches[,c(3,4)], 1,
                        calc_expected_hidden_nodes_per_dt,
                        lambda, mu) )

  draw_nodes <- function(x) {
    return(rpois(1, x))
  }

  obs_hidden <- sapply(exp_hidden, draw_nodes)
  return(list("lambda" = exp_hidden,
              "observed" = obs_hidden))
}

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
#' @param lambda birth rate
#' @param mu death rate
#' @return phyDat object
#' @export
sim_unlinked_tree <- function(phy,
                              Q1 = NULL,
                              Q2 = NULL,
                              rate1 = 1,
                              rate2 = 1,
                              l = 1000,
                              bf = NULL,
                              rootseq = NULL,
                              node_time = 1e-3,
                              lambda = NULL,
                              mu = NULL) {

  if (!is.null(rootseq) && length(rootseq) != l) {
    stop(
      "'rootseq' must have the same length as 'l'. \n",
      "length 'rootseq': ", length(rootseq), " \n",
      "value of 'l': ", l, " \n"
    )
  }

  if (length(geiger::is.extinct(phy)) > 0) {
    warning("this function should only be used with reconstructed trees, \n
            removing extinct branches")
    phy <- geiger::drop.extinct(phy)
  }

  if (is.null(lambda) || is.null(mu)) {
    stop('lambda and mu have to be provided')
  }

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

  hidden_nodes <- add_hidden_nodes(phy, lambda, mu )
  obs_hidden_nodes <- hidden_nodes$observed


  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parent[!match(parent, child, 0)][1])



  total_node_subs <- 0
  total_branch_subs <- 0

  phy_no_extinct <- geiger::drop.extinct(phy)

  P_n <- nodeSub::slow_matrix(eig_q2, node_time, rate2)    #get_p_matrix(node_time, eig_q2, rate2)
  if (any(P_n < 0)) P_n[P_n < 0] <- 0


  check_to_extinct_tip <- function(number) {
    if(number > length(phy$tip.label)) return(TRUE)
    if(number <= length(phy_no_extinct$tip.label)) return(TRUE)
    return(FALSE)
  }

  res <- matrix(NA, l, num_nodes)
  res[, root] <- rootseq
  tl <- phy$edge.length


  for (i in seq_along(tl)) {
    from <- parent[i]
    to <- child[i]
    node_subs <- 0
    branch_subs <- 0

    # first we do substitutions due to the node model:
    before_mut_seq <- res[, from]
    after_mut_seq <- before_mut_seq
    for (j in 1:m) {
      ind <- before_mut_seq == levels[j]
      after_mut_seq[ind] <- sample(levels, sum(ind), replace = TRUE,
                                   prob = P_n[, j])
    }
    res[, to] <- after_mut_seq
    node_subs <- node_subs + sum(res[, to] != res[, from])


    from <- to # the parent is now the individual again
    # then we add hidden node substitutions, if any
    if ( obs_hidden_nodes[i] > 0) {
      # we have to chop up the branch length, and mutate for each node
      total_bl <- tl[i]
      focal_t <- 0
      bl <- c()
      for(j in 1:obs_hidden_nodes[i]) {
        dt <- rexp(1, 1 / hidden_nodes$lambda[i])
        while(focal_t + dt >= total_bl) {
          dt <- rexp(1, 1 / hidden_nodes$lambda[i])
        }
        focal_t <- focal_t + dt
        bl[j] <- dt
      }
      bl <- c(bl, total_bl - sum(bl))
      assertthat::assert_that(abs(sum(bl) - total_bl) < 1e-9)

      # branch until first node:
      #P_b <- get_p_matrix(bl[1], eig_q1, rate1)
      P_b <-  nodeSub::slow_matrix(eig_q1, bl[1], rate1)
      # avoid numerical problems for larger P and small t
      if (any(P_b < 0)) P_b[P_b < 0] <- 0
      before_mut_seq <- res[, from]
      after_mut_seq <- before_mut_seq
      for (j in 1:m) {
        ind <- before_mut_seq == levels[j]
        after_mut_seq[ind] <- sample(levels, sum(ind), replace = TRUE,
                                     prob = P_b[, j])
      }
      branch_subs <- branch_subs + sum(after_mut_seq != before_mut_seq)
      res[, to] <- after_mut_seq

      bl <- bl[-1]

      for(h in 1:obs_hidden_nodes[i]) {
        before_mut_seq <- res[, from]
        after_mut_seq <- before_mut_seq
        for (j in 1:m) {
          ind <- before_mut_seq == levels[j]
          after_mut_seq[ind] <- sample(levels, sum(ind), replace = TRUE,
                                       prob = P_n[, j])
        }
        res[, to] <- after_mut_seq
        hidden_node_subs <- sum(res[, to] != before_mut_seq)
        node_subs <- node_subs + hidden_node_subs

        # and the subsequent branch
        #P_b <- get_p_matrix(bl[1], eig_q1, rate1)
        P_b <-  nodeSub::slow_matrix(eig_q1, bl[1], rate1)

        # avoid numerical problems for larger P and small t
        if (any(P_b < 0)) P_b[P_b < 0] <- 0

        before_mut_seq <- res[, from]
        after_mut_seq <- before_mut_seq
        for (j in 1:m) {
          ind <- before_mut_seq == levels[j]
          after_mut_seq[ind] <- sample(levels, sum(ind), replace = TRUE,
                                       prob = P_b[, j])
        }
        branch_subs <- branch_subs + sum(after_mut_seq != before_mut_seq)
        res[, to] <- after_mut_seq
        bl <- bl[-1]
      }
    } else {

      # "only" extra substitutions along the branch:
      #P <- get_p_matrix(tl[i], eig_q1, rate1)
      P <-  nodeSub::slow_matrix(eig_q1, tl[i], rate1)

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
      res[, to] <- after_mut_seq
    }

    if(check_to_extinct_tip(to)) {
      total_branch_subs <- total_branch_subs + branch_subs
      total_node_subs <- total_node_subs + node_subs
    }
  }

  k <- length(phy$tip.label)
  label <- c(phy$tip.label, as.character( (k + 1):num_nodes))
  colnames(res) <- label
  res <- res[, phy_no_extinct$tip.label, drop = FALSE]
  alignment_phydat <- phyDat.DNA( as.data.frame(res, stringsAsFactors = FALSE))

  output <- list("alignment" = alignment_phydat,
                 "root_seq" = rootseq,
                 "total_branch_substitutions" = total_branch_subs,
                 "total_node_substitutions" = total_node_subs)
  return(output)
}

#' simulate a sequence assuming node substitutions are shared among the offspring
#' optimize to obtain an equal amount of substitutions as a given alignment
#' @param input_tree tree for which to simulate sequences
#' @param focal_alignment alignment to match information content with
#' @param Q1 substitution matrix along the branches, default = JC
#' @param Q2 substitution matrix on the nodes, default = JC
#' @param rate1 mutation rate along the branch, default = 1
#' @param rate2 mutation rate on the node, default = 1
#' @param l number of base pairs to simulate
#' @param bf base frequencies, default = c(0.25, 0.25, 0.25, 0.25)
#' @param rootseq sequence at the root, simulated by default
#' @param node_time amount of time spent at the nodes
#' @param lambda birth rate
#' @param mu death rate
#' @param verbose should intermediate output be displayed?
#' @return phyDat object
#' @export
create_equal_alignment_tree <- function(input_tree,
                                       focal_alignment = NULL,
                                       Q1 = NULL,
                                       Q2 = NULL,
                                       rate1 = 1,
                                       rate2 = 1,
                                       l = 1000,
                                       bf = NULL,
                                       rootseq = NULL,
                                       fraction = 0.0,
                                       lambda = NULL,
                                       mu = NULL,
                                       verbose = FALSE) {

  if(length(geiger::is.extinct(input_tree)) > 0) {
    warning("found extinct tips, removing all extinct tips")
    input_tree <- geiger::drop.extinct(input_tree)
  }



  if(is.null(focal_alignment)) {
    warning("not found input alignment, simulating vanilla alignment")
    focal_alignment <- nodeSub::sim_normal(input_tree, l = l, Q = Q1, bf = bf, rootseq = rootseq, rate = rate1)
    rootseq <- focal_alignment$rootseq
    focal_alignment <- focal_alignment$alignment
  }

  num_emp_subs <- sum(calc_dist(focal_alignment, rootseq))

  adjusted_rate <-  rate1 + rate1 * fraction / (1 - fraction)

  node_time <- nodeSub::calc_required_node_time(input_tree, s = fraction)



  propose_alignments <- function(buffer, focal_rate) {
    proposed_alignment <- nodeSub::sim_unlinked_tree(input_tree, Q1 = Q1, Q2 = Q2,
                                                     rate1 = focal_rate, rate2 = focal_rate, l = l, bf = bf, rootseq = rootseq,
                                                     node_time = node_time, lambda = lambda, mu = mu)
    return(proposed_alignment$alignment)
  }

  calc_subs <- function(local_alignment) {
    sum(calc_dist(local_alignment, rootseq))
  }


  proposed_alignment <- propose_alignments(buffer = c(), focal_rate = adjusted_rate)

  proposed_subs <- calc_subs(proposed_alignment)
  cnt <- 1

  stored_factor <- 0
  while (proposed_subs != num_emp_subs) {

    alignments <- vector("list", 10)
    if (abs(stored_factor - 1) < 0.01) alignments <- vector("list", 100)
    alignments <- lapply(alignments, propose_alignments, adjusted_rate)
    all_subs <- unlist(lapply(alignments, calc_subs))

    num_matches <- length(which(all_subs == num_emp_subs))

    if (num_matches > 0) {
      a <- which(all_subs == num_emp_subs)[[1]]
      return(list("alignment" = alignments[[a]],
                  "rate" = adjusted_rate))
    } else {
      avg_sub <- mean(all_subs, na.rm = TRUE)
      factor <- num_emp_subs / avg_sub
      stored_factor <- factor
      adjusted_rate <- adjusted_rate * factor
    }

    cnt <- cnt + length(all_subs)
    if (verbose) cat(cnt, adjusted_rate, mean(all_subs, na.rm = TRUE),
                     num_emp_subs, factor, "\n")
  }

  return(list("alignment" = proposed_alignment,
              "rate" = adjusted_rate))
}
