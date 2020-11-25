#' @keywords internal
ab_n <- function(ab_nmin1, n) {
  if (n == 0) {
    ab_n <- c(1, 0)
  } else {
    ab_n <- c(3 * ab_nmin1[2] / n, (2 * ab_nmin1[2] + ab_nmin1[1]) / n)
  }
  return(ab_n)
}

#' @keywords internal
ab <- function(n) {
  ab_mat <- matrix(0,
                   nrow = n + 1,
                   ncol = 2)
  ab_mat[1, ] <- c(1, 0)
  for (i in 2:(n + 1)) {
    ab_mat[i, ] <- ab_n(ab_mat[i - 1, ], i - 1)
  }
  return(ab_mat)
}

#' @keywords internal
p_n <- function(n, mu, t) {
  ab_mat <- ab(n)
  p <- ab_mat * (mu * t) ^ (0:n) * exp(-3 * mu * t)
  return(p)
}

#' @keywords internal
mutate_seq_explicit <- function(local_sequence, pn) {
  bases <- c("a", "c", "g", "t")
  seq_before_mut <- local_sequence
  seq_after_mut  <- seq_before_mut
  num_mut <- 0
  for (j in 1:4) {
    ind <- which(seq_before_mut == bases[j])
    a <- sample(x = 1:length(pn),
                size = length(ind),
                replace = T,
                prob = pn)

    chosen_col <- rep(1, length(a))
    chosen_col[a > nrow(pn)] <- 2
    b <- which(chosen_col == 1)
    num_mut <- num_mut + sum(a[b] - 1)
    b <- which(chosen_col == 2)
    num_mut <- num_mut + sum((a[b] - nrow(pn) - 1), na.rm = TRUE)

    if (length(b) > 0) {
        other_bases <- bases[-which(bases == bases[j])]
        chosen_bases <- sample(other_bases, size = length(b), replace = TRUE)
        mutated_bases <- ind[b]
        seq_after_mut[mutated_bases] <- chosen_bases
    }
  }
  return(list("seq" = seq_after_mut,
              "num_mut" = num_mut))
}



#' Simulate sequences.
#'
#' Simulate sequences for a given evolutionary tree.
#'
#' \code{simSeq} is now a generic function to simulate sequence alignments to
#' along a phylogeny. It
#' is quite flexible and allows to generate DNA, RNA, amino acids, codon  or
#' binary sequences.  It is possible to give a \code{pml} object as input simSeq
#' return a \code{phyDat} from these model.  There is also a more low level
#' version, which lacks rate variation, but one can combine different
#' alignments having their own rate (see example). The rate parameter acts like
#' a scaler for the edge lengths.
#'
#' can be supplied.
#'
#' @param x a phylogenetic tree \code{tree}, i.e. an object of class
#' \code{phylo} or and object of class \code{pml}.
#' @param l length of the sequence to simulate.
#' @param Q the rate matrix.
#' @param bf base frequencies.
#' @param rootseq a vector of length l containing the root sequence, other root
#' sequence is randomly generated.
#' @param rate mutation rate or scaler for the edge length, a numerical value
#' greater than zero.
#' @return \code{simSeq} returns an object of class phyDat.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @keywords cluster
#' @export
sim_normal_explicit <- function(x,
                                l = 1000,
                                Q = NULL,   # nolint
                                bf = NULL,
                                rootseq = NULL,
                                rate = 1) {


  levels <- c("a", "c", "g", "t")

  lbf <- length(levels)

  if (is.null(bf)) bf <- rep(1 / lbf, lbf)
  if (is.null(Q)) {
    Q <- rep(1, lbf * (lbf - 1) / 2) # nolint
  }
  if (is.matrix(Q)) Q <- Q[lower.tri(Q)]  # nolint
  # capital Q is retained to conform to mathematical notation on wikipedia
  # and in the literature

  if (is.null(rootseq)) rootseq <- sample(levels, l, replace = TRUE, prob = bf)
  x <- stats::reorder(x)
  edge <- x$edge
  num_nodes <- max(edge)
  res <- matrix(NA, l, num_nodes)
  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parent[!match(parent, child, 0)][1])
  res[, root] <- rootseq
  tl <- x$edge.length

  total_branch_subs <- 0
  daughter_subs <- rep(0, length(parent))

  for (i in seq_along(tl)) {
    from <- parent[i]
    to <- child[i]
    P <-  p_n(100, rate, tl[i])  # get_p_matrix(tl[i], eig, rate)  # nolint
    # capital P is retained to conform to mathematical notation on wikipedia
    # and in the literature

    seq_before_mut <- res[, from]
    seq_after_mut <- mutate_seq_explicit(seq_before_mut, P)

    res[, to] <- seq_after_mut$seq

    branch_subs <- sum(res[, from] != res[, to])
    total_branch_subs <- total_branch_subs + branch_subs
    daughter_subs[i] <- seq_after_mut$num_mut
  }

  # now, given the daughter subs string, we need to calculate the total
  # accumulated divergence
  updated_subs <- calc_accumulated_substitutions(x, daughter_subs)

  phy_no_extinct <- geiger::drop.extinct(x)

  k <- length(x$tip.label)
  label <- c(x$tip.label, as.character((k + 1):num_nodes))
  colnames(res) <- label
  res <- res[, phy_no_extinct$tip.label, drop = FALSE]
  alignment_phydat <- phyDat.DNA(as.data.frame(res, stringsAsFactors = FALSE))

  total_inferred_substitutions <- sum(calc_dist(alignment_phydat, rootseq))

  output <- list("alignment" = alignment_phydat,
                 "root_seq" = rootseq,
                 "total_branch_substitutions" = updated_subs$total_branch_subs,
                 "total_node_substitutions" = updated_subs$total_node_subs,
                 "total_inferred_substitutions" = total_inferred_substitutions,
                 "total_accumulated_substitutions" =
                   updated_subs$total_accumulated_substitutions)
          #"incl_reverse_subs" = reverse_subs$total_accumulated_substitutions)

  return(output)
}
