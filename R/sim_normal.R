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
sim_normal <- function(x,
                       l = 1000,
                       Q = NULL,
                       bf = NULL,
                       rootseq = NULL,
                       rate = 1) {


  levels <- c("a", "c", "g", "t")

  lbf <- length(levels)

  if (is.null(bf)) bf <- rep(1 / lbf, lbf)
  if (is.null(Q)) {
    Q <- rep(1, lbf * (lbf - 1) / 2)
  }
  if (is.matrix(Q)) Q <- Q[lower.tri(Q)]

  eig <- phangorn::edQt(Q, bf)

  m <- length(levels)

  if (is.null(rootseq)) rootseq <- sample(levels, l, replace = TRUE, prob = bf)
  x <- stats::reorder(x)
  edge <- x$edge
  nNodes <- max(edge)
  res <- matrix(NA, l, nNodes)
  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parent[!match(parent, child, 0)][1])
  res[, root] <- rootseq
  tl <- x$edge.length

  total_branch_subs <- 0

  for (i in seq_along(tl)) {
    from <- parent[i]
    to <- child[i]
    P <- get_p_matrix(tl[i], eig, rate)
    # avoid numerical problems for larger P and small t
    if (any(P < 0)) P[P < 0] <- 0
    for (j in 1:m) {
      ind <- res[, from] == levels[j]
      res[ind, to] <- sample(levels, sum(ind), replace = TRUE, prob = P[, j])
    }

    branch_subs <- sum(res[,from] != res[,to])
    total_branch_subs <- total_branch_subs + branch_subs
  }
  phy_no_extinct <- geiger::drop.extinct(x)

  k <- length(x$tip.label)
  label <- c(x$tip.label, as.character( (k + 1):nNodes))
  colnames(res) <- label
  res <- res[, phy_no_extinct$tip.label, drop = FALSE]
  alignment_phydat <- phyDat.DNA(as.data.frame(res, stringsAsFactors = FALSE))

  output <- list("alignment" = alignment_phydat,
                "root_seq" = rootseq,
                "total_branch_substitutions" = total_branch_subs,
                "total_node_substitutions" = 0)
 return(output)
}
