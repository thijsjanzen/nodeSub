calc_expected_hidden_nodes_per_dt <- function(t0, t1, lambda, mu) {  # nolint
  numerator <- 1 - lambda / (mu) * exp((lambda - mu) * t0)
  denominator <- 1 - lambda / (mu) * exp((lambda - mu) * t1)
  mean_val <- 2 * lambda * (t0 - t1) - 2 * log(numerator / denominator)
  return(mean_val)
}

#' simulate substitutions along an arbitrary, imaginary, branch.
#' @param t0 starting time, assuming t1 > t0 > 0
#' @param t1 end time, assuming t1 > t0 > 0
#' @param lambda speciation rate
#' @param mu extinction rate
#' @param node_time time spent on the nodes
#' @param Q substitution matrix, default = JC
#' @param bf base frequencies, default = c(0.25, 0.25, 0.25, 0.25)
#' @param l number of base pairs to simulate
#' @param mut_rate mutation rate along the branch
#' @param node_rate mutation rate at the node
#' @return list with number of substitutions
#' @export
simulate_along_branch <- function(t0,
                                  t1,
                                  lambda,
                                  mu,
                                  node_time,
                                  Q = NULL,  # nolint,
                                  bf = NULL, # nolint
                                  l,
                                  mut_rate,
                                  node_rate) {
  levels <- c("a", "c", "g", "t")
  m <- length(levels) # always 4 (bases)

  if (is.null(bf)) bf <- rep(1 / m, m)
  if (is.null(Q)) {
    Q <- rep(1, m * (m - 1) / 2) # nolint
  }

  if (is.matrix(Q)) Q <- Q[lower.tri(Q)]  # nolint
  # capital Q is retained to conform to mathematical notation on wikipedia
  # and in the literature
  eig <- phangorn::edQt(Q, bf)

  bl <- t1 - t0
  root_seq <- sample(levels, l, replace = TRUE, prob = bf)

   P_b <- get_p_matrix(bl, eig, mut_rate)  # nolint
   P_n <- get_p_matrix(node_time, eig, node_rate)  # nolint

  # first we simulate substitutions regularly along the branch:
  old_seq <- root_seq
  new_seq <- root_seq

  for (j in 1:m) {
    ind <- root_seq == levels[j]
    new_seq[ind] <- sample(levels, sum(ind), replace = TRUE, prob = P_b[, j])
  }
  num_regular_sub <- sum(new_seq != old_seq)

  # then we do the same, using hidden nodes:
  lambda_nodes <- calc_expected_hidden_nodes_per_dt(t1, t0, lambda, mu)

  num_hidden_nodes <- 0
  if (!is.nan(lambda_nodes)) num_hidden_nodes <- stats::rpois(1, lambda_nodes)

  total_bl <- t1 - t0
  focal_t <- 0
  bl <- c()
  if (num_hidden_nodes > 0) {
    for (j in 1:num_hidden_nodes) {
      dt <- stats::rexp(1, 1 / lambda_nodes)
      while (focal_t + dt >= total_bl) {
        dt <- stats::rexp(1, 1 / lambda_nodes)
      }
      focal_t <- focal_t + dt
      bl[j] <- dt
    }
  }
  bl <- c(bl, total_bl - sum(bl))
  assertthat::assert_that(abs(sum(bl) - total_bl) < 1e-9)
  # allright, now that we have chopped everything up, we can start looping

  # branch until first node:
  P_b <- get_p_matrix(bl[1], eig, mut_rate)  # nolint
  # avoid numerical problems for larger P and small t
  before_mut_seq <- root_seq
  after_mut_seq <- root_seq
  for (j in 1:m) {
    ind <- before_mut_seq == levels[j]
    after_mut_seq[ind] <- sample(levels, sum(ind), replace = TRUE,
                                 prob = P_b[, j])
  }

  bl <- bl[-1]

  if (num_hidden_nodes > 0) {
    for (h in 1:num_hidden_nodes) {
      before_mut_seq <- after_mut_seq
      after_mut_seq <- before_mut_seq
      for (j in 1:m) {
        ind <- before_mut_seq == levels[j]
        after_mut_seq[ind] <- sample(levels, sum(ind), replace = TRUE,
                                     prob = P_n[, j])
      }

      # and the subsequent branch
      P_b <- get_p_matrix(bl[1], eig, mut_rate)  # nolint

      before_mut_seq <- after_mut_seq
      after_mut_seq <- before_mut_seq

      for (j in 1:m) {
        ind <- before_mut_seq == levels[j]
        after_mut_seq[ind] <- sample(levels, sum(ind), replace = TRUE,
                                     prob = P_b[, j])
      }
      bl <- bl[-1]
    }
  }

  num_substitutions_hidden_nodes <- sum(root_seq != after_mut_seq)

  return(list("normal_subs" = num_regular_sub,
              "node_subs" = num_substitutions_hidden_nodes))
}
