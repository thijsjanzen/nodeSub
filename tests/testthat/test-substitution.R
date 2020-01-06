context("mutations")

test_that("mutations use", {

  levels <- c("a", "c", "g", "t")
  lbf <- length(levels)

  # default is c(0.25, 0.25, 0.25, 0.25)
  bf <- rep(1 / lbf, lbf)
  Q <- rep(1, lbf * (lbf - 1) / 2) # default is JC69      # nolint

  # only extract the 6 important rates.
  if (is.matrix(Q)) Q <- Q[lower.tri(Q)] # nolint

  eig_q <- phangorn::edQt(Q, bf) # eigen values
  m <- length(levels) # always 4 (bases)

  node_transition_matrix <- make_transition_matrix(mut_double = 0)
  eigen_obj <- eigen(node_transition_matrix, FALSE)
  eigen_obj$inv <- solve.default(eigen_obj$vec)

  parentseq <- rep("a", 100)

  found <- c()
  for (repl in 1:100) {
    rate <- 0.3
    node_time <- 1

    P <- nodeSub::slow_matrix(eig = eig_q,  # nolint
                              node_time,
                              rate = rate)

    offspr1 <- parentseq
    offspr2 <- parentseq

    for (j in 1:m) {
      ind <- parentseq == levels[j]
      offspr1[ind] <- sample(levels, sum(ind), replace = TRUE, prob = P[, j])
      offspr2[ind] <- sample(levels, sum(ind), replace = TRUE, prob = P[, j])
    }
    offspr1
    offspr2
    # total substitutions
    sub1 <- sum(offspr1 != parentseq)
    sub2 <- sum(offspr2 != parentseq)

    # and now with the nodesub method:

    p_matrix <- slow_matrix(eig = eigen_obj,
                                     node_time,
                                     rate = rate)
    p_matrix_bl <- slow_matrix(eig = eig_q,
                                        node_time,
                                        rate = rate)

    offsprings <- get_mutated_sequences(parentseq, p_matrix)

    sub3 <- sum(offsprings[[1]] != parentseq)
    sub4 <- sum(offsprings[[2]] != parentseq)
    found <- rbind(found, cbind("regular", c(sub1, sub2)))
    found <- rbind(found, cbind("nodesub", c(sub3, sub4)))
  }

  colnames(found) <- c("method", "substitutions")
  found <- tibble::as_tibble(found)
  found$substitutions <- as.numeric(found$substitutions)

  vx <- found %>%
    dplyr::group_by(method) %>%
    dplyr::summarise("mean_sub" = mean(substitutions))

  testthat::expect_equal(vx$mean_sub[1], vx$mean_sub[2] / 2,
                         tolerance = 0.5)
})


test_that("matrices", {

  levels <- c("a", "c", "g", "t")
  lbf <- length(levels)

  # default is c(0.25, 0.25, 0.25, 0.25)
  bf <- rep(1 / lbf, lbf)
  Q <- rep(1, lbf * (lbf - 1) / 2) # nolint

  # only extract the 6 important rates.
  if (is.matrix(Q)) Q <- Q[lower.tri(Q)]  # nolint

  eig_q <- phangorn::edQt(Q, bf) # eigen values

  branch_length <- 1
  rate <- 0.1

  P1 <- nodeSub::slow_matrix(eig = eig_q,  # nolint
                             branch_length,
                             rate = rate)

  P2 <- get_p_matrix(branch_length, eig_q, rate) # nolint

  testthat::expect_equal(P1, P2)
})
