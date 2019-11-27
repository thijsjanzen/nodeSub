context("calc_required_node_time")

test_that("required_node_time use", {
  set.seed(1)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1,
                                     mu = 0,
                                     complete = FALSE)[[1]]
  node_time <- 0.5
  predicted_frac <- calc_expected_fraction(focal_tree,
                                           node_time = node_time,
                                           model = "linked",
                                           lambda = 1,
                                           mu = 0)

  found <- c()
  sub_rate <- 0.01
  for (repl in 1:10) {
    vx <- sim_dual_linked(focal_tree,
                          rate = sub_rate,
                          node_mut_rate_double = sub_rate,
                          l = 1000,
                          node_time = node_time)
    found[repl] <- vx$total_node_substitutions /
      (vx$total_node_substitutions + vx$total_branch_substitutions)

  }

  testthat::expect_equal(predicted_frac, mean(found), tolerance = 0.05)

  node_time <- 0.5
  predicted_frac <- calc_expected_fraction(focal_tree,
                                           node_time = node_time,
                                           model = "unlinked",
                                           lambda = 1,
                                           mu = 0)

  found <- c()
  sub_rate <- 0.01
  for (repl in 1:10) {
    vx <- sim_dual_independent(focal_tree,
                               rate1 = sub_rate,
                               rate2 = sub_rate,
                               l = 1000,
                               node_time = node_time)
    found[repl] <- vx$total_node_substitutions /
                  (vx$total_node_substitutions + vx$total_branch_substitutions)

  }

  testthat::expect_equal(predicted_frac, mean(found), tolerance = 0.05)
})
