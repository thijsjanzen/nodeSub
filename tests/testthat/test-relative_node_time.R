context("calc_required_node_time")

test_that("required_node_time use", {
  set.seed(1)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1,
                                     mu = 0,
                                     complete = FALSE)[[1]]

  chosen_fraction <- 0.1
  req_node_time <- nodeSub::calc_required_node_time(focal_tree,
                                                    s = chosen_fraction)

  obs_fraction <- nodeSub::calc_fraction(focal_tree, node_time = req_node_time)

  testthat::expect_equal(obs_fraction, 0.1)

  ali <- nodeSub::sim_dual_independent(focal_tree, rate1 = 0.001, rate2 = 0.001,
                                       l = 10000, node_time = req_node_time)

  obs_fraction <- ali$total_node_substitutions /
          (ali$total_node_substitutions + ali$total_branch_substitutions)

  testthat::expect_equal(obs_fraction, chosen_fraction, tol = 0.01)

})
