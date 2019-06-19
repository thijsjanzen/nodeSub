context("calc_required_node_time")

test_that("required_node_time use", {
  set.seed(666)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0,
                                     complete = FALSE)[[1]]

  # first calculate the estimated time spent at the node
  # for node time = 0.1
  node_time <- 0.1
  predicted_frac <- calc_time_spent_at_node(focal_tree,
                                            node_time = node_time,
                                            lambda = 1,
                                            mu = 0,
                                            is_birth_death = FALSE,
                                            model = "parent")

  predicted_node_time <- calc_required_node_time(focal_tree,
                                                 fraction = predicted_frac,
                                                 lambda = 1,
                                                 mu = 0,
                                                 is_birth_death = FALSE,
                                                 model = "parent")

  testthat::expect_equal(node_time, predicted_node_time)
})

test_that("required_node_time use", {
  set.seed(1)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1,
                                     mu = 0,
                                     complete = FALSE)[[1]]
  node_time <- 0.5
  predicted_frac <- calc_time_spent_at_node(focal_tree,
                                            node_time = node_time,
                                            lambda = 1,
                                            mu = 0,
                                            is_birth_death = FALSE,
                                            model = "parent")

  found <- c()
  sub_rate <- 0.001
  for(repl in 1:10) {
    vx <- sim_dual_parent(focal_tree,
                          rate1 = sub_rate,
                          rate2 = sub_rate,
                          l = 10000,
                          node_time = node_time)
    found[repl] <- vx$total_node_substitutions /
      (vx$total_node_substitutions + vx$total_branch_substitutions)

  }

  testthat::expect_equal(predicted_frac, mean(found), tolerance = 0.05)

  node_time <- 0.5
  predicted_frac <- calc_time_spent_at_node(focal_tree,
                                            node_time = node_time,
                                            lambda = 1,
                                            mu = 0,
                                            is_birth_death = FALSE,
                                            model = "independent")

  found <- c()
  sub_rate <- 0.001
  for(repl in 1:10) {
    vx <- sim_dual_independent(focal_tree,
                               rate1 = sub_rate,
                               rate2 = sub_rate,
                               l = 10000,
                               node_time = node_time)
    found[repl] <- vx$total_node_substitutions /
                  (vx$total_node_substitutions + vx$total_branch_substitutions)

  }

  testthat::expect_equal(predicted_frac, mean(found), tolerance = 0.05)

  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 10, numbsim = 1, lambda = 1, mu = 0,
                                     complete = FALSE)[[1]]
  node_time <- 0.5
  predicted_frac <- calc_time_spent_at_node(focal_tree,
                                            node_time = node_time,
                                            lambda = 1,
                                            mu = 0,
                                            is_birth_death = FALSE,
                                            model = "conditional")

  found <- c()
  sub_rate <- 0.001
  for(repl in 1:2) {
    vx <- sim_dual_linked(focal_tree, rate = sub_rate,
                          node_mut_rate_double = 0,
                          l = 10000, node_time = node_time)
    found[repl] <- vx$total_node_substitutions /
      (vx$total_node_substitutions + vx$total_branch_substitutions)

  }

  testthat::expect_equal(predicted_frac, mean(found), tolerance = 0.05)
})
