context("calc_required_node_time")

test_that("required_node_time use", {
  set.seed(666)
  yule_tree <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0,
                                     complete = FALSE)[[1]]

  # first calculate the estimated time spent at the node
  # for node time = 0.1
  node_time <- 0.1

  test_match <- function(model, lambda, mu, is_birth_death, focal_tree) {
    predicted_frac <- calc_time_spent_at_node(focal_tree,
                                              node_time = node_time,
                                              lambda,
                                              mu,
                                              is_birth_death,
                                              model)

    predicted_node_time <- calc_required_node_time(focal_tree,
                                                   fraction = predicted_frac,
                                                   lambda,
                                                   mu,
                                                   is_birth_death,
                                                   model)

    testthat::expect_equal(node_time, predicted_node_time)
  }


  # parent model, yule tree
  test_match("parent", lambda = 1, 0, is_birth_death = FALSE, yule_tree)

  # independent model, yule tree
  test_match("independent", lambda = 1, 0, is_birth_death = FALSE, yule_tree)

  # conditional model, yule tree
  test_match("conditional", lambda = 1, 0, is_birth_death = FALSE, yule_tree)



  # birth death tree
  bd_tree <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1,
                                  lambda = 1, mu = 0.3,
                                     complete = FALSE)[[1]]

  # parent model, yule tree
  test_match("parent", lambda = 1, mu = 0.3,
             is_birth_death = TRUE, bd_tree)

  # independent model, yule tree
  test_match("independent", lambda = 1, mu = 0.3,
             is_birth_death = TRUE, bd_tree)

  # conditional model, yule tree
  test_match("conditional", lambda = 1, mu = 0.3,
             is_birth_death = TRUE, bd_tree)
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
  for (repl in 1:10) {
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
  for (repl in 1:10) {
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
  for (repl in 1:2) {
    vx <- sim_dual_linked(focal_tree, rate = sub_rate,
                          node_mut_rate_double = 0,
                          l = 10000, node_time = node_time)
    found[repl] <- vx$total_node_substitutions /
      (vx$total_node_substitutions + vx$total_branch_substitutions)

  }

  testthat::expect_equal(predicted_frac, mean(found), tolerance = 0.05)
})


test_that("required_node_time use no rates", {
  test_match <- function(model, is_birth_death, focal_tree, node_time) {
    predicted_frac <-
      nodeSub::calc_time_spent_at_node(phy = focal_tree,
                                       node_time = node_time,
                                       is_birth_death = is_birth_death,
                                       model = model)

    predicted_node_time <-
      nodeSub::calc_required_node_time(phy = focal_tree,
                                       fraction = predicted_frac,
                                       is_birth_death = is_birth_death,
                                       model = model)

    testthat::expect_equal(node_time, predicted_node_time)
  }

  set.seed(666)
  yule_tree <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0,
                                    complete = FALSE)[[1]]

  # first calculate the estimated time spent at the node
  # for node time = 0.1
  node_time <- 0.1

  test_match(model = "parent", is_birth_death = FALSE,
             yule_tree, node_time)
  test_match(model = "independent", is_birth_death = FALSE,
             yule_tree, node_time)
  test_match(model = "conditional", is_birth_death = FALSE,
             yule_tree, node_time)


  bd_tree <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1,
                                  lambda = 1, mu = 0.3,
                                  complete = FALSE)[[1]]

  test_match(model = "parent", is_birth_death = TRUE,
             bd_tree, node_time)
  test_match(model = "independent", is_birth_death = TRUE,
             bd_tree, node_time)
  test_match(model = "conditional", is_birth_death = TRUE,
             bd_tree, node_time)
})

test_that("required_node_time abuse", {
  node_time <- 0.1
  yule_tree <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0,
                                    complete = FALSE)[[1]]

  testthat::expect_error(
    nodeSub::calc_time_spent_at_node(phy = yule_tree),
    "Have to provide a node time to calculate the expected fraction of \n
         time spent on the nodes."
  )
  testthat::expect_error(
    nodeSub::calc_required_node_time(phy = yule_tree),
    "Have to provide a fraction to calculate the expected \n
         time spent on the nodes."
  )

  testthat::expect_error(
    nodeSub::calc_time_spent_at_node(phy = 1, node_time = 0.1),
    "Did you forget to provide an input phylogeny?"
  )
  testthat::expect_error(
    nodeSub::calc_required_node_time(phy = 1, fraction = 0.1),
    "Did you forget to provide an input phylogeny?"
  )
})
