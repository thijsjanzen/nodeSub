context("count_hidden and reduce_tree")

test_that("count_hidden use", {

  for (mu in seq(0, 0.5, 0.1)) {
    test_tree <- TreeSim::sim.bd.taxa(n = 30, numbsim = 1,
                                    lambda = 1, mu,
                                    complete = TRUE)[[1]]

    reduced_tree <- nodeSub::reduce_tree(test_tree)

    num_hidden <- nodeSub::count_hidden(test_tree)

    num_extinct <- length(geiger::is.extinct(reduced_tree))
    testthat::expect_equal(num_hidden, num_extinct)
  }
})

test_that("count_hidden expectation", {
  set.seed(42)
  for (mu in seq(0, 0.5, 0.1)) {
    test_tree <- TreeSim::sim.bd.taxa(n = 50, numbsim = 10,
                                      lambda = 1, mu,
                                      complete = TRUE)




    num_hidden <- lapply(test_tree, nodeSub::count_hidden)
    avg_num_hidden <- mean(unlist(num_hidden))

    test_tree2 <- lapply(test_tree, geiger::drop.extinct)

    exp_hidden <- lapply(test_tree2, calc_expected_hidden_nodes,
                         lambda = 1, mu = mu)

    avg_exp_hidden <- mean(unlist(exp_hidden))

    testthat::expect_equal(avg_num_hidden, avg_exp_hidden, tolerance = 2)
  }
})
