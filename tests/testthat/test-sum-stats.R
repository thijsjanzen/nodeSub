context("calc_sum_stats")

test_that("calc_sum_stats", {

  phy1 <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0)[[1]]
  brts <- ape::branching.times(phy1)

  phy <- nodeSub::create_balanced_tree(brts)

  input <- list(phy, phy)
  class(input) <- "multiPhylo"

  testthat::expect_true(class(input) == "multiPhylo")
  testthat::expect_true(length(input) == 2)

  stats1 <- nodeSub::calc_sum_stats(input, phy)
  stats2 <- nodeSub::calc_sum_stats(input, phy, verbose = TRUE)

  testthat::expect_true(stats1$stats$beta[[1]] >= 9.9)
  testthat::expect_true(stats1$stats$beta[[2]] >= 9.9)

  testthat::expect_true(sum(stats1$differences, na.rm = TRUE) == 0)

  stats2 <- nodeSub::calc_sum_stats(input, phy1)
  testthat::expect_true(stats2$stats$beta[[1]] >= 9.9)
  testthat::expect_true(stats2$stats$beta[[2]] >= 9.9)

  testthat::expect_true(sum(stats2$differences, na.rm = TRUE) != 0)

  testthat::expect_true(all.equal(stats1$stats, stats2$stats))


  phy <- nodeSub::create_unbalanced_tree(brts)

  stats1 <- nodeSub::calc_sum_stats(phy, phy)

  testthat::expect_true(class(input) == "multiPhylo")
  testthat::expect_true(length(input) == 2)

  testthat::expect_true(stats1$stats$beta[[1]] < 0.0)

  testthat::expect_true(sum(stats1$differences, na.rm = TRUE) == 0)


  phy1 <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0.5)[[1]]
  phy2 <- geiger::drop.extinct(phy1)


  testthat::expect_warning(
    nodeSub::calc_sum_stats(phy2, phy1)
  )

})

test_that("calc_sum_stats abuse", {
  phy1 <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0)

  testthat::expect_error(nodeSub::calc_sum_stats(phy1, phy1))
})
