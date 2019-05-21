context("test_find_normal_branching_rate")

test_that("find_normal_branching_rate", {

  set.seed(42)
  phy <- TESS::tess.sim.taxa.age(n = 1, nTaxa = 100, age = 10, lambda = 1, mu = 0)[[1]]

  sim_regular <- function(phy, rate) {
    nodeSub::sim_normal(x = phy, l = 1000, rate = rate)
  }

  focal_rate <- 0.1

  focal_alignment <- sim_regular(phy, rate = focal_rate)


  vvv <- nodeSub::find_normal_branching_rate(phy,
                                             focal_alignment,
                                             sim_regular,
                                             replicates = 1)

  testthat::expect_equal(vvv, focal_rate, tolerance = 0.1, scale = focal_rate)

  focal_alignment <- nodeSub::sim_dual_parent(phy = phy,
                                              rate1 = focal_rate,
                                              rate2 = focal_rate * 3,
                                              node_time = 0.1)

  vvv <- nodeSub::find_normal_branching_rate(phy,
                                             focal_alignment,
                                             sim_regular,
                                             replicates = 1)

  testthat::expect_gt(vvv, focal_rate)
})

