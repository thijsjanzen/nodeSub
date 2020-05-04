context("calc_re_contr_node_time")

test_that("calc_fraction use", {
  phy <- TreeSim::sim.bd.taxa(n = 30, numbsim = 1, lambda = 1, mu = 0)[[1]]

  for (frac in seq(0, 0.9)) {
    req_node_time <- nodeSub::calc_required_node_time(phy, s = frac)

    estim_frac <- nodeSub::calc_fraction(phy, node_time = req_node_time)
    testthat::expect_equal(frac, estim_frac)

    req_node_time <- nodeSub::calc_required_node_time(phy, s = frac,
                                                      model = "linked")

    estim_frac <- nodeSub::calc_fraction(phy, node_time = req_node_time,
                                         model = "linked")
    testthat::expect_equal(frac, estim_frac)
  }
})
