context("test_marginal")

test_that("marginal", {
  skip("Not now, @thijsjanzen")
  set.seed(42)
  phy <- TreeSim::sim.bd.taxa(n = 10, numbsim = 1, lambda = 1, mu = 0)[[1]]

  seq_normal <- sim_normal(x = phy, l = 1000,  rate = 0.01)

  # skip("Not now, @thijsjanzen")
  skip_on_cran()
  phangorn::write.phyDat(seq_normal$alignment, file = "test.fasta", format = "fasta")
  ww <- nodeSub::estimate_marginal_models(fasta_filename = "test.fasta",
                                          use_yule_prior = TRUE)

  testthat::expect_true(length(ww$site_model_name) == 2)
  testthat::expect_true(sum(
      unique(ww$clock_model_name) %in% c("relaxed_log_normal", "strict")) == 2)

  testthat::expect_true(ww$weight[ww$clock_model_name == "strict"] >
                         ww$weight[ww$clock_model_name == "relaxed_log_normal"]  )
})
