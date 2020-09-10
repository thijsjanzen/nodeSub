context("test_marginal")

test_that("marginal", {
  set.seed(42)
  phy <- TreeSim::sim.bd.taxa(n = 20, numbsim = 1, lambda = 1, mu = 0)[[1]]

  seq_normal <- sim_normal(x = phy, l = 1000,  rate = 0.1)

 # skip("Takes too long")
  skip_on_cran()
  test_file_name <- "test.fasta"

  phangorn::write.phyDat(seq_normal$alignment,
                         file = test_file_name,
                         format = "fasta")
  ww <- nodeSub::estimate_marginal_models(fasta_filename = test_file_name,
                                          use_yule_prior = TRUE)

  testthat::expect_true(length(ww$site_model_name) == 2)
  testthat::expect_true(sum(
      unique(ww$clock_model_name) %in% c("relaxed_log_normal", "strict")) == 2)


  w1 <- ww$weight[ww$clock_model_name == "strict"]
  w2 <- ww$weight[ww$clock_model_name == "relaxed_log_normal"]

  testthat::expect_true( w1 > w2 )

  file.remove(test_file_name)
})
