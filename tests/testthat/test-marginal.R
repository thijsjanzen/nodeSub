context("test_marginal")

test_that("marginal", {
  set.seed(42)
  phy <- TreeSim::sim.bd.taxa(n = 10, numbsim = 1, lambda = 1, mu = 0)[[1]]

  seq_normal <- sim_normal(x = phy, l = 1000,  rate = 0.01)

  skip_on_cran()
  skip_on_ci()
  test_file_name <- "test.fasta"

  if (mauricer::is_beast2_ns_pkg_installed()) {

    phangorn::write.phyDat(seq_normal$alignment,
                         file = test_file_name,
                         format = "fasta")

    ww <- nodeSub::estimate_marginal_models(fasta_filename = test_file_name,
                                          use_yule_prior = TRUE,
                                          verbose = FALSE)

    testthat::expect_true(length(ww$site_model_name) == 2)
    testthat::expect_true(sum(
        unique(ww$clock_model_name) %in% c("relaxed_log_normal",
                                           "strict")) == 2)


    file.remove(test_file_name)
  }
})
