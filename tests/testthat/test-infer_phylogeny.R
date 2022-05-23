context("test_infer_phylogeny")

test_that("infer_phylogeny", {
  skip_on_cran()
  skip_on_ci()
  testthat::skip_on_os("mac")
  set.seed(42)
  phy  <- ape::read.tree(text = "(t1:10,(t3:2,t2:2):8);")

  seq_node_sub <- sim_normal(x = phy, l = 100,  rate = 0.1)


  testthat::expect_warning(
  all_trees <- infer_phylogeny(seq_node_sub$alignment,
                               treatment_name = "test",
                               chain_length = 1e5,
                               burnin = 0.1,
                               mcmc_seed = 42))
  mcc_tree <- all_trees$mcc_tree

  phangorn_dist <- phangorn::RF.dist(mcc_tree, phy)[[1]]
  testthat::expect_equal(phangorn_dist[[1]], 0)

  testthat::expect_equal(length(all_trees$all_trees), 1e5 / 5000)
  testthat::expect_true(class(all_trees$mcc_tree) == "phylo")
})
