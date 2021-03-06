context("test_infer_phylogeny")

test_that("infer_phylogeny", {
  set.seed(42)
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  seq_node_sub <- sim_normal(x = phy, l = 100,  rate = 0.1)

  skip_on_cran()
  skip_on_ci()
  all_trees <- infer_phylogeny(seq_node_sub$alignment,
                               treatment_name = "test",
                               chain_length = 1e5,
                               burnin = 0.1,
                               mcmc_seed = 42)
  mcc_tree <- all_trees$mcc_tree

  phangorn_dist <- phangorn::RF.dist(mcc_tree, phy)[[1]]
  testthat::expect_equal(phangorn_dist, 0)

  testthat::expect_equal(length(all_trees$all_trees), 3)
  testthat::expect_true(class(all_trees$mcc_tree) == "phylo")
  })
