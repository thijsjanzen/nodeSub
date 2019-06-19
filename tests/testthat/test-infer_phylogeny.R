context("test_infer_phylogeny")

test_that("infer_phylogeny", {
  set.seed(42)
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  seq_node_sub <- sim_normal(x = phy, l = 100,  rate = 0.1)

  all_trees <- infer_phylogeny(seq_node_sub$alignment,
                               treatment_name = "test",
                               burnin = 0.1)
  mcc_tree <- all_trees$mcc_tree

  testthat::expect_lt(nLTT::nltt_diff(mcc_tree, phy), 0.05)
  phangorn_dist <- phangorn::RF.dist(mcc_tree, phy)[[1]]
  testthat::expect_equal(phangorn_dist, 0)
})
