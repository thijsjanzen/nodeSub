context("test_infer_phylogeny")

test_that("infer_phylogeny", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  seq_nodeSub <- sim_normal(x = phy, l = 10000,  rate = 0.1)

  all_trees <- infer_phylogeny(seq_nodeSub,
                               treatment_name = "test",
                               burnin = 0.5,
                               chain_length = 1e6)
  mcc_tree <- all_trees$mcc_tree


  testthat::expect_lt(nLTT::nltt_diff(mcc_tree, phy), 0.01)
  testthat::expect_equal(phangorn::RF.dist(mcc_tree, phy), 0)
})



