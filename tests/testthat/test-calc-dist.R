context("calc_dist")

test_that("calc_dist use", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")
  seq_node_sub <- nodeSub::sim_normal(x = phy, l = 10000,  rate = 0.0)

  distances <- calc_dist(seq_node_sub$alignment, seq_node_sub$root_seq)

  testthat::expect_equal(sum(distances), 0)
})

test_that("calc_dist abuse", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")
  seq_node_sub <- nodeSub::sim_normal(x = phy, l = 10000,  rate = 0.0)

  seq_dnabin <- ape::as.DNAbin(seq_node_sub$alignment)
  testthat::expect_error(calc_dist(seq_dnabin, seq_node_sub$root_seq))

  testthat::expect_error(calc_dist(seq_node_sub$alignment))
})
