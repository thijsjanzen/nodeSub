context("calc_sum_stats")

test_that("calc_sum_stats", {
  skip_on_cran()
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  seq_node_sub <- sim_normal(x = phy, l = 100,  rate = 0.1)

  skip("Not now, @thijsjanzen")
  all_trees <- infer_phylogeny(seq_node_sub$alignment,
                               treatment_name = "test",
                               burnin = 0.1)

  stats <- nodeSub::calc_sum_stats(all_trees$all_trees,
                                   phy, verbose = FALSE)

  testthat::expect_equal(length(stats$stats$beta), length(all_trees$all_trees))
  testthat::expect_equal(length(colnames(stats$stats)), 4)
  testthat::expect_equal(length(colnames(stats$differences)), 5)
})
