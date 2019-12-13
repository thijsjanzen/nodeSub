context("dual_independent")


test_that("sim_dual_independent", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  sequences <- sim_dual_independent(phy)
  testthat::expect_true(class(sequences$alignment) == "phyDat")

  Q_JC <- matrix(1, nrow = 4, ncol = 4)
  sequences <- sim_dual_independent(phy, Q1 = Q_JC, Q2 = Q_JC)
  testthat::expect_true(class(sequences$alignment) == "phyDat")


  sequences <- sim_unlinked(phy)
  testthat::expect_true(class(sequences) == "phyDat")
  sequences <- sim_unlinked(phy, Q1 = Q_JC, Q2 = Q_JC)
  testthat::expect_true(class(sequences) == "phyDat")
})

test_that("zeros", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  sequences <- sim_dual_independent(phy, rate1 = 0, rate2 = 0)
  testthat::expect_true(class(sequences$alignment) == "phyDat")

  dist_node_sub <- phangorn::dist.ml(sequences$alignment)
  testthat::expect_equal(min(dist_node_sub), max(dist_node_sub))
  testthat::expect_equal(min(dist_node_sub), 0)



  sequences <- sim_unlinked(phy, rate1 = 0, rate2 = 0)
  testthat::expect_true(class(sequences) == "phyDat")

  dist_node_sub <- phangorn::dist.ml(sequences)
  testthat::expect_equal(min(dist_node_sub), max(dist_node_sub))
  testthat::expect_equal(min(dist_node_sub), 0)
})
