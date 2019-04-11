context("get coverage")

test_that("sim_dual_independent", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  sequences <- sim_dual_independent(phy)
  testthat::expect_true(class(sequences) == "phyDat")

  Q_JC <- matrix(1, nrow = 4, ncol = 4)
  sequences <- sim_dual_independent(phy, Q1 = Q_JC, Q2 = Q_JC)
  testthat::expect_true(class(sequences) == "phyDat")
})

test_that("sim_dual_parent", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  sequences <- sim_dual_parent(phy)
  testthat::expect_true(class(sequences) == "phyDat")

  Q_JC <- matrix(1, nrow = 4, ncol = 4)
  sequences <- sim_dual_parent(phy, Q1 = Q_JC, Q2 = Q_JC)
  testthat::expect_true(class(sequences) == "phyDat")

})

test_that("sim_dual_linked", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  sequences <- sim_dual_linked(phy)
  testthat::expect_true(class(sequences) == "phyDat")

  Q_JC <- matrix(1, nrow = 4, ncol = 4)

  sequences <- sim_dual_linked(phy, Q = Q_JC)
  testthat::expect_true(class(sequences) == "phyDat")
})





