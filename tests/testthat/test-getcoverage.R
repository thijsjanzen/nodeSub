context("get coverage")

test_that("sim_dual_independent", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  sequences <- sim_dual_independent(phy)
  testthat::expect_true(class(sequences) == "phyDat")
})

test_that("sim_dual_parent", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  sequences <- sim_dual_parent(phy)
  testthat::expect_true(class(sequences) == "phyDat")
})

test_that("sim_dual_linked", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  sequences <- sim_dual_parent(phy)
  testthat::expect_true(class(sequences) == "phyDat")
})





