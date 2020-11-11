context("sim_linked")

test_that("sim_linked", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  sequences <- nodeSub::sim_linked(phy)
  testthat::expect_true(class(sequences$alignment) == "phyDat")

  Q_JC <- matrix(1, nrow = 4, ncol = 4)  # nolint

  sequences <- sim_linked(phy, Q = Q_JC)
  testthat::expect_true(class(sequences$alignment) == "phyDat")

  sequences <- sim_linked(phy)
  testthat::expect_true(class(sequences$alignment) == "phyDat")
  sequences <- sim_linked(phy, Q = Q_JC)
  testthat::expect_true(class(sequences$alignment) == "phyDat")
})


test_that("zeros", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  sequences <- sim_linked(phy,
                          rate = 0,
                          node_time = 0)
  testthat::expect_true(class(sequences$alignment) == "phyDat")

  dist_node_sub <- phangorn::dist.ml(sequences$alignment)
  testthat::expect_equal(min(dist_node_sub), max(dist_node_sub))
  testthat::expect_equal(min(dist_node_sub), 0)


  sequences <- sim_linked(phy,
                          rate = 0,
                          node_time = 0)
  testthat::expect_true(class(sequences$alignment) == "phyDat")
  dist_node_sub <- phangorn::dist.ml(sequences$alignment)
  testthat::expect_equal(min(dist_node_sub), max(dist_node_sub))
  testthat::expect_equal(min(dist_node_sub), 0)
})

test_that("abuse", {

  root_sequence <- "acgt"

  testthat::expect_error(
    sim_linked(
      phy = ape::rcoal(3),
      rootseq = strsplit(root_sequence, split = "")[[1]],
      l = 12345
    ),
    "'rootseq' must have the same length as 'l'"
  )
})
