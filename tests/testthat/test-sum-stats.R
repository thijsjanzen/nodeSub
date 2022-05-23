context("calc_sum_stats")

test_that("calc_sum_stats", {
  testthat::skip_on_cran()
  #  skip on cran, RPANDA seems to cause BLAS errors

  if (requireNamespace("TreeSim")) {
    phy1 <- TreeSim::sim.bd.taxa(n = 100,
                                       numbsim = 1, lambda = 1, mu = 0)[[1]]
  } else {
    if (requireNamespace("ape")) {
      phy1 <- ape::rphylo(n = 100, birth = 1, death = 0.0)
    } else {
      stop("could not use TreeSim or ape to simulate tree")
    }
  }


  brts <- ape::branching.times(phy1)

  phy <- nodeSub::create_balanced_tree(brts)

  input <- list(phy, phy)
  class(input) <- "multiPhylo"

  testthat::expect_true(class(input) == "multiPhylo")
  testthat::expect_true(length(input) == 2)

  testthat::expect_warning(
    stats1 <- nodeSub::calc_sum_stats(input, phy)
  )
  testthat::expect_warning(
    stats2 <- nodeSub::calc_sum_stats(input, phy, verbose = TRUE)
  )

  testthat::expect_true(stats1$stats$beta[[1]] >= 9.9)
  testthat::expect_true(stats1$stats$beta[[2]] >= 9.9)

  testthat::expect_true(sum(stats1$differences, na.rm = TRUE) == 0)

  testthat::expect_warning(
    stats2 <- nodeSub::calc_sum_stats(input, phy1)
  )
  testthat::expect_true(stats2$stats$beta[[1]] >= 9.9)
  testthat::expect_true(stats2$stats$beta[[2]] >= 9.9)

  testthat::expect_true(sum(stats2$differences, na.rm = TRUE) != 0)

  testthat::expect_true(all.equal(stats1$stats, stats2$stats))


  phy <- nodeSub::create_unbalanced_tree(brts)

  testthat::expect_warning(
    stats1 <- nodeSub::calc_sum_stats(phy, phy)
  )

  testthat::expect_true(class(input) == "multiPhylo")
  testthat::expect_true(length(input) == 2)

  testthat::expect_true(stats1$stats$beta[[1]] < 0.0)

  testthat::expect_true(sum(stats1$differences, na.rm = TRUE) == 0)

  if (requireNamespace("TreeSim")) {
    phy1 <- TreeSim::sim.bd.taxa(n = 100,
                                 numbsim = 1, lambda = 1, mu = 0.5,
                                 complete = TRUE)[[1]]
  } else {
    if (requireNamespace("ape")) {
      phy1 <- ape::rphylo(n = 100, birth = 1, death = 0.5, fossils = TRUE)
    } else {
      stop("could not use TreeSim or ape to simulate tree")
    }
  }

  phy2 <- geiger::drop.extinct(phy1)


  testthat::expect_warning(
    nodeSub::calc_sum_stats(phy2, phy1)
  )

})

test_that("calc_sum_stats abuse", {
  if (requireNamespace("TreeSim")) {
    phy1 <- TreeSim::sim.bd.taxa(n = 100,
                                 numbsim = 1, lambda = 1, mu = 0)[[1]]
  } else {
    if (requireNamespace("ape")) {
      phy1 <- ape::rphylo(n = 100, birth = 1, death = 0.0)
    } else {
      stop("could not use TreeSim or ape to simulate tree")
    }
  }

  testthat::expect_warning(
    nodeSub::calc_sum_stats(phy1, phy1)
  )
})

