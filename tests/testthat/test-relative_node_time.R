context("calc_required_node_time")

test_that("required_node_time use", {
  set.seed(1)

  if (requireNamespace("TreeSim")) {
    focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                       numbsim = 1, lambda = 1, mu = 0)[[1]]
  } else {
    if (requireNamespace("ape")) {
      focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)
    } else {
      stop("could not use TreeSim or ape to simulate tree")
    }
  }

  chosen_fraction <- 0.1
  req_node_time <- nodeSub::calc_required_node_time(focal_tree,
                                                    s = chosen_fraction)

  obs_fraction <- nodeSub::calc_fraction(focal_tree, node_time = req_node_time)

  testthat::expect_equal(obs_fraction, 0.1)

  ali <- nodeSub::sim_unlinked(focal_tree, rate1 = 0.001, rate2 = 0.001,
                                       l = 10000, node_time = req_node_time)

  obs_fraction <- ali$total_node_substitutions /
          (ali$total_node_substitutions + ali$total_branch_substitutions)

  testthat::expect_equal(obs_fraction, chosen_fraction, tol = 0.01)
})

test_that("required_node_time abuse", {
  focal_tree <- 1


  chosen_fraction <- 0.1

  testthat::expect_error(
    req_node_time <- nodeSub::calc_required_node_time(focal_tree,
                                                      s = chosen_fraction),
   "phy needs to be a valid phylo object"
  )

  testthat::expect_error(
    obs_fraction <- nodeSub::calc_fraction(focal_tree,
                                           node_time = req_node_time),
    "phy needs to be a valid phylo object"
  )

  testthat::expect_error(
    nodeSub::calc_expected_hidden_nodes(focal_tree,
                                        lambda = 1,
                                        mu = 0),
    "requires a valid phylogeny as input"
  )
})

test_that("calc_expected_hidden_nodes, use", {
  set.seed(1)
  found <- c()
  for (repl in 1:10) {
    if (requireNamespace("TreeSim")) {
      focal_tree <- TreeSim::sim.bd.taxa(n = 10,
                                         numbsim = 1, lambda = 1, mu = 0.3,
                                         complete = TRUE)[[1]]
    } else {
      if (requireNamespace("ape")) {
        focal_tree <- ape::rphylo(n = 10, birth = 1, death = 0.3,
                                  fossils = TRUE)
      } else {
        stop("could not use TreeSim or ape to simulate tree")
      }
    }

    obs_hidden_nodes <- nodeSub::count_hidden(focal_tree)
    exp_hidden_nodes <-
       nodeSub::calc_expected_hidden_nodes(geiger::drop.extinct(focal_tree),
                                           lambda = 1,
                                           mu = 0.3)
    found <- rbind(found, c(obs_hidden_nodes, exp_hidden_nodes))
  }
  found <- colMeans(found)
  testthat::expect_equal(found[1], found[2], tolerance = 0.1)
})

test_that("calc_expected_hidden_nodes, abuse", {
  if (requireNamespace("TreeSim")) {
    focal_tree <- TreeSim::sim.bd.taxa(n = 10,
                                       numbsim = 1, lambda = 1, mu = 0.3,
                                       complete = TRUE)[[1]]
  } else {
    if (requireNamespace("ape")) {
      focal_tree <- ape::rphylo(n = 10, birth = 1, death = 0.3,
                                fossils = TRUE)
    } else {
      stop("could not use TreeSim or ape to simulate tree")
    }
  }

  testthat::expect_error(
    nodeSub::calc_expected_hidden_nodes(geiger::drop.extinct(focal_tree))
  )

  testthat::expect_warning(
    nodeSub::calc_expected_hidden_nodes(focal_tree, lambda = 1, mu = 0.3)
  )
})
