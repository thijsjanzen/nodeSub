context("nodeSub_sim_normal")

test_that("sim_normal", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")
  set.seed(43) # for reproducability
  seq_phangorn <- phangorn::simSeq(x = phy, l = 10000, rate = 0.01)
  dist_phangorn <- phangorn::dist.ml(seq_phangorn)

  seq_node_sub <- sim_normal(x = phy, l = 10000,  rate = 0.01)
  dist_node_sub <- phangorn::dist.ml(seq_node_sub$alignment)

  testthat::expect_equal(dist_phangorn[1], dist_node_sub[1], tolerance = 0.05)
  testthat::expect_equal(dist_phangorn[2], dist_node_sub[2], tolerance = 0.05)
  testthat::expect_equal(dist_phangorn[3], dist_node_sub[3], tolerance = 0.05)
})

test_that("JC", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  Q_JC <- matrix(1, nrow = 4, ncol = 4)  # nolint
  set.seed(42) # for reproducability
  seq_phangorn <- phangorn::simSeq(x = phy, Q = Q_JC, l = 10000, rate = 0.01)
  dist_phangorn <- phangorn::dist.ml(seq_phangorn)

  seq_node_sub <- sim_normal(x = phy, Q = Q_JC, l = 10000, rate = 0.01)
  dist_node_sub <- phangorn::dist.ml(seq_node_sub$alignment)

  testthat::expect_equal(dist_phangorn[1], dist_node_sub[1], tolerance = 0.05)
  testthat::expect_equal(dist_phangorn[2], dist_node_sub[2], tolerance = 0.05)
  testthat::expect_equal(dist_phangorn[3], dist_node_sub[3], tolerance = 0.05)
})

test_that("zero", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  set.seed(42) # for reproducability
  seq_phangorn <- phangorn::simSeq(x = phy, l = 10000, rate = 0.0)
  dist_phangorn <- phangorn::dist.ml(seq_phangorn)

  seq_node_sub <- sim_normal(x = phy, l = 10000, rate = 0.0)
  dist_node_sub <- phangorn::dist.ml(seq_node_sub$alignment)

  testthat::expect_equal(dist_phangorn[1], dist_node_sub[1], tolerance = 0.05)
  testthat::expect_equal(dist_phangorn[2], dist_node_sub[2], tolerance = 0.05)
  testthat::expect_equal(dist_phangorn[3], dist_node_sub[3], tolerance = 0.05)

  testthat::expect_equal(min(dist_phangorn), max(dist_phangorn))
  testthat::expect_equal(min(dist_phangorn), 0)

  testthat::expect_equal(min(dist_node_sub), max(dist_node_sub))
  testthat::expect_equal(min(dist_node_sub), 0)
})
