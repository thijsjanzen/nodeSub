context("nodeSub_sim_normal")

test_that("sim_normal", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")
  set.seed(43) # for reproducability
  seq_phangorn <- phangorn::simSeq(x = phy, l = 10000, rate = 0.01)
  dist_phangorn <- phangorn::dist.ml(seq_phangorn)

  seq_nodeSub <- sim_normal(x = phy, l = 10000,  rate = 0.01)
  dist_nodeSub <- phangorn::dist.ml(seq_nodeSub)

  testthat::expect_equal(dist_phangorn[1], dist_nodeSub[1], tolerance = 0.05)
  testthat::expect_equal(dist_phangorn[2], dist_nodeSub[2], tolerance = 0.05)
  testthat::expect_equal(dist_phangorn[3], dist_nodeSub[3], tolerance = 0.05)
})

test_that("JC", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  Q_JC <- matrix(1, nrow = 4, ncol = 4)
  set.seed(42) # for reproducability
  seq_phangorn <- phangorn::simSeq(x = phy, Q = Q_JC, l = 10000, rate = 0.01)
  dist_phangorn <- phangorn::dist.ml(seq_phangorn)

  seq_nodeSub <- sim_normal(x = phy, Q = Q_JC, l = 10000,  rate = 0.01)
  dist_nodeSub <- phangorn::dist.ml(seq_nodeSub)

  testthat::expect_equal(dist_phangorn[1], dist_nodeSub[1], tolerance = 0.05)
  testthat::expect_equal(dist_phangorn[2], dist_nodeSub[2], tolerance = 0.05)
  testthat::expect_equal(dist_phangorn[3], dist_nodeSub[3], tolerance = 0.05)
})



