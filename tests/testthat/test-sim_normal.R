context("nodeSub_sim_normal")

test_that("sim_normal", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  seq_phangorn <- phangorn::simSeq(phy, rate = 0.01, l = 10000)
  dist_phangorn <- phangorn::dist.ml(seq_phangorn)

  seq_nodeSub <- sim_normal(phy, rate = 0.01, l = 10000)
  dist_nodeSub <- phangorn::dist.ml(seq_nodeSub)

  testthat::expect_equal(dist_phangorn[1], dist_nodeSub[1], tolerance = 0.05)
  testthat::expect_equal(dist_phangorn[2], dist_nodeSub[2], tolerance = 0.05)
  testthat::expect_equal(dist_phangorn[3], dist_nodeSub[3], tolerance = 0.05)
})



