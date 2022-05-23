context("nodeSub_sim_normal_explicit")

test_that("sim_normal_explicit", {
  phy  <- ape::read.tree(text = "(t1:10,(t3:2,t2:2):8);")
  set.seed(43) # for reproducability
  seq_implicit <- sim_normal(x = phy, l = 10000,  rate = 0.01)
  seq_explicit <- sim_normal_explicit(x = phy, l = 10000, rate = 0.01,
                                      rootseq = seq_implicit$root_seq)

  a1 <- calc_dist(seq_implicit$alignment, seq_implicit$root_seq)
  a2 <- calc_dist(seq_explicit$alignment, seq_implicit$root_seq)

  testthat::expect_equal(length(a1), 3)
  testthat::expect_equal(length(a2), 3)
  testthat::expect_equal(sum(a1), sum(a2), tolerance = 0.05)
})
