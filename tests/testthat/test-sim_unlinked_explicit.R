context("sim_unlinked_explicit")

test_that("sim_unlinked_explicit", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")

  sequences <- nodeSub::sim_unlinked(phy, node_time = 0.1, rate1 = 0.01, rate2 = 0.01)
  testthat::expect_true(class(sequences$alignment) == "phyDat")

  sequences_2 <- nodeSub::sim_unlinked_explicit(phy = phy,
                                                rootseq = sequences$root_seq,
                                                node_time = 0.1,
                                                rate1 = 0.01,
                                                rate2 = 0.01)

  d1 <- calc_dist(sequences$alignment, root_sequence  = sequences$root_seq)
  d2 <- calc_dist(sequences_2$alignment, root_sequence  = sequences$root_seq)

  testthat::expect_equal(sum(d1), sum(d2), tolerance = 0.1)
})
