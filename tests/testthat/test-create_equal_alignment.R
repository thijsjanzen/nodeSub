context("create_equal_alignment")

test_that("create_equal_alignment", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")
  set.seed(666)

  # manual testing:
  # phy <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0)[[1]]


  sub_rate <- 1e-2
  seqlen <- 1000

  tau <- beautier::get_crown_age(phy) * 0.1
  seq_node_sub <- sim_unlinked(phy = phy, l = seqlen,
                               rate1 = sub_rate,
                               rate2 = sub_rate,
                               node_time = tau)

  sim_regular <- function(phy, rate, rootseq) {
    sim_normal(x = phy, l = seqlen, rootseq = rootseq, rate = rate)
  }

  seq_alt <- create_equal_alignment(input_tree = phy,
                                    alignment_result = seq_node_sub,
                                    sub_rate = sub_rate,
                                    verbose = TRUE)

  testthat::expect_equal(length(seq_alt$alignment),
                         length(seq_node_sub$alignment))

  subs_normal <- seq_node_sub$total_accumulated_substitutions
  subs_alt    <- seq_alt$total_accumulated_substitutions

  expected_rate <- sub_rate * (1 +  nodeSub::calc_fraction(phy = phy,
                                                           node_time = tau))

  testthat::expect_equal(expected_rate, seq_alt$adjusted_rate, tolerance = 0.01)

  testthat::expect_equal(sum(subs_normal), sum(subs_alt))
})
