context("create_equal_alignment")

test_that("create_equal_alignment", {
  # phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")
  set.seed(666)
  phy <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 1, mu = 0)[[1]]


  sub_rate <- 0.01
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
                                    alt_model = sim_regular,
                                    sub_rate = sub_rate,
                                    verbose = TRUE)

  testthat::expect_equal(length(seq_alt$alignment$alignment),
                         length(seq_node_sub$alignment))

  subs_normal <- seq_node_sub$total_accumulated_substitutions
  subs_alt    <- seq_alt$alignment$total_accumulated_substitutions


  testthat::expect_equal(sum(subs_normal), sum(subs_alt))
})
