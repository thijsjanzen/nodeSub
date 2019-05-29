context("create_equal_alignment")

test_that("create_equal_alignment", {
  phy  <- phytools::read.newick(text = "(t1:10,(t3:2,t2:2):8);")
  seq_node_sub <- sim_normal(x = phy, l = 1000,  rate = 0.1)

  sim_regular <- function(phy, rate, rootseq) {
    sim_normal(x = phy, l = 1000, rootseq = rootseq, rate = rate)
  }


  seq_alt <- create_equal_alignment(input_tree = phy,
                                    focal_alignment = seq_node_sub$alignment,
                                    root_sequence = seq_node_sub$root_seq,
                                    alt_model = sim_regular)

  testthat::expect_equal( length(seq_alt$alignment),
                          length(seq_node_sub$alignment) )

  num_diff_base <- calc_dist(seq_node_sub$alignment, seq_node_sub$root_seq)
  num_diff_alt  <- calc_dist(seq_alt$alignment, seq_node_sub$root_seq)

  testthat::expect_equal(sum(num_diff_base), sum(num_diff_alt))
})