test_that("abuse", {

  root_sequence <- "acgt"

  expect_error(
      sim_dual_independent(
        phy = ape::rcoal(3),
        rootseq = strsplit(root_sequence, split = "")[[1]],
        l = 1234500
      ),
    "'rootseq' must have the same length as 'l'"
  )
})
