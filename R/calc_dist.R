#' function to calculate the number of substitutions along a tree
#' @param alignment_phydat alignment object
#' @param root_sequence root sequence
#' @export
calc_dist <- function(alignment_phydat,
                      root_sequence = NULL) {

  if (is.null(root_sequence)) {
    stop("can not calculate distance from root sequence without root sequence")
  }

  testit::assert(class(alignment_phydat) == "phyDat")

  alignment_rawer <- phangorn::phyDat2alignment(alignment_phydat)

  alignment_matrix <- alignment_rawer$seq

  n_mutations <- rep(NA, length(alignment_matrix))
  for (i in seq_along(alignment_matrix)) {
    sequence_vector <- strsplit(alignment_matrix[i], split = "")[[1]]
    n_mutations[i] <- sum(root_sequence != sequence_vector)
  }
  return(n_mutations)
}
