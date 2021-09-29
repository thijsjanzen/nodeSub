#' function create an alignment with identical information content, using the
#' explicit method to simulate substitutions
#' @param input_tree phylogeny for which to generate alignment
#' @param sub_rate substitution rate used in the original phylogeny
#' @param alignment_result result of sim_normal, sim_linked or sim_unlinked
#' @param verbose provide intermediate output
#' @return list with four properties: 1) alignment: the alignment itself,
#' 2) adjusted rate: the substitution rate used to obtain identical information
#' content 3) total_accumulated_substitutions: the total number of
#' substitutions accumulated. 4) total_node_substitutions: total number of
#' substitutions accumulated on the nodes 5) total_branch_substitutions: total
#' number of substitutions accumulated on the branches.
#' @export
create_equal_alignment_explicit <- function(input_tree,
                                            sub_rate,
                                            alignment_result,
                                            verbose = FALSE) {

  num_emp_subs <- alignment_result$total_accumulated_substitutions

  adjusted_rate <- sub_rate +
    sub_rate * alignment_result$total_node_substitutions /
    alignment_result$total_branch_substitutions

  seqlen <- length(alignment_result$root_seq)

  proposed_alignment <- sim_normal_explicit(x = input_tree,
                                            l = seqlen,
                                            rootseq = alignment_result$root_seq,
                                            rate = adjusted_rate)

  proposed_subs <- proposed_alignment$total_accumulated_substitutions

  while (proposed_subs != num_emp_subs) {
    proposed_alignment <- sim_normal_explicit(x = input_tree,
                                              l = seqlen,
                                              rootseq =
                                                   alignment_result$root_seq,
                                              rate = adjusted_rate)

    proposed_subs <- proposed_alignment$total_accumulated_substitutions
    if (verbose) cat(proposed_subs, num_emp_subs, "\n")
  }
  proposed_alignment$adjusted_rate <- adjusted_rate

  return(proposed_alignment)
}
