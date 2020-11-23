#' function create an alignment with highly similar information content
#' @param input_tree phylogeny for which to generate alignment
#' @param sub_rate substitution rate used in the original phylogeny
#' @param alignment_result result of sim_normal, sim_linked or sim_unlinked
#' @param sim_function function that accepts a tree, sequence length,
#' rootsequence and substitution rate (in that order). Default is sim_normal
#' @param verbose provide intermediate output
#' @return new alignment, with added property "adjusted rate"
#' @export
create_equal_alignment <- function(input_tree,
                                   sub_rate,
                                   alignment_result,
                                   sim_function = NULL,
                                   verbose = FALSE) {

  num_emp_subs <- alignment_result$total_accumulated_substitutions

  adjusted_rate <- sub_rate +
              sub_rate * alignment_result$total_node_substitutions /
                         alignment_result$total_branch_substitutions

  seqlen <- length(alignment_result$root_seq)

  if (is.null(sim_function)) {
    sim_function <- function(input_tree, seqlen, rootseq, rate) {
      sim_normal(x = input_tree,
                 l = seqlen,
                 rootseq = rootseq,
                 rate = rate)
    }
  }

  proposed_alignment <- sim_function(input_tree, seqlen,
                                     alignment_result$root_seq, adjusted_rate)

  proposed_subs <- proposed_alignment$total_accumulated_substitutions

  while (proposed_subs != num_emp_subs) {
    proposed_alignment <- sim_function(input_tree, seqlen,
                                       alignment_result$root_seq, adjusted_rate)

    proposed_subs <- proposed_alignment$total_accumulated_substitutions
    if (verbose) cat(proposed_subs, "\n")
  }
  proposed_alignment$adjusted_rate <- adjusted_rate

  return(proposed_alignment)
}
