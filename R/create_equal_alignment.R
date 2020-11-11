#' function create an alignment with highly similar information content
#' @param input_tree phylogeny for which to generate alignment
#' @param sub_rate substitution rate used in the original phylogeny
#' @param alignment_result result of sim_normal, sim_linked or sim_unlinked
#' @param verbose provide intermediate output
#' @return new alignment, with added property "adjusted rate"
#' @export
create_equal_alignment <- function(input_tree,
                                   sub_rate,
                                   alignment_result,
                                   verbose = FALSE) {

  num_emp_subs <- alignment_result$total_accumulated_substitutions

  adjusted_rate <- sub_rate +
              sub_rate * alignment_result$total_node_substitutions /
                         alignment_result$total_branch_substitutions

  initial_guess <- adjusted_rate

  seqlen <- length(alignment_result$root_seq)

  proposed_alignment <- sim_normal(x = input_tree,
                                   l = seqlen,
                                   rootseq = alignment_result$root_seq,
                                   rate = adjusted_rate)

  propose_alignments <- function(buffer) {
    sim_normal(x = input_tree,
               l = seqlen,
               rootseq = alignment_result$root_seq,
               rate = adjusted_rate)
  }

  calc_subs <- function(local_alignment) {
     return(local_alignment$total_accumulated_substitutions)
  }

  proposed_subs <- proposed_alignment$total_accumulated_substitutions
  cnt <- 1
  while (proposed_subs != num_emp_subs) {

    alignments <- vector("list", 10)
    alignments <- lapply(alignments, propose_alignments)

    all_subs <- unlist(lapply(alignments, calc_subs))
    num_matches <- length(which(all_subs == num_emp_subs))

    if (num_matches > 0) {
      a <- which(all_subs == num_emp_subs)[[1]]
      output_alignment <- alignments[[a]]
      output_alignment$adjusted_rate <- adjusted_rate
      return(output_alignment)
    }
    cnt <- cnt + length(all_subs)
    if (verbose) cat(cnt, adjusted_rate, mean(all_subs, na.rm = TRUE), "\n")
  }

  proposed_alignment$adjusted_rate <- adjusted_rate

  return(proposed_alignment)
}
