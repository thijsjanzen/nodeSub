#' function create an alignment with highly similar information content
#' @param input_tree phylogeny for which to generate alignment
#' @param sub_rate substitution rate used in the original phylogeny
#' @param alignment_result result of sim_normal, sim_linked or sim_unlinked
#' @param verbose provide intermediate output
#' @return list with alignment and inferred rate
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

  proposed_subs <- proposed_alignment$total_accumulated_substitutions
  cnt <- 1

  propose_alignments <- function(buffer, focal_rate) {
    proposed_alignment <- sim_normal(x = input_tree,
                                     l = seqlen,
                                     rootseq = alignment_result$root_seq,
                                     rate = adjusted_rate)
    return(proposed_alignment)
  }

  calc_subs <- function(local_alignment) {
    return(local_alignment$total_accumulated_substitutions)
  }

  stored_factor <- 0
  while (proposed_subs != num_emp_subs) {

    alignments <- vector("list", 10)
    # if (abs(stored_factor - 1) < 0.01) alignments <- vector("list", 100)
    alignments <- lapply(alignments, propose_alignments, adjusted_rate)
    all_subs <- unlist(lapply(alignments, calc_subs))

    num_matches <- length(which(all_subs == num_emp_subs))

    if (num_matches > 0) {
      a <- which(all_subs == num_emp_subs)[[1]]
      return(list("alignment" = alignments[[a]],
                  "rate" = adjusted_rate))
    } else {
      avg_sub <- mean(all_subs, na.rm = TRUE)
      factor <- num_emp_subs / avg_sub
      stored_factor <- factor
      adjusted_rate <- adjusted_rate * factor
    }

    cnt <- cnt + length(all_subs)
    if (verbose) cat(cnt, adjusted_rate, mean(all_subs, na.rm = TRUE),
                     num_emp_subs, factor, "\n")
  }

  return(list("alignment" = proposed_alignment,
              "rate" = adjusted_rate,
              "initial_rate" = initial_guess))
}
