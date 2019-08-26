#' function create an alignment with highly similar information content
#' @param input_tree phylogeny for which to generate alignment
#' @param focal_alignment alignment to match information content with
#' @param alt_model alternative substitution model
#' @param root_sequence root sequence
#' @param verbose provide intermediate output
#' @return list with alignment and inferred rate
#' @export
create_equal_alignment <- function(input_tree,
                                   focal_alignment,
                                   root_sequence,
                                   alt_model,
                                   verbose = FALSE) {

  num_emp_subs <- sum(calc_dist(focal_alignment, root_sequence))

  t_mrca <- calc_tree_height(input_tree)

  # make an educated guess
  adjusted_rate <- num_emp_subs /
    (length(root_sequence) * t_mrca * length(input_tree$tip.label))

  proposed_alignment <- alt_model(input_tree,
                                  adjusted_rate,
                                  root_sequence)$alignment

  proposed_subs <- sum(calc_dist(proposed_alignment, root_sequence))
  cnt <- 1

  propose_alignments <- function(buffer, focal_rate) {
    proposed_alignment <- alt_model(phy = input_tree,
                                    rate = focal_rate,
                                    rootseq = root_sequence)$alignment
    return(proposed_alignment)
  }

  calc_subs <- function(local_alignment) {
    sum(calc_dist(local_alignment, root_sequence))
  }

  stored_factor <- 0
  while (proposed_subs != num_emp_subs) {

    alignments <- vector("list", 10)
    if(abs(stored_factor - 1) < 0.01) alignments <- vector("list", 100)
    alignments <- lapply(alignments, propose_alignments, adjusted_rate)
    all_subs <- unlist(lapply(alignments, calc_subs))

    if(sum(all_subs == num_emp_subs)) {
      a <- which(all_subs == num_emp_subs)
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
              "rate" = adjusted_rate))
}
