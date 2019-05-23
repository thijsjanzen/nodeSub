#' function create an alignment with highly similar information content
#' @param input_tree phylogeny for which to generate alignment
#' @param focal_alignment alignment to match information content with
#' @param alt_model alternative substitution model
#' @param root_sequence root sequence
#' @return list with alignment and inferred rate
#' @export
create_equal_alignment <- function(input_tree,
                                   focal_alignment,
                                   root_sequence,
                                   alt_model) {

  num_emp_subs <- sum(calc_dist(focal_alignment, root_sequence))

  t_mrca <- calc_tree_height(input_tree)

  # make an educated guess
  adjusted_rate <- num_emp_subs / (length(root_sequence) * t_mrca * length(input_tree$tip.label))

  proposed_alignment <- alt_model(input_tree, adjusted_rate, root_sequence)$alignment
  proposed_subs <- sum(calc_dist(proposed_alignment, root_sequence))
  rolling_avg <- rep(NA, 9)
  cnt <- 1
  while(proposed_subs != num_emp_subs) {
    proposed_alignment <- alt_model(phy = input_tree,
                                    rate = adjusted_rate,
                                    rootseq = root_sequence)$alignment
    proposed_subs <- sum(calc_dist(proposed_alignment, root_sequence))
    rolling_avg[cnt] <- proposed_subs
   # cat(
    factor <- 1
    if(cnt >= 10) {
      avg_sub <- mean(rolling_avg, na.rm=T)
      factor <- num_emp_subs / avg_sub
     # if(num_emp_subs < avg_sub) change <- 0.95
      adjusted_rate <- adjusted_rate * factor
      cnt <- 0
    }
    cat(cnt, adjusted_rate, proposed_subs, num_emp_subs, factor,"\n")
    cnt <- cnt + 1
  }
  return(list("alignment" = proposed_alignment,
              "rate" = adjusted_rate))
}
