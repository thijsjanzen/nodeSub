#' function create an alignment with identical information content
#' @param input_tree phylogeny for which to generate alignment
#' @param sub_rate substitution rate used in the original phylogeny
#' @param alignment_result result of sim_normal, sim_linked or sim_unlinked
#' @param sim_function function that accepts a tree, sequence length,
#' rootsequence and substitution rate (in that order). Default is sim_normal
#' @param verbose provide intermediate output
#' @param node_time node time
#' @param input_alignment_type was the input alignment simulated with a node
#' substitution model or a normal substitution model? Used to calculate the
#' twin mutation rate. Options are "nodesub" and "normal".
#' @return list with four properties: 1) alignment: the alignment itself,
#' 2) adjusted rate: the substitution rate used to obtain identical information
#' content 3) total_accumulated_substitutions: the total number of
#' substitutions accumulated. 4) total_node_substitutions: total number of
#' substitutions accumulated on the nodes 5) total_branch_substitutions: total
#' number of substitutions accumulated on the branches.
#' @export
create_equal_alignment <- function(input_tree,
                                   sub_rate,
                                   alignment_result,
                                   sim_function = NULL,
                                   verbose = FALSE,
                                   node_time = NULL,
                                   input_alignment_type = "nodesub") {

  num_emp_subs <- alignment_result$total_accumulated_substitutions

  adjusted_rate <- sub_rate +
              sub_rate * alignment_result$total_node_substitutions /
                         alignment_result$total_branch_substitutions

  if (input_alignment_type == "normal") {
    if (is.null(node_time)) {
      stop("Node time needs to be provided")
    }
    total_node_sub <- node_time * 2 * input_tree$Nnode
    total_branch_time <- sum(input_tree$edge.length)
    frac <- 1 + total_node_sub /
      total_branch_time
    adjusted_rate <- sub_rate / frac
  }

  if (input_alignment_type == "fix_sub_rate") {
    adjusted_rate <- sub_rate
  }

  seqlen <- length(alignment_result$root_seq)

  if (is.null(sim_function)) {
    sim_function <- function(input_tree, seqlen, rootseq, rate) {
      nodeSub::sim_normal(x = input_tree,
                 l = seqlen,
                 rootseq = rootseq,
                 rate = rate)
    }
  }

  proposed_alignment <- sim_function(input_tree,
                                     seqlen,
                                     alignment_result$root_seq,
                                     adjusted_rate)

  proposed_subs <- proposed_alignment$total_accumulated_substitutions

  while (proposed_subs != num_emp_subs) {
    proposed_alignment <- sim_function(input_tree,
                                       seqlen,
                                       alignment_result$root_seq,
                                       adjusted_rate)

    proposed_subs <- proposed_alignment$total_accumulated_substitutions
    if (verbose) cat(proposed_subs, " ",
                     num_emp_subs, " ",
                     sub_rate, " ",
                     adjusted_rate, "\n")
  }
  proposed_alignment$adjusted_rate <- adjusted_rate

  return(proposed_alignment)
}
