#' function create an alignment with highly similar information content
#' @param input_tree phylogeny for which to generate alignment
#' @param focal_alignment alignment to match information content with
#' @param alt_model alternative substitution model
#' @param max_rate maximum rate for inference
#' @return list with alignment and inferred rate
#' @export
create_equal_alignment <- function(input_tree,
                                   focal_alignment,
                                   alt_model,
                                   max_rate = 1) {

  emp_dist <- calc_dist(focal_alignment)

  local_env <- new.env()
  local_env$output_alignment  <- list()
  local_env$cnt <- 0
  local_env$fit <- c()
  local_env$rate <- c()

  output_alignment <- c()

  find_ll <- function(params) {
    vy <- sim_regular(input_tree, params)
    obs_dist <- calc_dist(vy)
    fit <- sum(abs(sort(obs_dist) - sort(emp_dist)))

    local_env$cnt <- local_env$cnt + 1
    local_env$output_alignment[[ local_env$cnt ]] <- vy
    local_env$fit[  local_env$cnt ] <- fit
    local_env$rate[ local_env$cnt ] <- params[[1]]
    return(fit)
  }

  vv1 <- stats::optimize(f = find_ll, interval = c(0, max_rate))

  best_fit <- which.min(local_env$fit)
  output_alignment <- local_env$output_alignment[[best_fit]]
  best_rate <- local_env$rate[best_fit]

  return(list("alignment" = output_alignment,
              "rate" = best_rate))
}
