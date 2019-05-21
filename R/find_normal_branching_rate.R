#' @keywords internal
calc_dist <- function(alignment) {
  num_seq <- length(alignment)
  dist_dist <- rep(NA, num_seq * (num_seq-1)/2)
  cnt <- 1
  for(i in 1:num_seq) {
    for(j in 1:i) {
      if(i != j) {
        a <- alignment[[i]]
        b <- alignment[[j]]
        local_diff <- length(a) - sum(a == b)
        dist_dist[cnt] <- local_diff
        cnt <- cnt + 1
      }
    }
  }
  dist_dist <- dist_dist[!is.na(dist_dist)]
  return(dist_dist)
}

#' find a matching substitution rate
#' @param phy tree for which to simulate sequences
#' @param alignment alignment of alternative model
#' @param sim_regular function to simulate an alignment under a normal model, taking two arguments: phy and rate
#' @param replicates number of replicate runs for which to estimate
#' @param max_rate maximum rate to be considered in optimization
#' @return rate
#' @export
find_normal_branching_rate <- function(phy,
                                       alignment,
                                       sim_regular,
                                       replicates = 1,
                                       max_rate = 1.0) {

  emp_dist <- calc_dist(alignment)

  found <- rep(NA, replicates)
  for(i in 1:replicates) {
    find_ll <- function(params) {
      vy <- sim_regular(phy, params)
      obs_dist <- calc_dist(vy)
      fit <- sum(abs(sort(obs_dist) - sort(emp_dist)))
      return(fit)
    }
    vv <- stats::optimize(f = find_ll, interval = c(0, max_rate))
    found[i] <- vv$minimum
  }
  return(found)
}
