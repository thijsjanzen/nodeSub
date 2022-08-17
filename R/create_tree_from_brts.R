#' create an unbalanced tree out of branching times
#' @param brts vector of branching times
#' @return phylo phylo object
#' @export
create_unbalanced_tree <- function(brts) {
  brts <- as.vector(sort(brts, decreasing = TRUE))

  ltab <- c()
  to_add1 <- c(brts[1], 0, -1, -1)
  to_add2 <- c(brts[1], -1, 2, -1)
  ltab <- rbind(ltab, to_add1, to_add2)

  leftcnt <- 2

  cnt <- 2
  for (i in 2:length(brts)) {
    focal_brt <- brts[i]
    cnt <- cnt + 1
    parent <- leftcnt
    daughter <- cnt * -1
    leftcnt <- daughter
    to_add <- c(focal_brt, parent, daughter, -1)
    ltab <- rbind(ltab, to_add)
  }

  phylo_tree <- DDD::L2phylo(ltab)
  return(phylo_tree)
}


#' create a balanced tree out of branching times
#' @param brts vector of branching times
#' @return phylo phylo object
#' @export
create_balanced_tree <- function(brts) {
  brts <- as.vector(sort(brts, decreasing = TRUE))

  ltab <- c()
  to_add1 <- c(brts[1], 0, -1, -1)
  to_add2 <- c(brts[1], -1, 2, -1)
  ltab <- rbind(ltab, to_add1, to_add2)

  cnt <- 3

  current_tips_left <- c(2)
  current_tips_right <- c(-1)
  brts <- brts[-1]

  while (length(brts) > 0) {

    num_to_do <- length(current_tips_left)
    for (i in 1:num_to_do) {
      if (i > length(brts)) {
        break
      }
      brt_left <- brts[i]
      parent <- current_tips_left[i]
      daughter <- cnt
      to_add <- c(brt_left, parent, daughter, -1)
      current_tips_left <- c(current_tips_left, daughter)
      ltab <- rbind(ltab, to_add)
      cnt <- cnt + 1
    }
    brts <- brts[-c(1:num_to_do)]

    num_to_do <- length(current_tips_right)

    for (i in 1:num_to_do) {
      if (i > length(brts)) {
        break
      }
      brt_right <- brts[i]
      parent <- current_tips_right[i]
      daughter <- cnt * -1
      to_add <- c(brt_right, parent, daughter, -1)
      current_tips_right <- c(current_tips_right, daughter)
      ltab <- rbind(ltab, to_add)
      cnt <- cnt + 1
    }
    brts <- brts[-c(1:num_to_do)]

  }
  phylo_tree <- DDD::L2phylo(ltab)
  return(phylo_tree)
}
