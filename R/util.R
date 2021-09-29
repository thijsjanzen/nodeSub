#' @keywords internal
#' this is an internal function from the phangorn package.
fast.table <- function(data) {  # nolint
  if (!is.data.frame(data))
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  da <- do.call("paste", c(data, sep = "\r"))
  ind <- !duplicated(da)
  levels <- da[ind]
  cat <- factor(da, levels = levels)
  nl <- length(levels(cat))
  bin <- (as.integer(cat) - 1)
  pd <- nl
  bin <- bin[!is.na(bin)]
  if (length(bin)) bin <- bin + 1
  y <- tabulate(bin, pd)
  result <- list(index = bin, weights = y, data = data[ind, ])
  result
}

#' @keywords internal
#' this is an internal function from the phangorn package.
phyDat.DNA <- function(data) {  # nolint
  nam <- names(data)

  data <- as.data.frame(data, stringsAsFactors = FALSE)

  data <- data.frame(tolower(as.matrix(data)), stringsAsFactors = FALSE)

  ac <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y",
          "k", "v", "h", "d", "b", "n", "?", "-")
  AC <- matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1), # nolint
                 c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1),
                 c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1)),
               18, 4, dimnames = list(NULL, c("a", "c", "g", "t")))

  ddd <- fast.table(data)
  data <- ddd$data
  weight <- ddd$weight
  index <- ddd$index
  p <- length(data[[1]])
  att <- attributes(data)

  data <- lapply(data, match, ac)
  attributes(data) <- att
  row.names(data) <- as.character(1:p)
  data <- stats::na.omit(data)
  rn <- as.numeric(rownames(data))

  aaa <- match(index, attr(data, "na.action"))
  index <- index[is.na(aaa)]
  index <- match(index, unique(index))
  rn <- as.numeric(rownames(data))
  attr(data, "na.action") <- NULL

  weight <- weight[rn]
  p <- dim(data)[1]
  names(data) <- nam
  attr(data, "row.names") <- NULL
  attr(data, "weight") <- weight
  attr(data, "nr") <- p
  attr(data, "nc") <- 4
  attr(data, "index") <- index
  attr(data, "levels") <- c("a", "c", "g", "t")
  attr(data, "allLevels") <- ac
  attr(data, "type") <- "DNA"
  attr(data, "contrast") <- AC
  class(data) <- "phyDat"
  data
}

#' calculate p matrix
#' @rawNamespace import(Rcpp)
#' @rawNamespace useDynLib(nodeSub)
#' @description calculates the p matrix
#' @param branch_length branch length
#' @param eig eigen object
#' @param rate rate
#' @return p matrix
get_p_matrix <- function(branch_length, eig = phangorn::edQt(), rate = 1.0) {
  res <- get_p_m_rcpp(eig, branch_length, rate)
  if (any(res < 0)) res[res < 0] <- 0
  return(res)
}


#' this function calculates the p matrix within R
#' this is slower than the C++ implementation in \code{get_p_matrix}
#' but provides a way to debug and verify
#' @param eig eigen object
#' @param branch_length branch length
#' @param rate substitution rate
#' @return p matrix
#' @export
slow_matrix <- function(eig,
                        branch_length,
                        rate) {

  eva <- eig$values
  ev <- eig$vectors
  evei <- eig$inv

  dim_size <- ncol(evei)

  P <- matrix(NA, nrow = dim_size, ncol = dim_size)  # nolint
  # capital P is retained to conform to mathematical notation on wikipedia
  # and in the literature

  if (branch_length == 0 || rate <= 0) {
    P <- matrix(0, nrow = dim_size, ncol = dim_size)  # nolint
    diag(P) <- 1
    return(P)
  }

  tmp <- exp(eva * rate * branch_length)

  for (i in 1:dim_size) {
    for (j in 1:dim_size) {
      res <- 0.0
      for (h in 1:dim_size) {
        res <- res + ev[i, h] * tmp[h] * evei[h, j]
      }
      P[i, j] <- res
    }
  }
  return(P)
}

#' @keywords internal
draw_bases <- function(focal_base, trans_matrix) {
  bases <- c("a", "t", "c", "g")

  focus <- which(focal_base == bases)
  picked_index <- sample(1:10, 1, prob = trans_matrix[focus, ], replace = TRUE)

  output <- switch(picked_index,
                   c("a", "a"),
                   c("t", "t"),
                   c("c", "c"),
                   c("g", "g"),
                   c("a", "c"),
                   c("a", "t"),
                   c("a", "g"),
                   c("t", "c"),
                   c("t", "g"),
                   c("c", "g"))

  if (stats::runif(1, 0, 1) < 0.5) return(rev(output))

  return(output)
}

#' @keywords internal
pick_rate <- function(matching_bases, double_rate, focal_target) {
  output <- NA

  if (matching_bases == 2) {
    output <- NA  # e.g. A -> AA, to be filled in later
  }
  if (matching_bases == 1) {
    output <- 1 # e.g. A -> AT, one mutation
  }
  if (matching_bases == 0) { # double mutation!
    a <- substr(focal_target, 1, 1)
    b <- substr(focal_target, 2, 2)
    if (a != b) {
      output <- double_rate  # e.g. A -> TG
    } else {
      output <- 0            # e.g. A -> TT, impossible by definition
    }
    # if a == b, use default output = 0
  }
  return(output)
}

#' @keywords internal
make_transition_matrix <- function(mut_double) {

  output <- matrix(NA, nrow = 4, ncol = 10)

  order_events <- c("aa", "tt", "cc", "gg",
                    "ac", "at", "ag",
                    "tc", "tg",
                    "cg")

  bases <- c("a", "t", "c", "g")
  for (j in 1:4) {
    focal <- bases[j]
    to_add <- c()
    for (i in seq_along(order_events)) {
      focal_target <- order_events[i]
      matching_bases <- stringr::str_count(focal_target, focal)

      to_add[i] <- pick_rate(matching_bases, mut_double, focal_target)
    }
    output[j, ] <- to_add
  }

  for (i in 1:4) {
    output[i, ] <- output[i, ] / sum(output[i, ], na.rm = TRUE)
  }

  for (i in 1:4) {
    output[i, i] <-  - sum(output[i, ], na.rm = TRUE)
  }
  for (j in 1:6) {
    output <- rbind(output, rep(0, length(output[1, ])))
  }

  return(output)
}

#' @keywords internal
get_mutated_sequences <- function(parent_seq, trans_matrix) {

  vx <- sapply(parent_seq, draw_bases, trans_matrix)

  child1_seq <- vx[1, ]
  child2_seq <- vx[2, ]

  return(list(child1_seq, child2_seq))
}

#' @keywords internal
calc_dist <- function(alignment_phydat,
                      root_sequence = NULL) {

  if (is.null(root_sequence)) {
    root_sequences = alignment_phydat$root_seq
  }

  if (class(alignment_phydat) != "phyDat") {
    stop("input alignment has to be of type phyDat")
  }

  alignment_rawer <- phangorn::phyDat2alignment(alignment_phydat)

  alignment_matrix <- alignment_rawer$seq

  n_mutations <- rep(NA, length(alignment_matrix))
  for (i in seq_along(alignment_matrix)) {
    sequence_vector <- strsplit(alignment_matrix[i], split = "")[[1]]
    n_mutations[i] <- sum(root_sequence != sequence_vector)
  }
  return(n_mutations)
}

#' @keywords internal
calc_accumulated_substitutions <- function(phy, branch_subs, node_subs = NULL) {
  if (is.null(node_subs)) {
    node_subs <- rep(0, length(branch_subs))
  }
  edge <- phy$edge
  no_extinct_phy <- geiger::drop.extinct(phy)
  tips_to_check <- no_extinct_phy$tip.label
  # not used:
  num_extant <- length(tips_to_check)
  root_parent <- edge[1, 1]
  edge <- cbind(edge, 0)

  # we start at each extant tip
  # and then traverse the tree to the root and mark each edge as being connected
  # to an extant tip. Only those should be taken into account.
  for (i in 1:num_extant) {
    parent_index <- which(edge[, 2] == i)
    parent <- edge[parent_index, 1]
    edge[parent_index, 3] <- 1
    while (parent != root_parent) {
      parent_index <- which(edge[, 2] == parent)
      edge[parent_index, 3] <- 1
      parent <- edge[parent_index, 1]
    }
  }

  total_branch_subs <- sum(branch_subs * edge[, 3])
  total_node_subs   <- sum(node_subs * edge[, 3])
  total_accum_substitutions <- sum(total_branch_subs, total_node_subs)

  return(list("total_branch_subs" = total_branch_subs,
              "total_node_subs" = total_node_subs,
              "total_accumulated_substitutions" =
                total_accum_substitutions))
}
