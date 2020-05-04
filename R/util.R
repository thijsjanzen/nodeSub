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
  if (is.matrix(data))
    nam <- row.names(data)
  else nam <- names(data)
  if (inherits(data, "DNAbin"))
    data <- as.character(data)
  if (inherits(data, "character")) data <- as.matrix(data)
  if (is.matrix(data))
    data <- as.data.frame(t(data), stringsAsFactors = FALSE)
  else data <- as.data.frame(data, stringsAsFactors = FALSE)

  data <- data.frame(tolower(as.matrix(data)), stringsAsFactors = FALSE)

  ac <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", "y",
          "k", "v", "h", "d", "b", "n", "?", "-")
  AC <- matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1), # nolint
                 c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1),
                 c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1)),
               18, 4, dimnames = list(NULL, c("a", "c", "g", "t")))

  compress <- TRUE
  if (length(data[[1]]) == 1) compress <- FALSE
  if (compress) {
    ddd <- fast.table(data)
    data <- ddd$data
    weight <- ddd$weight
    index <- ddd$index
  } else{
    p <- length(data[[1]])
    weight <- rep(1, p)
    index <- 1:p
  }
  p <- length(data[[1]])
  att <- attributes(data)

  data <- lapply(data, match, ac)
  attributes(data) <- att
  row.names(data) <- as.character(1:p)
  data <- stats::na.omit(data)
  rn <- as.numeric(rownames(data))

  if (!is.null(attr(data, "na.action"))) {
    warning("Found unknown characters. Deleted sites with with unknown states.")
  }

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
#' this is slower than the C++ implementation
#' but provides a way to debug
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

#' generate a tree conditional on n and t
#' @description generates a tree, conditional on number of tips and crown age
#' @param b birth rate (speciation rate)
#' @param m death rate (extinction rate)
#' @param focal_num_tips total number of tips to be conditioned on
#' @param focal_crown_age crown age to be conditioned on
#' @return phylo tree
#' @export
get_tree <- function(b,
                     m,
                     focal_num_tips,
                     focal_crown_age) {
  result_tree <-  phytools::pbtree(b = b, d = m,
                                   n = focal_num_tips,
                                   t = focal_crown_age,
                                   nsim = 1, quiet = FALSE)

  no_extinct_tree <- geiger::drop.extinct(result_tree)

  while (beautier::get_crown_age(no_extinct_tree) != focal_crown_age) {
    result_tree <- phytools::pbtree(b = b, d = m,
                                    n = focal_num_tips,
                                    t = focal_crown_age,
                                    nsim = 1, quiet = FALSE)
    no_extinct_tree <- geiger::drop.extinct(result_tree)
  }
  return(result_tree)
}




