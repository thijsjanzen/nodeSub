#' @keywords internal
fast.table <- function(data) {
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
phyDat.DNA <- function(data) {
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
  AC <- matrix(c(c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1),
                 c(0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1),
                 c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1)),
               18, 4, dimnames = list(NULL, c("a", "c", "g", "t")))

  compress <- TRUE
  if (length(data[[1]]) == 1) compress <- FALSE
  if (compress){
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
slow_matrix <- function(eig, branch_length, rate) {

  eva <- eig$values
  ev <- eig$vectors
  evei <- eig$inv

  dim_size <- ncol(evei)

  P <- matrix(NA, nrow = dim_size, ncol = dim_size)

  if (branch_length == 0 || rate <= 0) {
    for (i in 1:dim_size) {
      for (j in 1:dim_size) {
        if (i != j) P[i, j] <- 0
        if (i == j) P[i, j] <- 1
      }
    }
    return(P)
  }

  tmp <- rep(NA, dim_size)
  for (i in 1:dim_size) {
    tmp[i] <- exp(1.0 * eva[i] * rate * branch_length)
  }

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
