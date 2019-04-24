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

#' @keywords internal
#' @rawNamespace useDynLib(nodeSub)
get_p_matrix <- function(branch_length, eig = phangorn::edQt(), rate = 1.0) {
  res <- get_p_m_rcpp(eig, branch_length, rate)
  return(res)
}
