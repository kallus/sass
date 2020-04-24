glasso <- function(x, lambda) {
  d <- ncol(x)
  sol <- list()
  w <- NULL
  wi <- NULL
  S <- cor(x)
  for (i in 1:length(lambda)) {
    if (is.null(w)) {
      res <- glassoFast::glassoFast(S, lambda[i], start='cold')
    } else {
      res <- glassoFast::glassoFast(S, lambda[i], start='warm', w.init=w, wi.init=wi)
    }
    w <- res$w
    wi <- res$wi
    diag(res$wi) <- 0
    sol[[i]] <- Matrix::drop0(res$wi != 0)
  }
  return(sol)
}

glasso_strength <- function(x, lambda_max, net) {
  snet <- sum(net)
  d <- ncol(x)
  w <- NULL
  wi <- NULL
  S <- cor(x)
  lambda <- lambda_max
  repeat {
    if (is.null(w)) {
      res <- glassoFast::glassoFast(S, lambda, start='cold')
    } else {
      res <- glassoFast::glassoFast(S, lambda, start='warm', w.init=w, wi.init=wi)
    }
    w <- res$w
    wi <- res$wi
    diag(res$wi) <- 0
    sol <- Matrix::drop0(res$wi != 0)
    if (sum(sol) >= snet) {
      break
    }
    lambda <- 0.99 * lambda
  }
  return(list(w=w, wi=wi))
}

neigh_sel <- function(x, lambda) {
  d <- ncol(x)
  sol <- list()
  for (i in 1:length(lambda)) {
    sol[[i]] <- Matrix::Matrix(neighborhood_selection(x, lambda[i]))
  }
  return(sol)
}

read_grab_mask <- function(filename) {
  lines <- readLines(filename)
  blocksstr <- strsplit(lines, ' ')
  blocks <- lapply(blocksstr, function(s) as.numeric(s))
  d <- max(sapply(blocks, function(is) max(is)))
  mask <- matrix(F, d+1, d+1)
  for (i in 1:length(blocks)) {
    mask[blocks[[i]] + 1, blocks[[i]] + 1] <- T
  }
  return(mask)
}
