smooth_inner_minima <- function(x, minx, maxx) {
  B <- maxx-minx+1
  counts <- as.numeric(table(c(x, minx:maxx))-1)
  h <- parsearch(function(h) ndenscv(counts, h, minx, maxx), 0, 0.05*B)
  y <- ndens(counts, h, minx, maxx)
  inner_minima(minx:maxx, y)
}

# exactly one inner minimum or
# exactly one inner minimum and one minimum at right end or
# no inner minimum and one minimum at right end
# return the minimum with smallest y or NA otherwise
best_smooth_min <- function(x, minx, maxx) {
  B <- maxx-minx+1
  counts <- as.numeric(table(c(x, minx:maxx))-1)
  h <- parsearch(function(h) ndenscv(counts, h, minx, maxx), 0, 0.05*B)
  y <- ndens(counts, h, minx, maxx)
  inmin <- inner_minima(minx:maxx, y)
  if (length(inmin) != 1) {
    inmin <- NA
  }
  rightmin <- y[length(y)] < y[length(y)-1]
  if (!rightmin) {
    return(inmin)
  }
  if (is.na(inmin)) {
    return(length(y))
  }
  if (y[inmin] <= y[length(y)]) {
    return(inmin)
  } else {
    return(length(y))
  }
}

ndens <- function(x, h, minx, maxx) {
  ix <- minx:maxx-minx+1
  d <- logitdist(ix)
  K <- exp(-(d^2)/(2*h^2))
  K <- apply(K, 2, function(v) v/sum(v))
  K <- t(apply(K, 1, function(v) v/sum(v)))
  y <- K %*% x
  return(y/sum(y))
}

logitdist <- function(ix) {
  n <- length(ix)
  pos <- seq(0, 1, length=n+2)[-c(1, n+2)]
  pos <- -log(1/pos-1)
  abs(outer(pos, pos, '-'))
}

ndenscv <- function(x, h, minx, maxx) {
  ix <- minx:maxx-minx+1
  d <- logitdist(ix)
  cv <- rep(0, length(h))
  for (j in 1:length(h)) {
    if (h[j] <= 0) {
      cv[j] <- -Inf
      next
    }
    K <- exp(-(d^2)/(2*h[j]^2))
    K <- apply(K, 2, function(v) v/sum(v))
    K <- t(apply(K, 1, function(v) v/sum(v)))
    y <- K %*% x
    for (i in ix) {
      if (x[i] > 0) {
        yi <- y - K[, i] * (x[i]-1)
        yi <- yi/sum(yi)
        cv[j] <- cv[j] + x[i]*log(yi[i])
      }
    }
  }
  return(cv)
}

parsearch <- function(f, minx, maxx) {
  initial_minx <- minx
  initial_maxx <- maxx
  for (i in 1:3) {
    x <- seq(minx, maxx, length=10)
    cv <- mcsapply(x, f)
    imax <- which.max(cv)
    minx <- x[max(1, imax-1)]
    maxx <- x[min(imax+1, length(cv))]
  }
  if (x[imax] == initial_minx) {
    warning('maximum at lower boundary')
  }
  if (x[imax] == initial_maxx) {
    warning('maximum at upper boundary')
  }
  return(x[imax])
}

mcsapply <- function(i, f) {
  simplify2array(parallel::mclapply(i, f))
}

inner_minima <- function(x, y) {
  x[which(diff(sign(diff(y))) == 2) + 1]
}
