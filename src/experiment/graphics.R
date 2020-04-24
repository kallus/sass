hist_fit <- function(x, minx, maxx) {
  B <- maxx-minx+1
  counts <- as.numeric(table(c(x, minx:maxx))-1)
  h <- parsearch(function(h) ndenscv(counts, h, minx, maxx), 0, 0.05*B)
  ndens(counts, h, minx, maxx)
}

hist_fit_nologit <- function(x, minx, maxx) {
  B <- maxx-minx+1
  counts <- as.numeric(table(c(x, minx:maxx))-1)
  h <- parsearch(function(h) ndenscv_nologit(counts, h, minx, maxx), 0, 0.05*B)
  ndens_nologit(counts, h, minx, maxx)
}

ndens_nologit <- function(x, h, minx, maxx) {
  ix <- minx:maxx-minx+1
  d <- abs(outer(ix, ix, '-'))
  K <- exp(-(d^2)/(2*h^2))
  K <- apply(K, 2, function(v) v/sum(v))
  K <- t(apply(K, 1, function(v) v/sum(v)))
  y <- K %*% x
  return(y/sum(y))
}

ndenscv_nologit <- function(x, h, minx, maxx) {
  ix <- minx:maxx-minx+1
  d <- abs(outer(ix, ix, '-'))
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

