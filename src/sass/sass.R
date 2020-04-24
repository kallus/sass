library(Matrix)
# lambda: factor of lambda_max

#lambda <- 0.1; B <- 100; verbose <- T; allow_res_skip <- F

#' @export
sass <- function(x, lambda=0.1, B=100, verbose=T, allow_res_skip=T) {
  subsamples <- subsample_observations(nrow(x), B)
  lambda <- lambda * find_lambda_max(x)
  
  catv(verbose, 'Subsampling neighborhood selection...')
  time <- sec_elapsed(networks <- subsample_neighsel(x, lambda, B, subsamples))
  catv(verbose, 'done.\t', fmttime(time), '\n')
  
  catv(verbose, 'Consensus community detection...')
  time <- sec_elapsed(modules <- consensus_modules(networks, 2*B))
  catv(verbose, 'done.\t', fmttime(time), '\n')

  catv(verbose, 'Adaptively thresholding network...')
  time <- sec_elapsed({
    freq <- selection_freq(networks)
    threshold <- FDR_threshold(freq, FDR=0.1, length(networks))
    network <- threshold_frequencies(freq, threshold)
    factors <- hierarchy_factors(2*B*freq, 2*B, modules$module)
    freq <- adapt_frequencies(freq, factors, modules$module, B)
    network <- adaptively_threshold_frequencies(freq, network)
  })
  catv(verbose, 'done.\t', fmttime(time), '\n')
  
  rownames(network) <- colnames(network) <- colnames(x)
  modules$tree$labels <- colnames(x)

  out <- list(network=network, modules=modules$module, order=modules$order,
    marks=modules$module_marks, tree=modules$tree, factors=factors, freq=freq)
  class(out) <- 'sass'
  return(out)
}

catv <- function(verbose, ...) {
  if (verbose) {
    cat(...)
  }
}

sec_elapsed <- function(...) {
  system.time(...)[3]
}

fmttime <- function(seconds) {
  seconds <- as.integer(seconds)
  time <- rep(0, 3)
  if (seconds >= 3600) {
    time[1] <- floor(seconds/3600)
    seconds <- seconds - 3600 * time[1]
  }
  if (seconds >= 60) {
    time[2] <- floor(seconds/60)
    seconds <- seconds - 60 * time[2]
  }
  time[3] <- seconds
  sprintf('\tTime consumed: %02d:%02d:%02d', time[1], time[2], time[3])
}

threshold_frequencies_match_density <- function(freq, net) {
  ut <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)]
  }
  edges <- sum(ut(net) != 0)
  utfreq <- ut(freq)
  f <- function(th) {
    s <- sum(utfreq >= th)
    if (s <= edges) {
      return(edges-s)
    }
    return(.Machine$double.xmax)
  }
  thr <- optimize(f, c(0, 1))$minimum
  #thr <- binsearch(function(th) sum(utfreq >= th), 0, 1, edges, 1e-4)
  #thr <- Inf
  #values <- unique(round(uptri(as.matrix(freq)), 3))
  #utfreq <- ut(freq)
  #for (th in rev(sort(values))) {
  #  print(th)
  #  if (sum(utfreq >= th) <= edges) {
  #    thr <- th
  #  } else {
  #    break
  #  }
  #}
  Matrix::drop0(freq >= thr)
}

binsearch <- function(f, low, high, val, tol) {
  repeat {
    mid <- (low+high)/2
    midval <- f(mid)
    print(paste(mid, midval))
    if (midval == val) {
      return(mid)
    }
    if (midval < val) {
      high <- mid
    } else {
      low <- mid
    }
    if (high-low < tol) {
      return(mid)
    }
  }
}

solutions_match_density <- function(solutions, net) {
  ut <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)]
  }
  edges <- sum(ut(net) != 0)
  for (sol in solutions) {
    if (sum(ut(sol)) >= edges) {
      return(Matrix::drop0(sol))
    }
  }
  Matrix::drop0(solutions[[length(solutions)]])
}

uptri <- function(x) {
  x <- as.matrix(x)
  x[upper.tri(x)]
}

hierarchy_factors_mask2 <- function(freq, max_freq, mask) {
  factors <- rep(NA, 2)
  smins <- smooth_inner_minima(uptri(freq)[uptri(mask)], 0, max_freq)
  if (length(smins) == 0) smins <- max_freq
  factors[2] <- smins[length(smins)]
  smins <- smooth_inner_minima(uptri(freq)[uptri(!mask)], 0, max_freq)
  if (length(smins) == 0) smins <- max_freq
  factors[1] <- smins[length(smins)]
  factors <- rev(stats::isoreg(rev(factors))$yf)
  return(factors)
}

#' @export
hierarchy_factors_mask <- function(freq, max_freq, mask) {
  factors <- rep(NA, 2)
  factors[2] <- smoothmin(uptri(freq)[uptri(mask)], 0, max_freq)
  factors[1] <- smoothmin(uptri(freq)[uptri(!mask)], 0, max_freq)
  factors <- rev(stats::isoreg(rev(factors))$yf)
  return(factors)
}

#' @export
hierarchy_factors <- function(freq, max_freq, cl) {
  d <- ncol(freq)
  mask <- matrix(F, d, d)
  factors <- rep(NA, ncol(cl)+1)
  for (i in ncol(cl):1) {
    mask <- !mask & outer(cl[, i], cl[, i], '==')
    factors[i+1] <- smoothmin(uptri(freq)[uptri(mask)], 0, max_freq)
    mask <- outer(cl[, i], cl[, i], '==')
  }
  mask <- !mask
  factors[1] <- smoothmin(uptri(freq)[uptri(mask)], 0, max_freq)
  factors <- rev(stats::isoreg(rev(factors))$yf)
  return(factors)
}

smoothmin <- function(data, lower, upper) {
  bw <- 1
  while(T) {
    dens <- density(data, bw=bw, from=lower, to=upper)
    min_ix <- which(diff(sign(diff(dens$y))) == 2)
    if (length(min_ix) == 1) {
      min_x <- dens$x[min_ix]
      return(min_x)
    }
    if (dens$bw > (upper-lower)/2 || length(min_ix) == 0) {
      warning('Could not find smooth inner minima, falling back to non-smooth')
      counts <- hist(data, plot=F, breaks=(lower:(upper+1))-0.5)$counts
      min_x <- which(counts == min(counts)) + lower - 1
      return(max(min_x))
    }
    bw <- bw + 1
  }
}

#' @export
adapt_frequencies_mask2 <- function(freq, inside, between, mask, B) {
  suppressMessages({ # this maybe inefficient according to Matrix package
    freq[!mask] <- reshape(freq[!mask], between/2/B, inside/2/B)
    #freq[!mask] <- freq[!mask]^(log(inside/2/B)/log(between/2/B))
  })
  freq <- methods::as(freq, 'symmetricMatrix')
}

#' @export
adapt_frequencies_mask <- function(freq, factors, mask, B) {
  suppressMessages({ # this maybe inefficient according to Matrix package
    freq[mask] <- reshape(freq[mask], factors[2]/2/B, 0.5)
    freq[!mask] <- reshape(freq[!mask], factors[1]/2/B, 0.5)
  })
  freq <- as(freq, 'symmetricMatrix')
}

#' @export
adapt_frequencies <- function(freq, factors, cl, B) {
  d <- ncol(freq)
  mask <- matrix(F, d, d)
  for (i in ncol(cl):1) {
    mask <- !mask & outer(cl[, i], cl[, i], '==')
    freq[mask] <- reshape(freq[mask], factors[i+1]/2/B, 0.5)
    mask <- outer(cl[, i], cl[, i], '==')
  }
  mask <- !mask
  freq[mask] <- reshape(freq[mask], factors[1]/2/B, 0.5)
  freq <- as(freq, 'symmetricMatrix')
}

#' @export
adaptively_threshold_frequencies <- function(freq, network) {
  threshold <- freq@x[rank(-freq@x, ties='first') == sum(network@x)]
  Matrix::drop0(freq >= threshold, is.Csparse=T)
}

reshape <- function(data, x1, x2) {
  newdata <- data
  newdata[data <= x1] <- data[data <= x1] * x2/x1
  newdata[data > x1] <- x2 + (data[data > x1] - x1) * (1-x2)/(1-x1)
  return(newdata)
}

#' @S3method plot mrsics
plot.mrsics <- function(x, ...) {
  immbw(as.matrix(x$network[x$order, x$order]), ...)
  d <- ncol(x$network)
  for (i in 1:length(x$marks)) {
    marks <- x$marks[[i]]
    col <- grDevices::rgb(0.1, 0.1, 0.1, sqrt(1/length(marks)))
    for (j in 1:length(marks)) {
      graphics::lines(c(0, d), rep(marks[j], 2), lwd=i, lty=1, col=col)
      graphics::lines(rep(marks[j], 2), c(0, d), lwd=i, lty=1, col=col)
    }
  }
}

resolutions <- function(d) {
  if (d <= 200) {
    return(ceiling(d/20))
  }
  if (d <= 1000) {
    return(ceiling(d/c(20, 200)))
  }
  resrec <- function(d) {
    if (d <= 25) {
      return(ceiling(d/5))
    }
    return(c(ceiling(d/5), resrec(ceiling(d/5))))
  }
  return(c(ceiling(d/c(20, 200)), resrec(ceiling(d/200))))
}

#' @export
consensus_modules <- function(networks, B) {
  d <- ncol(networks[[1]])
  cls <- parallel::mclapply(1:B, function(i) {
    network <- threshold_networks(networks, FDR=0.5, subsample=T)
    g <- igraph::graph_from_adjacency_matrix(network, mode='undirected')
    igraph::cluster_louvain(g)$membership
  })
  check_mclapply_errors(cls)
  counts <- matrix(0, d, d)
  n_clusters <- rep(NA, B)
  for (i in 1:B) {
    n_clusters[i] <- max(cls[[i]])
    counts <- counts + outer(cls[[i]], cls[[i]], '==')
  }
  hcl <- stats::hclust(stats::as.dist(1-counts/B), method='ward.D2')
  counts <- counts[upper.tri(counts)]
  counts <- counts/B
  cl <- stats::cutree(hcl, k=round(mean(n_clusters)))
  k <- 1
  if (length(k) == 1) {
    cl <- matrix(cl, ncol=1)
  }
  ix <- do.call(order, lapply(1:length(k), function(i) cl[, i]))
  clmarks <- lapply(length(k):1, function(i) {
    which(diff(cl[ix, i]) != 0) + 0.5
  })
  return(list(module=cl, order=ix, module_marks=clmarks, tree=hcl, counts=counts))
}

#' @export
consensus_module_counts <- function(networks, B) {
  d <- ncol(networks[[1]])
  cls <- parallel::mclapply(1:B, function(i) {
    network <- threshold_networks(networks, FDR=0.5, subsample=T)
    g <- igraph::graph_from_adjacency_matrix(as.matrix(network), mode='undirected')
    igraph::cluster_louvain(g)$membership
  })
  check_mclapply_errors(cls)
  counts <- matrix(0, d, d)
  n_clusters <- rep(NA, B)
  for (i in 1:B) {
    n_clusters[i] <- max(cls[[i]])
    counts <- counts + outer(cls[[i]], cls[[i]], '==')
  }
  return(counts)
}

##' @export
#consensus_cuthillmckee_counts <- function(networks, B) {
#  d <- ncol(networks[[1]])
#  counts <- matrix(0, d, d)
#  for (i in 1:B) {
#    network <- threshold_networks(networks, FDR=0.5, subsample=T)
#    ix <- cuthill_mckee(network)
#    invix <- order(ix)
#    counts <- counts +
#      inside_cuthill_mckee(as.matrix(network[ix, ix]))[invix, invix]
#  }
#  return(counts)
#}

##' @export
#consensus_specseq_counts <- function(networks, B) {
#  d <- ncol(networks[[1]])
#  counts <- matrix(0, d, d)
#  centrality <- 1-outer(1:d, 1:d, function(i, j) abs(i-j))/(d-1)
#  for (i in 1:B) {
#    network <- threshold_networks(networks, FDR=0.5, subsample=T)
#    ix <- specseq(network)
#    invix <- order(ix)
#    counts <- counts + centrality[invix, invix]
#  }
#  return(round(counts))
#}

#' @export
hub_mask_counts <- function(networks, B) {
  d <- ncol(networks[[1]])
  edges_per_node <- parallel::mclapply(1:B, function(i) {
    network <- threshold_networks(networks, FDR=0.5, subsample=T)
    Matrix::rowSums(network)
  })
  check_mclapply_errors(edges_per_node)
  counts <- matrix(0, d, d)
  for (i in 1:B) {
    counts <- counts + outer(edges_per_node[[i]], edges_per_node[[i]], pmax)
  }
  counts[upper.tri(counts)]
}

#' @export
hub_mask <- function(networks, B) {
  d <- ncol(networks[[1]])
  edges_per_node <- parallel::mclapply(1:B, function(i) {
    network <- threshold_networks(networks, FDR=0.5, subsample=T)
    Matrix::rowSums(network)
  })
  check_mclapply_errors(edges_per_node)
  counts <- rep(0, d)
  for (i in 1:B) {
    counts <- counts + edges_per_node[[i]]
  }
  counts <- counts/B
  sm <- smoothmin(counts, 0, max(1, max(counts)*1.1))
  if (any(counts > sm)) {
    mask <- matrix(F, d, d)
    mask[counts > sm, ] <- T
    mask[, counts > sm] <- T
    return(mask)
  } else {
    return(matrix(F, d, d))
  }
}

gaussian_bic <- function(x) {
  v <- var(x)
  n <- length(x)
  n*log(2*pi*v) + sum((x - mean(x))^2)/v + log(n)*2
}

check_mclapply_errors <- function(results) {
  for (i in 1:length(results)) {
    if (class(results[[i]]) == 'try-error') {
      print(results[[i]])
      stop('Error while executing in parallel')
    }
  }
}

#' @export
selection_freq <- function(networks, subsample=F) {
  d <- ncol(networks[[1]])
  freq <- dsCMatrix(d)
  if (subsample) {
    ix <- sample(length(networks), floor(length(networks)/2))
  } else {
    ix <- 1:length(networks)
  }
  for (i in ix) {
    freq <- freq + networks[[i]]
  }
  freq / length(ix)
}

#' @export
threshold_frequencies <- function(freq, threshold) {
  Matrix::drop0(freq >= threshold, is.Csparse=T)
}

#' @export
threshold_networks <- function(networks, FDR, subsample=F) {
  freq <- selection_freq(networks, subsample)
  B <- ifelse(subsample, floor(length(networks)/2), length(networks))
  thr <- FDR_threshold(freq, FDR, B)
  threshold_frequencies(freq, thr)
}

# The following function is adapted from package stabs on CRAN
# https://cran.r-project.org/package=stabs
# License: GPL-2
# Authors: Benjamin Hofner, Torsten Hothorn
stabs_maxQ <- function (p, B) {
  if (B <= 1) {
    stop("B must be at least 2", call. = FALSE)
  }
  fact_1 <- 4 * B/p
  tmpfct <- function(q) ceiling(q * fact_1) + 1 - 2 * B
  res <- tmpfct(1:p)
  length(res[res < 0])
}

# The following function is adapted from package stabs on CRAN
# https://cran.r-project.org/package=stabs
# License: GPL-2
# Authors: Benjamin Hofner, Torsten Hothorn
stabs_minD <- function (q, p, pi, B, r=c(-1/2, -1/4)) {
  which <- ceiling(signif(pi/(1/(2 * B)), 10))
  min(c(1, stabs_D(q^2/p^2, which - B, B, r[1]),
    stabs_D(q/p, which, 2 * B, r[2])))
}

# The following function is adapted from package stabs on CRAN
# https://cran.r-project.org/package=stabs
# License: GPL-2
# Authors: Benjamin Hofner, Torsten Hothorn
stabs_D <- function (theta, which, B, r) {
  s <- 1/r
  thetaB <- theta * B
  k_start <- (ceiling(2 * thetaB) + 1)
  if (which < k_start) {
    return(1)
  }
  if (k_start > B) {
    stop("theta to large")
  }
  Find.a <- function(prev_a) stats::uniroot(Calc.a, lower=1e-05, upper=prev_a,
    tol=.Machine$double.eps^0.75)$root
  Calc.a <- function(a) {
    denom <- sum((a + 0:k)^s)
    num <- sum((0:k) * (a + 0:k)^s)
    num/denom - thetaB
  }
  OptimInt <- function(a, t, k, thetaB, s) {
    num <- (k + 1 - thetaB) * sum((a + 0:(t - 1))^s)
    denom <- sum((k + 1 - (0:k)) * (a + 0:k)^s)
    1 - num/denom
  }
  a_vec <- rep(1e+05, B)
  for (k in k_start:B) a_vec[k] <- Find.a(a_vec[k - 1])
  cur_optim <- rep(0, B)
  for (k in k_start:(B - 1)) cur_optim[k] <- stats::optimize(f = OptimInt,
    lower = a_vec[k + 1], upper = a_vec[k], t = which, k = k,
    thetaB = thetaB, s = s, maximum = TRUE)$objective
  return(max(cur_optim))
}

#' @export
FDR_threshold <- function(freq, FDR, B) {
  d <- ncol(freq)
  q <- sum(freq@x)
  p <- choose(d, 2)
  thresholds <- (B*max(freq)):ceiling(q/p*B) / B
  fdr <- 0
  for (i in 1:length(thresholds)) {
    if (q > stabs_maxQ(p, B/2)) {
      stop('CPSS assumptions violated. Try with higher lambda.')
    }
    ev <- stabs_minD(q, p, thresholds[i], B/2) * p
    n_selected <- sum(freq@x >= thresholds[i])
    fdr <- min(max(ev/n_selected, fdr), 1) # fdr is increasing and <= 1
    if (fdr > FDR) {
      return(thresholds[max(1, i-1)])
    }
  }
  warning('Threshold too high due to assumptions in CPSS.',
    'Try with lower lambda.')
  return(thresholds[i])
}

#' @export
subsample_observations <- function(n, B) {
  # if n is odd make sure observations are left out fairly
  ix <- matrix(F, n, 2*B)
  if (n %% 2 == 1) {
    leftout <- rep(sample(n), ceiling(B/n))
    subn <- (n-1)/2
    ix <- parallel::mclapply(1:B, function(i) {
      iix <- matrix(F, n, 2)
      iix[-leftout[i], 1] <- 1:(n-1) %in% sample(n-1, subn) 
      iix[-leftout[i], 2] <- !iix[-leftout[i], 1]
      iix
    })
    check_mclapply_errors(ix)
  } else {
    subn <- n/2
    ix <- parallel::mclapply(1:B, function(i) {
      iix <- matrix(F, n, 2)
      iix[, 1] <- 1:n %in% sample(n, subn) 
      iix[, 2] <- !iix[, 1]
      iix
    })
    check_mclapply_errors(ix)
  }
  return(do.call(cbind, ix))
}

#' @export
find_lambda_max <- function(x) {
  d <- ncol(x)
  S <- stats::cor(x)
  max(max(S - diag(d)), -min(S - diag(d)))
}

#' @export
subsample_neighsel <- function(x, lambda, B, subsamples) {
  networks <- parallel::mclapply(1:(2*B), function(i) {
    #n <- sum(subsamples[, i])
    #net <- huge::huge(stats::cor(x[subsamples[, i], ]), lambda, sym='or',
    #  verbose=F, scr=T)$path[[1]]
    #  #verbose=F, scr=T, scr.num=n/log(n))$path[[1]]
    #Matrix::Matrix(net != 0, forceCheck=T)
    neighborhood_selection(x[subsamples[, i], ], lambda)
  })
  check_mclapply_errors(networks)
  return(networks)
}

#' @export
subsample_glasso <- function(x, lambda, B, subsamples) {
  networks <- parallel::mclapply(1:(2*B), function(i) {
    res <- glassoFast::glassoFast(cor(x[subsamples[, i], ]), lambda, start='cold')
    diag(res$wi) <- 0
    Matrix::Matrix(res$wi != 0)
  })
  check_mclapply_errors(networks)
  return(networks)
}

neighborhood_selection <- function(x, lambda) {
  d <- ncol(x)
  n <- nrow(x)
  g <- lgCMatrix(d, d)
  for (i in 1:d) {
    corr <- abs(stats::cor(x[, i], x[, -i]))
    corr_limit <- -sort(-corr)[ceiling(n/log(n))]
    exclude <- which(corr < corr_limit)
    g[i, -i] <- 0 != glmnet::glmnet(x[, -i], x[, i], lambda=lambda,
      exclude=which(corr < corr_limit))$beta[, 1]
  }
  return(Matrix::Matrix(g | Matrix::t(g), forceCheck=T))
}

lgCMatrix <- function(nrow, ncol) {
  methods::as(Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0),
    dims=c(nrow, ncol), symmetric=F), 'lgCMatrix')
}

lsCMatrix <- function(d) {
  methods::as(dsCMatrix(d), 'lsCMatrix')
}

dsCMatrix <- function(d) {
  Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(d, d),
    symmetric=T)
}

sass_adaptive <- function(x, lambda_max, freq, net) {
  snet <- sum(net)
  d <- ncol(x)
  w <- NULL
  wi <- NULL
  S <- cor(x)
  lambda <- lambda_max
  freq0 <- as.matrix(freq == 0)
  repeat {
    rho <- ifelse(freq0, lambda_max, as.matrix(lambda/freq))
    if (is.null(w)) {
      res <- glassoFast::glassoFast(S, rho, start='cold')
    } else {
      res <- glassoFast::glassoFast(S, rho, start='warm', w.init=w, wi.init=wi)
    }
    w <- res$w
    wi <- res$wi
    diag(res$wi) <- 0
    sol <- Matrix::drop0(res$wi != 0)
    if (any(sol & freq0)) {
      warning('edge infered where lambda is infinity')
    }
    if (sum(sol) >= snet) {
      break
    }
    lambda <- 0.99 * lambda
  }
  return(sol)
}

sass_adaptive_strength <- function(x, lambda_max, freq, net) {
  snet <- sum(net)
  d <- ncol(x)
  w <- NULL
  wi <- NULL
  S <- cor(x)
  lambda <- lambda_max
  freq0 <- as.matrix(freq == 0)
  repeat {
    rho <- ifelse(freq0, lambda_max, as.matrix(lambda/freq))
    if (is.null(w)) {
      res <- glassoFast::glassoFast(S, rho, start='cold')
    } else {
      res <- glassoFast::glassoFast(S, rho, start='warm', w.init=w, wi.init=wi)
    }
    w <- res$w
    wi <- res$wi
    diag(res$wi) <- 0
    sol <- Matrix::drop0(res$wi != 0)
    if (any(sol & freq0)) {
      warning('edge infered where lambda is infinity')
    }
    if (sum(sol) >= snet) {
      break
    }
    lambda <- 0.99 * lambda
  }
  return(list(w=w, wi=wi))
}

#specseq <- function(x) {
#  g <- igraph::graph_from_adjacency_matrix(x, mode='undirected')
#  components <- igraph::clusters(g)
#  new_order <- order(components$membership)
#  for (i in 1:components$no) {
#    if (components$csize[i] > 2) {
#      ix <- components$membership == i
#      fiedler <- RSpectra::eigs(diag(Matrix::rowSums(x[ix, ix])) - x[ix, ix], 2, 'SM')$vectors[, 1]
#      ix <- ix[new_order]
#      new_order[ix] <- new_order[ix][order(fiedler)]
#    }
#  }
#  return(new_order)
#}

#cuthill_mckee <- function(x) {
#  n <- ncol(x)
#  result <- rep(NA, n)
#  rnd <- runif(n, 1e-8 - 0.5, 0.5 - 1e-8)      # randomize order for equal degree
#  by_degree <- order(Matrix::rowSums(x) + rnd) # order nodes by degree
#  current <- by_degree[1]                      # make lowest degree node current
#  result[1] <- current                         # add current node
#  queue <- c()
#  for (j in 2:n) {
#    nb <- which(x[current, ])                  # neighbors of current node
#    nb <- intersect(by_degree, nb)             # order neighbors by degree
#    queue <- union(queue, nb)                  # add neighbors to queue
#    queue <- setdiff(queue, result)            # remove nodes in result from queue
#    if (length(queue) > 0) {
#      current <- queue[1]                      # make node first in queue current
#      queue <- queue[-1]
#    } else {
#      current <- setdiff(by_degree, result)[1] # make lowest degree unadded node current
#    }
#    result[j] <- current                       # add current node
#  }
#  return(result)
#}

#inside_cuthill_mckee <- function(x) {
#  n <- ncol(x)
#  for (i in 1:n) {
#    w <- which(x[i, ])
#    j <- w[1]
#    if (!is.na(j) && j < i) {
#      x[i, j:i] <- T
#      x[j:i, i] <- T
#    }
#    j <- rev(w)[1]
#    if (!is.na(j) && j > i) {
#      x[i, i:j] <- T
#      x[i:j, i] <- T
#    }
#  }
#  return(x)
#}

# visualize matrix in same order as it is printed
immbw <- function(x, asp=TRUE, legend.outside=NA, xlab='Column index',
    ylab='Row index', ...) {
  if (asp) {
    asp <- 1
  } else {
    asp <- NA
  }
  d <- dim(x)+1
  colors <- c('white', 'dodgerblue4')
  graphics::image(1:d[2]-0.5, 1:d[1]-0.5, t(x), xlim=c(1, d[2])-0.5,
    ylim=c(d[1], 1)-0.5, xlab=xlab, ylab=ylab, asp=asp, axes=FALSE,
    col=colors, ...)
  atx <- pretty(1:d[1])
  atx <- atx[atx < d[1]]
  atx <- c(atx, d[1]-1)
  aty <- pretty(1:d[2])
  aty <- aty[aty < d[2]]
  aty <- c(aty, d[2]-1)
  graphics::axis(1, aty, pos=d[1]-0.5)
  graphics::axis(2, atx, pos=0.5)
  graphics::axis(3, aty, labels=F, pos=0.5)
  graphics::axis(4, atx, labels=F, pos=d[2]-0.5)
}

immcomp <- function(x, xhat, asp=TRUE, legend.outside=NA, xlab='Column index',
    ylab='Row index', ...) {
  if (asp) {
    asp <- 1
  } else {
    asp <- NA
  }
  d <- dim(x)+1
  colors <- c('white', 'black')
  graphics::image(1:d[2]-0.5, 1:d[1]-0.5, t(x & xhat), xlim=c(1, d[2])-0.5,
    ylim=c(d[1], 1)-0.5, xlab=xlab, ylab=ylab, asp=asp, axes=FALSE,
    col=colors, ...)
  colors <- c(rgb(0, 0, 0, 0), 'red')
  graphics::image(1:d[2]-0.5, 1:d[1]-0.5, t(xhat & !x), xlim=c(1, d[2])-0.5,
    ylim=c(d[1], 1)-0.5, xlab=xlab, ylab=ylab, asp=asp, axes=FALSE, add=T,
    col=colors, ...)
  colors <- c(rgb(0, 0, 0, 0), 'deepskyblue')
  graphics::image(1:d[2]-0.5, 1:d[1]-0.5, t(x & !xhat), xlim=c(1, d[2])-0.5,
    ylim=c(d[1], 1)-0.5, xlab=xlab, ylab=ylab, asp=asp, axes=FALSE, add=T,
    col=colors, ...)
  atx <- pretty(1:d[1])
  atx <- atx[atx < d[1]]
  atx <- c(atx, d[1]-1)
  aty <- pretty(1:d[2])
  aty <- aty[aty < d[2]]
  aty <- c(aty, d[2]-1)
  graphics::axis(1, aty, pos=d[1]-0.5)
  graphics::axis(2, atx, pos=0.5)
  graphics::axis(3, aty, labels=F, pos=0.5)
  graphics::axis(4, atx, labels=F, pos=d[2]-0.5)
}
