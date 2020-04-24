compute_fdr <- function(net, sol, mask) {
  ut <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)]
  }
  div0 <- function(num, denom) {
    if(denom == 0) {
      return(0)
    }
    return(num/denom)
  }
  sapply(rev(sort(unique(ut(sol)))), function(thr) c(
    div0(sum(ut(sol >= thr)), length(ut(sol))), #density
    div0(sum(ut(sol >= thr & !net)), sum(ut(sol >= thr))), #fdr
    div0(sum(ut(sol >= thr & !net & mask)), sum(ut(sol >= thr & mask))), #fdr inside
    div0(sum(ut(sol >= thr & !net & !mask)), sum(ut(sol >= thr & !mask))) #fdr between
  ))
}

compute_fdr1_nomask <- function(net, sol) {
  ut <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)]
  }
  div0 <- function(num, denom) {
    if(denom == 0) {
      return(0)
    }
    return(num/denom)
  }
  c(
    div0(sum(ut(sol)), length(ut(sol))), #density
    div0(sum(ut(sol & !net)), sum(ut(sol))), #fdr
    div0(sum(ut(sol & !net)), sum(ut(sol))), #fdr inside
    div0(sum(ut(sol & !net)), sum(ut(sol))) #fdr between
  )
}

compute_fdr1 <- function(net, sol, mask) {
  ut <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)]
  }
  div0 <- function(num, denom) {
    if(denom == 0) {
      return(0)
    }
    return(num/denom)
  }
  c(
    div0(sum(ut(sol)), length(ut(sol))), #density
    div0(sum(ut(sol & !net)), sum(ut(sol))), #fdr
    div0(sum(ut(sol & !net & mask)), sum(ut(sol & mask))), #fdr inside
    div0(sum(ut(sol & !net & !mask)), sum(ut(sol & !mask))) #fdr between
  )
}

compute_bias <- function(net, sol, mask) {
  ut <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)]
  }
  div0 <- function(num, denom) {
    if(denom == 0) {
      return(0)
    }
    return(num/denom)
  }
  #fnr_i <- div0(sum(ut(!sol & net & mask)), sum(ut(net & mask)))
  #fpr_o <- div0(sum(ut(sol & !net & !mask)), sum(ut(!net & !mask)))
  #(fnr_i + fpr_o) / 2
  fp_i <- sum(ut(sol & !net & mask))
  fp_o <- sum(ut(sol & !net & !mask))
  fn_i <- sum(ut(!sol & net & mask))
  fn_o <- sum(ut(!sol & net & !mask))
  p_i <- sum(ut(net & mask))
  p_o <- sum(ut(net & !mask))
  abs(fp_o-fn_o)/2/p_o + abs(fp_i-fn_i)/2/p_i
}

compute_biases <- function(net, sol, mask) {
  ut <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)]
  }
  div0 <- function(num, denom) {
    if(denom == 0) {
      return(0)
    }
    return(num/denom)
  }
  sapply(rev(sort(unique(ut(sol)))), function(thr) {
    solthr <- sol >= thr
    fp_i <- sum(ut(solthr & !net & mask))
    fp_o <- sum(ut(solthr & !net & !mask))
    fn_i <- sum(ut(!solthr & net & mask))
    fn_o <- sum(ut(!solthr & net & !mask))
    p_i <- sum(ut(net & mask))
    p_o <- sum(ut(net & !mask))
    abs(fp_o-fn_o)/2/p_o + abs(fp_i-fn_i)/2/p_i
  })
}

compute_roc <- function(net, sol) {
  ut <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)]
  }
  div0 <- function(num, denom) {
    if(denom == 0) {
      return(0)
    }
    return(num/denom)
  }
  snnet <- sum(!net)
  snet <- sum(net)
  thrnnet <- sol * (!net)
  thrnet <- sol * net
  sapply(rev(sort(unique(ut(sol)))), function(thr) c(
    div0(sum(thrnnet >= thr), snnet), #fpr
    div0(sum(thrnet >= thr), snet) #tpr
  ))
}

compute_roc_solutions <- function(net, solutions) {
  ut <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)]
  }
  div0 <- function(num, denom) {
    if(denom == 0) {
      return(0)
    }
    return(num/denom)
  }
  sapply(solutions, function(sol) c(
    div0(sum(ut(sol & !net)), sum(ut(!net))), #fpr
    div0(sum(ut(sol & net)), sum(ut(net))) #tpr
  ))
}

compute_roc1 <- function(net, sol) {
  ut <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)]
  }
  div0 <- function(num, denom) {
    if(denom == 0) {
      return(0)
    }
    return(num/denom)
  }
  c(
    div0(sum(ut(sol & !net)), sum(ut(!net))), #fpr
    div0(sum(ut(sol & net)), sum(ut(net))) #tpr
  )
}

compute_mcc <- function(net, sol) {
  ut <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)]
  }
  net <- ut(net)
  sol <- ut(sol)
  tp <- sum(sol != 0 & net != 0)
  fp <- sum(sol != 0 & net == 0)
  tn <- sum(sol == 0 & net == 0)
  fn <- sum(sol == 0 & net != 0)
  n <- length(sol)
  s <- (tp + fn)/n
  p <- (tp + fp)/n
  mcc <- (tp/n-s*p)/sqrt(s*p*(1-s)*(1-p))
  if (is.infinite(mcc) || is.na(mcc)) {
    mcc <- 0
  }
  return(mcc)
}

compute_mcc_mask <- function(net, sol, mask) {
  ut <- function(x) {
    x <- as.matrix(x)
    x[upper.tri(x)]
  }
  net <- ut(net)
  sol <- ut(sol)
  mask <- ut(mask)
  tp <- sum(sol != 0 & net != 0 & mask)
  fp <- sum(sol != 0 & net == 0 & mask)
  tn <- sum(sol == 0 & net == 0 & mask)
  fn <- sum(sol == 0 & net != 0 & mask)
  n <- sum(mask)
  s <- (tp + fn)/n
  p <- (tp + fp)/n
  mcc <- (tp/n-s*p)/sqrt(s*p*(1-s)*(1-p))
  if (is.infinite(mcc) || is.na(mcc)) {
    mcc <- 0
  }
  return(mcc)
}

compute_auc <- function(roc) {
  to <- c(0.01, 0.05, 0.1, 0.5, 1)
  sapply(to, function(to) {
    MESS::auc(c(0, roc[1, ]), c(0, roc[2, ]), to=to)/to
  })
}


compute_adjrand <- function(net) {
  g <- igraph::graph_from_adjacency_matrix(as.matrix(net), mode='undirected')
  cl <- igraph::cluster_louvain(g)$membership
  cl_true <- sort(rep(1:20, 20))
  mclust::adjustedRandIndex(cl, cl_true)
}

compute_adjrand_vsnet <- function(net, truenet) {
  g <- igraph::graph_from_adjacency_matrix(net, mode='undirected')
  cl <- igraph::cluster_louvain(g)$membership
  g <- igraph::graph_from_adjacency_matrix(truenet, mode='undirected')
  cl_true <- igraph::cluster_louvain(g)$membership
  mclust::adjustedRandIndex(cl, cl_true)
}
