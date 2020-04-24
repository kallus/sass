source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/real.R')
source('src/experiment/scores.R')

library(Matrix)

sif_to_adjmat <- function(sif, names) {
  d <- length(names)
  adj <- Matrix::Matrix(F, d, d, dimnames=list(names, names))
  for (i in 1:nrow(sif)) {
    cat(i / nrow(sif)); cat('\r')
    try(adj[sif[i, 1], sif[i, 3]] <- T, silent=T)
    try(adj[sif[i, 3], sif[i, 1]] <- T, silent=T)
  }
  return(adj)
}

step_adjmat <- function(adj, k) {
  ret <- list()
  m <- adj
  ret[[1]] <- adj != 0
  if (k >= 2) {
    for (i in 2:k) {
      m <- m %*% adj
      ret[[i]] <- (m + ret[[i-1]]) != 0
    }
  }
  return(ret)
}

sif_to_adjmats <- function(sif, keep_names, k) {
  names <- unique(c(sif[, 1], sif[, 3]))
  names <- union(names, keep_names)
  d <- length(names)
  adj <- Matrix::Matrix(F, d, d, dimnames=list(names, names))
  for (i in 1:nrow(sif)) {
    cat(i / nrow(sif)); cat('\r')
    adj[sif[i, 1], sif[i, 3]] <- T
    adj[sif[i, 3], sif[i, 1]] <- T
  }
  adjs <- step_adjmat(adj, k)
  lapply(adjs, function(adj) {
    adj[keep_names, keep_names]
  })
}

sass <- readRDS('cache/contender_output/sass_real')
genenames <- colnames(sass$truenet$data)

reactome <- read.table('data/PathwayCommons12.reactome.hgnc.sif.gz',
  stringsAsFactors=F)
#reactome_adjmat <- cached(c('sif_to_adjmat', 'reactome'), reactome, genenames)
max_k <- 1
#adjmats <- step_adjmat(reactome_adjmat, max_k)
adjmats <- cached(c('sif_to_adjmats', 'reactome'), reactome, genenames, max_k)

auc_diff <- function(a, b) {
  mean(a-b)
}

reps <- 100
B <- 100
lambda_ratio <- 0.1
datafile <- 'tcga_gbm_rnaseq'
diffs <- matrix(NA, reps, max_k)
stabsel_auroc <- matrix(NA, reps, max_k)
sass_auroc <- matrix(NA, reps, max_k)
for (i in 1:reps) {
  cat(i / reps); cat('\n')
  set.seed(i)
  net <- cached(c('real_data', i), T)
  n <- nrow(net$data)
  lambda_max <- cached(c('find_lambda_max', datafile, i), net$data)
  subsamples <- subsample_observations(n, B)
  networks <- cached(c('subsample_neighsel', datafile, i), net$data,
    lambda_ratio*lambda_max, B, subsamples)
  stabsel_freq <- cached(c('selection_freq', datafile, i), networks)
  comm_counts <- cached(c('consensus_module_counts', 'sass', datafile, i),
    networks, B)
  comm_threshold <- cached(c('smooth_inner_minima', 'sass', datafile, i),
    uptri(comm_counts), 0, B)
  if (length(comm_threshold) != 1) {
    comm_threshold <- Inf
  }
  cmask <- comm_counts > comm_threshold
  sass_freq <- stabsel_freq
  if (any(cmask)) {
    inside_min <- cached(
      c('smooth_inner_minima', 'sass', 'inside', datafile, i),
      2*B*uptri(stabsel_freq)[uptri(cmask)], 0, 2*B)
    if (length(inside_min) == 1) {
      between_min <- cached(
        c('best_smooth_min', 'sass', 'between', datafile, i),
        2*B*uptri(stabsel_freq)[uptri(!cmask)], 0, 2*B)
      if (!is.na(between_min) && inside_min < between_min) {
        sass_freq <- cached(c('adapt_frequencies_mask2', 'sass', datafile, i),
          stabsel_freq, inside_min, between_min, cmask, B)
      }
    }
  }
  parallel::mclapply(1:max_k, function(k) {
    stabsel_roc <- cached(c('compute_roc', 'stabsel', datafile, i, k),
      adjmats[[k]], stabsel_freq)
    sass_roc <- cached(c('compute_roc', 'sass', datafile, i, k),
      adjmats[[k]], sass_freq)
  })
  for (k in 1:max_k) {
    stabsel_roc <- cached(c('compute_roc', 'stabsel', datafile, i, k),
      adjmats[[k]], stabsel_freq)
    sass_roc <- cached(c('compute_roc', 'sass', datafile, i, k),
      adjmats[[k]], sass_freq)
    stabsel_roc <- stabsel_roc[, 1:(ncol(stabsel_roc)-1)]
    sass_roc <- sass_roc[, 1:(ncol(sass_roc)-1)]
    stabsel_roc <- cbind(c(0, 0), stabsel_roc, c(1, 1))
    sass_roc <- cbind(c(0, 0), sass_roc, c(1, 1))
    xout <- seq(0, 1, length=1000)
    stabsel_approx <- approx(stabsel_roc[1, ], stabsel_roc[2, ], xout)
    sass_approx <- approx(sass_roc[1, ], sass_roc[2, ], xout)
    stabsel_auroc[i, k] <- mean(stabsel_approx$y)
    sass_auroc[i, k] <- mean(sass_approx$y)
    diffs[i, k] <- cached(c('auc_diff', datafile, i, k), sass_approx$y,
      stabsel_approx$y)
  }
  print(diffs[i, ])
}

diffslist <- list()
stabsel_aurocs_list <- list()
sass_aurocs_list <- list()
for (k in 1:ncol(diffs)) {
  diffslist[[k]] <- diffs[, k]
  stabsel_aurocs_list[[k]] <- stabsel_auroc[, k]
  sass_aurocs_list[[k]] <- sass_auroc[, k]
}

tmp <- list()
tmp[['TCGA GBM RNASeq']] <- diffslist[[1]]
diffslist <- tmp

adjmat <- cached(c('sif_to_adjmats', 'reactome'), reactome, genenames, 1)[[1]]

lambda_ratio <- 0.1
datafile <- 'tcga_gbm_rnaseq_full'
net <- real_data()
n <- nrow(net$data)
lambda_max <- cached(c('find_lambda_max', datafile), net$data)
B <- 100
i <- 1
set.seed(i)
subsamples <- subsample_observations(n, B)
networks <- cached(c('subsample_neighsel', datafile, B), net$data,
  lambda_ratio*lambda_max, B, subsamples)
stabsel_freq <- cached(c('selection_freq', datafile, B), networks)
comm_counts <- cached(c('consensus_module_counts', 'sass', datafile, B),
  networks, B)
comm_threshold <- cached(c('smooth_inner_minima', 'sass', datafile, B),
  uptri(comm_counts), 0, B)
if (length(comm_threshold) != 1) {
  comm_threshold <- Inf
}
cmask <- comm_counts > comm_threshold
sass_freq <- stabsel_freq
if (any(cmask)) {
  inside_min <- cached(
    c('smooth_inner_minima', 'sass', 'inside', datafile, B),
    2*B*uptri(stabsel_freq)[uptri(cmask)], 0, 2*B)
  if (length(inside_min) == 1) {
    between_min <- cached(
      c('best_smooth_min', 'sass', 'between', datafile, B),
      2*B*uptri(stabsel_freq)[uptri(!cmask)], 0, 2*B)
    if (!is.na(between_min) && inside_min < between_min) {
      sass_freq <- cached(c('adapt_frequencies_mask2', 'sass', datafile, B),
        stabsel_freq, inside_min, between_min, cmask, B)
    }
  }
}
cpss <- cached(c('compute_roc', 'stabsel', datafile, B), adjmat, stabsel_freq)
sass <- cached(c('compute_roc', 'sass', datafile, B), adjmat,
  floor(2*B*sass_freq)/2/B)
cpss <- cbind(cpss[, 1:(ncol(cpss)-1)], c(1, 1))
sass <- cbind(sass[, 1:(ncol(sass)-1)], c(1, 1))

#pdf('fig/auc.pdf', width=1.19, height=1.88)
pdf('fig/auc.pdf', width=4, height=3)
layout(matrix(c(1, 1, 1, 2, 2), 1))
par(mar=c(6, 4, 4, 2))
plot(cpss[1, ], cpss[2, ], t='l', xlab='FPR', ylab='TPR', col='red', lty=2,
  xlim=c(0, 0.05), ylim=c(0, 0.1), mgp=c(2, 1, 0), lwd=2, main='A')
lines(sass[1, ], sass[2, ], lwd=2)
legend('bottomright', c('CPSS (NS)', 'SASS (NS)'), lwd=2, col=c(2, 1), lty=c(2, 1))
#layout(matrix(c(1, 2), 1))
#par(mar=c(1, 3, 1, 1), cex=0.65, mgp=c(2, 1, 0))
boxplot(diffslist[['TCGA GBM RNASeq']], main='B', vertical=T, outpch=NA)
title(ylab='AUC difference', line=2)
stripchart(diffslist, vertical=T, pch=20, add=T, method='jitter')
abline(h=0, lty=2)
dev.off()

for (name in names(diffslist)) {
  print(name)
  print(t.test(diffslist[[name]], alternative='greater'))
}

print(table(sign(diffslist[['TCGA GBM RNASeq']])))
