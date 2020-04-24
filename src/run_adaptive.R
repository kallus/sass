source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/real.R')
source('src/experiment/scores.R')

B <- 100
lambda_ratio <- 0.1

to_w <- function(theta, u=0.1, v=0.3) {
	omega = as.matrix(theta)*v
	diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
	stats::cov2cor(solve(omega))
}

pcor <- function(corr) {
  P <- solve(corr)
  Pii <- diag(P) %*% matrix(1, 1, ncol(P))
  -P / sqrt(Pii * t(Pii))
}

# run contending methods on data files
cor_rmse_sassa <- c()
cor_rmse_glasso <- c()
pcor_rmse_sassa <- c()
pcor_rmse_glasso <- c()
for (datafile in list.files(path='cache/data/')) {
  if (!startsWith(datafile, 'cluster')) {
    next
  }
  print(datafile)
  net <- readRDS(paste('cache/data', datafile, sep='/'))
  n <- nrow(net$data)
  lambda_max <- cached(c('find_lambda_max', datafile), net$data)
  set.seed(1)
  subsamples <- subsample_observations(n, B)
  networks <- cached(c('subsample_neighsel', datafile),
    net$data, lambda_ratio*lambda_max, B, subsamples)
  stabsel_freq <- cached(c('selection_freq', datafile), networks)
  threshold <- cached(c('FDR_threshold', 'stabsel', datafile),
    stabsel_freq, FDR=0.1, length(networks))
  stabsel_net <- cached(c('threshold_frequencies', 'stabsel', datafile),
    stabsel_freq, threshold)

  sassa_freq <- stabsel_freq
  sassa_net <- stabsel_net
  comm_counts <- cached(c('consensus_module_counts', 'sass', datafile), networks, B)
  comm_threshold <- cached(c('smooth_inner_minima', 'sass', datafile),
    uptri(comm_counts), 0, B)
  if (length(comm_threshold) != 1) {
    comm_threshold <- Inf
  }
  cmask <- comm_counts > comm_threshold
  if (any(cmask)) {
    inside_min <- cached(c('smooth_inner_minima', 'sass', 'inside', datafile),
      2*B*uptri(stabsel_freq)[uptri(cmask)], 0, 2*B)
    if (length(inside_min) == 1) {
      between_min <- cached(c('best_smooth_min', 'sass', 'between', datafile),
        2*B*uptri(stabsel_freq)[uptri(!cmask)], 0, 2*B)
      if (!is.na(between_min) && inside_min < between_min) {
        sassa_freq <- cached(c('adapt_frequencies_mask2', 'sass', datafile),
          stabsel_freq, inside_min, between_min, cmask, B)
      }
    }
  }
  sassa_wwi <- cached(c('sass_adaptive_strength', 'sassa', datafile),
    net$data, lambda_max, sassa_freq, stabsel_net)

  glasso_wwi <- cached(c('glasso_strength', datafile), net$data, lambda_max,
    stabsel_net)

  true_wwi <- list()
  true_wwi$w=to_w(net$theta)
  true_wwi$wi=solve(true_wwi$w)

  # compare partial correlations and/or correlations
  # RMSE
  cor_rmse_sassa <- c(cor_rmse_sassa,
    sqrt(mean(c(true_wwi$w - cov2cor(sassa_wwi$w))^2)))
  cor_rmse_glasso <- c(cor_rmse_glasso,
    sqrt(mean(c(true_wwi$w - cov2cor(glasso_wwi$w))^2)))
  pcor_rmse_sassa <- c(pcor_rmse_sassa,
    sqrt(mean(c(pcor(true_wwi$w) - pcor(cov2cor(sassa_wwi$w)))^2)))
  pcor_rmse_glasso <- c(pcor_rmse_glasso,
    sqrt(mean(c(pcor(true_wwi$w) - pcor(cov2cor(glasso_wwi$w)))^2)))
}

write.table(cbind(cor_rmse_sassa, cor_rmse_glasso, pcor_rmse_sassa,
  pcor_rmse_glasso), 'cache/glasso_sassa_cmp.tsv')
