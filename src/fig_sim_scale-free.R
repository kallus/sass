source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/scores.R')
source('src/experiment/graphics.R')

graph <- 'scale-free'
methods <- c('neighsel', 'stabsel', 'sass', 'glasso', 'stabselglasso',
  'sass-glasso', 'grab')
reps <- 100

mccs <- list()
mccs[[graph]] <- list()
for (method in methods) {
  mccs[[graph]][[method]] <- rep(NA, reps)
  for (repetition in 1:reps) {
    outputfile <- paste(method, '_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
    mccs[[graph]][[method]][repetition] <- mcc
  }
}

# print
for (method in methods) {
  mn <- mean(mccs[[graph]][[method]])
  print(paste('MCC', method, ':', round(mn, 3)))
  print(paste('MCC confint', method, ':', round(mn-confint(lm(mccs[[graph]][[method]]~1))[1], 3)))
}

# remove
for (method in c('grab')) {
  mccs[[graph]][[method]] <- NULL
}

tmp <- mccs
mccs[[graph]] <- list()
mccs[[graph]][['NS']] <- tmp[[graph]][['neighsel']]
mccs[[graph]][['CPSS (NS)']] <- tmp[[graph]][['stabsel']]
mccs[[graph]][['SASS (NS)']] <- tmp[[graph]][['sass']]
mccs[[graph]][['GL']] <- tmp[[graph]][['glasso']]
mccs[[graph]][['CPSS (GL)']] <- tmp[[graph]][['stabselglasso']]
mccs[[graph]][['SASS (GL)']] <- tmp[[graph]][['sass-glasso']]

mccs_impr <- list()
mccs_impr[[graph]] <- list()
mccs_impr[[graph]][['NS']] <- rep(NA, reps)
for (repetition in 1:reps) {
  outputfile <- paste('stabsel_', graph, '_300_', repetition, sep='')
  res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
  mcc_denom <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
  outputfile <- paste('sass_', graph, '_300_', repetition, sep='')
  res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
  mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
  mccs_impr[[graph]][['NS']][repetition] <- mcc/mcc_denom
}
mccs_impr[[graph]][['GL']] <- rep(NA, reps)
for (repetition in 1:reps) {
  outputfile <- paste('stabselglasso_', graph, '_300_', repetition, sep='')
  res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
  mcc_denom <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
  outputfile <- paste('sass-glasso_', graph, '_300_', repetition, sep='')
  res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
  mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
  mccs_impr[[graph]][['GL']][repetition] <- mcc/mcc_denom
}

pdf('fig/sim_scale-free.pdf', width=3.9, height=3)
layout(matrix(c(1, 1, 1, 1, 1, 2, 2, 2), 1))
par(mar=c(6, 4, 4, 2))
boxplot(mccs[[graph]], main='A', ylab='MCC', vertical=T, las=2,
  outpch=NA)
stripchart(mccs[[graph]], vertical=T, pch=20, add=T, method='jitter')
boxplot(mccs_impr[[graph]], main='B', ylab='SASS MCC / CPSS MCC', vertical=T, las=2,
  outpch=NA)
stripchart(mccs_impr[[graph]], vertical=T, add=T, pch=20, method='jitter')
abline(h=1, lty=2)
dev.off()
