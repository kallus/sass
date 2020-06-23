source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/scores.R')
source('src/experiment/graphics.R')

graph <- 'cluster'
methods <- c('neighsel', 'stabsel', 'sass', 'glasso', 'stabselglasso',
  'sass-glasso', 'sasskncl', 'grab')
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

mccs_inside <- list()
mccs_inside[[graph]] <- list()
for (method in methods) {
  mccs_inside[[graph]][[method]] <- rep(NA, reps)
  for (repetition in 1:reps) {
    outputfile <- paste(method, '_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc <- cached(c('compute_mcc_mask', outputfile, 'inside'), res$truenet$theta, res$estnet,
      res$truenet$cmask)
    mccs_inside[[graph]][[method]][repetition] <- mcc
  }
}

mccs_outside <- list()
mccs_outside[[graph]] <- list()
for (method in methods) {
  mccs_outside[[graph]][[method]] <- rep(NA, reps)
  for (repetition in 1:reps) {
    outputfile <- paste(method, '_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc <- cached(c('compute_mcc_mask', outputfile, 'outside'), res$truenet$theta, res$estnet,
      !res$truenet$cmask)
    mccs_outside[[graph]][[method]][repetition] <- mcc
  }
}

mccs_diff <- list()
mccs_diff[[graph]] <- list()
for (method in methods) {
  mccs_diff[[graph]][[method]] <- rep(NA, reps)
  for (repetition in 1:reps) {
    mccs_diff[[graph]][[method]][repetition] <- 
      mccs_inside[[graph]][[method]][repetition] /
      mccs_outside[[graph]][[method]][repetition]
  }
}

# print
for (method in methods) {
  mn <- mean(mccs[[graph]][[method]])
  print(paste('MCC', method, ':', round(mn, 3)))
  print(paste('MCC confint', method, ':', round(mn-confint(lm(mccs[[graph]][[method]]~1))[1], 3)))
  mn <- mean(mccs_diff[[graph]][[method]])
  print(paste('diff', method, ':', round(mn, 3)))
  print(paste('diff confint', method, ':', round(mn-confint(lm(mccs_diff[[graph]][[method]]~1))[1], 3)))
}

# remove
for (method in c('sasskncl', 'grab')) {
  mccs[[graph]][[method]] <- NULL
  mccs_diff[[graph]][[method]] <- NULL
}

tmp <- mccs
mccs[[graph]] <- list()
mccs[[graph]][['NS']] <- tmp[[graph]][['neighsel']]
mccs[[graph]][['CPSS (NS)']] <- tmp[[graph]][['stabsel']]
mccs[[graph]][['SASS (NS)']] <- tmp[[graph]][['sass']]
mccs[[graph]][['GL']] <- tmp[[graph]][['glasso']]
mccs[[graph]][['CPSS (GL)']] <- tmp[[graph]][['stabselglasso']]
mccs[[graph]][['SASS (GL)']] <- tmp[[graph]][['sass-glasso']]
tmp_diff <- mccs_diff
mccs_diff[[graph]] <- list()
mccs_diff[[graph]][['NS']] <- tmp_diff[[graph]][['neighsel']]
mccs_diff[[graph]][['CPSS (NS)']] <- tmp_diff[[graph]][['stabsel']]
mccs_diff[[graph]][['SASS (NS)']] <- tmp_diff[[graph]][['sass']]
mccs_diff[[graph]][['GL']] <- tmp_diff[[graph]][['glasso']]
mccs_diff[[graph]][['CPSS (GL)']] <- tmp_diff[[graph]][['stabselglasso']]
mccs_diff[[graph]][['SASS (GL)']] <- tmp_diff[[graph]][['sass-glasso']]

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

pdf('fig/sim_cluster.pdf', width=6.5, height=3)
layout(matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3), 1))
par(mar=c(6, 4, 4, 2))
boxplot(mccs[[graph]], main='A', ylab='MCC', vertical=T, las=2,
  outpch=NA)
stripchart(mccs[[graph]], vertical=T, pch=20, add=T, method='jitter')
boxplot(mccs_diff[[graph]], main='B', ylab='MCC within / MCC between',
  vertical=T, las=2, outpch=NA)
stripchart(mccs_diff[[graph]], vertical=T, pch=20, add=T, method='jitter')
abline(h=1, lty=2)
boxplot(mccs_impr[[graph]], main='C', ylab='SASS MCC / CPSS MCC', vertical=T, las=2,
  outpch=NA)
stripchart(mccs_impr[[graph]], vertical=T, add=T, pch=20, method='jitter')
abline(h=1, lty=2)
dev.off()

for (name in names(mccs_impr[[graph]])) {
  print(name)
  print(t.test(mccs_impr[[graph]][[name]]-1, alternative='greater'))
}
