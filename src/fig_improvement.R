source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/scores.R')
source('src/experiment/graphics.R')

graphs <- c('cluster', 'overlapping', 'scale-free')
reps <- 100

mccs <- list()
for (graph in graphs) {
  mccs[[graph]] <- list()
  mccs[[graph]][['neighsel']] <- rep(NA, reps)
  for (repetition in 1:reps) {
    outputfile <- paste('stabsel_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc_denom <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
    outputfile <- paste('sass_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
    mccs[[graph]][['neighsel']][repetition] <- mcc/mcc_denom
  }
  mccs[[graph]][['glasso']] <- rep(NA, reps)
  for (repetition in 1:reps) {
    outputfile <- paste('stabselglasso_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc_denom <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
    outputfile <- paste('sass-glasso_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
    mccs[[graph]][['glasso']][repetition] <- mcc/mcc_denom
  }
  mccs[[graph]][['sassadaptive']] <- rep(NA, reps)
  for (repetition in 1:reps) {
    outputfile <- paste('sassadaptive_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc_denom <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
    outputfile <- paste('sass_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
    mccs[[graph]][['sassadaptive']][repetition] <- mcc/mcc_denom
  }
}

mccs_inside <- list()
for (graph in graphs) {
  if (graph == 'scale-free') {
    next
  }
  mccs_inside[[graph]] <- list()
  mccs_inside[[graph]][['neighsel']] <- rep(NA, reps)
  for (repetition in 1:reps) {
    outputfile <- paste('stabsel_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc_denom <- cached(c('compute_mcc_mask', outputfile, 'inside'), res$truenet$theta, res$estnet,
      res$truenet$cmask)
    outputfile <- paste('sass_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc <- cached(c('compute_mcc_mask', outputfile, 'inside'), res$truenet$theta, res$estnet,
      res$truenet$cmask)
    mccs_inside[[graph]][['neighsel']][repetition] <- mcc/mcc_denom
  }
  mccs_inside[[graph]][['glasso']] <- rep(NA, reps)
  for (repetition in 1:reps) {
    outputfile <- paste('stabselglasso_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc_denom <- cached(c('compute_mcc_mask', outputfile, 'inside'), res$truenet$theta, res$estnet,
      res$truenet$cmask)
    outputfile <- paste('sass-glasso_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc <- cached(c('compute_mcc_mask', outputfile, 'inside'), res$truenet$theta, res$estnet,
      res$truenet$cmask)
    mccs_inside[[graph]][['glasso']][repetition] <- mcc/mcc_denom
  }
  mccs_inside[[graph]][['sassadaptive']] <- rep(NA, reps)
  for (repetition in 1:reps) {
    outputfile <- paste('sassadaptive_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc_denom <- cached(c('compute_mcc_mask', outputfile, 'inside'), res$truenet$theta, res$estnet,
      res$truenet$cmask)
    outputfile <- paste('sass_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc <- cached(c('compute_mcc_mask', outputfile, 'inside'), res$truenet$theta, res$estnet,
      res$truenet$cmask)
    mccs_inside[[graph]][['sassadaptive']][repetition] <- mcc/mcc_denom
  }
}

mccs_outside <- list()
for (graph in graphs) {
  if (graph == 'scale-free') {
    next
  }
  mccs_outside[[graph]] <- list()
  mccs_outside[[graph]][['neighsel']] <- rep(NA, reps)
  for (repetition in 1:reps) {
    outputfile <- paste('stabsel_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc_denom <- cached(c('compute_mcc_mask', outputfile, 'outside'), res$truenet$theta, res$estnet,
      !res$truenet$cmask)
    outputfile <- paste('sass_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc <- cached(c('compute_mcc_mask', outputfile, 'outside'), res$truenet$theta, res$estnet,
      !res$truenet$cmask)
    mccs_outside[[graph]][['neighsel']][repetition] <- mcc/mcc_denom
  }
  mccs_outside[[graph]][['glasso']] <- rep(NA, reps)
  for (repetition in 1:reps) {
    outputfile <- paste('stabselglasso_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc_denom <- cached(c('compute_mcc_mask', outputfile, 'outside'), res$truenet$theta, res$estnet,
      !res$truenet$cmask)
    outputfile <- paste('sass-glasso_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc <- cached(c('compute_mcc_mask', outputfile, 'outside'), res$truenet$theta, res$estnet,
      !res$truenet$cmask)
    mccs_outside[[graph]][['glasso']][repetition] <- mcc/mcc_denom
  }
  mccs_outside[[graph]][['sassadaptive']] <- rep(NA, reps)
  for (repetition in 1:reps) {
    outputfile <- paste('sassadaptive_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc_denom <- cached(c('compute_mcc_mask', outputfile, 'outside'), res$truenet$theta, res$estnet,
      !res$truenet$cmask)
    outputfile <- paste('sass_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    mcc <- cached(c('compute_mcc_mask', outputfile, 'outside'), res$truenet$theta, res$estnet,
      !res$truenet$cmask)
    mccs_outside[[graph]][['sassadaptive']][repetition] <- mcc/mcc_denom
  }
}

pdf('fig/improvement.pdf', width=8, height=8)

par(mfrow=c(3, 3))
for (graph in graphs) {
  boxplot(mccs[[graph]], main=graph, ylab='SASS MCC / MCC', vertical=T, las=2,
    outpch=NA)
  mtext('All', side=4)
  stripchart(mccs[[graph]], vertical=T, add=T, pch=20, method='jitter')
  abline(h=1, lty=2)
}
for (graph in graphs) {
  if (graph == 'scale-free') {
    plot.new()
    next
  }
  boxplot(mccs_inside[[graph]], main=graph, ylab='SASS MCC / MCC', vertical=T,
    las=2, outpch=NA)
  mtext('Inside', side=4)
  stripchart(mccs_inside[[graph]], vertical=T, add=T, pch=20, method='jitter')
  abline(h=1, lty=2)
}
for (graph in graphs) {
  if (graph == 'scale-free') {
    plot.new()
    next
  }
  boxplot(mccs_outside[[graph]], main=graph, ylab='SASS MCC / MCC', vertical=T,
    las=2, outpch=NA)
  mtext('Between', side=4)
  stripchart(mccs_outside[[graph]], vertical=T, add=T, pch=20, method='jitter')
  abline(h=1, lty=2)
}

dev.off()
