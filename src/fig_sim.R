source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/scores.R')
source('src/experiment/graphics.R')

graphs <- c('cluster', 'overlapping', 'scale-free')
methods <- c('neighsel', 'stabsel', 'sass', 'sasskncl', 'glasso', 'stabselglasso',
  'sass-glasso', 'sassadaptive')#, 'grab')
reps <- 100

mccs <- list()
for (graph in graphs) {
  mccs[[graph]] <- list()
  for (method in methods) {
    if (method %in% c('grab', 'sasskncl') && graph == 'scale-free') {
      next
    }
    mccs[[graph]][[method]] <- rep(NA, reps)
    for (repetition in 1:reps) {
      outputfile <- paste(method, '_', graph, '_300_', repetition, sep='')
      res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
      mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
      mccs[[graph]][[method]][repetition] <- mcc
    }
  }
}

mccs_inside <- list()
for (graph in graphs) {
  if (graph == 'scale-free') {
    next
  }
  mccs_inside[[graph]] <- list()
  for (method in methods) {
    if (method == 'grab' && graph == 'scale-free') {
      next
    }
    mccs_inside[[graph]][[method]] <- rep(NA, reps)
    for (repetition in 1:reps) {
      outputfile <- paste(method, '_', graph, '_300_', repetition, sep='')
      res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
      mcc <- cached(c('compute_mcc_mask', outputfile, 'inside'), res$truenet$theta, res$estnet,
        res$truenet$cmask)
      mccs_inside[[graph]][[method]][repetition] <- mcc
    }
  }
}

mccs_outside <- list()
for (graph in graphs) {
  if (graph == 'scale-free') {
    next
  }
  mccs_outside[[graph]] <- list()
  for (method in methods) {
    if (method == 'grab' && graph == 'scale-free') {
      next
    }
    mccs_outside[[graph]][[method]] <- rep(NA, reps)
    for (repetition in 1:reps) {
      outputfile <- paste(method, '_', graph, '_300_', repetition, sep='')
      res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
      mcc <- cached(c('compute_mcc_mask', outputfile, 'outside'), res$truenet$theta, res$estnet,
        !res$truenet$cmask)
      mccs_outside[[graph]][[method]][repetition] <- mcc
    }
  }
}

pdf('fig/sim.pdf', width=8, height=8)
par(mfrow=c(3, 3))
for (graph in graphs) {
  boxplot(mccs[[graph]], main=graph, ylab='MCC', vertical=T, las=2,
    outpch=NA)
  mtext('All', side=4)
  stripchart(mccs[[graph]], vertical=T, pch=20, add=T, method='jitter')
}
for (graph in graphs) {
  if (graph == 'scale-free') {
    plot.new()
    next
  }
  boxplot(mccs_inside[[graph]], main=graph, ylab='MCC', vertical=T, las=2,
    outpch=NA)
  mtext('Inside', side=4)
  stripchart(mccs_inside[[graph]], vertical=T, pch=20, add=T, method='jitter')
}
for (graph in graphs) {
  if (graph == 'scale-free') {
    plot.new()
    next
  }
  boxplot(mccs_outside[[graph]], main=graph, ylab='MCC', vertical=T, las=2,
    outpch=NA)
  mtext('Between', side=4)
  stripchart(mccs_outside[[graph]], vertical=T, pch=20, add=T, method='jitter')
}
dev.off()
