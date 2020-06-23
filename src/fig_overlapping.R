source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/scores.R')
source('src/experiment/graphics.R')

graph <- 'overlapping'

plotroc <- function(graph, repetition, show_legend=F) {
  col <- 0
  colors <- c(2, 1, 4)
  # count methods
  methods <- c('stabsel', 'sass')
  for (method in methods) {
    col <- col + 1
    outputfile <- paste(method, '_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    roc <- cached(c('compute_roc', outputfile), res$truenet$theta, res$freq)
    roc1 <- cached(c('compute_roc1', outputfile), res$truenet$theta, res$estnet)
    lines(roc[1, ], roc[2, ], col=colors[col], lty=colors[col], lwd=2)
    points(roc1[1], roc1[2], col=1, lwd=2)
  }

  # grab
  col <- col + 1
  outputfile <- paste('grab_', graph, '_300_', repetition, sep='')
  res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
  roc1 <- cached(c('compute_roc1', outputfile), res$truenet$theta, res$estnet)
  points(roc1[1], roc1[2], col=1, lwd=2)
  lambdas <- seq(0.48, 0.1, -0.02)
  roccache <- paste('cache/cache/grab_roc_curve', graph, repetition, sep='_')
  if (file.exists(roccache)) {
    roc <- readRDS(roccache)
  } else {
    roc <- matrix(NA, 2, length(lambdas))
    for (i in 1:length(lambdas)) {
      lambda <- lambdas[i]
      filename <- paste('cache/grab/output/net_20/', graph, '_300_', repetition, '_', lambda, sep='')
      if (file.exists(filename)) {
        freq <- as.matrix(read.table(filename))
        diag(freq) <- 0
        tmp <- compute_roc1(res$truenet$theta, freq != 0)
        roc[1, i] <- tmp[1]
        roc[2, i] <- tmp[2]
      }
    }
    saveRDS(roc, roccache)
  }
  lines(roc[1, ], roc[2, ], col=colors[col], lty=colors[col], lwd=2)
  methods <- c(methods, 'grab')

  if (show_legend) {
    legend('bottomright', legend=methods, col=1:length(methods),
      lty=1:length(methods), lwd=2)
  }
}

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
  ix <- is.finite(mccs_diff[[graph]][[method]])
  mn <- mean(mccs_diff[[graph]][[method]][ix])
  print(paste('diff', method, ':', round(mn, 3)))
  print(paste('diff confint', method, ':', round(mn-confint(lm(mccs_diff[[graph]][[method]][ix]~1))[1], 3)))
}

# remove
for (method in c('sasskncl', 'neighsel', 'glasso', 'stabselglasso',
    'sass-glasso')) {
  mccs[[graph]][[method]] <- NULL
}

tmp <- mccs
mccs[[graph]] <- list()
mccs[[graph]][['CPSS (NS)']] <- tmp[[graph]][['stabsel']]
mccs[[graph]][['SASS (NS)']] <- tmp[[graph]][['sass']]
mccs[[graph]][['GRAB']] <- tmp[[graph]][['grab']]

pdf('fig/overlapping.pdf', width=4, height=3)
layout(matrix(c(1, 1, 1, 2, 2), 1))
par(mar=c(6, 4, 4, 2))
reps <- 1
for (repetition in 1:reps) {
  # zoom out
  plot(-1, -1, xlab='FPR', ylab='TPR', xlim=c(0, 0.1), ylim=c(0, 1), main='A',
    mgp=c(2, 1, 0))
  plotroc(graph, repetition, show_legend=F)
  legend('bottomright', c('CPSS (NS)', 'SASS (NS)', 'GRAB', 'CPSS density'),
    lwd=2, col=c(2, 1, 4, 1), lty=c(2, 1, 4, NA), pch=c(NA, NA, NA, 1))
}
boxplot(mccs[[graph]], main='B', ylab='MCC', vertical=T, las=2,
  outpch=NA)
stripchart(mccs[[graph]], vertical=T, pch=20, add=T, method='jitter')
dev.off()
