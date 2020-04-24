source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/scores.R')
source('src/experiment/graphics.R')

plotroc <- function(graph, methods, repetition, show_legend=F) {
  col <- 0
  colors <- c(3, 2, 1, 4)
  # count methods
  for (method in methods) {
    if (method %in% c('grab')) {
      next
    }
    col <- col + 1
    outputfile <- paste(method, '_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    if (method %in% c('neighsel', 'glasso')) {
      roc <- res$roc
    } else {
      roc <- cached(c('compute_roc', outputfile), res$truenet$theta, res$freq)
    }
    roc1 <- cached(c('compute_roc1', outputfile), res$truenet$theta, res$estnet)
    lines(roc[1, ], roc[2, ], col=colors[col], lty=colors[col], lwd=2)
    points(roc1[1], roc1[2], col=1, lwd=2)
  }

  if ('grab' %in% methods) {
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
  }

  if (show_legend) {
    legend('bottomright', legend=methods, col=1:length(methods),
      lty=1:length(methods), lwd=2)
  }
}

graphs <- c('cluster', 'scale-free', 'overlapping')

ylimmin <- 0.2
xlimmax <- 0.05

pdf('fig/rocs.pdf', width=5.0, height=6.75)
par(mfrow=c(3, 2), mar=c(3, 4, 2, 1))

plot(-1, -1, xlab='FPR', ylab='TPR', xlim=c(0, xlimmax), ylim=c(ylimmin, 1), main='A',
  mgp=c(2, 1, 0))
plotroc('cluster', c('neighsel', 'stabsel', 'sass', 'sasskncl'), 1, show_legend=F)
legend('bottomright', c('NS', 'CPSS (NS)', 'SASS (NS)', 'SASS* (NS)', 'CPSS density'),
  lwd=2, col=c(3, 2, 1, 4, 1), lty=c(3, 2, 1, 4, NA), pch=c(NA, NA, NA, NA, 1))

plot(-1, -1, xlab='FPR', ylab='TPR', xlim=c(0, xlimmax), ylim=c(ylimmin, 1), main='B',
  mgp=c(2, 1, 0))
plotroc('cluster', c('glasso', 'stabselglasso', 'sass-glasso', 'grab'), 1, show_legend=F)
legend('bottomright', c('GL', 'CPSS (GL)', 'SASS (GL)', 'GRAB', 'CPSS density'),
  lwd=2, col=c(3, 2, 1, 4, 1), lty=c(3, 2, 1, 4, NA), pch=c(NA, NA, NA, NA, 1))

plot(-1, -1, xlab='FPR', ylab='TPR', xlim=c(0, xlimmax), ylim=c(ylimmin, 1), main='C',
  mgp=c(2, 1, 0))
plotroc('scale-free', c('neighsel', 'stabsel', 'sass'), 1, show_legend=F)
legend('bottomright', c('NS', 'CPSS (NS)', 'SASS (NS)', 'CPSS density'),
  lwd=2, col=c(3, 2, 1, 1), lty=c(3, 2, 1, NA), pch=c(NA, NA, NA, 1))

plot(-1, -1, xlab='FPR', ylab='TPR', xlim=c(0, xlimmax), ylim=c(ylimmin, 1), main='D',
  mgp=c(2, 1, 0))
plotroc('scale-free', c('glasso', 'stabselglasso', 'sass-glasso', 'grab'), 1, show_legend=F)
legend('bottomright', c('GL', 'CPSS (GL)', 'SASS (GL)', 'GRAB', 'CPSS density'),
  lwd=2, col=c(3, 2, 1, 4, 1), lty=c(3, 2, 1, 4, NA), pch=c(NA, NA, NA, NA, 1))

plot(-1, -1, xlab='FPR', ylab='TPR', xlim=c(0, xlimmax), ylim=c(ylimmin, 1), main='E',
  mgp=c(2, 1, 0))
plotroc('overlapping', c('neighsel', 'stabsel', 'sass', 'sasskncl'), 1, show_legend=F)
legend('bottomright', c('NS', 'CPSS (NS)', 'SASS (NS)', 'SASS* (NS)', 'CPSS density'),
  lwd=2, col=c(3, 2, 1, 4, 1), lty=c(3, 2, 1, 4, NA), pch=c(NA, NA, NA, NA, 1))

plot(-1, -1, xlab='FPR', ylab='TPR', xlim=c(0, xlimmax), ylim=c(ylimmin, 1), main='F',
  mgp=c(2, 1, 0))
plotroc('overlapping', c('glasso', 'stabselglasso', 'sass-glasso', 'grab'), 1, show_legend=F)
legend('bottomright', c('GL', 'CPSS (GL)', 'SASS (GL)', 'GRAB', 'CPSS density'),
  lwd=2, col=c(3, 2, 1, 4, 1), lty=c(3, 2, 1, 4, NA), pch=c(NA, NA, NA, NA, 1))

dev.off()
