source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/scores.R')
source('src/experiment/graphics.R')

hist <- function(data, main, ...) {
  x <- sort(unique(data))
  y <- rep(NA, length(x))
  for (i in 1:length(x)) {
    y[i] <- sum(data > x[i]-1e-4 & data < x[i]+1e-4)
  }
  plot(x, y, pch=20, ylab='Frequency', mgp=c(2, 1, 0), ...)
  title(main, line=0.5)
  return(length(data))
}

pdf('fig/hist.pdf', width=6.5, height=6.5)

B <- 100
outputfile <- 'sass_cluster_300_1'
res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))

par(mfrow=c(3, 2), mar=c(3, 3, 2, 1))

ylimedge <- 50
len <- hist(res$counts$all, main='A) All node pairs', xlab='Count', ylim=c(0, ylimedge))
stabselthreshold <- cached(c('FDR_threshold', 'stabsel', 'cluster_300_1'))
abline(v=2*B*stabselthreshold, lty=2, lwd=2)
legend('bottomleft', c('Histogram', 'CPSS threshold'), pch=c(20, NA),
  lty=c(NA, 2), lwd=2, bg='white')

len <- hist(res$counts$community, main='B) Community co-assignment',
  xlab='Count', ylim=c(0, 500))
lines(0:B, len*cached(c('hist_fit_nologit', outputfile, 'figedge'),
  res$counts$community, 0, B), col='blue', lty=6, lwd=2)
lines(0:B, len*cached(c('hist_fit', outputfile, 'figedge'),
  res$counts$community, 0, B), col='red', lwd=2)
minpoint <- cached(c('smooth_inner_minima', outputfile, 'figedge'), res$counts$community, 0, B)
points(minpoint, min(len*cached(c('hist_fit', outputfile, 'figedge'))),
  col='red', cex=2, lwd=3, pch=4)
legend('topright', c('Histogram', 'Gauss. smooth.',
  'logit Gauss. sm.', 'Est. minimum'), pch=c(20, NA, NA, 4), inset=c(0.1, 0),
  lty=c(NA, 6, 1, NA), col=c('black', 'blue', 'red', 'red'), lwd=c(1, 2, 2, 3),
  pt.cex=c(1, 1, 1, 2), bg='white')

len <- hist(res$counts$low, main='C) Between-community node pairs',
  xlab='Count', ylim=c(0, ylimedge))
lines(0:(2*B), len*cached(c('hist_fit_nologit', outputfile, 'figlow'),
  res$counts$low, 0, 2*B), col='blue', lty=6, lwd=2)
lines(0:(2*B), len*cached(c('hist_fit', outputfile, 'figlow'), res$counts$low,
  0, 2*B), col='red', lwd=2)
minpoint <- cached(c('smooth_inner_minima', outputfile, 'figlow'), res$counts$low, 0, 2*B)
abline(v=2*B*stabselthreshold, lty=2, lwd=2)
points(minpoint, min(len*cached(c('hist_fit', outputfile, 'figlow'))),
  col='red', cex=2, lwd=3, pch=4)
legend('bottomleft', c('Histogram', 'CPSS threshold', 'Gauss. smooth.',
  'logit Gauss. sm.', 'Est. minimum'), pch=c(20, NA, NA, NA, 4),
  lty=c(NA, 2, 6, 1, NA), col=c('black', 'black', 'blue', 'red', 'red'),
  lwd=c(1, 2, 2, 2, 3),
  pt.cex=c(1, 1, 1, 1, 2), bg='white')


len <- hist(res$counts$high, main='D) Within-community node pairs',
  xlab='Count', ylim=c(0, ylimedge))
lines(0:(2*B), len*cached(c('hist_fit_nologit', outputfile, 'fighigh'),
  res$counts$high, 0, 2*B), col='blue', lty=6, lwd=2)
lines(0:(2*B), len*cached(c('hist_fit', outputfile, 'fighigh'), res$counts$high,
  0, 2*B), col='red', lwd=2)
minpoint <- cached(c('smooth_inner_minima', outputfile, 'fighigh'), res$counts$high, 0, 2*B)
abline(v=2*B*stabselthreshold, lty=2, lwd=2)
points(minpoint, min(len*cached(c('hist_fit', outputfile, 'fighigh'))),
  col='red', cex=2, lwd=3, pch=4)
legend('top', c('Histogram', 'CPSS threshold', 'Gauss. smooth.',
  'logit Gauss. sm.', 'Est. minimum'), pch=c(20, NA, NA, NA, 4),
  lty=c(NA, 2, 6, 1, NA), col=c('black', 'black', 'blue', 'red', 'red'),
  lwd=c(1, 2, 2, 2, 3),
  pt.cex=c(1, 1, 1, 1, 2), bg='white')

len <- hist(res$counts$lowrescaled, main='E) Between-community node pairs',
  xlab='Pseudo count', ylim=c(0, ylimedge))
sassthreshold <- min(res$freq[res$estnet])
abline(v=2*B*sassthreshold, lwd=2)
points(minpoint, min(len*cached(c('hist_fit', outputfile, 'figlow'))),
  col='red', cex=2, lwd=3, pch=4)
legend('bottomleft', c('Histogram', 'SASS thresh.', 'Est. minimum'),
  pch=c(20, NA, 4),, pt.cex=c(1, 1, 2), col=c('black', 'black', 'red'),
  lty=c(NA, 1, NA), lwd=c(2, 2, 3), bg='white')

len <- hist(res$counts$highrescaled, main='F) Within-community node pairs',
  xlab='Pseudo count', ylim=c(0, ylimedge))
abline(v=2*B*sassthreshold, lwd=2)
points(minpoint, min(len*cached(c('hist_fit', outputfile, 'fighigh'))),
  col='red', cex=2, lwd=3, pch=4)
legend('topleft', c('Histogram', 'SASS threshold', 'Est. minimum'), inset=c(0.2, 0),
  pch=c(20, NA, 4),, pt.cex=c(1, 1, 2), col=c('black', 'black', 'red'),
  lty=c(NA, 1, NA), lwd=c(2, 2, 3), bg='white')


dev.off()
