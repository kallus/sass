source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/scores.R')
source('src/experiment/graphics.R')

immcomp <- function(x, xhat, asp=TRUE, legend.outside=NA, xlab='Column index',
    ylab='Row index', ...) {
  if (asp) {
    asp <- 1
  } else {
    asp <- NA
  }
  d <- dim(x)+1
  colors <- c('white', 'black')
  graphics::image(1:d[2]-0.5, 1:d[1]-0.5, t(x & xhat), xlim=c(1, d[2])-0.5,
    ylim=c(d[1], 1)-0.5, xlab=xlab, ylab=ylab, asp=asp, axes=FALSE,
    col=colors, ...)
  colors <- c(rgb(0, 0, 0, 0), 'red')
  graphics::image(1:d[2]-0.5, 1:d[1]-0.5, t(xhat & !x), xlim=c(1, d[2])-0.5,
    ylim=c(d[1], 1)-0.5, xlab=xlab, ylab=ylab, asp=asp, axes=FALSE, add=T,
    col=colors, ...)
  colors <- c(rgb(0, 0, 0, 0), 'deepskyblue')
  graphics::image(1:d[2]-0.5, 1:d[1]-0.5, t(x & !xhat), xlim=c(1, d[2])-0.5,
    ylim=c(d[1], 1)-0.5, xlab=xlab, ylab=ylab, asp=asp, axes=FALSE, add=T,
    col=colors, ...)
  atx <- pretty(1:d[1])
  atx <- atx[atx < d[1]]
  atx <- c(atx, d[1]-1)
  aty <- pretty(1:d[2])
  aty <- aty[aty < d[2]]
  aty <- c(aty, d[2]-1)
  lines(c(0.5, d[2]-0.5), c(0.5, 0.5))
  lines(c(0.5, d[2]-0.5), c(d[1]-0.5, d[1]-0.5))
  lines(c(0.5, 0.5), c(0.5, d[1]-0.5))
  lines(c(d[2]-0.5, d[2]-0.5), c(0.5, d[1]-0.5))
  #abline(h=0.5)
  #abline(h=d[1]-0.5)
  #abline(v=0.5)
  #abline(v=d[2]-0.5)
}

graph <- 'cluster'
outputfile <- paste('stabsel', graph, 300, 1, sep='_')
res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
truenet <- res$truenet$theta
stabsel <- res$estnet
outputfile <- paste('sass', graph, 300, 1, sep='_')
res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
sass <- res$estnet

ix <- which(sapply(1:ncol(stabsel),
  function(i) !all(stabsel[i, ] == sass[i, ])))

pdf('fig/adjdiff.pdf', width=6.5, height=3.5)
par(mfrow=c(1, 2), mar=c(1, 1, 1, 1))
immcomp(as.matrix(truenet)[ix, ix], as.matrix(stabsel)[ix, ix], main='A')
legend('bottomleft', c('True positive', 'False positive', 'False negative'),
  pch='.', col=c('black', 'red', 'deepskyblue'), pt.cex=1, cex=0.65, bg='white')
immcomp(as.matrix(truenet)[ix, ix], as.matrix(sass)[ix, ix], main='B')
legend('bottomleft', c('True positive', 'False positive', 'False negative'),
  pch='.', col=c('black', 'red', 'deepskyblue'), pt.cex=1, cex=0.65, bg='white')
dev.off()
