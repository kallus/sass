source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/scores.R')
source('src/experiment/graphics.R')

filename <- 'cache/contender_output/sass_real'
res <- readRDS(filename)

comm_counts <- cached(c('consensus_module_counts', 'sass', 'real'))

resolutions <- function(d) {
  if (d <= 200) {
    return(ceiling(d/20))
  }
  if (d <= 1000) {
    return(ceiling(d/c(20, 200)))
  }
  resrec <- function(d) {
    if (d <= 25) {
      return(ceiling(d/5))
    }
    return(c(ceiling(d/5), resrec(ceiling(d/5))))
  }
  return(c(ceiling(d/c(20, 200)), resrec(ceiling(d/200))))
}

# visualize matrix in same order as it is printed
immbw <- function(x, asp=TRUE, legend.outside=NA, xlab='Column index',
    ylab='Row index', ...) {
  if (asp) {
    asp <- 1
  } else {
    asp <- NA
  }
  d <- dim(x)+1
  colors <- c('white', 'dodgerblue4')
  graphics::image(1:d[2]-0.5, 1:d[1]-0.5, t(x), xlim=c(1, d[2])-0.5,
    ylim=c(d[1], 1)-0.5, xlab=xlab, ylab=ylab, asp=asp, axes=FALSE,
    col=colors, ...)
}

plot.mrsics <- function(res) {
  immbw(as.matrix(res$network[res$order, res$order][1:855+1582, 1:855+1582]), ann=F)
  d <- ncol(res$network)
  for (i in 1:2){#length(res$marks)) {
    marks <- res$marks[[i]]
    col <- rgb(0.1, 0.1, 0.1, sqrt(1/length(marks)))
    col <- 1
    for (j in 1:length(marks)) {
      lines(c(0, d), rep(marks[j]-1582, 2), lwd=i^3/5, lty=1, col=3-i)
      lines(rep(marks[j]-1582, 2), c(0, d), lwd=i^3/5, lty=1, col=3-i)
    }
  }
}

B <- 100
d <- ncol(comm_counts)
k <- resolutions(d)
hcl <- hclust(as.dist(B-comm_counts), method='ward.D2')
cl <- cutree(hcl, rev(k))
if (length(k) == 1) {
  cl <- matrix(cl, ncol=1)
}
for (i in 1:ncol(cl)) {
  map <- order(order(-table(cl[, i])))
  cl[, i] <- map[cl[, i]]
}
ix <- do.call(order, lapply(1:length(k), function(i) cl[, i]))
clmarks <- lapply(length(k):1, function(i) {
  which(diff(cl[ix, i]) != 0) + 0.5
})

for (i in 1:ncol(cl)) {
  print(table(cl[, i]))
}

png('fig/hierarchic.png', antialias='none', res=1200, width=6.5, height=6.5,
  units='in')
par(mar=c(0, 0, 0, 0))
plot.mrsics(list(network=res$estnet, order=ix, marks=clmarks))
legend('bottomleft', c('Edge', 'Community', 'Super community'),
  pch=c('.', NA, NA), lty=c(NA, 1, 1), lwd=c(1, 1/5, 8/5), col=c('dodgerblue4',
  'red', 'black'), pt.cex=1, cex=0.65)
dev.off()

# hist <- function(data) {
#   td <- table(data)
#   barplot(td)
#   abline(h=mean(td), lty=2)
# }
# par(mfrow=c(1, 5))
# for (i in 1:4) {
#   hist(cl[, i])
# }
# densities <- rep(NA, 5)
# mask <- outer(cl[, 1], cl[, 1], '==')
# densities[1] <- mean(res$estnet[!mask])
# densities[2] <- mean(res$estnet[mask])
# mask <- outer(cl[, 2], cl[, 2], '==')
# densities[3] <- mean(res$estnet[mask])
# mask <- outer(cl[, 3], cl[, 3], '==')
# densities[4] <- mean(res$estnet[mask])
# mask <- outer(cl[, 4], cl[, 4], '==')
# densities[5] <- mean(res$estnet[mask])
# barplot(densities)
