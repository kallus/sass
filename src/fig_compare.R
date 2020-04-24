source('sass/sass.R')
source('sass/smoothing.R')
source('experiment/generator.R')
source('experiment/caching.R')
source('experiment/contenders.R')
source('experiment/scores.R')
source('experiment/graphics.R')

methods <- c('stabsel', 'sass', 'neighsel', 'stabselglasso', 'sass-glasso', 'glasso', 'grab')
graphs <- c('cluster', 'overlapping', 'hub', 'scale-free', 'erdos-renyi')

pdf('compare.pdf')

# roc curves
for (graph in graphs) {
  par(mfrow=c(1, 1))
  plot(-1, -1, xlab='FPR', ylab='TPR', xlim=c(0, 0.0075), ylim=c(0, 1), main=graph)
  col <- 0
  for (method in methods) {
    pattern <- paste('^', method, '_', graph, '_', sep='')
    col <- col + 1
    for (outputfile in list.files(path='contender_output/', pattern=pattern)) {
      res <- readRDS(paste('contender_output', outputfile, sep='/'))
      if (method %in% c('neighsel', 'glasso')) {
        roc <- res$roc
      } else {
        roc <- cached(c('compute_roc', outputfile), res$truenet$theta, res$freq)
      }
      roc1 <- cached(c('compute_roc1', outputfile), res$truenet$theta, res$estnet)
      lines(roc[1, ], roc[2, ], col=col)
      points(roc1[1], roc1[2], col=col)
      break
    }
  }
  legend('bottomright', legend=methods, col=1:col, lty=1)
}

# fdr curves
for (graph in intersect(graphs, c('cluster', 'overlapping'))) {
  par(mfrow=c(1, 1))
  plot(-1, -1, xlab='Network density', ylab='FDR and bias', xlim=c(0.005, 0.01), ylim=c(0, 1), main=graph)
  col <- 0
  for (method in methods) {
    if (method %in% c('neighsel', 'glasso')) {
      next
    }
    pattern <- paste('^', method, '_', graph, '_', sep='')
    col <- col + 1
    for (outputfile in list.files(path='contender_output/', pattern=pattern)) {
      res <- readRDS(paste('contender_output', outputfile, sep='/'))
      fdr <- cached(c('compute_fdr', outputfile), res$truenet$theta, res$freq, res$truenet$cmask)
      bias <- cached(c('compute_biases', outputfile), res$truenet$theta, res$freq, res$truenet$cmask)
      lines(fdr[1, ], fdr[2, ], col=col)
      lines(fdr[1, ], bias, col=col, lty=2)
      break
    }
  }
  legend('topleft', legend=methods, col=1:col, lty=1)
}

# boxplots mcc
for (graph in graphs) {
  mccs <- list()
  for (method in methods) {
    mccs[[method]] <- c()
    pattern <- paste('^', method, '_', graph, '_', sep='')
    for (outputfile in list.files(path='contender_output/', pattern=pattern)) {
      res <- readRDS(paste('contender_output', outputfile, sep='/'))
      mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
      mccs[[method]] <- c(mccs[[method]], mcc)
    }
  }
  par(mfrow=c(1, 1))
  stripchart(mccs, main=graph, ylab='MCC', method='jitter',
    vertical=T, pch=20, las=2)
}

# boxplots fdr
for (graph in graphs) {
  fdrs <- list()
  for (method in methods) {
    fdrs[[method]] <- c()
    pattern <- paste('^', method, '_', graph, '_', sep='')
    for (outputfile in list.files(path='contender_output/', pattern=pattern)) {
      res <- readRDS(paste('contender_output', outputfile, sep='/'))
      fdr <- cached(c('compute_fdr1_nomask', outputfile), res$truenet$theta, res$estnet)
      fdrs[[method]] <- c(fdrs[[method]], fdr)
    }
  }
  par(mfrow=c(1, 1))
  stripchart(fdrs, main=graph, ylab='FDR', method='jitter',
    vertical=T, pch=20, las=2)
}

# boxplots bias
for (graph in intersect(graphs, c('cluster', 'overlapping'))) {
  biases <- list()
  for (method in methods) {
    biases[[method]] <- c()
    pattern <- paste('^', method, '_', graph, '_', sep='')
    for (outputfile in list.files(path='contender_output/', pattern=pattern)) {
      res <- readRDS(paste('contender_output', outputfile, sep='/'))
      bias <- cached(c('compute_bias', outputfile), res$truenet$theta, res$estnet, res$truenet$cmask)
      biases[[method]] <- c(biases[[method]], bias)
    }
  }
  par(mfrow=c(1, 1))
  stripchart(biases, main=graph, ylab='Bias', method='jitter',
    vertical=T, pch=20, las=2)
  abline(h=0, lty=2)
}

# scatter mcc vs bias
for (graph in intersect(graphs, c('cluster', 'overlapping'))) {
  mccs <- list()
  biases <- list()
  par(mfrow=c(1, 1))
  plot(-1, -1, xlim=c(0, 1), ylim=c(0, 2), main=graph, xlab='MCC', ylab='Bias')
  i <- 0
  for (method in methods) {
    mccs[[method]] <- c()
    biases[[method]] <- c()
    pattern <- paste('^', method, '_', graph, '_', sep='')
    for (outputfile in list.files(path='contender_output/', pattern=pattern)) {
      res <- readRDS(paste('contender_output', outputfile, sep='/'))
      mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
      mccs[[method]] <- c(mccs[[method]], mcc)
      bias <- cached(c('compute_bias', outputfile), res$truenet$theta, res$estnet, res$truenet$cmask)
      biases[[method]] <- c(biases[[method]], bias)
    }
    i <- i + 1
    points(mccs[[method]], biases[[method]], pch=i, col=i)
  }
  legend('topleft', legend=methods, pch=1:i, col=1:i)
}

# scatter mcc vs fdr difference, zoom 1
for (graph in intersect(graphs, c('cluster', 'overlapping'))) {
  mccs <- list()
  biases <- list()
  zoom_methods <- c('sass', 'stabsel', 'neighsel')
  xlim <- c(1, 0)
  ylim <- c(1, 0)
  for (method in zoom_methods) {
    mccs[[method]] <- c()
    biases[[method]] <- c()
    pattern <- paste('^', method, '_', graph, '_', sep='')
    for (outputfile in list.files(path='contender_output/', pattern=pattern)) {
      res <- readRDS(paste('contender_output', outputfile, sep='/'))
      mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
      if (mcc < xlim[1]) xlim[1] <- mcc
      if (mcc > xlim[2]) xlim[2] <- mcc
      mccs[[method]] <- c(mccs[[method]], mcc)
      bias <- cached(c('compute_bias', outputfile), res$truenet$theta, res$estnet, res$truenet$cmask)
      biases[[method]] <- c(biases[[method]], bias)
      if (bias < ylim[1]) ylim[1] <- bias
      if (bias > ylim[2]) ylim[2] <- bias
    }
  }
  par(mfrow=c(1, 1))
  xlim[1] <- xlim[1] - 0.01*(xlim[2]-xlim[1])
  xlim[2] <- xlim[2] + 0.01*(xlim[2]-xlim[1])
  ylim[1] <- ylim[1] - 0.01*(ylim[2]-ylim[1])
  ylim[2] <- ylim[2] + 0.01*(ylim[2]-ylim[1])
  plot(-i, -1, xlim=xlim, ylim=ylim, main=graph, xlab='MCC', ylab='Bias')
  i <- 0
  for (method in zoom_methods) {
    i <- i + 1
    points(mccs[[method]], biases[[method]], pch=i, col=i)
  }
  legend('topleft', legend=zoom_methods, pch=1:i, col=1:i)
  for (i in 1:1000) {
    mccs <- list()
    biases <- list()
    for (method in zoom_methods) {
      outputfile <- paste(method, '_', graph, '_300_', i, sep='')
      if (file.exists(paste('contender_output', outputfile, sep='/'))) {
        res <- readRDS(paste('contender_output', outputfile, sep='/'))
        mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
        mccs[[method]] <- c(mccs[[method]], mcc)
        bias <- cached(c('compute_bias', outputfile), res$truenet$theta, res$estnet, res$truenet$cmask)
        biases[[method]] <- c(biases[[method]], bias)
      }
    }
    add_lines <- function(from, to) {
      lines(c(mccs[[from]], mccs[[to]]),
        c(biases[[from]], biases[[to]]), col=rgb(0, 0, 0, 0.25))
    }
    add_lines('stabsel', 'sass')
    add_lines('neighsel', 'sass')
  }
}

# scatter mcc vs fdr difference, zoom 2
for (graph in intersect(graphs, c('cluster', 'overlapping'))) {
  mccs <- list()
  fdrd <- list()
  zoom_methods <- c('sass-glasso', 'stabselglasso', 'glasso')
  xlim <- c(1, 0)
  ylim <- c(1, 0)
  for (method in zoom_methods) {
    mccs[[method]] <- c()
    fdrd[[method]] <- c()
    pattern <- paste('^', method, '_', graph, '_', sep='')
    for (outputfile in list.files(path='contender_output/', pattern=pattern)) {
      res <- readRDS(paste('contender_output', outputfile, sep='/'))
      mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
      if (mcc < xlim[1]) xlim[1] <- mcc
      if (mcc > xlim[2]) xlim[2] <- mcc
      mccs[[method]] <- c(mccs[[method]], mcc)
      fdr1 <- cached(c('compute_fdr1', outputfile), res$truenet$theta, res$estnet, res$truenet$cmask)
      fdrd[[method]] <- c(fdrd[[method]], abs(fdr1[4]-fdr1[3]))
      fdrdi <- abs(fdr1[4]-fdr1[3])
      if (fdrdi < ylim[1]) ylim[1] <- fdrdi
      if (fdrdi > ylim[2]) ylim[2] <- fdrdi
    }
  }
  par(mfrow=c(1, 1))
  xlim[1] <- xlim[1] - 0.01*(xlim[2]-xlim[1])
  xlim[2] <- xlim[2] + 0.01*(xlim[2]-xlim[1])
  ylim[1] <- ylim[1] - 0.01*(ylim[2]-ylim[1])
  ylim[2] <- ylim[2] + 0.01*(ylim[2]-ylim[1])
  plot(-i, -1, xlim=xlim, ylim=ylim, main=graph, xlab='MCC', ylab='Bias')
  i <- 0
  for (method in zoom_methods) {
    i <- i + 1
    points(mccs[[method]], fdrd[[method]], pch=i, col=i)
  }
  legend('topleft', legend=zoom_methods, pch=1:i, col=1:i)
  for (i in 1:1000) {
    mccs <- list()
    fdrd <- list()
    for (method in zoom_methods) {
      outputfile <- paste(method, '_', graph, '_300_', i, sep='')
      if (file.exists(paste('contender_output', outputfile, sep='/'))) {
        res <- readRDS(paste('contender_output', outputfile, sep='/'))
        mcc <- cached(c('compute_mcc', outputfile), res$truenet$theta, res$estnet)
        mccs[[method]] <- c(mccs[[method]], mcc)
        fdr1 <- cached(c('compute_fdr1', outputfile), res$truenet$theta, res$estnet, res$truenet$cmask)
        fdrd[[method]] <- c(fdrd[[method]], abs(fdr1[4]-fdr1[3]))
      }
    }
    add_lines <- function(from, to) {
      lines(c(mccs[[from]], mccs[[to]]),
        c(fdrd[[from]], fdrd[[to]]), col=rgb(0, 0, 0, 0.25))
    }
    add_lines('glasso', 'sass-glasso')
    add_lines('stabselglasso', 'sass-glasso')
  }
}

# boxplots community detection
for (graph in intersect(graphs, c('cluster'))) {
  arand <- list()
  for (method in methods) {
    arand[[method]] <- c()
    pattern <- paste('^', method, '_', graph, '_', sep='')
    for (outputfile in list.files(path='contender_output/', pattern=pattern)) {
      res <- readRDS(paste('contender_output', outputfile, sep='/'))
      adjrandindex <- cached(c('compute_adjrand', outputfile), res$estnet)
      arand[[method]] <- c(arand[[method]], adjrandindex)
    }
  }
  par(mfrow=c(1, 1))
  stripchart(arand, main=graph, ylab='Adj. Rand index', method='jitter',
    vertical=T, pch=20, las=2)
  abline(h=0, lty=2)
}

# adjacency matrices
for (graph in graphs) {
  for (method in methods) {
    pattern <- paste('^', method, '_', graph, '_', sep='')
    for (outputfile in list.files(path='contender_output/', pattern=pattern)) {
      res <- readRDS(paste('contender_output', outputfile, sep='/'))
      immcomp(as.matrix(res$truenet$theta), as.matrix(res$estnet),
        main=paste(graph, method))
      break
    }
  }
}

# histograms
B <- 100
for (graph in union(graphs, 'real')) {
  pattern <- paste('^sass_', graph, sep='')
  for (outputfile in list.files(path='contender_output/', pattern=pattern)) {
    res <- readRDS(paste('contender_output', outputfile, sep='/'))
    par(mfrow=c(3, 2))
    bins <- 0:201-0.5
    hm <- max(hist(res$counts$all, bins, plot=F)$density)
    hist(res$counts$all, bins, freq=F, main=graph, xlab='All edge counts', ylim=c(0, min(0.01, hm)))
    bins <- 0:101-0.5
    hm <- max(hist(res$counts$community, bins, plot=F)$density)
    hist(res$counts$community, bins, freq=F, main=graph, xlab='Community counts', ylim=c(0, min(0.01, hm)))
    lines(0:B, cached(c('hist_fit_nologit', outputfile, 'figedge'), res$counts$community, 0, B), col='blue')
    lines(0:B, cached(c('hist_fit', outputfile, 'figedge'), res$counts$community, 0, B), col='red')
    minpoint <- cached(c('smooth_inner_minima', outputfile, 'figedge'), res$counts$community, 0, B)
    abline(v=minpoint, col='red')
    bins <- 0:201-0.5
    hm <- max(hist(res$counts$low, bins, plot=F)$density)
    hist(res$counts$low, bins, freq=F, main=graph, xlab='Outside community edge counts', ylim=c(0, min(0.01, hm)))
    lines(0:(2*B), cached(c('hist_fit_nologit', outputfile, 'figlow'), res$counts$low, 0, 2*B), col='blue')
    lines(0:(2*B), cached(c('hist_fit', outputfile, 'figlow'), res$counts$low, 0, 2*B), col='red')
    minpoint <- cached(c('smooth_inner_minima', outputfile, 'figlow'), res$counts$low, 0, 2*B)
    abline(v=minpoint, col='red')
    hm <- max(hist(res$counts$high, bins, plot=F)$density)
    if (is.na(hm)) {
      hm <- 0.01
    }
    hist(res$counts$high, bins, freq=F, main=graph, xlab='Inside community edge counts', ylim=c(0, min(0.01, hm)))
    lines(0:(2*B), cached(c('hist_fit_nologit', outputfile, 'fighigh'), res$counts$high, 0, 2*B), col='blue')
    lines(0:(2*B), cached(c('hist_fit', outputfile, 'fighigh'), res$counts$high, 0, 2*B), col='red')
    minpoint <- cached(c('smooth_inner_minima', outputfile, 'fighigh'), res$counts$high, 0, 2*B)
    abline(v=minpoint, col='red')
    hm <- max(hist(res$counts$lowrescaled, bins, plot=F)$density)
    hist(res$counts$lowrescaled, bins, freq=F, main=graph, xlab='Outside community edge counts rescaled',
      ylim=c(0, min(0.01, hm)))
    abline(v=minpoint, col='red')
    hm <- max(hist(res$counts$highrescaled, bins, plot=F)$density)
    if (is.na(hm)) {
      hm <- 0.01
    }
    hist(res$counts$highrescaled, bins, freq=F, main=graph, xlab='Inside community edge counts rescaled',
      ylim=c(0, min(0.01, hm)))
    abline(v=minpoint, col='red')
  }
}

dev.off()
