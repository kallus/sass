source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/scores.R')
source('src/experiment/graphics.R')

graphs <- c('cluster', 'overlapping', 'hub', 'scale-free', 'erdos-renyi')
graphs <- c('cluster', 'overlapping', 'scale-free')
plotroc <- function(graph, repetition, show_legend=F) {
  col <- 0
  # count methods
  methods <- c('stabsel', 'sass', 'stabselglasso', 'sass-glasso')
  methods <- c('sass')
  for (method in methods) {
    col <- col + 1
    outputfile <- paste(method, '_', graph, '_300_', repetition, sep='')
    res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
    roc <- cached(c('compute_roc', outputfile), res$truenet$theta, res$freq)
    roc1 <- cached(c('compute_roc1', outputfile), res$truenet$theta, res$estnet)
    lines(roc[1, ], roc[2, ], col=col, lty=col)
    points(roc1[1], roc1[2], col=col)
  }

  # grab
  col <- col + 1
  outputfile <- paste('grab_', graph, '_300_', repetition, sep='')
  res <- readRDS(paste('cache/contender_output', outputfile, sep='/'))
  roc1 <- cached(c('compute_roc1', outputfile), res$truenet$theta, res$estnet)
  points(roc1[1], roc1[2], col=col)
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
  lines(roc[1, ], roc[2, ], col=col, lty=col)
  methods <- c(methods, 'grab')

  ## neighsel
  #col <- col + 1
  #outputfile <- paste('neighsel_', graph, '_300_', repetition, sep='')
  #if (file.exists(paste('contender_output', outputfile, sep='/'))) {
  #  res <- readRDS(paste('contender_output', outputfile, sep='/'))
  #  roc1 <- cached(c('compute_roc1', outputfile), res$truenet$theta, res$estnet)
  #  points(roc1[1], roc1[2], col=col)
  #  datafile <- paste(graph, '_300_', repetition, sep='')
  #  roc <- cached(c('compute_roc_solutions', datafile, 'neighsel'), res$true$theta, res$solutions)
  #  lines(roc[1, ], roc[2, ], col=col)
  #  methods <- c(methods, 'neighsel')
  #}

  ## glasso
  #col <- col + 1
  #outputfile <- paste('glasso_', graph, '_300_', repetition, sep='')
  #if (file.exists(paste('contender_output', outputfile, sep='/'))) {
  #  res <- readRDS(paste('contender_output', outputfile, sep='/'))
  #  roc1 <- cached(c('compute_roc1', outputfile), res$truenet$theta, res$estnet)
  #  points(roc1[1], roc1[2], col=col)
  #  datafile <- paste(graph, '_300_', repetition, sep='')
  #  roc <- cached(c('compute_roc_solutions', datafile, 'glasso'), res$true$theta, res$solutions)
  #  lines(roc[1, ], roc[2, ], col=col)
  #  methods <- c(methods, 'glasso')
  #}

  if (show_legend) {
    legend('bottomright', legend=methods, col=1:length(methods),
      lty=1:length(methods))
  }
}

pdf('fig/roc.pdf', width=5, height=8)
par(mfrow=c(length(graphs), 2))

for (graph in graphs) {
  reps <- 1

  for (repetition in 1:reps) {
    # zoom out
    plot(-1, -1, xlab='FPR', ylab='TPR', xlim=c(0, 0.1), ylim=c(0, 1), main=graph)
    plotroc(graph, repetition, show_legend=T)

    # zoom in
    plot(-1, -1, xlab='FPR', ylab='TPR', xlim=c(0, 0.01), ylim=c(0, 1), main=graph)
    plotroc(graph, repetition, show_legend=T)
  }
}

dev.off()
