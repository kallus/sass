source('src/experiment/caching.R')
source('src/experiment/scores.R')

library(Matrix)

#list_connections <- function(adj, names) {
#  cons <- which(adj & upper.tri(adj), arr.ind=T)
#  data.frame(V1=names[cons[, 1]], V2=names[cons[, 2]], stringsAsFactors=F)
#}
#
#list_connection_diffs <- function(adj1, adj2, names) {
#  added <- which(!adj1 & adj2 & upper.tri(adj1), arr.ind=T)
#  removed <- which(!adj2 & adj1 & upper.tri(adj1), arr.ind=T)
#  list(not_in_adj1=data.frame(V1=names[added[, 1]], V2=names[added[, 2]],
#      stringsAsFactors=F),
#    not_in_adj2=data.frame(V1=names[removed[, 1]], V2=names[removed[, 2]],
#      stringsAsFactors=F))
#}
#
#prepare_db <- function(tbl) {
#  env <- new.env()
#  for (i in 1:nrow(tbl)) {
#    cat(i / nrow(tbl)); cat('\r')
#    n1 <- tbl$V1[i]
#    n2 <- tbl$V3[i]
#    env[[paste(n1, n2, sep='.')]] <- T
#    env[[paste(n2, n1, sep='.')]] <- T
#  }
#  return(env)
#}
#
#in_db <- function(cons, env) {
#  ret <- rep(F, nrow(cons))
#  if (nrow(cons) > 0) {
#    for (i in 1:nrow(cons)) {
#      cat(i / nrow(cons)); cat('\r')
#      hits <- exists(paste(cons[i, 1], cons[i, 2], sep='.'), envir=env)
#      if (hits) {
#        ret[i] <- T
#      }
#    }
#  }
#  return(ret)
#}

sif_to_adjmat <- function(sif, names) {
  d <- length(names)
  adj <- Matrix::Matrix(F, d, d, dimnames=list(names, names))
  for (i in 1:nrow(sif)) {
    cat(i / nrow(sif)); cat('\r')
    try(adj[sif[i, 1], sif[i, 3]] <- T, silent=T)
    try(adj[sif[i, 3], sif[i, 1]] <- T, silent=T)
  }
  return(adj)
}

sass <- readRDS('cache/contender_output/sass_real')
#stabsel <- readRDS('contender_output/stabsel_real')
#
genenames <- colnames(sass$truenet$data)
#
#sass$connections <- list_connections(sass$estnet, genenames)
#stabsel$connections <- list_connections(stabsel$estnet, genenames)


reactome <- read.table('data/PathwayCommons12.reactome.hgnc.sif.gz',
  stringsAsFactors=F)
#
#db <- prepare_db(reactome)
#sass$connections$in_db <- in_db(sass$connections, db)
#stabsel$connections$in_db <- in_db(stabsel$connections, db)
#
#difference <- list_connection_diffs(sass$estnet, stabsel$estnet, genenames)
#difference$not_in_adj1$in_db <- in_db(difference$not_in_adj1, db)
#difference$not_in_adj2$in_db <- in_db(difference$not_in_adj2, db)
#
#contingency <- matrix(NA, ncol=2, nrow=2,
#  dimnames=list(c('in db', 'not in db'), c('stabsel', 'removed')))
#contingency['in db', 'stabsel'] <- sum(stabsel$connections$in_db)
#contingency['not in db', 'stabsel'] <- sum(!stabsel$connections$in_db)
#contingency['in db', 'removed'] <- sum(difference$not_in_adj1$in_db)
#contingency['not in db', 'removed'] <- sum(!difference$not_in_adj1$in_db)

datafile <- 'real'
stabsel_freq <- cached(c('selection_freq', datafile))
sass_freq <- cached(c('adapt_frequencies_mask2', 'sass', datafile))

reactome_adjmat <- cached(c('sif_to_adjmat', 'reactome'), reactome, genenames)

stabsel_roc <- cached(c('compute_roc', 'realstabselpathwayfig'), reactome_adjmat, stabsel_freq)

sass_roc <- cached(c('compute_roc', 'realsasspathwayfig'), reactome_adjmat, sass_freq)

stabsel_roc <- stabsel_roc[, 1:(ncol(stabsel_roc)-1)]
sass_roc <- sass_roc[, 1:(ncol(sass_roc)-1)]

xout <- seq(max(stabsel_roc[1, 1], sass_roc[1, 1]),
  min(stabsel_roc[1, ncol(stabsel_roc)], sass_roc[1, ncol(sass_roc)]),
  length=1000)
stabsel_approx <- approx(stabsel_roc[1, ], stabsel_roc[2, ], xout)
sass_approx <- approx(sass_roc[1, ], sass_roc[2, ], xout)

pdf('fig/pathway.pdf', width=8, height=5)
par(mfrow=c(1, 2))
plot(stabsel_roc[1, ], stabsel_roc[2, ], xlab='FPR', ylab='TPR', t='l',
  main='ROC')
lines(sass_roc[1, ], sass_roc[2, ], col=2)
legend('topleft', c('stabsel', 'sass'), col=1:2, lty=1)
plot(sass_approx$x, cumsum(sass_approx$y-stabsel_approx$y), t='l',
  main='AUROC difference (sass - stabsel)', xlab='FPR',
  ylab='Integral of ROC difference')
dev.off()
