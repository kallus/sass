real_data <- function(subsample=F) {
  # utility methods
  row_apply <- function(x, f) {
    res <- apply(x, 1, f)
    if (is.matrix(res)) {
      return(t(res))
    }
    return(res)
  }
  col_apply <- function(x, f) apply(x, 2, f)

  top_k_mad_columns <- function(x, k) {
    x[, rank(-col_apply(x, function(x) mad(x, na.rm=TRUE)),
      ties.method='first') <= k]
  }

  # methods for normalization
  quantnormalize <- function(x.na) {
    x <- x.na[!is.na(x.na)]
    x <- rank(x, ties.method='average')
    x <- qnorm(x/(length(x)+1))
    x <- (x-mean(x))/sd(x)
    x.na[!is.na(x.na)] <- x
    return(x.na)
  }

  center_scale_columns <- function(x) scale(x, center=TRUE, scale=TRUE)
  center_columns <- function(x) scale(x, center=TRUE, scale=FALSE)
  center_rows <- function(x) t(scale(t(x), center=TRUE, scale=FALSE))
  quant_normalize_rows <- function(x) row_apply(x, quantnormalize)

  # load data
  x <- readRDS('data/tcga_gbm_rnaseq.rds.gz')
  clinical <- readRDS('data/tcga_gbm_clinical.rds.gz')

  if (subsample) {
    ix <- sample(nrow(x), round(0.9*nrow(x)))
    x <- x[ix, ]
    clinical <- clinical[ix, ]
  }

  x <- quant_normalize_rows(x)

  # remove 0 MAD columns
  x <- x[, which(col_apply(x, function(x) mad(x, na.rm=TRUE)) > 0)]

  # remove the 75% of genes with lowest variance
  qv <- quantile(col_apply(x, function(x) var(x, na.rm=TRUE)), 0.75)
  x <- x[, which(col_apply(x, function(x) var(x, na.rm=TRUE)) > qv)]

  # normalize
  x <- center_scale_columns(x)

  # remove normal tissue
  ix <- clinical$sampletype != 'Solid Tissue Normal'
  x <- x[ix, ]
  clinical <- clinical[ix, ]

  return(list(data=x))
}
