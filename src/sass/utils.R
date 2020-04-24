uptri <- function(x) {
  x <- as.matrix(x)
  x[upper.tri(x)]
}
