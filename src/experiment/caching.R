cached <- function(identifier, ...) {
  filename <- paste('cache/cache/', paste(identifier, collapse='_'), sep='')
  if (file.exists(filename)) {
    result <- try(readRDS(filename), silent=T)
    if (!inherits(result, 'try-error')) {
      return(result)
    }
  }
  result <- do.call(identifier[1], list(...))
  saveRDS(result, filename)
  return(result)
}

save_exists <- function(dir, ...) {
  name <- paste(..., sep='_')
  filename <- paste(dir, name, sep='/')
  file.exists(filename)
}

save <- function(object, dir, ...) {
  name <- paste(..., sep='_')
  filename <- paste(dir, name, sep='/')
  saveRDS(object, filename)
}
