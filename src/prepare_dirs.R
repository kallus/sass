dirs <- c('cache', 'fig', 'cache/cache', 'cache/contender_output', 'cache/data',
  'cache/grab', 'cache/grab/data', 'cache/grab/output')
for (dir in dirs) {
  if (!file.exists(dir)) {
    dir.create(dir)
  }
}
