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
comm_threshold <- cached(c('smooth_inner_minima', 'sass', 'real'))
cmask <- comm_counts > comm_threshold

g <- igraph::graph_from_adjacency_matrix(cmask, mode='undirected')
louvain <- igraph::cluster_louvain(g)
ix <- order(louvain$membership)

png('fig/clusters.png', antialias='none', res=300, width=4*1024, height=4*1024)
image(Matrix::drop0(cmask[ix, ix]))
dev.off()
