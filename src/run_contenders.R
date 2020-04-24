source('src/sass/sass.R')
source('src/sass/smoothing.R')
source('src/experiment/generator.R')
source('src/experiment/caching.R')
source('src/experiment/contenders.R')
source('src/experiment/real.R')
source('src/experiment/scores.R')

# generate data
#graphs <- c('hub', 'scale-free', 'erdos-renyi', 'cluster', 'overlapping')
graphs <- c('scale-free', 'cluster', 'overlapping')
seeds <- 1:100
n <- 300
for (seed in seeds) {
  for (graph in graphs) {
    if (!save_exists('cache/data', graph, n, seed)) {
      net <- generator(n, graph, seed)
      save(net, 'cache/data', graph, n, seed)
    }
  }
}

# prepare real data
if (!save_exists('cache/data', 'real')) {
  net <- real_data()
  save(net, 'cache/data', 'real')
}

# export data to Python for GRAB
for (datafile in list.files(path='cache/data/')) {
  net <- readRDS(paste('cache/data', datafile, sep='/'))
  filename <- paste('cache/grab/data/', datafile, sep='')
  if (!file.exists(filename)) {
    write.table(net$data, file=filename, row.names=F, col.names=F)
  }
}

B <- 100
lambda_ratio <- 0.1
lambda_ratio_glasso <- 1/3 # increased due to threshold=1

# run contending methods on data files
for (datafile in list.files(path='cache/data/')) {
  print(datafile)
  net <- readRDS(paste('cache/data', datafile, sep='/'))
  n <- nrow(net$data)
  lambda_max <- cached(c('find_lambda_max', datafile), net$data)
  set.seed(1)
  subsamples <- subsample_observations(n, B)
  networks <- cached(c('subsample_neighsel', datafile),
    net$data, lambda_ratio*lambda_max, B, subsamples)
  stabsel_freq <- cached(c('selection_freq', datafile), networks)
  threshold <- cached(c('FDR_threshold', 'stabsel', datafile),
    stabsel_freq, FDR=0.1, length(networks))
  stabsel_net <- cached(c('threshold_frequencies', 'stabsel', datafile),
    stabsel_freq, threshold)
#
  if (!save_exists('cache/contender_output', 'stabsel', datafile)) {
    save(list(freq=stabsel_freq, estnet=stabsel_net, truenet=net),
      'cache/contender_output', 'stabsel', datafile)
  }
#
  if (!save_exists('cache/contender_output', 'sassknowncl', datafile)) {
    if (!is.null(net$cmask)) {
      sasskncl_freq <- stabsel_freq
      sasskncl_net <- stabsel_net
      cmask <- net$cmask
      if (any(cmask)) {
        inside_min <- cached(c('smooth_inner_minima', 'sasskncl', 'inside', datafile),
          2*B*uptri(stabsel_freq)[uptri(cmask)], 0, 2*B)
        if (length(inside_min) == 1) {
          between_min <- cached(c('best_smooth_min', 'sasskncl', 'between', datafile),
            2*B*uptri(stabsel_freq)[uptri(!cmask)], 0, 2*B)
          if (!is.na(between_min) && inside_min < between_min) {
            sasskncl_freq <- cached(c('adapt_frequencies_mask2', 'sasskncl', datafile),
              stabsel_freq, inside_min, between_min, cmask, B)
            sasskncl_net <- cached(c('threshold_frequencies_match_density', 'sasskncl', datafile),
              sasskncl_freq, stabsel_net)
          }
        }
      }
      save(list(freq=sasskncl_freq, estnet=sasskncl_net, truenet=net,
        counts=list(all=2*B*uptri(stabsel_freq),
          highrescaled=2*B*uptri(sasskncl_freq)[uptri(cmask)],
          lowrescaled=2*B*uptri(sasskncl_freq)[uptri(!cmask)],
          high=2*B*uptri(stabsel_freq)[uptri(cmask)],
          low=2*B*uptri(stabsel_freq)[uptri(!cmask)])),
        'cache/contender_output', 'sasskncl', datafile)
    }
  }
#
  if (!save_exists('cache/contender_output', 'sassadaptive', datafile)) {
    sassadaptive_freq <- stabsel_freq
    sassadaptive_net <- stabsel_net
    comm_counts <- cached(c('consensus_module_counts', 'sass', datafile), networks, B)
    comm_threshold <- cached(c('smooth_inner_minima', 'sass', datafile),
      uptri(comm_counts), 0, B)
    if (length(comm_threshold) != 1) {
      comm_threshold <- Inf
    }
    cmask <- comm_counts > comm_threshold
    if (any(cmask)) {
      inside_min <- cached(c('smooth_inner_minima', 'sass', 'inside', datafile),
        2*B*uptri(stabsel_freq)[uptri(cmask)], 0, 2*B)
      if (length(inside_min) == 1) {
        between_min <- cached(c('best_smooth_min', 'sass', 'between', datafile),
          2*B*uptri(stabsel_freq)[uptri(!cmask)], 0, 2*B)
        if (!is.na(between_min) && inside_min < between_min) {
          sassadaptive_freq <- cached(c('adapt_frequencies_mask2', 'sass', datafile),
            stabsel_freq, inside_min, between_min, cmask, B)
        }
      }
    }
    sassadaptive_net <- cached(c('sass_adaptive', 'sassadaptive', datafile),
      net$data, lambda_max, sassadaptive_freq, stabsel_net)
    save(list(freq=sassadaptive_freq, estnet=sassadaptive_net, truenet=net,
      counts=list(all=2*B*uptri(stabsel_freq), community=uptri(comm_counts),
        highrescaled=2*B*uptri(sassadaptive_freq)[uptri(cmask)],
        lowrescaled=2*B*uptri(sassadaptive_freq)[uptri(!cmask)],
        high=2*B*uptri(stabsel_freq)[uptri(cmask)],
        low=2*B*uptri(stabsel_freq)[uptri(!cmask)])),
      'cache/contender_output', 'sassadaptive', datafile)
  }
#
  if (!save_exists('cache/contender_output', 'sass', datafile)) {
    sass_freq <- stabsel_freq
    sass_net <- stabsel_net
    comm_counts <- cached(c('consensus_module_counts', 'sass', datafile), networks, B)
    comm_threshold <- cached(c('smooth_inner_minima', 'sass', datafile),
      uptri(comm_counts), 0, B)
    if (length(comm_threshold) != 1) {
      comm_threshold <- Inf
    }
    cmask <- comm_counts > comm_threshold
    if (any(cmask)) {
      inside_min <- cached(c('smooth_inner_minima', 'sass', 'inside', datafile),
        2*B*uptri(stabsel_freq)[uptri(cmask)], 0, 2*B)
      if (length(inside_min) == 1) {
        between_min <- cached(c('best_smooth_min', 'sass', 'between', datafile),
          2*B*uptri(stabsel_freq)[uptri(!cmask)], 0, 2*B)
        if (!is.na(between_min) && inside_min < between_min) {
          sass_freq <- cached(c('adapt_frequencies_mask2', 'sass', datafile),
            stabsel_freq, inside_min, between_min, cmask, B)
          sass_net <- cached(c('threshold_frequencies_match_density', 'sass', datafile),
            sass_freq, stabsel_net)
        }
      }
    }
    save(list(freq=sass_freq, estnet=sass_net, truenet=net,
      counts=list(all=2*B*uptri(stabsel_freq), community=uptri(comm_counts),
        highrescaled=2*B*uptri(sass_freq)[uptri(cmask)],
        lowrescaled=2*B*uptri(sass_freq)[uptri(!cmask)],
        high=2*B*uptri(stabsel_freq)[uptri(cmask)],
        low=2*B*uptri(stabsel_freq)[uptri(!cmask)])),
      'cache/contender_output', 'sass', datafile)
  }
#
#  if (!save_exists('contender_output', 'sass-cuthill', datafile)) {
#    sass_freq <- stabsel_freq
#    sass_net <- stabsel_net
#    comm_counts <- cached(c('consensus_cuthillmckee_counts', 'sass-cuthill', datafile), networks, B)
#    comm_threshold <- cached(c('smooth_inner_minima', 'sass-cuthill', datafile),
#      uptri(comm_counts), 0, B)
#    if (length(comm_threshold) != 1) {
#      comm_threshold <- Inf
#    }
#    cmask <- comm_counts > comm_threshold
#    if (any(cmask)) {
#      inside_min <- cached(c('smooth_inner_minima', 'sass-cuthill', 'inside', datafile),
#        2*B*uptri(stabsel_freq)[uptri(cmask)], 0, 2*B)
#      if (length(inside_min) == 1) {
#        between_min <- cached(c('best_smooth_min', 'sass-cuthill', 'between', datafile),
#          2*B*uptri(stabsel_freq)[uptri(!cmask)], 0, 2*B)
#        if (!is.na(between_min) && inside_min < between_min) {
#          sass_freq <- cached(c('adapt_frequencies_mask2', 'sass-cuthill', datafile),
#            stabsel_freq, inside_min, between_min, cmask, B)
#          sass_net <- cached(c('threshold_frequencies_match_density', 'sass-cuthill', datafile),
#            sass_freq, stabsel_net)
#        }
#      }
#    }
#    save(list(freq=sass_freq, estnet=sass_net, truenet=net,
#      counts=list(all=2*B*uptri(stabsel_freq), community=uptri(comm_counts),
#        highrescaled=2*B*uptri(sass_freq)[uptri(cmask)],
#        lowrescaled=2*B*uptri(sass_freq)[uptri(!cmask)],
#        high=2*B*uptri(stabsel_freq)[uptri(cmask)],
#        low=2*B*uptri(stabsel_freq)[uptri(!cmask)])),
#      'contender_output', 'sass-cuthill', datafile)
#  }
#
  networksglasso <- cached(c('subsample_glasso', datafile),
    net$data, lambda_ratio_glasso*lambda_max, B, subsamples)
  stabselglasso_freq <- cached(c('selection_freq', 'stabselglasso', datafile), networksglasso)
  threshold <- cached(c('FDR_threshold', 'stabselglasso', datafile),
    stabselglasso_freq, FDR=0.1, length(networksglasso))
  stabselglasso_net <- cached(c('threshold_frequencies', 'stabselglasso', datafile),
    stabselglasso_freq, threshold)
  if (!save_exists('cache/contender_output', 'stabselglasso', datafile)) {
    save(list(freq=stabselglasso_freq, estnet=stabselglasso_net, truenet=net),
      'cache/contender_output', 'stabselglasso', datafile)
  }
#
  if (!save_exists('cache/contender_output', 'sass-glasso', datafile)) {
    sass_freq <- stabselglasso_freq
    sass_net <- stabselglasso_net
    comm_counts <- cached(c('consensus_module_counts', 'sass-glasso', datafile), networksglasso, B)
    comm_threshold <- cached(c('smooth_inner_minima', 'sass-glasso', datafile),
      uptri(comm_counts), 0, B)
    if (length(comm_threshold) != 1) {
      comm_threshold <- Inf
    }
    cmask <- comm_counts > comm_threshold
    if (any(cmask)) {
      inside_min <- cached(c('smooth_inner_minima', 'sass-glasso', 'inside', datafile),
        2*B*uptri(stabselglasso_freq)[uptri(cmask)], 0, 2*B)
      if (length(inside_min) == 1) {
        between_min <- cached(c('best_smooth_min', 'sass-glasso', 'between', datafile),
          2*B*uptri(stabselglasso_freq)[uptri(!cmask)], 0, 2*B)
        if (!is.na(between_min) && inside_min < between_min) {
          sass_freq <- cached(c('adapt_frequencies_mask2', 'sass-glasso', datafile),
            stabselglasso_freq, inside_min, between_min, cmask, B)
          sass_net <- cached(c('threshold_frequencies_match_density', 'sass-glasso', datafile),
            sass_freq, stabselglasso_net)
        }
      }
    }
    save(list(freq=sass_freq, estnet=sass_net, truenet=net,
      counts=list(all=2*B*uptri(stabselglasso_freq), community=uptri(comm_counts),
        highrescaled=2*B*uptri(sass_freq)[uptri(cmask)],
        lowrescaled=2*B*uptri(sass_freq)[uptri(!cmask)],
        high=2*B*uptri(stabselglasso_freq)[uptri(cmask)],
        low=2*B*uptri(stabselglasso_freq)[uptri(!cmask)])),
      'cache/contender_output', 'sass-glasso', datafile)
  }
#
  if (!save_exists('cache/contender_output', 'sass-grab', datafile)) {
    sass_freq <- stabsel_freq
    sass_net <- stabsel_net
    filename <- paste('cache/grab/output/net_20/', datafile, sep='')
    if (file.exists(filename)) {
      freq <- as.matrix(read.table(filename))
      diag(freq) <- 0
      grab_net <- threshold_frequencies_match_density(freq, stabsel_net)
      lambda <- min(freq[as.matrix(grab_net)])
      filename <- paste('cache/grab/output/blocks_20/', datafile, '_', lambda, sep='')
      if (file.exists(filename)) {
        cmask <- read_grab_mask(filename)
        if (any(cmask)) {
          inside_min <- cached(c('smooth_inner_minima', 'sass-grab', 'inside', datafile),
            2*B*uptri(stabsel_freq)[uptri(cmask)], 0, 2*B)
          if (length(inside_min) == 1) {
            between_min <- cached(c('best_smooth_min', 'sass-grab', 'between', datafile),
              2*B*uptri(stabsel_freq)[uptri(!cmask)], 0, 2*B)
            if (!is.na(between_min) && inside_min < between_min) {
              sass_freq <- cached(c('adapt_frequencies_mask2', 'sass-grab', datafile),
                stabsel_freq, inside_min, between_min, cmask, B)
              sass_net <- cached(c('threshold_frequencies_match_density', 'sass-grab', datafile),
                sass_freq, stabsel_net)
            }
          }
        }
        save(list(freq=sass_freq, estnet=sass_net, truenet=net,
          counts=list(all=2*B*uptri(stabsel_freq),
            highrescaled=2*B*uptri(sass_freq)[uptri(cmask)],
            lowrescaled=2*B*uptri(sass_freq)[uptri(!cmask)],
            high=2*B*uptri(stabsel_freq)[uptri(cmask)],
            low=2*B*uptri(stabsel_freq)[uptri(!cmask)])),
          'cache/contender_output', 'sass-grab', datafile)
      }
    }
  }
#
  if (datafile == 'real') next
  if (!save_exists('cache/contender_output', 'neighsel', datafile)) {
    solutions <- cached(c('neigh_sel', datafile), net$data, seq(1, 0.1, length=100) * lambda_max)
    neighsel_net <- cached(c('solutions_match_density', datafile, 'neighsel'), solutions, stabsel_net)
    roc <- cached(c('compute_roc_solutions', datafile, 'neighsel'), net$theta, solutions)
    save(list(estnet=neighsel_net, truenet=net, roc=roc, solutions=solutions),
      'cache/contender_output', 'neighsel', datafile)
  }
  if (!save_exists('cache/contender_output', 'glasso', datafile)) {
    solutions <- cached(c('glasso', datafile), net$data, seq(1, 0.1, length=100) * lambda_max)
    glasso_net <- cached(c('solutions_match_density', datafile, 'glasso'), solutions, stabselglasso_net)
    roc <- cached(c('compute_roc_solutions', datafile, 'glasso'), net$theta, solutions)
    save(list(estnet=glasso_net, truenet=net, roc=roc, solutions=solutions),
      'cache/contender_output', 'glasso', datafile)
  }
  if (!save_exists('cache/contender_output', 'grab', datafile)) {
    for (lambda in seq(0.48, 0.1, -0.02)) {
      filename <- paste('cache/grab/output/net_20/', datafile, '_', lambda, sep='')
      if (file.exists(filename)) {
        freq <- as.matrix(read.table(filename))
        diag(freq) <- 0
        if (sum(freq) >= sum(as.matrix(stabsel_net))) {
          save(list(freq=freq, estnet=freq, truenet=net),
            'cache/contender_output', 'grab', datafile)
          break
        }
      }
    }
  }
  #if (!save_exists('contender_output', 'grab18', datafile)) {
  #  filename <- paste('experiment/grab/output/net_18/', datafile, sep='')
  #  if (file.exists(filename)) {
  #    freq <- as.matrix(read.table(filename))
  #    diag(freq) <- 0
  #    grab_net <- threshold_frequencies_match_density(freq, stabsel_net)
  #    save(list(freq=freq, estnet=grab_net, truenet=net),
  #      'contender_output', 'grab18', datafile)
  #  }
  #}
  #if (!save_exists('contender_output', 'grab22', datafile)) {
  #  filename <- paste('experiment/grab/output/net_22/', datafile, sep='')
  #  if (file.exists(filename)) {
  #    freq <- as.matrix(read.table(filename))
  #    diag(freq) <- 0
  #    grab_net <- threshold_frequencies_match_density(freq, stabsel_net)
  #    save(list(freq=freq, estnet=grab_net, truenet=net),
  #      'contender_output', 'grab22', datafile)
  #  }
  #}
}

file.create('cache/contender_output/contenders_run')
