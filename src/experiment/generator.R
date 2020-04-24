#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.generator(): Data generator                                      #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tzhao5@jhu.edu> and <hanliu@cs.jhu.edu>                      #
# Date: July 15th 2011                                                  #
# Version: 1.1.0				                				                        #
# License: GPL-2                                                        #
#-----------------------------------------------------------------------#
# Edited by Jonatan Kallus

g <- 20
prob <- 0.3
prob2 <- 0.01
u = 0.1
v = 0.3
     
cluster <- function(d) {
  theta <- matrix(0, d, d)
  tmp <- matrix(runif(d^2, 0, 0.5), d, d)
  tmp <- tmp + t(tmp)
  theta[tmp < prob2] <- 1
	g.large = d%%g
	g.small = g - g.large
	n.small = floor(d/g)
	n.large = n.small+1
	g.list = c(rep(n.small,g.small),rep(n.large,g.large))
	g.ind = rep(c(1:g),g.list)
  for(i in 1:g){
    tmp = which(g.ind==i)
    tmp2 = matrix(stats::runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
    tmp2 = tmp2 + t(tmp2)		 	
    theta[tmp,tmp][tmp2<prob] = 1
    rm(tmp,tmp2)
    gc()
  }
  return(theta)
}

hub <- function(d) {
  theta <- matrix(0, d, d)
  tmp <- matrix(runif(d^2, 0, 4), d, d)
  mask <- matrix(0, d, d)
  mask[1:g, ] <- 1
  theta[mask == 1 & tmp < prob] <- 1
  theta[mask == 0 & tmp < prob2] <- 1
  theta <- theta + t(theta)
  theta[theta == 2] <- 1
  return(theta)
}

scalefree <- function(d) {
  as.matrix(huge::huge.generator(n=100, d=d, graph='scale-free', verbose=F)$theta)
}

erdosrenyi <- function(d) {
  theta <- matrix(0, d, d)
  tmp <- matrix(runif(d^2, 0, 0.5), d, d)
  tmp <- tmp + t(tmp)
  theta[tmp < sqrt(prod(prob, prob2))] <- 1
  return(theta)
}

overlapping <- function(d) {
  mask <- matrix(0, d, d)
  for (i in 1:g) {
    ix <- (i-1)*20 + 1:20
    if (i > 1) {
      ix <- c((i-1)*20 - 0:4, ix)
    }
    mask[ix, ix] <- 1
  }
  theta <- matrix(0, d, d)
  tmp <- matrix(runif(d^2, 0, 0.5), d, d)
  tmp <- tmp + t(tmp)
  theta[tmp < prob2] <- 1
  theta[tmp < prob & mask == 1] <- 1
  return(theta)
}

generator <- function(n, graph, seed, permute=0, d=400) {	
  set.seed(1)
  graph <- gsub('-', '', graph)
  theta <- do.call(graph, list(d))
  cmask <- NULL
  if (graph == 'cluster') {
    g.large = d%%g
    g.small = g - g.large
    n.small = floor(d/g)
    n.large = n.small+1
    g.list = c(rep(n.small,g.small),rep(n.large,g.large))
    g.ind = rep(c(1:g),g.list)
    cmask <- outer(g.ind, g.ind, '==')
  }
  if (graph == 'overlapping') {
    cmask <- matrix(F, d, d)
    for (i in 1:g) {
      ix <- (i-1)*20 + 1:20
      if (i > 1) {
        ix <- c((i-1)*20 - 0:4, ix)
      }
      cmask[ix, ix] <- T
    }
  }
	diag(theta) = 0
	omega = theta*v
	
	# make omega positive definite and standardized
	diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
	sigma = stats::cov2cor(solve(omega))
	omega = solve(sigma)
	
	# generate multivariate normal data
  set.seed(seed)
	x = MASS::mvrnorm(n,rep(0,d),sigma)
  if (permute > 0) {
    for (i in 1:ceiling(n*permute)) {
      x[i, ] <- x[i, sample(d)]
    }
  }

  return(list(data=x, theta=Matrix::Matrix(theta, sparse=T), graph=graph,
    seed=seed, cmask=cmask))
	sim = list(data = x, sigma = sigma, sigmahat = sigmahat, omega = omega,
    theta = Matrix::Matrix(theta,sparse = TRUE), sparsity= sum(theta)/(d*(d-1)),
    graph.type=graph, masks=masks, cl=cl)
	class(sim) = "sim" 
	return(sim)
	if(is.null(g)){
		g = 1
		if(graph == "hub" || graph == "cluster"){
			if(d > 40)	g = ceiling(d/20)
			if(d <= 40) g = 2
		}
    if(graph == "clustersofclusters") {
      if(d > 80)  g = ceiling(d/20)
      if(d <= 80) g = 4
    }
	}
	
	if(graph == "random"){
		if(is.null(prob))	prob = min(1, 3/d)
		prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
	}
	
	if(graph == "cluster"){
		if(is.null(prob)){
			if(d/g > 30)	prob = 0.3
			if(d/g <= 30)	prob = min(1,6*g/d)
		}
		prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
	}  
	 	
	if(graph == "clustersofclusters"){
		if(is.null(prob)){
			if(d/g > 60)	prob = 0.3
			if(d/g <= 60)	prob = min(1,6*g/d)
		}
		prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
	}  
	 	
	
	# parition variables into groups
	g.large = d%%g
	g.small = g - g.large
	n.small = floor(d/g)
	n.large = n.small+1
	g.list = c(rep(n.small,g.small),rep(n.large,g.large))
	g.ind = rep(c(1:g),g.list)
	rm(g.large,g.small,n.small,n.large,g.list)
	gc()
	
	# build the graph structure
	theta = matrix(0,d,d);
	if(graph == "band"){
		if(is.null(u)) u = 0.1
		if(is.null(v)) v = 0.3
		for(i in 1:g){
			diag(theta[1:(d-i),(1+i):d]) = 1
			diag(theta[(1+i):d,1:(d-1)]) = 1
		}	
	}
	if(graph == "cluster"){
		if(is.null(u)) u = 0.1
		if(is.null(v)) v = 0.3
    cl <- g.ind
		for(i in 1:g){
		 	tmp = which(g.ind==i)
		 	tmp2 = matrix(stats::runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
		 	tmp2 = tmp2 + t(tmp2)		 	
		 	theta[tmp,tmp][tmp2<prob] = 1
		 	rm(tmp,tmp2)
		 	gc()
		}
	}
  masks <- list()
	if(graph == "overlapping"){
    if (d != 390) {
      stop('set d to 390')
    }
    tmp <- matrix(runif(d^2, 0, 0.5), d, d)
    tmp <- tmp + t(tmp)
    cl <- rep(0, d)
    clust <- matrix(F, d, d)
    for (i in 0:29) {
      ix <- (i*13):min(d, (i*13 + 14))
      clust[ix, ix] <- T
    }
    masks[[1]] <- clust
    tmp[!clust] <- 1
    theta[tmp < prob] <- 1
		if(is.null(u)) u = 0.1
		if(is.null(v)) v = 0.3
	}
	if(graph == "clustersofclusters"){
		if(is.null(u)) u = 0.1
		if(is.null(v)) v = 0.3
    theta = matrix(stats::runif(d^2,0,0.5), d, d)
    theta = 1 * (theta + t(theta) < prob^4)
    tmp = which(g.ind<=round(g/2))
    cl <- matrix(c(g.ind, 1 + g.ind > round(g/2)), length(g.ind), 2)
    tmp2 = matrix(stats::runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
    tmp2 = tmp2 + t(tmp2)		 	
    theta[tmp,tmp] = 1 * (tmp2<prob^2.5)
    masks$between <- theta > 1
    masks$between[tmp, tmp] <- TRUE
    rm(tmp,tmp2)
    gc()
    tmp = which(g.ind>round(g/2))
    tmp2 = matrix(stats::runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
    tmp2 = tmp2 + t(tmp2)		 	
    theta[tmp,tmp] = 1 * (tmp2<prob^2.5)
    masks$between[tmp, tmp] <- TRUE
    rm(tmp,tmp2)
    gc()
    masks$inside <- theta > 1
		for(i in 1:g){
		 	tmp = which(g.ind==i)
		 	tmp2 = matrix(stats::runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
		 	tmp2 = tmp2 + t(tmp2)		 	
		 	theta[tmp,tmp] = 1 * (tmp2<prob)
      masks$inside[tmp, tmp] <- TRUE
		 	rm(tmp,tmp2)
		 	gc()
		}
    masks$between <- masks$between & !masks$inside
    masks$outside <- !masks$inside & !masks$between
    diag(masks$inside) <- FALSE
	}
	if(graph == "hub"){
		if(is.null(u)) u = 0.1
		if(is.null(v)) v = 0.3
		for(i in 1:g){
		 	tmp = which(g.ind==i)
		 	theta[tmp[1],tmp] = 1
		 	theta[tmp,tmp[1]] = 1
		 	rm(tmp)
		 	gc()
		}
	}
	if(graph == "random"){
		if(is.null(u)) u = 0.1
		if(is.null(v)) v = 0.3
		
		tmp = matrix(stats::runif(d^2,0,0.5),d,d)
		tmp = tmp + t(tmp)
		theta[tmp < prob] = 1
		#theta[tmp >= tprob] = 0
		rm(tmp)
		gc()
	}
	
	if(graph == "scale-free"){
		if(is.null(u)) u = 0.1
		if(is.null(v)) v = 0.3
    stop('Not implemented')
		#out = .C("SFGen",dd0=as.integer(2),dd=as.integer(d),G=as.integer(theta),PACKAGE="huge")
		#theta = matrix(as.numeric(out$G),d,d)
	}
	diag(theta) = 0
	omega = theta*v
	
	# make omega positive definite and standardized
	diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
	sigma = stats::cov2cor(solve(omega))
	omega = solve(sigma)
	
	# generate multivariate normal data
	x = MASS::mvrnorm(n,rep(0,d),sigma)
	sigmahat = stats::cor(x)
	
	# graph and covariance visulization
	if(vis == TRUE){
    stop('Not implemented')
	}
	if(verbose) cat("done.\n")
	rm(vis,verbose)
	gc()
	
	sim = list(data = x, sigma = sigma, sigmahat = sigmahat, omega = omega,
    theta = Matrix::Matrix(theta,sparse = TRUE), sparsity= sum(theta)/(d*(d-1)),
    graph.type=graph, masks=masks, cl=cl)
	class(sim) = "sim" 
	return(sim)
}
