#' @title construct H function
#' @description H is constructed based on a vector that designates the clusters
#' @param clustering (numeric) A vector of cluster assignments
#' @return (matrix) H based on the cluster assignments
construct.H <- function(clustering, cluster.max = max(clustering)) {
  clusters <- seq(max(cluster.max))
  
  if (min(clustering) < 1) {
    stop(simpleError('cluster indexing starts at 1'))
  }
  
  # find |A_k|
  cluster.sizes <- sapply(clusters, function(i) {
    length(clustering[clustering == i])
  })
  
  # construct H
  H <- sapply(clustering, function(i) {
    h <- rep(0, length(clusters))
    h[i] <- 1 / sqrt(cluster.sizes[i])
    return(h)
  }) %>% 
    t()
  
  return(H)
}

kernel.kmeans.obj <- function(K, clustering, cluster.max = max(clustering)) {
  H <- construct.H(clustering, cluster.max)
  psych::tr(t(H) %*% K %*% H)
}

ratio.cut.obj <- function(L, clustering, cluster.max = max(clustering)) {
  kernel.kmeans.obj(L, clustering, cluster.max)
}

#' @title k-means clustering with an initial clustering
#' @param X (numeric) A data matrix
#' @param clust.init (numeric) A vector of cluster assignments
#' @return The object returned by stats::kmeans
kmeans.initialized <- function(X, clust.init) {
  # dimensionality checks
  # assertthat::assert_that(nrow(X) == length(clust.init))
  assertthat::assert_that(min(clust.init) == 1)
  assertthat::assert_that(length(unique(clust.init)) == max(clust.init))
  
  # find initial centers
  if (is.matrix(X)) {
    centers.init <- sapply(unique(clust.init), function(i) {
      clust.ind <- which(clust.init == i)
      apply(X[clust.ind, ], 2, mean)
    }) %>% 
      t()
  } else {
    centers.init <- sapply(unique(clust.init), function(i) {
      clust.ind <- which(clust.init == i)
      mean(X[clust.ind])
    })
  }
  
  # apply stats::kmeans
  kmeans(X, centers.init)
}

wss.1d <- function(v, clust) {
  sapply(unique(clust), function(i) {
    v.clust <- v[clust == i]
    sum((v.clust - mean(v.clust)) ** 2)
  }) %>% 
    sum()
}

fiedler.vector.partitions <- function(L, parallel = TRUE) {
  # extract fiedler vector
  n <- nrow(L)
  L.eigen <- eigen(L)
  f <- L.eigen$vectors[, n - 1]
  
  # keep track of order
  f.order <- order(f)
  
  plyr::ldply(seq(n - 1), function(i) {
    # construct clustering
    clust <- rep(1, n)
    clust[-f.order[seq(i)]] <- 2
    clust.string <- paste0('(', paste(f.order[seq(i)], collapse = ', '), 
                           '), (', 
                           paste(f.order[-seq(i)], collapse = ', '), 
                           ')')
    
    dplyr::data_frame(clustering = clust.string, 
                      R = ratio.cut.obj(L, clust), 
                      W = wss.1d(f, clust))
  }, .parallel = parallel)
  return()
}

all.2way.partitions <- function(L, parallel = TRUE) {
  n <- nrow(L)
  K <- MASS::ginv(L)
  
  plyr::ldply(seq(floor(n / 2)), function(i) {
    partitions <- partitions::listParts(c(i, n - i))
    plyr::ldply(partitions, function(partition) {
      clustering <- rep(NA, n)
      clustering[partition[[1]]] <- 1
      clustering[partition[[2]]] <- 2
      
      dplyr::data_frame(partition = paste(as.character(partition), 
                                          collapse = ', '), 
                        W = kernel.kmeans.obj(K, clustering), 
                        R = ratio.cut.obj(L, clustering))
    })
  }, .parallel = parallel)
}

diag.matrix <- function(m) {
  # diagonal matrix of a matrix
  as.matrix(Matrix::Diagonal(nrow(m), diag(m)))
}

ECT <- function(W) {
  # find the matrix of expected commute times for weight matrix W
  
  # size of W (i.e., number of vertices)
  n <- nrow(W)
  
  # 1 vector
  e <- rep(1, n)
  
  # volume of graph
  w.tilde <- sum(W)
  
  # transition probabilities
  pi. <- (W %*% e) / w.tilde
  
  # total weight (and inverse)
  T <- Matrix::Diagonal(n, W %*% e)
  T.inv <- Matrix::Diagonal(n, (W %*% e) ** -1)
  
  # matrix of transition probabilities
  P <- T.inv %*% W
  
  # fundamental matrix of kemeny and snell
  Z <- solve(diag(n) - P + e %*% t(pi.))
  D <- w.tilde * T.inv
  
  # matrix of first passage times
  M <- as.matrix((diag(n) - Z + e %*% t(e) %*% diag.matrix(Z)) %*% D)
  
  # matrix of expected commute times
  C. <- (M - diag.matrix(M)) + t(M - diag.matrix(M))
  
  return(C.)
}

mds.kappa <- function(C) {
  #  Critchley's kappa operator on centered nxn matrices.
  #  This is the inverse of tau.
  
  n <- nrow(C)
  H <- matrix(1,nrow=n,ncol=n)
  d <- diag(C)
  H <- diag(d) %*% H
  H <- H + t(H) - 2 * C
  d <- seq(1, n^2, n + 1)
  H[d] <- 0
  return(H)
}

mds.tau <- function(H) {
  #  This function returns the double centering of the inputted matrix.
  #  See Critchley for details.
  n <- nrow(H)
  P <- diag(n) - 1/n
  return(-0.5 * P %*% H %*% P)
}

degenerate.graph.laplacian <- function(epsilon) {
  L <- rbind(c(eps, -eps, 0, 0), 
             c(-eps, 1 + eps, -1, 0),
             c(0, -1, 2, -1), 
             c(0, 0, -1, 1))
  return(L)
}

load.double.spiral <- function() {
  # file path is hard-coded :(
  W <- read.table('~/dev/ratio-cut/data/spiral-graph.txt') %>% 
    as.matrix() %>% 
    return()
}
