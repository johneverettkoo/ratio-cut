rc.sdp <- function(L, k = 2, 
                   init.clust = NA, 
                   verbose = FALSE, parallel = FALSE, cores = 4) {
  if (parallel) {
    doMC::registerDoMC(cores)
  }
  
  n <- nrow(L)
  
  Z <- CVXR::Variable(n, n, PSD = TRUE)
  objective <- CVXR::Minimize(CVXR::matrix_trace(L %*% Z))
  constraints <- list(CVXR::matrix_trace(Z) == k,
                      Z %*% rep(1, n) == rep(1, n))
  positive.constraints <- (Z >= 0)
  constraints %<>% 
    c(positive.constraints) %>%
    # c(CVXR::multiply(Z, Z) == Z) %>% 
    unname()
  problem <- CVXR::Problem(objective, constraints)
  out <- CVXR::solve(problem, parallel = parallel, verbose = verbose)
  Z.hat <- out$getValue(Z)
  
  Z.eigen <- eigen(Z.hat, symmetric = TRUE)
  embedding <- Z.eigen$vectors %>% 
    sweep(2, Z.eigen$values, '*')
  if (is.na(init.clust)) {
    clustering <- kmeans(embedding, 2)$cluster
  } else {
    clustering <- kmeans.initialized(embedding, init.clust)$cluster
  }
  return(list(Z = Z.hat, clustering = clustering))
}

kmeans.sdp <- function(K, k = 2, 
                       init.clust = NA, 
                       verbose = FALSE, parallel = FALSE, cores = 4) {
  if (parallel) {
    doMC::registerDoMC(cores)
  }
  
  n <- nrow(K)
  
  Z <- CVXR::Variable(n, n, PSD = TRUE)
  objective <- CVXR::Maximize(CVXR::matrix_trace(K %*% Z))
  constraints <- list(CVXR::matrix_trace(Z) == k,
                      Z %*% rep(1, n) == rep(1, n))
  positive.constraints <- (Z >= 0)
  constraints %<>% 
    c(positive.constraints) %>%
    # c(Z %*% Z == Z) %>% 
    unname()
  problem <- CVXR::Problem(objective, constraints)
  out <- CVXR::solve(problem, parallel = parallel, verbose = verbose)
  Z.hat <- out$getValue(Z)
  
  Z.eigen <- eigen(Z.hat, symmetric = TRUE)
  embedding <- Z.eigen$vectors %>% 
    sweep(2, Z.eigen$values, '*')
  if (is.na(init.clust)) {
    clustering <- kmeans(embedding, 2)$cluster
  } else {
    clustering <- kmeans.initialized(embedding, init.clust)$cluster
  }
  return(list(Z = Z.hat, clustering = clustering))
}
