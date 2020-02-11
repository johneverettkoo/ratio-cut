rc.sdp <- function(L, k = 2) {
  n <- nrow(L)
  
  Z <- Variable(n, n, PSD = TRUE)
  objective <- Minimize(matrix_trace(L %*% Z))
  constraints <- list(matrix_trace(Z) == k,
                      Z %*% rep(1, n) == rep(1, n))
  positive.constraints <- plyr::alply(expand.grid(i = seq(n), j = seq(n)), 1,
                                      function(row) {
                                        if (row$i <= row$j) {
                                          Z[row$i, row$j] >= 0
                                        }
                                      }) %>% 
    plyr::compact()
  constraints <- c(constraints, positive.constraints) %>% 
    unname()
  problem <- Problem(objective, constraints)
  out <- solve(problem, parallel = FALSE)
  Z.hat <- out$getValue(Z)
  
  Z.eigen <- eigen(Z.hat, symmetric = TRUE)
  embedding <- Z.eigen$vectors %>% 
    sweep(2, Z.eigen$values, '*')
  clustering <- kmeans(embedding, 2)$cluster
  return(clustering)
}

kmeans.sdp <- function(K, k = 2) {
  n <- nrow(K)
  
  Z <- Variable(n, n, PSD = TRUE)
  objective <- Maximize(matrix_trace(K %*% Z))
  constraints <- list(matrix_trace(Z) == k,
                      Z %*% rep(1, n) == rep(1, n))
  positive.constraints <- plyr::alply(expand.grid(i = seq(n), j = seq(n)), 1,
                                      function(row) {
                                        if (row$i <= row$j) {
                                          Z[row$i, row$j] >= 0
                                        }
                                      }) %>% 
    plyr::compact()
  constraints <- c(constraints, positive.constraints) %>% 
    unname()
  problem <- Problem(objective, constraints)
  out <- solve(problem, parallel = FALSE)
  Z.hat <- out$getValue(Z)
  
  Z.eigen <- eigen(Z.hat, symmetric = TRUE)
  embedding <- Z.eigen$vectors %>% 
    sweep(2, Z.eigen$values, '*')
  clustering <- kmeans(embedding, 2)$cluster
  return(clustering)
}
